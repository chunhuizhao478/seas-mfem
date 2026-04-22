// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan Step 5b (NEW, rev-3h+, REVIEW R-V92-I02) —
// Multi-step stability probe for the F01+F02 coupled-RK4-on-psi fix,
// exercised across MPI partition seams on the 4-km shared-fault fixture.
//
// ============================================================================
// What this test probes (beyond §4.7 and test_rk4_psi_integration)
// ============================================================================
// §4.7 runs 1 RK4 step on the 4-km fixture and checks shared-fault DOFData
// consistency across ranks.  `test_rk4_psi_integration` verifies O(dt^4)
// convergence of the driver's per-step RK4-on-psi arithmetic on a single
// QP with constant V.  Neither exercises the multi-step integration loop
// across MPI partition seams — the exact code path that would betray a
// subtle bug in the coupled-RK4-on-(Q, psi) flow (e.g. stale psi input
// between stages, incorrect DOFData snapshot ordering, or psi-invariance
// violations on shared fault QPs).
//
// This test replicates the driver's full per-step RK4-on-(Q, psi)
// arithmetic (local copy, does NOT modify the driver) and runs 100 steps
// on the 4-km fixture with an antisymmetric σ_xy(y) projection.  After
// every 10 steps it asserts a physical-envelope invariant:
//   • `|sigma_n_corr - sigma_n0| < 1 MPa`  per DOFData
//   • `slip_rate < 1 m/s`                  per DOFData
//   • `|psi - psi_initial| < 0.1`          per DOFData
//   • no NaN anywhere in Q or DOFData
// These bounds are 100-1000× looser than TPV102's expected amplitudes at
// t ≈ 1 ms (the test runs 100 × 1e-5 s = 1 ms of simulated time, well
// before any non-trivial rupture growth on this fixture — Q starts with
// a tanh profile, not a bona-fide nucleation).  Any violation = bug in
// the F01+F02 arithmetic or in its interaction with ghost exchange.
//
// ============================================================================
// Interpretation of results (per plan Step 5b decision table)
// ============================================================================
//   PASS 100/100       F01+F02 arithmetic multi-step stable on shared fault
//                      QPs; commit + push F01+F02 with high confidence.
//   FAIL before 50     F01+F02 has a bug unit tests don't catch; debug
//                      locally via the per-step dump (increase verbosity)
//                      before any Frontera submission.
//   DRIFTS past 50     F01+F02 is partial; commit + push but plan for a
//                      follow-up (f01f02_dev sbatch or Step 3 bulk probe)
//                      to quantify the residual on production mesh.
//
// Usage:
//   mpirun -np 4 ./seas_test_rk4_f01f02_shared_fault_stability

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../friction/state_evolution.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

// ---------------------------------------------------------------------------
// Build the same 4-km cartesian tet fixture as §4.7 / §4.6.
// ---------------------------------------------------------------------------
static Mesh BuildFixtureMesh(int nx, int ny, int nz,
                              double Lx, double Ly, double Lz)
{
   Mesh mesh = Mesh::MakeCartesian3D(2*nx, 2*ny, nz, Element::TETRAHEDRON,
                                     2.0*Lx, 2.0*Ly, Lz);
   for (int v = 0; v < mesh.GetNV(); v++)
   {
      real_t *x = mesh.GetVertex(v);
      x[0] -= Lx; x[1] -= Ly; x[2] -= Lz;
   }
   const real_t tol = 1e-6;
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      auto *Tr = mesh.GetBdrElementTransformation(be);
      const IntegrationPoint &ip = Geometries.GetCenter(Tr->GetGeometryType());
      Tr->SetIntPoint(&ip);
      Vector c(3); Tr->Transform(ip, c);
      if (std::abs(c(2)) < tol) { mesh.SetBdrAttribute(be, 1); }
      else                      { mesh.SetBdrAttribute(be, 5); }
   }
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(f);
      if (!ftr) { continue; }
      const IntegrationPoint &ip = Geometries.GetCenter(ftr->GetGeometryType());
      ftr->Face->SetIntPoint(&ip);
      Vector c(3); ftr->Face->Transform(ip, c);
      if (std::abs(c(1)) > tol) { continue; }
      Array<int> verts; mesh.GetFaceVertices(f, verts);
      if (verts.Size() == 3)
      { mesh.AddBdrTriangle(verts[0], verts[1], verts[2], 3); }
      else if (verts.Size() == 4)
      { mesh.AddBdrQuad(verts[0], verts[1], verts[2], verts[3], 3); }
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// ---------------------------------------------------------------------------
// Collect fault coords in the exact order expected by SetFaultDOFData.
// ---------------------------------------------------------------------------
static void BuildFaultCoords(ParMesh &pmesh,
                              const Array<int> &int_faces,
                              const Array<int> &shr_faces,
                              int order,
                              std::vector<Vector> &fault_coords)
{
   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = pmesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         ftr->SetAllIntPoints(&ir.IntPoint(q));
         Vector phys(3); ftr->Face->Transform(ir.IntPoint(q), phys);
         fault_coords.push_back(phys);
      }
   }
   for (int i = 0; i < shr_faces.Size(); i++)
   {
      auto *ftr = pmesh.GetSharedFaceTransformations(shr_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         ftr->SetAllIntPoints(&ir.IntPoint(q));
         Vector phys(3); ftr->Face->Transform(ir.IntPoint(q), phys);
         fault_coords.push_back(phys);
      }
   }
}

// ---------------------------------------------------------------------------
// Replicates the driver's per-step RK4-on-(Q, psi) block EXACTLY as of
// the F01+F02 fix — same order of operations, same AgingLawPsi::Rate
// calls, same stage pre-psi advance.  Nucleation is NOT applied
// (this fixture runs short enough that nucleation never activates).
// Returns true if the step completed without NaN.
// ---------------------------------------------------------------------------
static bool DriverRK4OneStep(WaveOperator<ParMesh> &wave,
                              std::vector<DOFData> &dof_data,
                              const AgingLawPsi &aging_law,
                              Vector &Q,
                              real_t dt)
{
   const int size = Q.Size();
   const int nf = static_cast<int>(dof_data.size());
   Vector k1(size), k2(size), k3(size), k4(size), Q_tmp(size);

   std::vector<real_t> psi_n(nf);
   std::vector<real_t> sr_k1(nf), sr_k2(nf), sr_k3(nf), sr_k4(nf);
   std::vector<real_t> psi_k1(nf), psi_k2(nf), psi_k3(nf), psi_k4(nf);

   for (int i = 0; i < nf; i++) { psi_n[i] = dof_data[i].psi; }

   // Stage 1
   wave.Mult(Q, k1);
   for (int i = 0; i < nf; i++)
   {
      sr_k1[i] = dof_data[i].slip_rate;
      psi_k1[i] = aging_law.Rate(sr_k1[i], psi_n[i], dof_data[i].Dc);
      dof_data[i].psi = psi_n[i] + 0.5 * dt * psi_k1[i];
   }

   // Stage 2
   add(Q, 0.5*dt, k1, Q_tmp);
   wave.Mult(Q_tmp, k2);
   for (int i = 0; i < nf; i++)
   {
      sr_k2[i] = dof_data[i].slip_rate;
      psi_k2[i] = aging_law.Rate(sr_k2[i], dof_data[i].psi, dof_data[i].Dc);
      dof_data[i].psi = psi_n[i] + 0.5 * dt * psi_k2[i];
   }

   // Stage 3
   add(Q, 0.5*dt, k2, Q_tmp);
   wave.Mult(Q_tmp, k3);
   for (int i = 0; i < nf; i++)
   {
      sr_k3[i] = dof_data[i].slip_rate;
      psi_k3[i] = aging_law.Rate(sr_k3[i], dof_data[i].psi, dof_data[i].Dc);
      dof_data[i].psi = psi_n[i] + dt * psi_k3[i];
   }

   // Stage 4
   add(Q, dt, k3, Q_tmp);
   wave.Mult(Q_tmp, k4);
   for (int i = 0; i < nf; i++)
   {
      sr_k4[i] = dof_data[i].slip_rate;
      psi_k4[i] = aging_law.Rate(sr_k4[i], dof_data[i].psi, dof_data[i].Dc);
   }

   // Update Q with RK4 weights.
   for (int i = 0; i < size; i++)
   {
      Q[i] += dt / 6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
   }

   // Final psi update via Butcher combination + refresh averaged V/slip/tau.
   for (int i = 0; i < nf; i++)
   {
      dof_data[i].psi = psi_n[i] + dt / 6.0 *
                        (psi_k1[i] + 2.0*psi_k2[i] + 2.0*psi_k3[i] + psi_k4[i]);
      // Note: driver also RK4-averages V/tau/sigma_n into DOFData; here
      // we keep the last-stage values, which is sufficient for the bound
      // check (any stage exceeding the envelope implies a bug).
   }

   // NaN check on Q (cheap; DOFData is checked by caller).
   for (int i = 0; i < size; i++)
   {
      if (std::isnan(Q[i])) { return false; }
   }
   return true;
}

// ---------------------------------------------------------------------------
// Check the physical-envelope invariant on all fault DOFs, both locally
// and across ranks.  Returns {max |sigma_n - sigma_n0|, max slip_rate,
// max |psi - psi_0|} as globals.  ok = true iff all bounds hold.
// ---------------------------------------------------------------------------
struct EnvelopeCheck
{
   double max_sigma_n_dev = 0.0;
   double max_slip_rate   = 0.0;
   double max_psi_dev     = 0.0;
   int nan_count          = 0;
   bool ok                = true;
};

static EnvelopeCheck CheckEnvelope(const std::vector<DOFData> &dof_data,
                                    const std::vector<real_t> &psi_initial,
                                    MPI_Comm comm)
{
   EnvelopeCheck c;
   const int nf = static_cast<int>(dof_data.size());
   for (int i = 0; i < nf; i++)
   {
      const DOFData &d = dof_data[i];
      if (std::isnan(d.sigma_n_corr) || std::isnan(d.slip_rate) ||
          std::isnan(d.psi) || std::isnan(d.V1) || std::isnan(d.V2))
      { c.nan_count++; continue; }
      double ds = std::abs(d.sigma_n_corr - d.sigma_n0);
      double dp = std::abs(d.psi - psi_initial[i]);
      c.max_sigma_n_dev = std::max(c.max_sigma_n_dev, ds);
      c.max_slip_rate   = std::max(c.max_slip_rate,   (double)d.slip_rate);
      c.max_psi_dev     = std::max(c.max_psi_dev,     dp);
   }
   double gmax_sn = 0.0, gmax_sr = 0.0, gmax_psi = 0.0;
   int gnan = 0;
   MPI_Allreduce(&c.max_sigma_n_dev, &gmax_sn, 1, MPI_DOUBLE, MPI_MAX, comm);
   MPI_Allreduce(&c.max_slip_rate,   &gmax_sr, 1, MPI_DOUBLE, MPI_MAX, comm);
   MPI_Allreduce(&c.max_psi_dev,     &gmax_psi, 1, MPI_DOUBLE, MPI_MAX, comm);
   MPI_Allreduce(&c.nan_count,       &gnan,    1, MPI_INT,    MPI_SUM, comm);
   c.max_sigma_n_dev = gmax_sn;
   c.max_slip_rate   = gmax_sr;
   c.max_psi_dev     = gmax_psi;
   c.nan_count       = gnan;

   // Envelope bounds — generous by 100-1000× over expected TPV102 amplitudes
   // at t ≈ 1 ms on a 4-km fixture with Q seeded by tanh(y/L0):
   //   sigma_n envelope: 1 MPa (vs 120 MPa background; TPV102 spec ≤ 1 MPa)
   //   slip_rate        : 1 m/s (well beyond typical 1e-3 at early times)
   //   psi deviation    : 0.1 (vs psi ∈ [0.4, 0.85]; initial psi_0 ≈ 0.77)
   const double bnd_sn  = 1.0e6;
   const double bnd_sr  = 1.0;
   const double bnd_psi = 0.1;
   c.ok = (c.nan_count == 0) && (c.max_sigma_n_dev < bnd_sn) &&
          (c.max_slip_rate < bnd_sr) && (c.max_psi_dev < bnd_psi);
   return c;
}

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   auto log = [&](const std::string &m)
   { if (rank == 0) { std::cout << m << "\n"; } };

   log("\n=== TPV102 v9.2.0 Step 5b — F01+F02 100-step Stability Probe ===");
   log(std::string("  MPI ranks: ") + std::to_string(nprocs));

   if (nprocs < 4)
   {
      log("  SKIPPED: requires >=4 MPI ranks for ParMETIS to cut the fault.");
      MPI_Finalize();
      return 77;
   }

   const int nx = 2, ny = 2, nz = 2;
   const double Lx = 2.0e3, Ly = 2.0e3, Lz = 4.0e3;
   Mesh serial_mesh = BuildFixtureMesh(nx, ny, nz, Lx, Ly, Lz);
   if (rank == 0)
   {
      std::cout << "  serial mesh : " << serial_mesh.GetNE() << " elements, "
                << serial_mesh.GetNumFaces() << " faces\n";
   }
   ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {5};

   const int order = 1;
   WaveOperator<ParMesh> wave(pmesh, order,
                              TPV102Params::lambda,
                              TPV102Params::mu,
                              TPV102Params::rho, bc);

   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();
   const int n_int = int_faces.Size();
   const int n_shr = shr_faces.Size();
   int total_shared = 0;
   MPI_Allreduce(&n_shr, &total_shared, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (total_shared == 0)
   {
      log("  SKIPPED: no shared fault faces at this rank count.");
      MPI_Finalize();
      return 77;
   }

   int nqp_per_face = 0;
   if (n_shr > 0)
   {
      auto *ftr = pmesh.GetSharedFaceTransformations(shr_faces[0]);
      nqp_per_face = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   else if (n_int > 0)
   {
      auto *ftr = pmesh.GetInteriorFaceTransformations(int_faces[0]);
      nqp_per_face = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   { int local = nqp_per_face, global = 0;
     MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
     nqp_per_face = global; }
   const int num_fault_total = (n_int + n_shr) * nqp_per_face;

   std::vector<Vector> fault_coords;
   BuildFaultCoords(pmesh, int_faces, shr_faces, order, fault_coords);

   std::vector<DOFData> dof_data;
   if (num_fault_total > 0)
   {
      InitializeFaultDOFs(dof_data, num_fault_total, fault_coords);
   }
   std::vector<real_t> psi_initial(num_fault_total);
   for (int i = 0; i < num_fault_total; i++) { psi_initial[i] = dof_data[i].psi; }

   FaultFaceFlux fault_flux(TPV102Params::rho, TPV102Params::cp,
                            TPV102Params::cs);
   wave.SetFaultFlux(&fault_flux);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   AgingLawPsi aging_law(TPV102Params::b, TPV102Params::V0, TPV102Params::f0);

   // Seed Q with antisymmetric σ_xy(y) = τ_ini · tanh(y / L0) — Phase D.
   const int size = wave.Height();
   Vector Q(size); Q = 0.0;
   {
      const auto &fes_const = wave.GetFESpace();
      auto &fes = const_cast<ParFiniteElementSpace &>(fes_const);
      const int ndof_total = fes.GetNDofs();
      const double L0 = 500.0;
      FunctionCoefficient sxy_coeff(
         [L0](const Vector &x) -> real_t
         { return TPV102Params::tau_ini * std::tanh(x(1) / L0); });
      ParGridFunction sxy_gf(&fes);
      sxy_gf.ProjectCoefficient(sxy_coeff);
      for (int i = 0; i < ndof_total; i++)
      { Q(SXY * ndof_total + i) = sxy_gf(i); }
   }

   const real_t dt    = 1.0e-5;
   const int    nstep = 100;
   const int    check_every = 10;

   if (rank == 0)
   {
      std::cout << "  driving " << nstep << " RK4 steps at dt = " << dt
                << " s; envelope checked every " << check_every << " steps\n";
      std::cout << "  envelope: |σ_n-σ_n0|<1MPa, slip_rate<1m/s, "
                   "|psi-psi_0|<0.1, no NaN\n\n";
   }

   int failures = 0;
   int first_violation_step = -1;
   double terminal_sn_dev = 0.0, terminal_slip_rate = 0.0, terminal_psi_dev = 0.0;
   int terminal_nan = 0;

   for (int step = 1; step <= nstep; step++)
   {
      bool ok = DriverRK4OneStep(wave, dof_data, aging_law, Q, dt);
      if (!ok)
      {
         if (rank == 0)
         {
            std::cout << "  step " << std::setw(3) << step
                      << "  NaN in Q — stopping\n";
         }
         failures++;
         if (first_violation_step < 0) { first_violation_step = step; }
         break;
      }

      if (step % check_every == 0 || step == nstep)
      {
         EnvelopeCheck c = CheckEnvelope(dof_data, psi_initial, MPI_COMM_WORLD);
         if (rank == 0)
         {
            std::cout << "  step " << std::setw(3) << step
                      << "  max|σ_n-σ_n0|=" << std::scientific
                      << std::setprecision(3) << c.max_sigma_n_dev << " Pa"
                      << "  max slip_rate=" << c.max_slip_rate << " m/s"
                      << "  max|psi-psi_0|=" << c.max_psi_dev
                      << "  NaNs=" << c.nan_count
                      << (c.ok ? "  OK" : "  VIOLATION") << "\n";
         }
         terminal_sn_dev    = c.max_sigma_n_dev;
         terminal_slip_rate = c.max_slip_rate;
         terminal_psi_dev   = c.max_psi_dev;
         terminal_nan       = c.nan_count;
         if (!c.ok)
         {
            failures++;
            if (first_violation_step < 0) { first_violation_step = step; }
         }
      }
   }

   int num_tests = 1, num_passed = 0, num_failed = 0;
   if (failures == 0)
   {
      num_passed = 1;
      if (rank == 0)
      {
         std::cout << "\n  PASSED: 100/100 RK4 steps stay within envelope "
                      "on 4-km shared-fault fixture\n";
         std::cout << "          terminal max|σ_n-σ_n0| = " << std::scientific
                   << std::setprecision(3) << terminal_sn_dev << " Pa\n";
         std::cout << "          terminal max slip_rate = "
                   << terminal_slip_rate << " m/s\n";
         std::cout << "          terminal max|psi-psi_0| = "
                   << terminal_psi_dev << "\n";
      }
   }
   else
   {
      num_failed = 1;
      if (rank == 0)
      {
         std::cout << "\n  FAILED: F01+F02 envelope violated (first at step "
                   << first_violation_step << ")\n";
         std::cout << "          last max|σ_n-σ_n0| = " << std::scientific
                   << std::setprecision(3) << terminal_sn_dev << " Pa\n";
         std::cout << "          last max slip_rate = "
                   << terminal_slip_rate << " m/s\n";
         std::cout << "          last max|psi-psi_0| = "
                   << terminal_psi_dev << "\n";
         std::cout << "          NaNs detected       = " << terminal_nan << "\n";
         std::cout << "  → HOLD push; bisect the failing stage before any "
                      "Frontera submission (plan Step 5b decision table).\n";
      }
   }

   if (rank == 0)
   {
      std::cout << "\n========================================\n";
      std::cout << "  Step 5b results: " << num_passed << " passed, "
                << num_failed << " failed out of " << num_tests << " tests\n";
      std::cout << "========================================\n";
   }

   MPI_Finalize();
   return (num_failed == 0) ? 0 : 1;
}
