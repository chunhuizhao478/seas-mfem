// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/unit/test_lsw_rk_shared_fault_mpi.cpp — Phase 14 RK + LSW, test L8.
// Plan: document/mixed_flux_dev/PLAN_rk45_lsw_mixed_flux_2026-05-29.md.
//
// L8 — MPI shared-fault LSW dispatch parity.  The serial unit tests L1–L7 only
// exercise the INTERIOR-fault Mult dispatch (wave_operator.inl :2617); the
// SHARED-fault dispatch (:3260) is a parallel-only path that a wrong/missing
// LSW arm there would leave untested (it would compile, pass `make test`, and
// fail only at np>1 — every production run).  This test forces the fault onto a
// rank seam (a 2-tet mesh, one tet per rank) so the shared-fault EvaluateLSW arm
// fires, and asserts it produces the SAME fault observables as the serial
// interior-fault path.
//
//   Setup (both serial + parallel): a TPV205 LSW nucleation-patch QP
//   (tau_nuc = 81.6 MPa, sigma_n = 120 MPa, mu_s = 0.677), fluctuation-Q at rest
//   (Q = 0).  The instantaneous EvaluateLSW closed form gives, on the FIRST
//   Mult, V = (tau_nuc − mu_s·sigma_n)/eta_s ≈ 0.0779 m/s at every fault QP.
//
//   Assertions:
//     1. Serial interior-fault Mult: max fault slip_rate == V_analytic.
//     2. Parallel shared-fault Mult: global-max (MPI_MAX) slip_rate == V_analytic
//        AND == the serial value to round-off (the :3260 LSW arm gives the same
//        physics as :2617).
//     3. The parallel mesh actually has shared fault QPs (so :3260 is exercised).
//     4. All dQdt + slip_rate finite & non-negative on every rank.
//
// Usage:  mpirun -np 2 ./seas_test_lsw_rk_shared_fault_mpi   (returns 77 if np<2)

#include "mfem.hpp"

#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv205_friction.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../config/tpv205_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

#include <mpi.h>

using namespace mfem;
using namespace mfem::seas;

namespace
{

// Two-tet mesh with a single interior fault face at y = 0 (one tet each side),
// fault attr = 3, free = 1.  Verbatim from test_interior_vs_shared_branch_live.
Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
      {0.0,  1.0, 0.0}, {0.0, -1.0, 0.0},
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);   // tet 0: − side
   mesh.AddTet(0, 1, 2, 3, 1);   // tet 1: + side
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy) < 1e-10) { mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3); }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Overwrite a fault QP to a TPV205 LSW nucleation-patch QP, preserving the
// impedances InitializeFaultDOFs set (TPV102 == TPV205 material).
void MakeQPLSW(DOFData &d)
{
   d.sigma_n0 = TPV205Params::sigma_n;     // 120 MPa
   d.tau1_0   = 0.0;
   d.tau2_0   = TPV205Params::tau_nuc;     // 81.6 MPa
   d.sigma_n_nuc = 0.0; d.tau1_nuc = 0.0; d.tau2_nuc = 0.0;
   d.lsw_mu_s = TPV205Params::mu_s; d.lsw_mu_d = TPV205Params::mu_d;
   d.lsw_d_c  = TPV205Params::d_c;
   d.a = 0.0; d.psi = 0.0; d.Dc = 0.0;
   d.slip1 = 0.0; d.slip2 = 0.0; d.slip_rate = 0.0; d.V1 = 0.0; d.V2 = 0.0;
}

// Templated fault setup (interior for Mesh; interior + shared for ParMesh),
// then overwrite every QP to the LSW nucleation patch (fluctuation-Q: prestress
// lives in DOFData, Q stays 0).  `ff` must outlive the wave.Mult call.
template <typename MeshT>
int SetupLSWFaultMPI(WaveOperator<MeshT> &wave, MeshT &mesh, int order,
                     std::vector<DOFData> &dof_data, FaultFaceFlux &ff)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();

   int nqp = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      MFEM_VERIFY(ftr, "interior fault face FTR null");
      nqp = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>)
   {
      if (nqp == 0 && shr_faces.Size() > 0)
      {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[0]);
         MFEM_VERIFY(ftr, "shared fault face FTR null");
         nqp = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
      }
   }
#endif

   std::vector<Vector> fault_coords;
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>)
   {
      for (int i = 0; i < shr_faces.Size(); i++)
      {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[i]);
         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);
            Vector phys(3); ftr->Face->Transform(ip, phys);
            fault_coords.push_back(phys);
         }
      }
   }
#endif

   const int n_fault = (int_faces.Size() + shr_faces.Size()) * nqp;
   if (n_fault > 0)
   {
      InitializeFaultDOFs(dof_data, n_fault, fault_coords);
      for (auto &d : dof_data) { MakeQPLSW(d); }
   }
   wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);
   return n_fault;
}

bool AllFinite(const Vector &v)
{
   for (int i = 0; i < v.Size(); i++) { if (!std::isfinite(v(i))) { return false; } }
   return true;
}

// Add a VX perturbation to the −y-side element ONLY (the side whose tet is on
// rank 0 under the explicit partition).  This creates a fault-normal velocity
// jump across the seam, so the shared-fault Q± exchange (the thing the :3260
// LSW dispatch consumes) is genuinely exercised — unlike the at-rest (Q=0) case
// where the trial traction is identically 0 on both sides.  Templated so the
// serial Mesh and the ParMesh apply the SAME one-sided kick (on ParMesh each
// rank sees only its local elements, so only rank 0's −side tet is kicked).
template <typename MeshT>
void KickMinusSideVX(MeshT &mesh, const FiniteElementSpace &fes, Vector &Q,
                     real_t val)
{
   const int ndof = fes.GetNDofs();
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      real_t cy = 0; Array<int> ev; mesh.GetElementVertices(e, ev);
      for (int v = 0; v < ev.Size(); v++) { cy += mesh.GetVertex(ev[v])[1]; }
      if (cy / ev.Size() >= 0.0) { continue; }   // −y side only
      Array<int> ed; fes.GetElementDofs(e, ed);
      for (int j = 0; j < ed.Size(); j++) { Q(VX * ndof + ed[j]) += val; }
   }
}

} // namespace

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if (rank == 0)
   {
      std::cout << "\n=== L8: RK+LSW shared-fault Mult dispatch parity (np="
                << nprocs << ") ===\n";
   }
   if (nprocs < 2)
   {
      if (rank == 0) { std::cout << "  SKIPPED: requires np >= 2\n"; }
      MPI_Finalize();
      return 77;
   }

   const int order = 1;
   const real_t V_analytic =
      (TPV205Params::tau_nuc - TPV205Params::mu_s * TPV205Params::sigma_n)
      / TPV205Params::eta_s;

   int num_failed = 0;
   auto check = [&](bool cond, const char *msg)
   {
      if (rank == 0)
      {
         std::cout << (cond ? "  PASSED: " : "  FAILED: ") << msg << "\n";
      }
      if (!cond) { num_failed++; }
   };

   // ---- Serial interior-fault reference (rank 0 only) ----
   const real_t V_kick = 1.0e-5;   // one-sided VX perturbation (cross-rank jump)
   real_t V_serial = 0.0, V_serial_kick = 0.0;
   bool serial_finite = true;
   if (rank == 0)
   {
      std::cout << "\n-- Serial (interior-fault branch, :2617) --\n";
      Mesh mesh = BuildTwoTetFaultMesh();
      BoundaryConfig bc; bc.natural_attrs = {1}; bc.fault_attr = 3;
      bc.absorbing_attrs = {};
      WaveOperator<Mesh> wave(mesh, order, TPV205Params::lambda,
                              TPV205Params::mu, TPV205Params::rho, bc);
      real_t zero_bg[NUM_STATE] = {0};
      wave.SetAbsorbingBackground(zero_bg);
      FaultFaceFlux ff(TPV205Params::rho, TPV205Params::cp, TPV205Params::cs);
      std::vector<DOFData> dof_data;
      const int n_fault = SetupLSWFaultMPI<Mesh>(wave, mesh, order, dof_data, ff);
      MFEM_VERIFY(n_fault > 0, "serial: expected interior fault QPs");

      Vector Q(wave.Height()); Q = 0.0;
      Vector k(wave.Height());
      wave.Mult(Q, k);                       // (a) at-rest
      serial_finite = AllFinite(k);
      for (const auto &d : dof_data)
      {
         if (!std::isfinite(d.slip_rate) || d.slip_rate < 0.0) { serial_finite = false; }
         V_serial = std::max(V_serial, d.slip_rate);
      }
      // (b) cross-rank jump: VX kick on the −side only, then re-solve.  Mult is
      // slip-stateless, so dof_data.slip stays 0 and a fresh solve is valid.
      KickMinusSideVX(mesh, wave.GetFESpace(), Q, V_kick);
      wave.Mult(Q, k);
      if (!AllFinite(k)) { serial_finite = false; }
      for (const auto &d : dof_data)
      {
         if (!std::isfinite(d.slip_rate) || d.slip_rate < 0.0) { serial_finite = false; }
         V_serial_kick = std::max(V_serial_kick, d.slip_rate);
      }
      std::cout << "    serial max slip_rate (at-rest) = " << std::scientific
                << std::setprecision(6) << V_serial
                << "  (analytic " << V_analytic << ")\n"
                << "    serial max slip_rate (kicked)  = " << V_serial_kick << "\n";
   }
   MPI_Bcast(&V_serial, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&V_serial_kick, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   // ---- Parallel shared-fault path (np=2, explicit partition: fault on seam) ----
   if (rank == 0)
   {
      std::cout << "\n-- Parallel (shared-fault branch, :3260, explicit partition) --\n";
   }
   Mesh serial_mesh = BuildTwoTetFaultMesh();
   std::vector<int> partitioning(serial_mesh.GetNE());
   for (int e = 0; e < serial_mesh.GetNE(); e++)
   {
      real_t cy = 0; Array<int> ev; serial_mesh.GetElementVertices(e, ev);
      for (int v = 0; v < ev.Size(); v++) { cy += serial_mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      partitioning[e] = (cy < 0) ? 0 : 1;   // − side → rank 0, + side → rank 1
   }
   ParMesh pmesh(MPI_COMM_WORLD, serial_mesh, partitioning.data());
   BoundaryConfig bc; bc.natural_attrs = {1}; bc.fault_attr = 3;
   bc.absorbing_attrs = {};
   WaveOperator<ParMesh> wave(pmesh, order, TPV205Params::lambda,
                              TPV205Params::mu, TPV205Params::rho, bc);
   real_t zero_bg[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(zero_bg);
   FaultFaceFlux ff(TPV205Params::rho, TPV205Params::cp, TPV205Params::cs);
   std::vector<DOFData> dof_data;
   SetupLSWFaultMPI<ParMesh>(wave, pmesh, order, dof_data, ff);

   // Confirm the fault is on the rank seam (shared QPs exist globally).
   int local_shared = wave.GetNumSharedFaultQPs();
   int global_shared = 0;
   MPI_Allreduce(&local_shared, &global_shared, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   Vector Q(wave.Height()); Q = 0.0;
   Vector k(wave.Height());
   wave.Mult(Q, k);                          // (a) at-rest
   bool par_finite = AllFinite(k);
   real_t V_local = 0.0;
   for (const auto &d : dof_data)
   {
      if (!std::isfinite(d.slip_rate) || d.slip_rate < 0.0) { par_finite = false; }
      V_local = std::max(V_local, d.slip_rate);
   }
   // (b) cross-rank jump: same one-sided VX kick (only rank 0's −side tet is
   // local to it), then re-solve through the SHARED-fault dispatch (:3260).
   KickMinusSideVX(pmesh, wave.GetFESpace(), Q, V_kick);
   wave.Mult(Q, k);
   if (!AllFinite(k)) { par_finite = false; }
   real_t V_local_kick = 0.0;
   for (const auto &d : dof_data)
   {
      if (!std::isfinite(d.slip_rate) || d.slip_rate < 0.0) { par_finite = false; }
      V_local_kick = std::max(V_local_kick, d.slip_rate);
   }
   real_t V_parallel = 0.0, V_parallel_kick = 0.0;
   MPI_Allreduce(&V_local, &V_parallel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&V_local_kick, &V_parallel_kick, 1, MPI_DOUBLE, MPI_MAX,
                 MPI_COMM_WORLD);
   int all_finite = (par_finite ? 1 : 0), global_finite = 1;
   MPI_Allreduce(&all_finite, &global_finite, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

   if (rank == 0)
   {
      std::cout << "    global shared fault QPs = " << global_shared << "\n"
                << "    parallel max slip_rate (at-rest) = " << std::scientific
                << std::setprecision(6) << V_parallel << "\n"
                << "    parallel max slip_rate (kicked)  = " << V_parallel_kick
                << "  (serial kicked " << V_serial_kick << ")\n";
   }

   // ---- Assertions (reported on rank 0) ----
   check(serial_finite, "serial interior Mult: dQdt + slip_rate finite, non-neg");
   check(std::abs(V_serial - V_analytic) <= 1e-3 * V_analytic,
         "serial interior slip_rate == analytic LSW V (0.0779 m/s)");
   check(global_shared > 0,
         "parallel mesh has shared fault QPs (the :3260 path is exercised)");
   check(global_finite == 1,
         "parallel shared Mult: dQdt + slip_rate finite, non-neg on all ranks");
   check(std::abs(V_parallel - V_analytic) <= 1e-3 * V_analytic,
         "parallel shared slip_rate == analytic LSW V");
   check(std::abs(V_parallel - V_serial) <= 1e-9 * std::max(V_serial, real_t(1)),
         "parallel shared slip_rate == serial interior slip_rate (:3260 == :2617)");
   // R-001: the one-sided VX kick creates a real cross-rank Q± jump, so this
   // parity now exercises the shared-fault Q-exchange (NOT just the at-rest
   // prestress).  A swapped +/− side or stale ghost in the :3260 pairing would
   // make V_parallel_kick != V_serial_kick.
   check(V_serial_kick > 0.0 &&
         std::abs(V_serial_kick - V_analytic) > 1e-6 * V_analytic,
         "VX kick actually perturbed V (cross-rank jump is non-trivial, not a no-op)");
   check(std::abs(V_parallel_kick - V_serial_kick)
         <= 1e-9 * std::max(V_serial_kick, real_t(1)),
         "KICKED parallel shared slip_rate == serial interior (:3260 Q-exchange == :2617)");

   int global_failed = 0;
   MPI_Allreduce(&num_failed, &global_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   if (rank == 0)
   {
      std::cout << "\n========================================\n"
                << (global_failed == 0 ? "  L8 PASSED" : "  L8 FAILED")
                << " (" << global_failed << " failures)\n"
                << "========================================\n";
   }
   MPI_Finalize();
   return (global_failed == 0) ? 0 : 1;
}
