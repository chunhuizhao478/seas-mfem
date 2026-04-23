// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper-bug: MULTI-STEP rupture-regime serial-vs-parallel
// cross-check.
//
// ============================================================================
// Motivation
// ============================================================================
// `test_interior_vs_shared_branch_live` proves that the interior-fault
// and shared-fault assembly branches produce bit-identical rhs in a
// SINGLE wave.Mult at pre-rupture amplitudes (V ~ 1e-5 m/s).  The
// pepper pattern observed in the 400-rank Frontera run 7672897 only
// appears INSIDE the rupture area (V ~ 1 m/s) and fills the whole
// ruptured disk at sub-element scale, not just partition boundaries.
//
// Two hypotheses remain untested by the existing pepper suite:
//   HR1 — multi-step accumulation: per-step divergence below the
//         single-step tolerance but accumulating over many steps.
//   HR2 — persistent-nucleation path: DOFData.tau2_nuc read during
//         EvaluateTotal may interact differently with the interior vs
//         shared branch when the friction solver is in the rupture
//         regime (V_abs ~ 0.1 - 1 m/s).
//
// This test exercises both:
//   1. Build the same 2-tet fault fixture SERIAL and np=2 PARALLEL.
//   2. Set DOFData.tau2_nuc = nuc_dtau at every fault QP so the
//      friction solver sees the full 25 MPa rupture drive.
//   3. Run N=20 ADER-2 steps on both runs.
//   4. After EACH step, compare bulk Q at the rank-0-local element's
//      fault-adjacent DOFs AND DOFData.slip_rate at matching fault
//      QPs between serial and np=2.
//   5. Assert max relative divergence < 1e-10 on every step and every
//      component.
//
// Interpretation of outcomes:
//   - PASS: serial and np=2 stay bit-identical through rupture, so
//     the pepper cannot be caused by the interior/shared dispatch or
//     by MPI ghost exchange on this fixture.  Pepper origin is
//     elsewhere (e.g. per-QP basis computation at rupture amplitudes,
//     triangle-orientation quadrature differences, etc.).
//   - FAIL at step K, component C: we've localized the divergence to
//     step K, rupture regime, component C.  The failing component
//     names whether the issue is in tangential (SXY/SXZ), normal
//     (SXX), or velocity (VX/VY/VZ) coupling.  Cuts the search space
//     by ~6×.
//
// ============================================================================
// Usage:  mpirun -np 2 ./seas_test_rupture_multistep_serial_vs_parallel
//         Requires MFEM built with MPI.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   double vv = (v), tt = (tol); \
   if (vv <= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", tol " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; } \
} while (0)

namespace {

constexpr real_t kL = 1000.0;       // 1 km cube — ADER-CFL-friendly for dt = 1e-4 s
constexpr real_t kDt = 1.0e-4;      // matches production sub-step for 200 m mesh
constexpr int    kNSteps = 20;
constexpr int    kOrder = 1;

const char *StateName(int c)
{
   static const char *kNames[NUM_STATE] =
   {"SXX", "SYY", "SZZ", "SXY", "SYZ", "SXZ", "VX", "VY", "VZ"};
   return (c >= 0 && c < NUM_STATE) ? kNames[c] : "?";
}

Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {kL, 0.0, 0.0}, {0.0, 0.0, kL},
      {0.0,  kL, 0.0},   // V3 on +y side
      {0.0, -kL, 0.0},   // V4 on -y side
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
         if (std::abs(cy) < 1e-8)
         { mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3); }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Populate fault_coords + DOFData from wave operator's fault-face lists.
// Returns number of fault QPs local to this rank.
template <typename MeshT>
int SetupFault(WaveOperator<MeshT> &wave, MeshT &mesh, int order,
               std::vector<DOFData> &dof_data, FaultFaceFlux &ff,
               std::vector<Vector> &fault_coords)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();

   int nqp = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      MFEM_VERIFY(ftr, "interior fault FTR null");
      nqp = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>)
   {
      if (nqp == 0 && shr_faces.Size() > 0)
      {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[0]);
         MFEM_VERIFY(ftr, "shared fault FTR null");
         nqp = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
      }
   }
#endif

   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                2*order);
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
         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                   2*order);
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

   const int n_fault =
      (int_faces.Size() + shr_faces.Size()) * nqp;
   if (n_fault > 0)
   {
      InitializeFaultDOFs(dof_data, n_fault, fault_coords);
      ZeroDOFDataPreStressTotal(dof_data, n_fault);
      // Persistent nucleation: set every fault QP's tau2_nuc to nuc_dtau.
      // Intentionally uniform across QPs (not a spatial NucleationPerturbation
      // call) so both serial and parallel see the IDENTICAL forcing — the
      // only remaining degree of freedom is then the branch dispatch.
      for (int i = 0; i < n_fault; i++)
      {
         dof_data[i].tau2_nuc = TPV102Params::nuc_dtau;
      }
   }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);
   return n_fault;
}

// Snapshot DOFData.slip_rate into a vector for comparison.
void SnapshotSlipRates(const std::vector<DOFData> &dof_data,
                       std::vector<real_t> &out)
{
   out.resize(dof_data.size());
   for (size_t i = 0; i < dof_data.size(); i++) { out[i] = dof_data[i].slip_rate; }
}

// Find the rank-0 element whose centroid has y < 0 ("minus side").
// Returns -1 if none.
template <typename MeshT>
int FindMinusSideElem(MeshT &mesh)
{
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      Array<int> ev; mesh.GetElementVertices(e, ev);
      real_t cy = 0;
      for (int v = 0; v < ev.Size(); v++) { cy += mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      if (cy < 0.0) { return e; }
   }
   return -1;
}

} // anonymous

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if (rank == 0)
   {
      std::cout << "\n=== TPV102 pepper-bug MULTI-STEP rupture-regime "
                << "serial-vs-parallel cross-check ===\n"
                << "  MPI ranks: " << nprocs << "\n"
                << "  fixture L = " << kL << " m, dt = " << std::scientific
                << std::setprecision(3) << kDt << " s, N = " << kNSteps
                << " ADER-2 steps\n"
                << "  tau2_nuc = " << TPV102Params::nuc_dtau
                << " Pa (full nucleation at every fault QP)\n";
   }

   if (nprocs != 2)
   {
      if (rank == 0)
      { std::cout << "  SKIPPED: requires np == 2\n"; }
      MPI_Finalize();
      return 77;
   }

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {};
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;

   // -----------------------------------------------------------------
   // Serial run (rank 0 only).  Record Q per step at the rank-0 element
   // and dof-0, plus DOFData.slip_rate[0] snapshot.
   // -----------------------------------------------------------------
   std::vector<std::vector<real_t>> Q_serial_per_step(kNSteps);
   std::vector<real_t> slip_serial_per_step(kNSteps, 0.0);
   if (rank == 0)
   {
      std::cout << "\n-- Serial run --\n";
      Mesh mesh = BuildTwoTetFaultMesh();
      WaveOperator<Mesh> wave(mesh, kOrder, TPV102Params::lambda,
                              TPV102Params::mu, TPV102Params::rho, bc);
      wave.SetAbsorbingBackground(bulk_bg);

      FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
      std::vector<DOFData> dof_data;
      std::vector<Vector> fault_coords;
      const int n_fault = SetupFault(wave, mesh, kOrder, dof_data, ff,
                                      fault_coords);
      MFEM_VERIFY(n_fault > 0, "serial: no fault QPs");

      const auto &fes = wave.GetFESpace();
      const int ndof_total = fes.GetNDofs();
      const int size = wave.Height();

      const int e_minus = FindMinusSideElem(mesh);
      MFEM_VERIFY(e_minus >= 0, "serial: no -side element");
      Array<int> edofs; fes.GetElementDofs(e_minus, edofs);
      const int dof0 = edofs[0];

      Vector Q(size);
      InitializeStateTotal(Q, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);

      Vector Q_new(size);
      for (int step = 0; step < kNSteps; step++)
      {
         wave.AdvanceADER(Q, kDt, /*ader_order=*/2, Q_new);
         Q.Swap(Q_new);
         Q_serial_per_step[step].resize(NUM_STATE);
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_serial_per_step[step][c] = Q(c * ndof_total + dof0);
         }
         // DOFData slip_rate at fault QP 0 (any QP — all are uniform here).
         slip_serial_per_step[step] = dof_data[0].slip_rate;
      }
      std::cout << "    final-step Q[SXY, minus-dof0] = " << std::scientific
                << std::setprecision(6)
                << Q_serial_per_step[kNSteps-1][SXY] << " Pa\n";
      std::cout << "    final-step slip_rate[0]       = "
                << slip_serial_per_step[kNSteps-1] << " m/s\n";
   }

   // Broadcast serial snapshots to all ranks.
   int dims[2] = {kNSteps, NUM_STATE};
   MPI_Bcast(dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
   std::vector<real_t> Q_serial_flat(kNSteps * NUM_STATE, 0.0);
   if (rank == 0)
   {
      for (int s = 0; s < kNSteps; s++)
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_serial_flat[s * NUM_STATE + c] = Q_serial_per_step[s][c];
         }
   }
   MPI_Bcast(Q_serial_flat.data(), kNSteps * NUM_STATE, MPI_DOUBLE, 0,
             MPI_COMM_WORLD);
   MPI_Bcast(slip_serial_per_step.data(), kNSteps, MPI_DOUBLE, 0,
             MPI_COMM_WORLD);

   // -----------------------------------------------------------------
   // Parallel run (np=2, explicit partitioning one tet per rank).
   // -----------------------------------------------------------------
   if (rank == 0) { std::cout << "\n-- Parallel run (np=2) --\n"; }

   Mesh serial_mesh = BuildTwoTetFaultMesh();
   std::vector<int> partitioning(serial_mesh.GetNE(), 0);
   for (int e = 0; e < serial_mesh.GetNE(); e++)
   {
      real_t cy = 0;
      Array<int> ev; serial_mesh.GetElementVertices(e, ev);
      for (int v = 0; v < ev.Size(); v++) { cy += serial_mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      partitioning[e] = (cy < 0.0) ? 0 : 1;  // -side to rank 0, +side to rank 1
   }
   ParMesh pmesh(MPI_COMM_WORLD, serial_mesh, partitioning.data());

   WaveOperator<ParMesh> wave_p(pmesh, kOrder, TPV102Params::lambda,
                                 TPV102Params::mu, TPV102Params::rho, bc);
   wave_p.SetAbsorbingBackground(bulk_bg);

   FaultFaceFlux ff_p(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   std::vector<DOFData> dof_p;
   std::vector<Vector> fault_coords_p;
   const int n_fault_p = SetupFault(wave_p, pmesh, kOrder, dof_p, ff_p,
                                     fault_coords_p);

   const auto &fes_p = wave_p.GetFESpace();
   const int ndof_total_p = fes_p.GetNDofs();
   const int size_p = wave_p.Height();

   const int e_minus_p = FindMinusSideElem(pmesh);
   const int dof0_p = [&]() -> int {
      if (e_minus_p < 0) { return -1; }
      Array<int> edofs; fes_p.GetElementDofs(e_minus_p, edofs);
      return edofs.Size() > 0 ? edofs[0] : -1;
   }();

   Vector Q_p(size_p);
   InitializeStateTotal(Q_p, ndof_total_p, TPV102Params::sigma_n,
                        TPV102Params::tau_ini);
   Vector Q_p_new(size_p);

   std::vector<std::vector<real_t>> Q_parallel_per_step(kNSteps,
                                                         std::vector<real_t>(NUM_STATE, 0.0));
   std::vector<real_t> slip_parallel_per_step(kNSteps, 0.0);

   for (int step = 0; step < kNSteps; step++)
   {
      wave_p.AdvanceADER(Q_p, kDt, /*ader_order=*/2, Q_p_new);
      Q_p.Swap(Q_p_new);
      if (rank == 0 && dof0_p >= 0)
      {
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_parallel_per_step[step][c] = Q_p(c * ndof_total_p + dof0_p);
         }
         if (n_fault_p > 0)
         {
            slip_parallel_per_step[step] = dof_p[0].slip_rate;
         }
      }
   }

   // -----------------------------------------------------------------
   // Compare serial vs parallel at each step, each component.
   // -----------------------------------------------------------------
   if (rank == 0)
   {
      std::cout << "\n-- Step-by-step comparison (rank 0 -side element dof 0) --\n";
      std::cout << std::setw(5) << "step"
                << std::setw(14) << "worst rel"
                << std::setw(6) << "comp"
                << std::setw(16) << "|slip_ser|"
                << std::setw(16) << "|slip_par|"
                << "\n";

      real_t worst_rel_overall = 0.0;
      int worst_step = -1, worst_comp = -1;
      real_t worst_slip_diff = 0.0;

      for (int step = 0; step < kNSteps; step++)
      {
         real_t worst_rel_step = 0.0;
         int comp_worst = -1;
         for (int c = 0; c < NUM_STATE; c++)
         {
            const real_t a = Q_serial_flat[step * NUM_STATE + c];
            const real_t b = Q_parallel_per_step[step][c];
            const real_t denom = std::max({std::abs(a), std::abs(b), real_t(1.0)});
            const real_t rel = std::abs(a - b) / denom;
            if (rel > worst_rel_step) { worst_rel_step = rel; comp_worst = c; }
         }
         if (worst_rel_step > worst_rel_overall)
         {
            worst_rel_overall = worst_rel_step;
            worst_step = step;
            worst_comp = comp_worst;
         }
         const real_t slip_diff = std::abs(slip_serial_per_step[step]
                                            - slip_parallel_per_step[step]);
         worst_slip_diff = std::max(worst_slip_diff, slip_diff);

         std::cout << std::setw(5) << step
                   << "  " << std::scientific << std::setprecision(3)
                   << worst_rel_step
                   << "  " << std::setw(3) << StateName(comp_worst)
                   << "  " << std::scientific << std::setprecision(3)
                   << std::abs(slip_serial_per_step[step])
                   << "  " << std::abs(slip_parallel_per_step[step])
                   << "\n";
      }

      std::cout << "\n  worst_rel over all steps = " << std::scientific
                << std::setprecision(6) << worst_rel_overall
                << "  at step " << worst_step
                << ", component " << StateName(worst_comp) << "\n"
                << "  worst |slip_serial - slip_parallel| = "
                << worst_slip_diff << " m/s\n";

      TEST_LE(worst_rel_overall, 1.0e-10,
              "rupture-regime serial-vs-parallel bulk Q matches to 1e-10");
      TEST_LE(worst_slip_diff, 1.0e-10,
              "rupture-regime serial-vs-parallel slip_rate matches to 1e-10");

      std::cout << "\n========================================\n"
                << "  Results: " << num_passed << " passed, "
                << num_failed << " failed out of " << num_tests << " tests\n"
                << "========================================\n";
   }

   int exit_code = (num_failed == 0) ? 0 : 1;
   MPI_Bcast(&exit_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Finalize();
   return exit_code;
}
