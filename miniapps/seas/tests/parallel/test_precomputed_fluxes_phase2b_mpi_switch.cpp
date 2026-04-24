// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 "Topology-Based Precomputed Face-Rotation" plan 2026-04-23
// Phase 2b acceptance gates (§7.3) — MPI shared-face dispatch switch.
//
//   P_SWITCH_ON_MPI_LINEAR_EQ — on a homogeneous Cartesian plane-wave
//     ParMesh fixture, AdvanceADER with the flag ON must match the
//     runtime-path output to 1e-10 (the plan relaxes from 1e-12 to
//     absorb MPI reduction ordering noise).
//
//   P_SWITCH_OFF_MPI_REGRESSION — same fixture with flag OFF: default
//     path unchanged (bit-identical re-runs of AdvanceADER).
//
// Also exercises the Phase 2b Init shared-face pass:
//   * `InitSharedFaces` must push exactly ONE FaceEntry per shared
//     non-fault face per rank, keyed by (mesh_face_idx, Elem1No).
//   * On a no-fault fixture, every shared face must have an entry.
//
// Run with: mpirun -np N ./seas_test_precomputed_fluxes_phase2b_mpi_switch
// (plan §7.3 targets np=4).

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/precomputed_face_fluxes.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define PAR_TEST_LE(v, tol, msg) do { \
   num_tests++; \
   const double vv = (v), tt = (tol); \
   if (vv <= tt) { num_passed++; \
      if (rank == 0) std::cout << "  PASSED: " << msg \
                               << "  (got " << std::scientific \
                               << std::setprecision(3) << vv \
                               << ", tol " << tt << ")\n"; } \
   else { num_failed++; \
      if (rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                               << "  (got " << std::scientific \
                               << std::setprecision(3) << vv \
                               << ", expected <= " << tt << ")\n"; } \
} while (0)

#define PAR_TEST_TRUE(expr, msg) do { \
   num_tests++; \
   if (expr) { num_passed++; \
      if (rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      if (rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

namespace
{

constexpr real_t kRho    = 2670.0;
constexpr real_t kLambda = 32.04e9;
constexpr real_t kMu     = 32.04e9;

struct PWave
{
   real_t Q0[NUM_STATE];
   real_t cp;
   real_t k;
};

PWave MakePWave(real_t k)
{
   PWave w;
   const real_t lp = kLambda + 2.0 * kMu;
   w.cp = std::sqrt(lp / kRho);
   w.k  = k;
   for (int c = 0; c < NUM_STATE; c++) { w.Q0[c] = 0.0; }
   w.Q0[SXX] = -lp;
   w.Q0[SYY] = -kLambda;
   w.Q0[SZZ] = -kLambda;
   w.Q0[VX]  =  w.cp;
   return w;
}

void FillPlaneWave(const FiniteElementSpace &fes, const PWave &w,
                   real_t t, Vector &Q)
{
   const int ndof_total = fes.GetNDofs();
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
   for (int e = 0; e < fes.GetNE(); e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      const IntegrationRule &nodes = fe->GetNodes();
      Array<int> edofs; fes.GetElementDofs(e, edofs);
      for (int j = 0; j < fe->GetDof(); j++)
      {
         const IntegrationPoint &ip = nodes.IntPoint(j);
         Vector x(3);
         Tr->Transform(ip, x);
         const real_t phase = w.k * (x(0) - w.cp * t);
         const real_t s = std::sin(phase);
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q[c * ndof_total + edofs[j]] = w.Q0[c] * s;
         }
      }
   }
}

real_t InfNormDiff(const Vector &a, const Vector &b)
{
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      m = std::max(m, std::abs(a(i) - b(i)));
   }
   return m;
}

real_t Normlinf(const Vector &a)
{
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++) { m = std::max(m, std::abs(a(i))); }
   return m;
}

real_t GlobalMax(real_t local, MPI_Comm comm)
{
   real_t g = 0.0;
   MPI_Allreduce(&local, &g, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   return g;
}

} // anonymous

// ==========================================================================
// Phase 2b gates on a homogeneous Cartesian plane-wave ParMesh fixture.
// ==========================================================================
void TestMPIPlaneWaveSwitch(int rank, int nprocs, MPI_Comm comm)
{
   if (rank == 0)
   {
      std::cout << "\n[TestMPIPlaneWaveSwitch]  nprocs=" << nprocs << "\n";
   }

   const real_t L = 1000.0;
   // 4x2x4 = 32 hexes (M_ref): matches plan §7.3 fixture.
   Mesh serial = Mesh::MakeCartesian3D(4, 2, 4, Element::TETRAHEDRON,
                                       L, L, L);
   ParMesh pmesh(comm, serial);

   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;

   const int p = 1;   // Phase 1 p=1 precondition.
   WaveOperator<ParMesh> wave(pmesh, p, kLambda, kMu, kRho, bc);
   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   const auto &fes = wave.GetFESpace();

   const real_t lambda_wave = 8.0 * L;
   const real_t k = 2.0 * M_PI / lambda_wave;
   PWave w = MakePWave(k);

   const real_t cp = w.cp;
   const real_t h = L / 4.0;
   const real_t dt_cfl = h / cp;
   const real_t dt = 0.01 * dt_cfl;

   Vector Q; FillPlaneWave(fes, w, 0.0, Q);

   // ---- Runtime path (flag off, default) ----
   PAR_TEST_TRUE(!wave.UsingPrecomputedFaceFluxes(),
                 "flag defaults to off");
   Vector Q_runtime;
   wave.AdvanceADER(Q, dt, /*order=*/2, Q_runtime);
   const real_t q_inf_local = Normlinf(Q_runtime);
   const real_t q_inf = GlobalMax(q_inf_local, comm);
   PAR_TEST_TRUE(q_inf > 0.0,
                 "runtime ADER step produces nonzero output on ≥1 rank");

   // ---- Precomputed path (flag on) ----
   wave.UsePrecomputedFaceFluxes(true);
   PAR_TEST_TRUE(wave.UsingPrecomputedFaceFluxes(),
                 "UsePrecomputedFaceFluxes(true) sets the flag");
   PAR_TEST_TRUE(wave.GetPrecomputedFaceFluxes().IsInitialized(),
                 "Init has run");

   // Phase 2b Init scope check: on a no-fault ParMesh fixture the
   // local shared non-fault face count must equal the number of
   // (face_idx, Elem1No) keys added by InitSharedFaces.  Elem1No is
   // always >= 0 for shared faces; `precomputed_face_fluxes_.HasEntry(
   // shared_mesh_face_idx, Elem1No)` must be true for every one.
   {
      int n_shared = pmesh.GetNSharedFaces();
      int n_shared_with_entry = 0;
      for (int sf = 0; sf < n_shared; sf++)
      {
         FaceElementTransformations *ftr =
            pmesh.GetSharedFaceTransformations(sf);
         if (!ftr) { continue; }
         int mesh_face_idx = pmesh.GetSharedFace(sf);
         int e1 = ftr->Elem1No;
         if (wave.GetPrecomputedFaceFluxes().HasEntry(mesh_face_idx, e1))
         {
            n_shared_with_entry++;
         }
      }
      int local_match = (n_shared_with_entry == n_shared) ? 1 : 0;
      int global_match = 0;
      MPI_Allreduce(&local_match, &global_match, 1, MPI_INT, MPI_MIN, comm);
      PAR_TEST_TRUE(global_match == 1,
                    "InitSharedFaces (§7.2): every shared non-fault face "
                    "has a FaceEntry keyed by (mesh_face_idx, Elem1No) "
                    "on every rank");
   }

   Vector Q_precomputed;
   wave.AdvanceADER(Q, dt, /*order=*/2, Q_precomputed);

   const real_t abs_diff_local = InfNormDiff(Q_precomputed, Q_runtime);
   const real_t abs_diff = GlobalMax(abs_diff_local, comm);

   if (rank == 0)
   {
      std::cout << "  global ||Q_precomputed - Q_runtime||_inf = "
                << std::scientific << abs_diff
                << "  (||Q_runtime||_inf = " << q_inf << ")\n";
   }
   // P_SWITCH_ON_MPI_LINEAR_EQ: plan §7.3 tol 1e-10 (relaxed from
   // 1e-12 to absorb MPI reduction-ordering noise).
   PAR_TEST_LE(abs_diff / std::max(q_inf, real_t(1.0)), 1e-10,
               "P_SWITCH_ON_MPI_LINEAR_EQ: ADER step matches runtime "
               "to 1e-10 rel on ParMesh");

   // ---- P_SWITCH_OFF_MPI_REGRESSION ----
   wave.UsePrecomputedFaceFluxes(false);
   PAR_TEST_TRUE(!wave.UsingPrecomputedFaceFluxes(),
                 "UsePrecomputedFaceFluxes(false) clears the flag");
   Vector Q_runtime_again;
   wave.AdvanceADER(Q, dt, /*order=*/2, Q_runtime_again);
   const real_t rerun_diff_local = InfNormDiff(Q_runtime_again, Q_runtime);
   const real_t rerun_diff = GlobalMax(rerun_diff_local, comm);
   PAR_TEST_LE(rerun_diff / std::max(q_inf, real_t(1.0)), 1e-15,
               "P_SWITCH_OFF_MPI_REGRESSION: flag-off ADER step "
               "bit-identical before/after flag toggling");
}

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank, nprocs;
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);

   if (rank == 0)
   {
      std::cout << "\n=== TPV102 Phase 2b MPI dispatch switch tests ===\n";
      std::cout << "  nprocs=" << nprocs << "\n";
   }

   TestMPIPlaneWaveSwitch(rank, nprocs, comm);

   // Reduce pass/fail counts across ranks so rank 0's summary reflects
   // the whole job.
   int t_local = num_tests, p_local = num_passed, f_local = num_failed;
   int t_sum = 0, p_sum = 0, f_sum = 0;
   MPI_Reduce(&t_local, &t_sum, 1, MPI_INT, MPI_SUM, 0, comm);
   MPI_Reduce(&p_local, &p_sum, 1, MPI_INT, MPI_SUM, 0, comm);
   MPI_Reduce(&f_local, &f_sum, 1, MPI_INT, MPI_SUM, 0, comm);

   if (rank == 0)
   {
      std::cout << "\n===================================================\n";
      std::cout << "Total:  " << t_sum << "  (per-rank: " << t_local << ")\n";
      std::cout << "Passed: " << p_sum << "\n";
      std::cout << "Failed: " << f_sum << "\n";
      std::cout << "===================================================\n";
   }

   int rc = (num_failed == 0) ? 0 : 1;
   int any_fail = 0;
   MPI_Allreduce(&rc, &any_fail, 1, MPI_INT, MPI_MAX, comm);
   MPI_Finalize();
   return any_fail;
}
