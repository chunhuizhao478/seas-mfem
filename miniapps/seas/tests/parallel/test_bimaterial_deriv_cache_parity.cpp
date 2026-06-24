// Parallel (ParMesh) BimaterialWaveOperator DerivMode + shared-CK parity test.
//
// GATE for the 2026-06-24 spatial-dyn optimization (document/code_optimization_dev/
// action_plan_2026-06-24.md, Phase 0):
//
//   - 0B (--deriv-cache): the existing scalar test_wave_operator_cached_parallel
//     proves Cached==OnTheFly only for the SCALAR `WaveOperator`.  The SAFS
//     production run uses `interior_flux="matrix"` -> `BimaterialWaveOperator`
//     with per-element A_d^e from the flux pool, which the scalar test does NOT
//     cover.  This test closes that gap: on a HOMOGENEOUS-Coefficient bimaterial
//     operator, the cached `D_d^e`/`S_d^e` mat-vec must match the OnTheFly
//     per-QP quadrature to <= 1e-12 (machine-eps, NOT bit-exact — R-002 round-off
//     re-association) for ApplySpatialDerivative, ComputeVolumeRHS, AND the full
//     AdvanceADER predictor-corrector (the production path).
//
//   - 0D / shared-CK (--shared-ck-recursion): the scalar TestSharedCKParity proves
//     fused==separate bit-exact only for `WaveOperator`.  This test adds the
//     bimaterial proof: ComputeADERSubStepStatesAndIntegral + AdvanceADER(...,&I)
//     == ComputeADERSubStepStates + AdvanceADER(...) to EXACTLY 0.0, for both
//     DerivModes, on the matrix path.
//
// The CK recursion + the cached kernels are element-LOCAL (no MPI), so on a
// partitioned ParMesh each rank runs them on its local elements; AdvanceADER
// additionally exercises the local + shared face correctors (and thus the
// R-1505 one-time bulk-bg consensus hoist) across `n_steps` macro-steps.
//
// Standalone, fault-free (bc.fault_attr=0).  TWO-material Coefficient fixture
// (A in x<0.5, B in x>=0.5) so per-element flux-pool divergence is exercised
// (REVIEW R-003); FE order swept over {1,2} so the PRODUCTION P1 deriv-cache
// regime is covered (REVIEW R-004); the consensus re-arm is exercised by a
// second SetAbsorbingBackground call (REVIEW R-006).
//
// Run with: mpirun -np 2 ./seas_test_bimaterial_deriv_cache_parity
//           mpirun -np 4 ./seas_test_bimaterial_deriv_cache_parity

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/bimaterial_wave_operator.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../domain/boundary_config.hpp"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

#ifdef MFEM_USE_MPI
// Per-rank max relative difference of two equal-length vectors.
static double RelDiff(const Vector &a, const Vector &b)
{
   double mr = 0.0, md = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      mr = std::max(mr, std::abs(a[i]));
      md = std::max(md, std::abs(a[i] - b[i]));
   }
   return md / (mr + 1e-300);
}

// Deterministic, rank-shifted fill so ranks carry distinct local state.
static void FillQ(Vector &Q, int rank, unsigned salt)
{
   unsigned s = salt + 7919u * static_cast<unsigned>(rank);
   for (int i = 0; i < Q.Size(); i++)
   {
      s = 1664525u * s + 1013904223u;
      Q[i] = (static_cast<double>(s) / 4294967296.0) - 0.5;
   }
}
#endif

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);

   if (rank == 0)
   {
      std::cout << "=== Parallel BimaterialWaveOperator DerivMode + shared-CK "
                << "parity (np=" << nprocs << ") ===\n";
   }

   Mesh serial_mesh = Mesh::MakeCartesian3D(4, 4, 4, Element::TETRAHEDRON,
                                            1.0, 1.0, 1.0);
   for (int b = 0; b < serial_mesh.GetNBE(); b++)
   {
      serial_mesh.SetBdrAttribute(b, 5);
   }
   serial_mesh.SetAttributes();

   ParMesh pmesh(comm, serial_mesh);
   BoundaryConfig bc;
   bc.absorbing_attrs = {5};
   bc.fault_attr = 0;   // fault-free: exercises the matrix path without a fault

   // TWO-material fixture (REVIEW R-003 / plan §0B): material A in x<0.5, B in
   // x>=0.5.  With distinct A_d^e per element, every access routes through the
   // per-element flux pool, so a pool index/divergence bug (or a leak of the
   // (1,1,1) `flux_` sentinel) produces a non-zero Cached==OnTheFly /
   // fused==separate delta — a homogeneous fixture is blind to it.
   const real_t lamA = 32.04e9, muA = 32.04e9, rhoA = 2670.0;
   const real_t lamB = 48.06e9, muB = 40.00e9, rhoB = 3000.0;
   FunctionCoefficient lam_c([=](const Vector &x){ return x[0] < 0.5 ? lamA : lamB; });
   FunctionCoefficient mu_c ([=](const Vector &x){ return x[0] < 0.5 ? muA  : muB;  });
   FunctionCoefficient rho_c([=](const Vector &x){ return x[0] < 0.5 ? rhoA : rhoB; });

   double local_max_cached = 0.0;   // Cached vs OnTheFly (tol 1e-12)
   double local_max_ck     = 0.0;   // fused vs separate (tol 0.0, bit-exact)

   // REVIEW R-004: sweep the FE order.  P1 (ndof=4, affine tets) is the
   // PRODUCTION stride for --deriv-cache — where OnTheFly inverts the per-QP
   // Jacobian and the cached matvec re-associates round-off (the reason the tol
   // is 1e-12, not 0.0); P2 keeps higher-order coverage.
   for (int fe_order : {1, 2})
   {
      BimaterialWaveOperator<ParMesh> wave(
         pmesh, fe_order, MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c), bc);

      // Total-Q background (0 is valid under fluctuation-Q dispatch); AdvanceADER
      // and the shared corrector require it (R-1505).
      real_t Q_bg[NUM_STATE];
      for (int c = 0; c < NUM_STATE; c++) { Q_bg[c] = 0.0; }
      wave.SetAbsorbingBackground(Q_bg);

      const int N = wave.Height();

      // -- Part A: ApplySpatialDerivative + ComputeVolumeRHS, Cached vs OnTheFly --
      {
         Vector Q(N);
         FillQ(Q, rank, 12345u);

         for (int d = 0; d < 3; d++)
         {
            wave.SetDerivMode(DerivMode::OnTheFly);
            Vector r1; wave.ApplySpatialDerivative(d, Q, r1);
            wave.SetDerivMode(DerivMode::Cached);
            Vector r2; wave.ApplySpatialDerivative(d, Q, r2);
            local_max_cached = std::max(local_max_cached, RelDiff(r1, r2));
         }

         wave.SetDerivMode(DerivMode::OnTheFly);
         Vector v1(N); v1 = 0.0; wave.ComputeVolumeRHS_ForTest(Q, v1);
         wave.SetDerivMode(DerivMode::Cached);
         Vector v2(N); v2 = 0.0; wave.ComputeVolumeRHS_ForTest(Q, v2);
         local_max_cached = std::max(local_max_cached, RelDiff(v1, v2));
      }

      // -- Part B: full AdvanceADER predictor-corrector, Cached vs OnTheFly,
      //    over several macro-steps (exercises the local + shared corrector and
      //    the R-1505 one-time bulk-bg consensus hoist on steps 2..n). --
      const real_t dt = 1e-4;
      const int ader_order = 2;
      {
         const int n_steps = 3;
         Vector Q0(N); FillQ(Q0, rank, 271828u);

         wave.SetDerivMode(DerivMode::OnTheFly);
         Vector q_otf(Q0);
         for (int s = 0; s < n_steps; s++)
         {
            Vector qn(N); wave.AdvanceADER(q_otf, dt, ader_order, qn); q_otf = qn;
         }

         wave.SetDerivMode(DerivMode::Cached);
         Vector q_cac(Q0);
         for (int s = 0; s < n_steps; s++)
         {
            Vector qn(N); wave.AdvanceADER(q_cac, dt, ader_order, qn); q_cac = qn;
         }
         wave.SetDerivMode(DerivMode::OnTheFly);
         local_max_cached = std::max(local_max_cached, RelDiff(q_otf, q_cac));

         // REVIEW R-006: re-arm coverage.  Call SetAbsorbingBackground again
         // (a second, still-valid background) — this resets bulk_bg_consensus_done_
         // so the NEXT corrector re-runs the one-time consensus.  Output must
         // still match Cached vs OnTheFly; the run must not hang.
         real_t Q_bg2[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++) { Q_bg2[c] = (c < 6) ? 1.0e3 * (c + 1) : 0.0; }
         wave.SetAbsorbingBackground(Q_bg2);
         wave.SetDerivMode(DerivMode::OnTheFly);
         Vector r_otf(N); wave.AdvanceADER(q_otf, dt, ader_order, r_otf);
         wave.SetDerivMode(DerivMode::Cached);
         Vector r_cac(N); wave.AdvanceADER(q_cac, dt, ader_order, r_cac);
         wave.SetDerivMode(DerivMode::OnTheFly);
         local_max_cached = std::max(local_max_cached, RelDiff(r_otf, r_cac));
         wave.SetAbsorbingBackground(Q_bg);   // restore zero bg
      }

      // -- Part C: shared-CK fusion bit-exact (fused==separate, tol 0.0), both
      //    DerivModes, on the bimaterial (matrix) path. --
      for (int order_c = 2; order_c <= 4; order_c++)
      {
         std::vector<real_t> tau_nodes(order_c);
         for (int o = 0; o < order_c; o++)
         {
            tau_nodes[o] = dt * (static_cast<real_t>(o) + 0.5)
                         / static_cast<real_t>(order_c);
         }

         for (int mode_i = 0; mode_i < 2; mode_i++)
         {
            wave.SetDerivMode(mode_i == 0 ? DerivMode::OnTheFly
                                          : DerivMode::Cached);
            Vector Q(N);
            FillQ(Q, rank, 31415u + static_cast<unsigned>(order_c)
                                 + 100u * static_cast<unsigned>(mode_i));

            std::vector<Vector> qpn_sep;
            wave.ComputeADERSubStepStates(Q, dt, order_c, tau_nodes, qpn_sep);
            Vector qn_sep; wave.AdvanceADER(Q, dt, order_c, qn_sep);

            std::vector<Vector> qpn_m; Vector I_m;
            wave.ComputeADERSubStepStatesAndIntegral(Q, dt, order_c, tau_nodes,
                                                     qpn_m, I_m);
            Vector qn_m; wave.AdvanceADER(Q, dt, order_c, qn_m, &I_m);

            double dmax = 0.0;
            for (std::size_t o = 0; o < qpn_sep.size(); o++)
            {
               for (int i = 0; i < N; i++)
               {
                  dmax = std::max(dmax, std::abs(qpn_sep[o][i] - qpn_m[o][i]));
               }
            }
            for (int i = 0; i < N; i++)
            {
               dmax = std::max(dmax, std::abs(qn_sep[i] - qn_m[i]));
            }
            local_max_ck = std::max(local_max_ck, dmax);
         }
      }
      wave.SetDerivMode(DerivMode::OnTheFly);   // restore default
   }   // fe_order sweep

   double global_max_cached = 0.0, global_max_ck = 0.0;
   MPI_Allreduce(&local_max_cached, &global_max_cached, 1, MPI_DOUBLE,
                 MPI_MAX, comm);
   MPI_Allreduce(&local_max_ck, &global_max_ck, 1, MPI_DOUBLE, MPI_MAX, comm);

   const int failed_cached = (global_max_cached <= 1e-12) ? 0 : 1;
   const int failed_ck     = (global_max_ck     <= 0.0)   ? 0 : 1;
   const int failed = failed_cached | failed_ck;

   if (rank == 0)
   {
      std::cout << (failed_cached ? "  FAILED" : "  PASSED")
                << ": bimaterial Cached==OnTheFly "
                << "(ApplySpatialDerivative + ComputeVolumeRHS + AdvanceADER, "
                << "max rel " << global_max_cached << ", tol 1e-12)\n";
      std::cout << (failed_ck ? "  FAILED" : "  PASSED")
                << ": bimaterial shared-CK fused==separate "
                << "(bit-exact, max abs " << global_max_ck << ", tol 0)\n";
      std::cout << (failed ? "=== FAILED ===\n" : "=== PASSED ===\n");
   }

   MPI_Finalize();
   return failed;
#else
   (void)argc; (void)argv;
   std::cout << "Parallel test requires MFEM_USE_MPI.\n";
   return 0;
#endif
}
