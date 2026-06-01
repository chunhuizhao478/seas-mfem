// =====================================================================
// bench_ader_hotpath.cpp — local micro-benchmark for the ADER hot kernels
// identified in document/caliper_perfgraph_dev/ader_hotpath_optimization_plan.
//
// Measures wall time of the three dominant per-macro-step operations on a
// SMALL serial tet mesh (NOT a full-mesh run — see memory
// `feedback-no-local-mesh-runs`; this is a unit-scale timing harness):
//
//   * ApplySpatialDerivative   — the 67.6%-of-runtime kernel (called 6*(O-1)x/step)
//   * ComputeADERSubStepStates — CK recursion #1 (the substep/friction predictor)
//   * AdvanceADER              — CK recursion #2 (predictor) + corrector
//
// It also reports the "macro-step" aggregate (SubStepStates + AdvanceADER),
// which is the quantity whose before/after ratio approximates the wall-time
// speedup of Levers 1 & 2 on a fault-OWNING rank.
//
// USAGE:
//   ./seas_bench_ader_hotpath [fe_order] [ader_order] [nx] [reps]
//   defaults: fe_order=2  ader_order=3  nx=12  reps=auto (>= ~150 ms/block)
//   (fe_order=2, ader_order=3 mirrors the tpv31_p2_aderO3 job.)
//
// HOW TO MEASURE A SPEEDUP:
//   1. On the BASELINE source, build + run; record the "mean ns/call" lines.
//   2. Apply Lever 1 / Lever 2; rebuild (make clean first — template/no-header-
//      deps hazard) + run with the SAME args.
//   3. speedup(region) = baseline_ns / optimized_ns.  The "macro-step
//      aggregate" ratio is the headline number.
//
// This is a TIMING tool, not a pass/fail test: it always exits 0 (unless a
// guard/return value is wrong).  Numerical correctness is gated by
// seas_test_wave_operator_spatial_derivative and the parity tests, NOT here.
// =====================================================================

#include "mfem.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../domain/boundary_config.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{
using Clock = std::chrono::steady_clock;

double SecondsSince(const Clock::time_point &t0)
{
   return std::chrono::duration<double>(Clock::now() - t0).count();
}

// Deterministic pseudo-random fill (no Math.random / time dependence): a
// cheap LCG so the benchmark is reproducible run-to-run and across branches.
void FillDeterministic(Vector &v, unsigned seed)
{
   unsigned s = seed ? seed : 1u;
   for (int i = 0; i < v.Size(); i++)
   {
      s = 1664525u * s + 1013904223u;          // Numerical Recipes LCG
      v[i] = (static_cast<double>(s) / 4294967296.0) - 0.5;   // [-0.5, 0.5)
   }
}

// Time `body` over `reps` calls; auto-grow reps until the block exceeds
// `min_seconds` so the per-call number is not quantization-limited.
// Returns mean nanoseconds per call; writes the chosen rep count to *used.
template <typename F>
double TimeBlock(F &&body, int reps, double min_seconds, long *used)
{
   long r = reps > 0 ? reps : 1;
   for (;;)
   {
      const auto t0 = Clock::now();
      for (long k = 0; k < r; k++) { body(); }
      const double secs = SecondsSince(t0);
      if (secs >= min_seconds || reps > 0)
      {
         if (used) { *used = r; }
         return (secs / static_cast<double>(r)) * 1e9;   // ns/call
      }
      r *= 2;                                             // grow and retry
   }
}
} // namespace

int main(int argc, char *argv[])
{
   const int  fe_order   = (argc > 1) ? std::atoi(argv[1]) : 2;
   const int  ader_order = (argc > 2) ? std::atoi(argv[2]) : 3;
   const int  nx         = (argc > 3) ? std::atoi(argv[3]) : 12;
   const int  reps_arg   = (argc > 4) ? std::atoi(argv[4]) : 0;   // 0 => auto
   const double kMinSec  = 0.15;

   MFEM_VERIFY(fe_order >= 1 && fe_order <= 4, "fe_order must be in {1..4}");
   MFEM_VERIFY(ader_order >= 2 && ader_order <= 4,
               "ader_order must be in {2,3,4}");

   // Small straight-tet mesh on the unit cube (matches TPV31's affine tets).
   Mesh mesh = Mesh::MakeCartesian3D(nx, nx, nx, Element::TETRAHEDRON,
                                     1.0, 1.0, 1.0);

   const real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;                          // no fault: pure bulk-wave cost

   WaveOperator<Mesh> wave(mesh, fe_order, lambda, mu, rho, bc);

   // AdvanceADER's corrector (ComputeADERFaceFluxRHS) requires the absorbing
   // background to have been set (has_bulk_bg_ guard); Q_bg = 0 is valid.
   real_t Q_bg[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++) { Q_bg[c] = 0.0; }
   wave.SetAbsorbingBackground(Q_bg);

   const FiniteElementSpace &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int N = NUM_STATE * ndof_total;

   Vector Q(N), dQ, Q_new(N);
   FillDeterministic(Q, 0xC0FFEEu);

   const real_t dt = 1e-4;
   std::vector<real_t> tau_nodes(ader_order);
   for (int o = 0; o < ader_order; o++)
   {
      tau_nodes[o] = dt * (static_cast<real_t>(o) + 0.5)
                   / static_cast<real_t>(ader_order);   // midpoints on [0,dt]
   }
   std::vector<Vector> Q_per_node;

   const int ndof_per_el = fes.GetFE(0)->GetDof();
   std::printf("\n=== ADER hot-path micro-benchmark ===\n");
   std::printf("fe_order=%d  ader_order=%d  mesh=%d^3 tets (ne=%d)  "
               "ndof/elem=%d  ndof_total=%d  state-dofs=%d\n",
               fe_order, ader_order, nx, mesh.GetNE(), ndof_per_el,
               ndof_total, N);

   // Time the three kernels in a given DerivMode (warmup + timed).  Cached
   // builds the D_d^e cache on SetDerivMode; the warmup also primes lazy
   // member buffers.
   // asd/sss/adv: the per-region times.  two_rec = sss+adv (the default
   // two-recursion macro-step).  merged = the Lever-2 single-recursion
   // macro-step (ComputeADERSubStepStatesAndIntegral + AdvanceADER(...,&I)).
   // vol = ComputeVolumeRHS alone (Lever 3 — part of the corrector).
   struct Timing { double asd, sss, adv, merged, vol; };
   Vector I_pre;
   Vector rhs_vol(N);
   auto measure = [&](DerivMode mode) -> Timing
   {
      wave.SetDerivMode(mode);
      // R-005: ComputeVolumeRHS accumulates (+=); zero per measure so it does
      // not grow unboundedly across reps / both modes.
      rhs_vol = 0.0;
      for (int w = 0; w < 3; w++)
      {
         wave.ApplySpatialDerivative(0, Q, dQ);
         wave.ComputeADERSubStepStates(Q, dt, ader_order, tau_nodes, Q_per_node);
         wave.AdvanceADER(Q, dt, ader_order, Q_new);
         wave.ComputeADERSubStepStatesAndIntegral(Q, dt, ader_order, tau_nodes,
                                                  Q_per_node, I_pre);
      }
      long used = 0;
      Timing t;
      // x3 dirs so the number reflects one CK-recursion level.
      t.asd = TimeBlock(
         [&]{ for (int d = 0; d < 3; d++) { wave.ApplySpatialDerivative(d, Q, dQ); } },
         reps_arg, kMinSec, &used);
      t.sss = TimeBlock(
         [&]{ wave.ComputeADERSubStepStates(Q, dt, ader_order, tau_nodes,
                                            Q_per_node); },
         reps_arg, kMinSec, &used);
      t.adv = TimeBlock(
         [&]{ wave.AdvanceADER(Q, dt, ader_order, Q_new); },
         reps_arg, kMinSec, &used);
      t.merged = TimeBlock(
         [&]{ wave.ComputeADERSubStepStatesAndIntegral(Q, dt, ader_order,
                                                       tau_nodes, Q_per_node, I_pre);
              wave.AdvanceADER(Q, dt, ader_order, Q_new, &I_pre); },
         reps_arg, kMinSec, &used);
      t.vol = TimeBlock(
         [&]{ wave.ComputeVolumeRHS_ForTest(Q, rhs_vol); },
         reps_arg, kMinSec, &used);
      return t;
   };

   const Timing otf = measure(DerivMode::OnTheFly);
   const Timing cac = measure(DerivMode::Cached);
   wave.SetDerivMode(DerivMode::OnTheFly);   // restore default

   auto row = [](const char *name, double a, double b)
   {
      std::printf("%-48s %12.1f %12.1f  %6.2fx\n", name, a, b,
                  (b > 0.0) ? a / b : 0.0);
   };
   std::printf("\n%-48s %12s %12s  %7s\n",
               "region (ns/call)", "OnTheFly", "Cached", "speedup");
   std::printf("---------------------------------------------------------------"
               "------------------\n");
   row("ApplySpatialDerivative (x3 dirs)",            otf.asd, cac.asd);
   row("ComputeADERSubStepStates (CK rec#1)",         otf.sss, cac.sss);
   row("ComputeVolumeRHS (corrector, Lever 3)",       otf.vol, cac.vol);
   row("AdvanceADER (CK rec#2 + corrector)",          otf.adv, cac.adv);
   std::printf("---------------------------------------------------------------"
               "------------------\n");
   row("macro-step: 2 recursions (default)",          otf.sss + otf.adv,
       cac.sss + cac.adv);
   row("macro-step: MERGED 1 recursion (Lever 2)",    otf.merged, cac.merged);
   std::printf("---------------------------------------------------------------"
               "------------------\n");
   std::printf("  Lever 1 (cache) on the 2-recursion macro-step : %6.2fx\n",
               (cac.sss + cac.adv) > 0.0
               ? (otf.sss + otf.adv) / (cac.sss + cac.adv) : 0.0);
   std::printf("  Lever 2 (merge) at fixed mode (OnTheFly/Cached): "
               "%6.2fx / %6.2fx\n",
               otf.merged > 0.0 ? (otf.sss + otf.adv) / otf.merged : 0.0,
               cac.merged > 0.0 ? (cac.sss + cac.adv) / cac.merged : 0.0);
   std::printf("  COMBINED (OnTheFly 2-rec baseline -> Cached merged): %6.2fx\n\n",
               cac.merged > 0.0 ? (otf.sss + otf.adv) / cac.merged : 0.0);

   // Touch results so the optimizer can't elide the timed work.
   real_t sink = dQ.Size() ? dQ[0] : 0.0;
   sink += Q_new.Size() ? Q_new[Q_new.Size() - 1] : 0.0;
   if (Q_per_node.size() && Q_per_node[0].Size()) { sink += Q_per_node[0][0]; }
   sink += rhs_vol.Size() ? rhs_vol[0] : 0.0;   // R-005: keep ComputeVolumeRHS live
   std::printf("(checksum %.3e — ignore)\n", static_cast<double>(sink));
   return 0;
}
