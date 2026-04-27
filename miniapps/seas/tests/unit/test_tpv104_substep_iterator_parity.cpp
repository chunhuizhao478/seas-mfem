// Round-7 R-602/R-603: Tpv104SubStepIterator parity contract at O=1.
//
// Contract:
//   At sub-step quadrature O = 1 with deltaT = {dt_macro},
//   time_weights = {1.0}, the new `AdvanceWithSubStepStates` method
//   (per-sub-step pointwise Q) and the legacy `Advance` method
//   (time-integrated I) must produce IDENTICAL DOFData updates and
//   I_imp_flat outputs when:
//      Q_pointwise[0]  ==  I_flat / dt_macro
//   on every QP.  This is the SeisSol-equivalent reduction at O = 1
//   — both APIs reduce to a single ComputeStageState call on Q̄.
//
//   We verify this by:
//     1. Picking a representative DOFData state (TPV102 locked-fault
//        snapshot at rest + some perturbation in V_ini).
//     2. Running Advance on a fixed I_flat.
//     3. Running AdvanceWithSubStepStates with Q_pointwise[0] = I_flat / dt
//        on a FRESH dof_data copy.
//     4. Asserting both DOFData snapshots and both I_imp outputs match.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv104_substep_iterator.hpp"
#include "../../friction/slip_law_srw_psi.hpp"
#include "../../config/tpv104_params.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

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

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

namespace
{

// Configure a representative TPV104-style DOFData at rest with a small
// perturbation so the friction solve is exercised non-trivially.
DOFData MakeTPV104DOF()
{
   DOFData d;
   d.Zp_plus  = TPV104Params::Zp;
   d.Zp_minus = TPV104Params::Zp;
   d.Zs_plus  = TPV104Params::Zs;
   d.Zs_minus = TPV104Params::Zs;
   d.eta_p    = TPV104Params::eta_p;
   d.eta_s    = TPV104Params::eta_s;
   d.sigma_n0 = TPV104Params::sigma_n;
   d.tau1_0   = 0.0;                      // pure strike-slip
   d.tau2_0   = 29.38e6;                  // representative pre-stress
   d.sigma_n_nuc = 0.0;
   d.tau1_nuc    = 0.0;
   d.tau2_nuc    = 0.0;                   // nucleation off for parity test
   d.a   = TPV104Params::a_in;
   d.Dc  = TPV104Params::L;
   d.psi = 0.564;                          // representative steady-state ψ
   d.slip_rate = 0.0;
   d.V1 = 0.0; d.V2 = 0.0;
   d.slip1 = 0.0; d.slip2 = 0.0;
   d.tau1_corr = 0.0; d.tau2_corr = 0.0;
   d.sigma_n_corr = 0.0;
   return d;
}

real_t MaxAbsDiff(const std::vector<real_t> &a,
                  const std::vector<real_t> &b)
{
   const size_t n = std::min(a.size(), b.size());
   real_t m = 0.0;
   for (size_t i = 0; i < n; i++)
   {
      real_t d = std::abs(a[i] - b[i]);
      if (d > m) { m = d; }
   }
   return m;
}

real_t MaxDOFFieldDiff(const DOFData &a, const DOFData &b)
{
   real_t m = 0.0;
   m = std::max(m, std::abs(a.psi          - b.psi));
   m = std::max(m, std::abs(a.slip_rate    - b.slip_rate));
   m = std::max(m, std::abs(a.V1           - b.V1));
   m = std::max(m, std::abs(a.V2           - b.V2));
   m = std::max(m, std::abs(a.slip1        - b.slip1));
   m = std::max(m, std::abs(a.slip2        - b.slip2));
   m = std::max(m, std::abs(a.tau1_corr    - b.tau1_corr));
   m = std::max(m, std::abs(a.tau2_corr    - b.tau2_corr));
   m = std::max(m, std::abs(a.sigma_n_corr - b.sigma_n_corr));
   return m;
}

} // anonymous

int main()
{
   std::cout << "\n=== Round-7 R-602/R-603: SubStep iterator parity at O=1 ===\n";

   FaultFaceFlux flux(TPV104Params::rho, TPV104Params::cp, TPV104Params::cs);

   // SlipLawSRWPsi for global friction scalars (b, V0, f0, muW); a_scalar
   // and V_w_default are not consulted in production-mode SRW free-function
   // dispatch (per-DOF a / V_w shadow them), so any nonzero placeholder
   // value is fine.
   SlipLawSRWPsi state_evo(
      TPV104Params::a_in, TPV104Params::b, TPV104Params::V0,
      TPV104Params::f0, TPV104Params::muW, TPV104Params::V_w_out);
   state_evo.SetProductionMode();

   Tpv104SubStepIterator iter(flux, state_evo);
   const real_t dt_macro = 1e-4;
   iter.SetSubSteps({dt_macro}, {1.0});

   // Single-DOF fixture.
   const int n = 1;
   std::vector<DOFData> dof_data_legacy = { MakeTPV104DOF() };
   std::vector<DOFData> dof_data_substep = { MakeTPV104DOF() };
   std::vector<Vector> fault_coords(n, Vector(3));
   fault_coords[0] = 0.0;
   fault_coords[0](2) = -7500.0;            // hypocenter z
   std::vector<real_t> V_w(n, TPV104Params::V_w_in);

   // Construct a representative fault-local-frame Q̄ at the fault QP.
   // Pure strike-slip rest state with a small VZ kick.
   std::vector<real_t> Q_avg_plus(NUM_STATE, 0.0);
   std::vector<real_t> Q_avg_minus(NUM_STATE, 0.0);
   Q_avg_plus[VZ]  = +0.001;     // small + side strike velocity (m/s)
   Q_avg_minus[VZ] = -0.001;     // small − side strike velocity (m/s)

   // Legacy I_flat = dt * Q̄.  AdvanceWithSubStepStates expects the
   // pointwise Q̄ directly (already at sub-step time, no integration scaling).
   std::vector<real_t> I_plus_flat(NUM_STATE * n, 0.0);
   std::vector<real_t> I_minus_flat(NUM_STATE * n, 0.0);
   for (int c = 0; c < NUM_STATE; c++)
   {
      I_plus_flat [c]  = dt_macro * Q_avg_plus[c];
      I_minus_flat[c]  = dt_macro * Q_avg_minus[c];
   }

   // Outputs.
   std::vector<real_t> I_imp_plus_legacy (NUM_STATE * n, 0.0);
   std::vector<real_t> I_imp_minus_legacy(NUM_STATE * n, 0.0);
   std::vector<real_t> I_imp_plus_substep (NUM_STATE * n, 0.0);
   std::vector<real_t> I_imp_minus_substep(NUM_STATE * n, 0.0);

   // (1) Legacy Advance — operates on I (time-integrated) and divides
   //     internally by dt_macro to get Q̄.
   iter.Advance(dof_data_legacy, fault_coords, V_w,
                I_plus_flat.data(), I_minus_flat.data(),
                dt_macro, /*t_macro_start=*/0.0,
                I_imp_plus_legacy.data(), I_imp_minus_legacy.data(),
                FrictionSolver::Method::NewtonRaphsonStable);

   // (2) AdvanceWithSubStepStates — operates on per-sub-step pointwise Q̄
   //     directly.  At O = 1, single sub-step with the SAME Q̄ as the
   //     legacy path's I/dt_macro reduction.
   std::vector<std::vector<real_t>> Q_plus_per({Q_avg_plus});
   std::vector<std::vector<real_t>> Q_minus_per({Q_avg_minus});
   iter.AdvanceWithSubStepStates(dof_data_substep, fault_coords, V_w,
                                 Q_plus_per, Q_minus_per,
                                 dt_macro, /*t_macro_start=*/0.0,
                                 I_imp_plus_substep.data(),
                                 I_imp_minus_substep.data(),
                                 FrictionSolver::Method::NewtonRaphsonStable);

   // ------------------------------------------------------------------
   // Gate 1: I_imp outputs match.
   // ------------------------------------------------------------------
   std::cout << "\n-- Gate 1: I_imp outputs at O=1 --\n";
   const real_t err_I_plus  =
      MaxAbsDiff(I_imp_plus_legacy,  I_imp_plus_substep);
   const real_t err_I_minus =
      MaxAbsDiff(I_imp_minus_legacy, I_imp_minus_substep);
   std::cout << "  max |I_imp_plus_legacy − I_imp_plus_substep|   = "
             << std::scientific << err_I_plus  << "\n";
   std::cout << "  max |I_imp_minus_legacy − I_imp_minus_substep| = "
             << std::scientific << err_I_minus << "\n";
   TEST_LE(err_I_plus,  0.0,
           "I_imp_plus parity at O=1 (bit-equal)");
   TEST_LE(err_I_minus, 0.0,
           "I_imp_minus parity at O=1 (bit-equal)");

   // ------------------------------------------------------------------
   // Gate 2: DOFData updates match.
   // ------------------------------------------------------------------
   std::cout << "\n-- Gate 2: DOFData updates at O=1 --\n";
   const real_t err_dof =
      MaxDOFFieldDiff(dof_data_legacy[0], dof_data_substep[0]);
   std::cout << "  max DOFData field diff = " << std::scientific
             << err_dof << "\n";
   std::cout << "    psi(legacy/substep)        = "
             << dof_data_legacy[0].psi << " / "
             << dof_data_substep[0].psi << "\n";
   std::cout << "    slip_rate(legacy/substep)  = "
             << dof_data_legacy[0].slip_rate << " / "
             << dof_data_substep[0].slip_rate << "\n";
   std::cout << "    V2(legacy/substep)         = "
             << dof_data_legacy[0].V2 << " / "
             << dof_data_substep[0].V2 << "\n";
   TEST_LE(err_dof, 0.0,
           "DOFData fields parity at O=1 (bit-equal)");

   // ------------------------------------------------------------------
   // Gate 3: Both APIs called with reasonable-magnitude inputs do not
   //         produce NaN/Inf in any output.  Defensive smoke.
   // ------------------------------------------------------------------
   std::cout << "\n-- Gate 3: outputs finite --\n";
   bool all_finite = true;
   for (real_t v : I_imp_plus_substep)  { if (!std::isfinite(v)) { all_finite = false; } }
   for (real_t v : I_imp_minus_substep) { if (!std::isfinite(v)) { all_finite = false; } }
   if (!std::isfinite(dof_data_substep[0].slip_rate) ||
       !std::isfinite(dof_data_substep[0].V1) ||
       !std::isfinite(dof_data_substep[0].V2) ||
       !std::isfinite(dof_data_substep[0].psi))
   {
      all_finite = false;
   }
   TEST_ASSERT(all_finite,
               "AdvanceWithSubStepStates produces finite outputs");

   // ------------------------------------------------------------------
   // R-1003 Test 5: ADER-O ∈ {2, 3, 4} stability + cross-order consistency.
   //
   // Pre-R-1003 the `*_parity.cpp` test family only exercised
   // `ader_order = 2` (existing test_tpv104_substep_iterator.cpp:192) and
   // implicitly the O=1 reduction above.  R-1003's MPI extension touches
   // the per-sub-step predictor and the canonical-frame routing on
   // shared-fault QPs; even at np=1 it pays to add explicit O=3 / O=4
   // coverage of the iterator path so latent CK-recursion or cumulative-
   // accumulator bugs surface here rather than at np>1 production runs.
   //
   // Contract on a *constant-Q* fixture (Q̄ identical across sub-steps,
   // no nucleation, no time-dependent perturbation):
   //   - For every O ∈ {2, 3, 4}, AdvanceWithSubStepStates produces finite
   //     I_imp_*_flat and finite DOFData (no NaN/Inf).
   //   - Writing the same Q̄ at every sub-step is mathematically equivalent
   //     to a single sub-step of size dt_macro: the per-sub-step ψ /
   //     slip-rate trajectory traverses the same friction curve in O steps
   //     of dt_macro/O.  Output V1/V2 (last-sub-step writeback) and
   //     accumulated I_imp at every O must match the O=1 baseline within
   //     the per-sub-step ψ-update rounding budget (analytic update is
   //     ≤ O(eps) in floating point per call, accumulating across O calls).
   //   - The cross-order |V2| spread is bounded — cheap regression net
   //     for a CK predictor that mismatches sub-step time nodes.
   // ------------------------------------------------------------------
   std::cout << "\n-- R-1003 Test 5: ADER orders 2, 3, 4 (np=1 stability) --\n";

   // Re-use the same Q̄ fixture from above (Q_avg_plus, Q_avg_minus,
   // dt_macro).  O=1 baseline already in dof_data_substep / I_imp_*_substep.
   const std::vector<int> orders_to_test = {2, 3, 4};
   const real_t order_consistency_tol = 1e-7;  // ψ-update analytic rounding
                                               // accumulates O calls × eps,
                                               // ≈ O · 1e-15 at most for
                                               // representative TPV104 ψ;
                                               // 1e-7 is a generous net
                                               // that still catches O(1)
                                               // CK / accumulator regressions.

   for (int O : orders_to_test)
   {
      std::cout << "\n  ADER-O" << O << ":\n";

      // Equal-width sub-steps with equal weights — same configuration the
      // driver wires at L1404 for the production substep dispatch path.
      std::vector<real_t> deltaT(O, dt_macro / static_cast<real_t>(O));
      std::vector<real_t> weights(O, 1.0 / static_cast<real_t>(O));
      iter.SetSubSteps(deltaT, weights);

      // Identical Q̄ at every sub-step.  This is the constant-Q path: the
      // CK predictor's per-sub-step Q reduces to the same time-average
      // ANY ADER order produces from a constant input, by inspection of
      // the Taylor recursion (D^{(k)} = 0 for k ≥ 1 on a spatially
      // constant Q, so Q(τ) = D^{(0)} = Q for every τ).  The driver does
      // not exercise this path — production Q has nontrivial spatial
      // variation — but the iterator's API contract is independent of
      // the predictor: feed it identical Q̄ across O sub-steps and the
      // accumulator must equal the O=1 result.
      std::vector<std::vector<real_t>> Q_plus_O(O, Q_avg_plus);
      std::vector<std::vector<real_t>> Q_minus_O(O, Q_avg_minus);

      std::vector<DOFData> dof_data_O = { MakeTPV104DOF() };
      std::vector<real_t> I_imp_plus_O (NUM_STATE * n, 0.0);
      std::vector<real_t> I_imp_minus_O(NUM_STATE * n, 0.0);

      iter.AdvanceWithSubStepStates(
         dof_data_O, fault_coords, V_w,
         Q_plus_O, Q_minus_O,
         dt_macro, /*t_macro_start=*/0.0,
         I_imp_plus_O.data(), I_imp_minus_O.data(),
         FrictionSolver::Method::NewtonRaphsonStable);

      // Gate 5a: outputs finite for this O.
      bool finite_O = true;
      for (real_t v : I_imp_plus_O)  { if (!std::isfinite(v)) { finite_O = false; } }
      for (real_t v : I_imp_minus_O) { if (!std::isfinite(v)) { finite_O = false; } }
      if (!std::isfinite(dof_data_O[0].psi)       ||
          !std::isfinite(dof_data_O[0].slip_rate) ||
          !std::isfinite(dof_data_O[0].V1)        ||
          !std::isfinite(dof_data_O[0].V2))
      { finite_O = false; }
      TEST_ASSERT(finite_O,
                  "ADER-O" + std::to_string(O) +
                  ": iterator produces finite outputs");

      // Gate 5b: I_imp accumulator matches O=1 baseline on constant Q̄.
      // (Σ_o weights[o] · dt_macro · Q_imp = dt_macro · Q_imp for any O
      // when Q_imp is the same at every sub-step.)
      const real_t err_I_plus_O  =
         MaxAbsDiff(I_imp_plus_O,  I_imp_plus_substep);
      const real_t err_I_minus_O =
         MaxAbsDiff(I_imp_minus_O, I_imp_minus_substep);
      std::cout << "    max |I_imp_plus(O="  << O << ") − I_imp_plus(O=1)|  = "
                << std::scientific << err_I_plus_O  << "\n";
      std::cout << "    max |I_imp_minus(O=" << O << ") − I_imp_minus(O=1)| = "
                << std::scientific << err_I_minus_O << "\n";
      TEST_LE(err_I_plus_O,  order_consistency_tol,
              "ADER-O" + std::to_string(O) +
              ": I_imp_plus accumulator matches O=1 on constant Q̄");
      TEST_LE(err_I_minus_O, order_consistency_tol,
              "ADER-O" + std::to_string(O) +
              ": I_imp_minus accumulator matches O=1 on constant Q̄");

      // Gate 5c: WriteBackState V1/V2 match O=1 baseline.  WriteBackState
      // commits the LAST sub-step's V1/V2/τ*_corr/σ_n_corr — at constant
      // Q̄ across O sub-steps, the last-step result should match the O=1
      // single-step result within the analytic ψ-update budget.
      const real_t v1_diff = std::abs(dof_data_O[0].V1
                                      - dof_data_substep[0].V1);
      const real_t v2_diff = std::abs(dof_data_O[0].V2
                                      - dof_data_substep[0].V2);
      const real_t psi_diff = std::abs(dof_data_O[0].psi
                                       - dof_data_substep[0].psi);
      std::cout << "    |ΔV1| = " << std::scientific << v1_diff
                << ",  |ΔV2| = " << v2_diff
                << ",  |Δψ| = "  << psi_diff << "\n";
      TEST_LE(v2_diff, order_consistency_tol,
              "ADER-O" + std::to_string(O) +
              ": V2 (last-sub-step writeback) matches O=1 on constant Q̄");
      // V1 should remain ≈ 0 (pure strike-slip); use absolute bound.
      TEST_LE(std::abs(dof_data_O[0].V1), order_consistency_tol,
              "ADER-O" + std::to_string(O) +
              ": V1 stays ≈ 0 for pure strike-slip (no dip excitation)");
      TEST_LE(psi_diff, order_consistency_tol,
              "ADER-O" + std::to_string(O) +
              ": ψ converges to O=1 on constant Q̄ within analytic budget");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
