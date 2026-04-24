// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Unit tests for dynamic/tpv104_nucleation.hpp (§4.10 Step 6 gates
// T_TPV104_NUC_1..5).

#include "test_macros.hpp"
#include "../../dynamic/tpv104_nucleation.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv104_params.hpp"
#include "../../config/tpv102_params.hpp"

#include <cstdlib>
#include <cmath>
#include <random>

using namespace mfem;
using namespace mfem::seas;

// ---------------------------------------------------------------------------
// Standalone inline reference (R-003 pattern): straight-line arithmetic for
// smoothStep / smoothStepIncrement, no calls into production.
// ---------------------------------------------------------------------------
static real_t reference_smooth_step(real_t t, real_t t0)
{
   if (t <= 0.0) { return 0.0; }
   if (t < t0)
   {
      const real_t tau = t - t0;
      return std::exp(tau * tau / (t * (t - 2.0 * t0)));
   }
   return 1.0;
}

static real_t reference_smooth_step_increment(real_t t, real_t dt, real_t t0)
{
   return reference_smooth_step(t, t0) - reference_smooth_step(t - dt, t0);
}

// ----------------------------------------------------------------------------
// T_TPV104_NUC_1 — SmoothStepIncrement matches standalone reference to 1e-15.
// ----------------------------------------------------------------------------
void TestSmoothStepIncrementByteMatch()
{
   std::cout << "\n[T_TPV104_NUC_1] SmoothStepIncrement byte-match\n";

   std::mt19937 rng(20260424u);
   std::uniform_real_distribution<real_t> t_rng(-0.5, 2.0);
   std::uniform_real_distribution<real_t> dt_rng(1e-6, 0.5);
   std::uniform_real_distribution<real_t> t0_rng(0.5, 2.0);

   const int N = 10000;
   int ok = 0;
   real_t max_rel = 0.0;
   for (int i = 0; i < N; ++i)
   {
      const real_t t   = t_rng(rng);
      const real_t dt  = dt_rng(rng);
      const real_t t0  = t0_rng(rng);
      const real_t ref = reference_smooth_step_increment(t, dt, t0);
      const real_t got = SmoothStepIncrement_TPV104(t, dt, t0);

      if (ref == 0.0 && got == 0.0) { ++ok; continue; }
      const real_t denom = std::max<real_t>(std::abs(ref), 1.0);
      const real_t rel = std::abs(got - ref) / denom;
      max_rel = std::max(max_rel, rel);
      if (rel < 1e-15) { ++ok; }
   }
   std::cout << "  max rel err = " << max_rel << " / " << N << " samples\n";
   TEST_ASSERT(ok == N,
               "All 10000 SmoothStepIncrement samples < 1e-15 rel");
}

// ----------------------------------------------------------------------------
// T_TPV104_NUC_2 — accumulator telescopes to the full perturbation.
// Run dt=1e-3, T_nuc=1.0, 1000 sub-steps; final tau2_nuc[hypo] ≈ Δτ₀·F(r).
// ----------------------------------------------------------------------------
void TestAccumulatorTelescopes()
{
   std::cout << "\n[T_TPV104_NUC_2] accumulator telescopes to full perturbation\n";

   // Single-QP fixture at the hypocenter.
   std::vector<DOFData> dof_data(1);
   std::vector<Vector> fault_coords(1);
   fault_coords[0].SetSize(3);
   fault_coords[0](0) = TPV104Params::hypo_along_strike;
   fault_coords[0](1) = 0.0;
   // fault_coords uses z<0 = depth; hypo_down_dip is the magnitude.
   fault_coords[0](2) = -TPV104Params::hypo_down_dip;

   ResetNucleationAccumulator_TPV104(dof_data);
   TEST_NEAR(dof_data[0].tau2_nuc, 0.0, 1e-30, "accumulator starts at 0");

   // R3-002 (review round 3): derive dt from an integer step count so
   // the loop deterministically terminates at t_end = nuc_T.  The
   // previous form `int(nuc_T / 1e-3)` truncates to 999 under IEEE
   // double rounding on some toolchains (1.0 / 0.001000... = 999.999...),
   // placing the telescoping sum precisely at the 1e-6 tolerance
   // boundary — a one-ULP flip can flip PASS↔FAIL.
   constexpr int n_steps = 1000;
   const real_t dt = TPV104Params::nuc_T / n_steps;

   for (int k = 1; k <= n_steps; ++k)
   {
      const real_t t_end = k * dt;
      ApplyNucleationIncremental_TPV104(dof_data, fault_coords, t_end, dt);
   }

   // At the hypocenter: r = 0, F(0) = exp(0/(0-R²)) = exp(0) = 1.
   const real_t expected = TPV104Params::nuc_dtau * 1.0;
   const real_t got = dof_data[0].tau2_nuc;
   const real_t rel = std::abs(got - expected) / std::abs(expected);
   std::cout << "  tau2_nuc(hypo) after T_nuc = " << got
             << " (expected " << expected << ", rel err " << rel << ")\n";
   TEST_ASSERT(rel < 1e-6,
               "tau2_nuc[hypo] ≈ Δτ₀·F(0) = Δτ₀·1 to 1e-6 rel");

   // R3-001 (review round 3): TPV104 is pure strike-slip — tau1_nuc and
   // sigma_n_nuc must remain at 0 after every accumulator call.  Any
   // copy-paste bug in a Step-7 refactor that routes the increment to
   // the wrong channel would pass tau2_nuc tests and silently break
   // Phase-3 probes; this invariant catches it immediately.
   TEST_NEAR(dof_data[0].tau1_nuc, 0.0, 1e-30,
             "tau1_nuc remains at 0 (no dip-shear nucleation)");
   TEST_NEAR(dof_data[0].sigma_n_nuc, 0.0, 1e-30,
             "sigma_n_nuc remains at 0 (no normal-stress nucleation)");
}

// ----------------------------------------------------------------------------
// T_TPV104_NUC_3 — spatial factor F(r) matches TPV102's NucleationSpatial.
// Both return 0 at r ≥ R and exp(r²/(r²−R²)) for r < R.
// ----------------------------------------------------------------------------
void TestSpatialFactorMatchesTPV102()
{
   std::cout << "\n[T_TPV104_NUC_3] NucleationSpatial_TPV104 ≡ TPV102's form\n";

   const std::vector<real_t> r_samples = {0.0, 500.0, 1000.0, 1500.0,
                                          2000.0, 2500.0, 2999.0,
                                          3000.0, 3001.0, 5000.0};
   int ok = 0;
   for (real_t r : r_samples)
   {
      const real_t f104 = NucleationSpatial_TPV104(r);
      const real_t f102 = NucleationSpatial(r);   // TPV102 version
      // Both use R = 3 km — identical formula.
      if (f104 == f102) { ++ok; }
      else
      {
         std::cerr << "  FAIL r=" << r
                   << " f104=" << f104 << " f102=" << f102 << "\n";
      }
   }
   TEST_ASSERT(ok == (int)r_samples.size(),
               "F(r) bit-identical between TPV104 and TPV102");
}

// ----------------------------------------------------------------------------
// T_TPV104_NUC_4 — TPV102 baseline unaffected.  Running
// ApplyNucleationPrestress (TPV102 overwrite pattern) after Step 6 lands
// produces the same tau2_nuc as calling it pre-Step-6.
// Concretely: Step 6 adds new symbols; it does not mutate any TPV102
// code path.  We exercise TPV102's overwrite pattern to confirm it is
// deterministic and un-tangled from the new accumulator.
// ----------------------------------------------------------------------------
void TestTPV102Unaffected()
{
   std::cout << "\n[T_TPV104_NUC_4] TPV102 overwrite pattern unaffected\n";

   std::vector<DOFData> dof_data(1);
   std::vector<Vector> fault_coords(1);
   fault_coords[0].SetSize(3);
   fault_coords[0](0) = TPV102Params::hypo_along_strike;
   fault_coords[0](1) = 0.0;
   fault_coords[0](2) = -TPV102Params::hypo_down_dip;

   // Overwrite pattern (TPV102): tau2_nuc is replaced each call.
   ApplyNucleationPrestress(dof_data, fault_coords, /*t=*/0.5);
   const real_t val_at_0_5 = dof_data[0].tau2_nuc;

   ApplyNucleationPrestress(dof_data, fault_coords, /*t=*/1.0);
   const real_t val_at_1_0 = dof_data[0].tau2_nuc;

   // Overwrite semantics: val_at_1_0 is Δτ_TPV102 · F(0) · G(1) = 25 MPa,
   // NOT an accumulation of the 0.5 s call.
   const real_t expected_full = TPV102Params::nuc_dtau * 1.0 * 1.0;
   TEST_NEAR(val_at_1_0, expected_full, 1e-3,
             "TPV102 ApplyNucleationPrestress(t=T_nuc) = Δτ·F·G = 25 MPa");
   TEST_ASSERT(val_at_0_5 != val_at_1_0,
               "TPV102 overwrite pattern replaces (not accumulates)");

   // Confirm TPV104 and TPV102 nucleation amplitudes still differ as
   // expected — this guards against an accidental constant-sharing bug.
   TEST_ASSERT(TPV104Params::nuc_dtau != TPV102Params::nuc_dtau,
               "TPV104 Δτ0 (45 MPa) ≠ TPV102 Δτ0 (25 MPa)");
}

// ----------------------------------------------------------------------------
// T_TPV104_NUC_5 — guard bounds.
//   t ≤ 0:        increment = 0
//   t ≥ T_nuc:    increment = 0
//   sum over [0, T_nuc] == Δτ · F(r)
// ----------------------------------------------------------------------------
void TestGuardBounds()
{
   std::cout << "\n[T_TPV104_NUC_5] smoothStepIncrement guard bounds\n";

   const real_t t0 = TPV104Params::nuc_T;  // = 1.0

   // t ≤ 0 branch: smoothStep(t, t0) = 0 for both endpoints.
   const real_t inc_neg = SmoothStepIncrement_TPV104(-0.1, 0.01, t0);
   TEST_NEAR(inc_neg, 0.0, 1e-30,
             "increment = 0 for t ≤ 0 (both endpoints < 0)");

   // Post-ramp branch: smoothStep(t, t0) = 1 for both endpoints.
   const real_t inc_post = SmoothStepIncrement_TPV104(1.5, 0.01, t0);
   TEST_NEAR(inc_post, 0.0, 1e-30,
             "increment = 0 for t > T_nuc (both endpoints ≥ t0)");

   // Crossing t = 0 from negative: first call spans [t-dt, t] partially
   // across 0.  Increment = smoothStep(t, t0) − 0 = smoothStep(t, t0).
   const real_t t_cross = 0.005, dt_cross = 0.01;
   const real_t inc_cross = SmoothStepIncrement_TPV104(t_cross, dt_cross, t0);
   const real_t ss_t_cross = SmoothStep_TPV104(t_cross, t0);
   TEST_NEAR(inc_cross, ss_t_cross, 1e-15,
             "crossing t=0 → increment = smoothStep(t, t0)");

   // Crossing t = T_nuc from below: last call spans [t-dt, t] across t0.
   //   increment = 1 − smoothStep(t - dt, t0)  (nonzero ramp-residual).
   const real_t t_above = t0 + 0.005, dt_above = 0.01;
   const real_t inc_above = SmoothStepIncrement_TPV104(t_above, dt_above, t0);
   const real_t residual = 1.0 - SmoothStep_TPV104(t_above - dt_above, t0);
   TEST_NEAR(inc_above, residual, 1e-15,
             "crossing t=T_nuc → increment = 1 − smoothStep(t-dt, t0)");

   // Telescoping sum over [0, T_nuc]: accumulator sums to exactly
   // smoothStep(T_nuc, T_nuc) - smoothStep(0, T_nuc) = 1 - 0 = 1.
   // Then at the hypocenter (F=1), tau2_nuc ends at Δτ₀.
   std::vector<DOFData> dof_data(1);
   std::vector<Vector> fault_coords(1);
   fault_coords[0].SetSize(3);
   fault_coords[0](0) = TPV104Params::hypo_along_strike;
   fault_coords[0](1) = 0.0;
   fault_coords[0](2) = -TPV104Params::hypo_down_dip;
   ResetNucleationAccumulator_TPV104(dof_data);

   const int n = 500;
   const real_t dt = t0 / n;
   for (int k = 1; k <= n; ++k)
   {
      ApplyNucleationIncremental_TPV104(dof_data, fault_coords, k * dt, dt);
   }
   const real_t expected = TPV104Params::nuc_dtau;
   const real_t rel = std::abs(dof_data[0].tau2_nuc - expected) /
                      std::abs(expected);
   TEST_ASSERT(rel < 1e-6,
               "telescoping sum over [0, T_nuc] = Δτ₀·F(0) to 1e-6 rel");

   // After T_nuc, further sub-steps must not grow the accumulator.
   const real_t frozen = dof_data[0].tau2_nuc;
   ApplyNucleationIncremental_TPV104(dof_data, fault_coords, t0 + 0.5, 0.1);
   ApplyNucleationIncremental_TPV104(dof_data, fault_coords, t0 + 1.0, 0.1);
   TEST_NEAR(dof_data[0].tau2_nuc, frozen, 1e-30,
             "accumulator frozen after T_nuc (no spurious growth)");

   // R3-001 invariant (also guards the pre-ramp and post-ramp regimes):
   TEST_NEAR(dof_data[0].tau1_nuc, 0.0, 1e-30,
             "tau1_nuc remains at 0 across ramp boundaries");
   TEST_NEAR(dof_data[0].sigma_n_nuc, 0.0, 1e-30,
             "sigma_n_nuc remains at 0 across ramp boundaries");
}

// ----------------------------------------------------------------------------
// R3-001 dedicated sentinel test — the TPV104 strike-slip invariant.
// Pre-populate `tau1_nuc` and `sigma_n_nuc` with sentinel values that
// would be overwritten if the accumulator ever touched them.  Run a
// full ramp (200 sub-steps over 2·T_nuc) and confirm the sentinels are
// preserved bit-exactly.  A copy-paste bug that routed the increment to
// tau1_nuc would overwrite the sentinel and fire this test; a bug that
// added to sigma_n_nuc would do the same.
// ----------------------------------------------------------------------------
void TestNoDipOrNormalStressNucleation()
{
   std::cout << "\n[R3-001] strike-slip invariant (sentinel check)\n";

   std::vector<DOFData> dof_data(1);
   std::vector<Vector> fault_coords(1);
   fault_coords[0].SetSize(3);
   fault_coords[0](0) = 0.0;                                // at hypocenter
   fault_coords[0](1) = 0.0;
   fault_coords[0](2) = -TPV104Params::hypo_down_dip;

   // Sentinel values: deliberately extreme so any overwrite is obvious.
   const real_t sent_t1 = 1.0e99;
   const real_t sent_sn = 2.0e99;
   dof_data[0].tau1_nuc    = sent_t1;
   dof_data[0].sigma_n_nuc = sent_sn;
   dof_data[0].tau2_nuc    = 0.0;

   // Run a full ramp and then past the end (t > T_nuc) to cover all
   // three smoothStep regimes: ramp-entry, ramp-interior, ramp-exit.
   const real_t dt = 0.01;
   for (int k = 1; k <= 200; ++k)
   {
      ApplyNucleationIncremental_TPV104(dof_data, fault_coords,
                                        k * dt, dt);
   }

   // Accumulator writes only to tau2_nuc.  The sentinels must be
   // BIT-EXACTLY preserved (tolerance 0.0).
   TEST_ASSERT(dof_data[0].tau1_nuc == sent_t1,
               "tau1_nuc sentinel preserved (no dip-shear write)");
   TEST_ASSERT(dof_data[0].sigma_n_nuc == sent_sn,
               "sigma_n_nuc sentinel preserved (no normal-stress write)");

   // Sanity: tau2_nuc did receive the expected ramp.
   TEST_ASSERT(dof_data[0].tau2_nuc > 0.0,
               "tau2_nuc received the ramp (accumulator is live)");
   const real_t expected = TPV104Params::nuc_dtau;
   const real_t rel = std::abs(dof_data[0].tau2_nuc - expected)
                      / std::abs(expected);
   TEST_ASSERT(rel < 1e-6,
               "tau2_nuc reaches Δτ₀ at t = 2·T_nuc");
}

int main(int argc, char *argv[])
{
   TestSmoothStepIncrementByteMatch();
   TestAccumulatorTelescopes();
   TestSpatialFactorMatchesTPV102();
   TestTPV102Unaffected();
   TestGuardBounds();
   TestNoDipOrNormalStressNucleation();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
