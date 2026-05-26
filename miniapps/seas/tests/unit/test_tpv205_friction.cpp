// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Unit tests for dynamic/tpv205_friction.hpp covering REVIEW.md
// findings R-002 (strength-barrier slip-weakening collapse) and R-003
// (tensile σ_n unlocking the barrier).
//
// Both findings are CRITICAL: a barrier QP that ever slips for any
// reason (R-003 enables it via tensile σ_n; numerical asymmetry can
// also seed it) silently weakens to μ_d = 0.525 (R-002), at which point
// the rupture front escapes the 30 km × 15 km rupture area and the
// SCEC TPV5 benchmark comparison is meaningless.

#include "test_macros.hpp"
#include "../../dynamic/tpv205_friction.hpp"
#include "../../config/tpv205_params.hpp"

#include <cmath>

using namespace mfem;
using namespace mfem::seas;

// ---------------------------------------------------------------------------
// R-002: strength barrier MUST NOT slip-weaken when δ ≥ d_c.
// Pre-fix the LSWFrictionCoefficient_TPV205 returned μ_d (0.525) once
// δ exceeded d_c, regardless of whether the QP was inside the rupture
// area or in the strength barrier (μ_s = 10000).
// ---------------------------------------------------------------------------
static void test_R002_BarrierDoesNotCollapseAtLargeSlip()
{
   const real_t mu_s_barrier = TPV205Params::mu_s_barrier;  // 10000
   const real_t mu_d         = TPV205Params::mu_d;          // 0.525
   const real_t d_c          = TPV205Params::d_c;           // 0.40 m

   // δ = 0 in the barrier zone — should return μ_s_barrier.
   {
      const real_t mu = LSWFrictionCoefficient_TPV205(
                           /*delta=*/0.0, mu_s_barrier, mu_d, d_c);
      TEST_NEAR(mu, mu_s_barrier, 1e-12,
                "R-002: barrier μ at δ=0 must equal μ_s_barrier");
   }

   // δ = 0.5 d_c (mid-weakening regime in the rupture area) — barrier
   // QP must STAY at μ_s_barrier, not interpolate.
   {
      const real_t mu = LSWFrictionCoefficient_TPV205(
                           /*delta=*/0.2, mu_s_barrier, mu_d, d_c);
      TEST_ASSERT(mu >= 0.5 * mu_s_barrier,
                  "R-002: barrier μ at δ=0.2 m must stay barrier-sized");
   }

   // δ = 2 d_c (well past d_c; pre-fix this branch returned μ_d).
   {
      const real_t mu = LSWFrictionCoefficient_TPV205(
                           /*delta=*/0.8, mu_s_barrier, mu_d, d_c);
      TEST_ASSERT(mu >= 0.5 * mu_s_barrier,
                  "R-002: barrier μ at δ=0.8 m must NOT collapse to μ_d");
   }

   // δ = 1 m (the reviewer's exact test value).
   {
      const real_t mu = LSWFrictionCoefficient_TPV205(
                           /*delta=*/1.0, mu_s_barrier, mu_d, d_c);
      TEST_ASSERT(mu >= 0.5 * mu_s_barrier,
                  "R-002: barrier μ at δ=1.0 m must stay barrier-sized");
   }

   // INSIDE the rupture area (μ_s = 0.677), the standard LSW law
   // applies — μ at δ ≥ d_c MUST equal μ_d.  This sanity check makes
   // sure the R-002 fix did not over-clamp the rupture-area path.
   {
      const real_t mu_s_in = TPV205Params::mu_s;  // 0.677
      const real_t mu_at_dc =
         LSWFrictionCoefficient_TPV205(d_c, mu_s_in, mu_d, d_c);
      TEST_NEAR(mu_at_dc, mu_d, 1e-12,
                "R-002 sanity: rupture-area μ at δ=d_c must equal μ_d");

      const real_t mu_at_2dc =
         LSWFrictionCoefficient_TPV205(2.0 * d_c, mu_s_in, mu_d, d_c);
      TEST_NEAR(mu_at_2dc, mu_d, 1e-12,
                "R-002 sanity: rupture-area μ at δ=2 d_c must equal μ_d");

      const real_t mu_mid =
         LSWFrictionCoefficient_TPV205(0.5 * d_c, mu_s_in, mu_d, d_c);
      const real_t expected = mu_s_in - 0.5 * (mu_s_in - mu_d);
      TEST_NEAR(mu_mid, expected, 1e-12,
                "R-002 sanity: rupture-area linear weakening at δ=0.5 d_c");
   }
}

// ---------------------------------------------------------------------------
// R-003: tensile σ_n MUST NOT unlock the strength barrier.
// Pre-fix, max(σ_n, 0) clamping zeroed τ_strength under tensile
// transients in the barrier zone, allowing V ≈ 15 m/s; combined with
// R-002 this permanently collapsed the barrier within ~26 ADER steps.
// ---------------------------------------------------------------------------
static void test_R003_BarrierLockedUnderTensileSigmaN()
{
   const real_t eta_s = TPV205Params::eta_s;
   const real_t mu_s_barrier = TPV205Params::mu_s_barrier;

   // Barrier QP with tensile σ_n_total = -150 MPa, |τ_total| = 70 MPa.
   // Pre-fix: τ_strength = 10000 × max(-150e6, 0) = 0, so V = 70 MPa /
   // η_s ≈ 15 m/s.  Post-fix: τ_strength = 10000 × |-150e6| = 1.5e12 Pa,
   // so V = 0.
   real_t V_abs, V1, V2, t1c, t2c;
   SolveLSW_TPV205(
      /*tau1_trial=*/0.0,    /*tau2_trial=*/70.0e6,
      /*tau1_total=*/0.0,    /*tau2_total=*/70.0e6,
      /*sigma_n_total=*/-150.0e6, eta_s,
      /*mu_eff=*/mu_s_barrier,
      V_abs, V1, V2, t1c, t2c);
   TEST_NEAR(V_abs, 0.0, 1e-12,
             "R-003: barrier QP under tensile σ_n must have V_abs = 0");
   TEST_NEAR(V1, 0.0, 1e-12,
             "R-003: barrier QP V1 must be 0 under tensile σ_n");
   TEST_NEAR(V2, 0.0, 1e-12,
             "R-003: barrier QP V2 must be 0 under tensile σ_n");

   // Same QP but with compressive σ_n_total = +120 MPa — the
   // standard barrier compression case.  Should also have V = 0.
   SolveLSW_TPV205(
      /*tau1_trial=*/0.0,    /*tau2_trial=*/70.0e6,
      /*tau1_total=*/0.0,    /*tau2_total=*/70.0e6,
      /*sigma_n_total=*/120.0e6, eta_s,
      /*mu_eff=*/mu_s_barrier,
      V_abs, V1, V2, t1c, t2c);
   TEST_NEAR(V_abs, 0.0, 1e-12,
             "R-003 sanity: barrier QP under compressive σ_n has V_abs=0");

   // INSIDE the rupture area (μ_eff ~ 0.677), tensile σ_n SHOULD open
   // the fault — V_abs = |τ| / η_s.  This sanity check makes sure the
   // R-003 fix didn't break the standard "fault opens under tension"
   // semantic in the rupture area.
   {
      const real_t mu_eff_in = TPV205Params::mu_s;  // 0.677
      SolveLSW_TPV205(
         /*tau1_trial=*/0.0,    /*tau2_trial=*/70.0e6,
         /*tau1_total=*/0.0,    /*tau2_total=*/70.0e6,
         /*sigma_n_total=*/-50.0e6, eta_s,
         /*mu_eff=*/mu_eff_in,
         V_abs, V1, V2, t1c, t2c);
      const real_t expected_V = 70.0e6 / eta_s;
      TEST_REL_NEAR(V_abs, expected_V, 1e-10,
                    "R-003 sanity: rupture-area QP opens under tensile σ_n");
   }
}

// ---------------------------------------------------------------------------
// Combined regression: under R-001 + R-002 + R-003 the original failure
// chain is "tensile transient sets V > 0 in barrier QP → δ accumulates
// past d_c → μ collapses to μ_d → barrier permanently broken".  Verify
// the post-fix path keeps V = 0 in the barrier even when arbitrary
// (δ, σ_n) pairs are imposed.
// ---------------------------------------------------------------------------
static void test_R002_R003_combined_BarrierStaysLockedThroughRupture()
{
   const real_t eta_s = TPV205Params::eta_s;
   const real_t mu_s_barrier = TPV205Params::mu_s_barrier;
   const real_t mu_d         = TPV205Params::mu_d;
   const real_t d_c          = TPV205Params::d_c;

   // Sweep δ from 0 to 2·d_c and σ_n from -150 MPa to +150 MPa.  At
   // every (δ, σ_n) pair the barrier QP must yield V_abs = 0.
   for (real_t delta_test = 0.0; delta_test <= 2.0 * d_c; delta_test += 0.1)
   {
      for (real_t sigma_n_test = -150.0e6; sigma_n_test <= 150.0e6;
           sigma_n_test += 50.0e6)
      {
         const real_t mu_eff = LSWFrictionCoefficient_TPV205(
                                  delta_test, mu_s_barrier, mu_d, d_c);
         real_t V_abs, V1, V2, t1c, t2c;
         SolveLSW_TPV205(
            /*tau1_trial=*/0.0,    /*tau2_trial=*/70.0e6,
            /*tau1_total=*/0.0,    /*tau2_total=*/70.0e6,
            sigma_n_test, eta_s,
            mu_eff,
            V_abs, V1, V2, t1c, t2c);
         char msg[256];
         std::snprintf(msg, sizeof(msg),
                       "R-002+R-003 combined: barrier V=0 at "
                       "δ=%.2f m, σ_n=%.0f MPa",
                       delta_test, sigma_n_test / 1e6);
         TEST_NEAR(V_abs, 0.0, 1e-12, msg);
      }
   }
}

// ---------------------------------------------------------------------------
// σ_n strength floor (sliver-blowup plan 2026-05-26 §Phase 2): the trailing
// `sigma_n_floor` argument of SolveLSW_TPV205 floors the σ_n that enters the
// SHEAR STRENGTH only.  The default 0.0 reproduces the historical "fault
// opens under tension" free-slip (`max(σ_n, 0)`) byte-exactly; a positive
// floor saturates the strength at μ_eff·floor so a tensile σ_n excursion can
// no longer collapse the strength to 0 and free-slide.  Numbers verbatim from
// the plan acceptance: σ_n=−2.9e9, μ_eff=0.3, η_s=4.6e6, |τ|=1.7e9.
// ---------------------------------------------------------------------------
static void test_lsw_strength_floor()
{
   const real_t sigma_n = -2.9e9;   // tensile (compression is +)
   const real_t mu_eff  = 0.3;      // < 0.5·mu_s_barrier ⇒ NOT a barrier QP
   const real_t eta_s   = 4.6e6;
   const real_t tau_mag = 1.7e9;    // |τ| carried entirely in the strike chan
   const real_t tau_abs = tau_mag;

   // --- Floor disabled (0.0): standard free slip, τ_strength = 0. ---
   real_t V_free, V1f, V2f, t1c_f, t2c_f;
   SolveLSW_TPV205(/*tau1_trial=*/0.0, /*tau2_trial=*/tau_mag,
                   /*tau1_total=*/0.0, /*tau2_total=*/tau_mag,
                   sigma_n, eta_s, mu_eff,
                   V_free, V1f, V2f, t1c_f, t2c_f,
                   /*sigma_n_floor=*/0.0);
   const real_t V_free_expected = tau_abs / eta_s;
   TEST_REL_NEAR(V_free, V_free_expected, 1e-12,
                 "floor=0: tensile sigma_n free-slides, V_abs = |tau|/eta_s");
   // τ_strength backed out of the radiation balance: τ_strength = |τ| − η_s·V.
   const real_t strength_free = tau_abs - eta_s * V_free;
   TEST_NEAR(strength_free, 0.0, 1e-3,
             "floor=0: tau_strength = 0 (free slip)");

   // --- Default-argument byte-exactness: omitting the arg == floor 0.0. ---
   real_t V_def, V1d, V2d, t1c_d, t2c_d;
   SolveLSW_TPV205(0.0, tau_mag, 0.0, tau_mag,
                   sigma_n, eta_s, mu_eff,
                   V_def, V1d, V2d, t1c_d, t2c_d);  // default sigma_n_floor=0
   TEST_NEAR(V_def, V_free, 0.0,
             "default sigma_n_floor arg is byte-exact with explicit 0.0");

   // --- Floor 10 MPa: strength saturates at μ_eff·floor = 3 MPa. ---
   const real_t floor = 10.0e6;
   real_t V_fl, V1l, V2l, t1c_l, t2c_l;
   SolveLSW_TPV205(0.0, tau_mag, 0.0, tau_mag,
                   sigma_n, eta_s, mu_eff,
                   V_fl, V1l, V2l, t1c_l, t2c_l,
                   floor);
   const real_t strength_floor   = mu_eff * floor;             // 3e6
   const real_t V_floor_expected = (tau_abs - strength_floor) / eta_s;
   TEST_REL_NEAR(V_fl, V_floor_expected, 1e-12,
                 "floor=10MPa: V_abs = (|tau| - mu_eff*floor)/eta_s (bounded)");
   // τ_strength backed out == μ_eff·σ_n_floor.
   const real_t strength_floored = tau_abs - eta_s * V_fl;
   TEST_REL_NEAR(strength_floored, strength_floor, 1e-9,
                 "floor=10MPa: tau_strength == mu_eff*sigma_n_floor (= 3 MPa)");
   // The floored slip rate is strictly below the free-slip rate.
   TEST_ASSERT(V_fl < V_free,
               "floor=10MPa: floored V_abs < free-slip V_abs");

   // --- Compressive σ_n is unaffected by a floor below it: max(σ_n,floor) ---
   // == σ_n when σ_n > floor, so floor=10MPa and floor=0 agree at σ_n=120MPa.
   real_t Vc0, Vc1, Vc2, tc1, tc2, Vc0b, Vc1b, Vc2b, tc1b, tc2b;
   SolveLSW_TPV205(0.0, tau_mag, 0.0, tau_mag, /*sigma_n=*/120.0e6, eta_s,
                   mu_eff, Vc0, Vc1, Vc2, tc1, tc2, /*floor=*/0.0);
   SolveLSW_TPV205(0.0, tau_mag, 0.0, tau_mag, /*sigma_n=*/120.0e6, eta_s,
                   mu_eff, Vc0b, Vc1b, Vc2b, tc1b, tc2b, /*floor=*/10.0e6);
   TEST_NEAR(Vc0, Vc0b, 0.0,
             "compressive sigma_n above floor: floor=10MPa == floor=0");
}

int main(int /*argc*/, char ** /*argv*/)
{
   std::cout << "Running TPV205 friction unit tests "
             << "(R-002, R-003 regression + sigma_n strength floor)...\n\n";

   std::cout << "[test_R002_BarrierDoesNotCollapseAtLargeSlip]\n";
   test_R002_BarrierDoesNotCollapseAtLargeSlip();

   std::cout << "\n[test_R003_BarrierLockedUnderTensileSigmaN]\n";
   test_R003_BarrierLockedUnderTensileSigmaN();

   std::cout << "\n[test_R002_R003_combined_BarrierStaysLockedThroughRupture]\n";
   test_R002_R003_combined_BarrierStaysLockedThroughRupture();

   std::cout << "\n[test_lsw_strength_floor]\n";
   test_lsw_strength_floor();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
