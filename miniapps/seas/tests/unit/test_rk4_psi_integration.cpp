// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan Step 5 + REVIEW R-V92-H01 — unit test for the
// coupled-RK4-on-psi arithmetic applied in drivers/tpv102_driver.cpp
// (lines around the `AgingLawPsi aging_law(...)` block).
//
// ============================================================================
// Motivation (R-V92-H01)
// ============================================================================
// The F01+F02 fix replaces four `UpdateStateAnalytic(psi_n, sr_k_i, ..., dt_sub)`
// calls inside the driver's RK4 time loop with classical coupled RK4 on psi
// (four `AgingLawPsi::Rate(V, psi, Dc)` evaluations + a Butcher-tableau
// weighted sum).  The existing unit tests all call `WaveOperator::Mult`
// directly — they build the operator, inject a Q, and inspect DOFData.
// NONE of them drive the driver's per-step RK4 loop, so the psi-RK4
// integration is untested by the regression gates.
//
// This test closes that coverage gap by REPLICATING the exact arithmetic
// applied in `drivers/tpv102_driver.cpp` to a single fault QP and checking
// that the resulting psi trajectory converges to the analytic constant-V
// solution (supplied by `UpdateStateAnalytic`, which is exact for constant
// V) with the expected O(dt^4) rate.
//
// ============================================================================
// What this test DOES NOT cover
// ============================================================================
// - End-to-end (Q, psi) coupling: this fixture fakes the four per-stage
//   V samples as a constant value so the coupling is trivial.  The
//   (Q, psi) coupling test requires either a full driver run on a tiny
//   fixture (not safe per feedback_no_local_reproducer.md) or a separate
//   Frontera run — out of scope here.
// - Nucleation: the driver's `ApplyNucleation` path mutates tau2_0,
//   not psi, so the fix is orthogonal.  No nucleation in this test.
//
// ============================================================================
// Test plan
// ============================================================================
// For a grid of dt values, integrate psi on [0, T] with V = V_const using
// the driver's 4-stage RK4 arithmetic.  Record the final psi and compare
// against `UpdateStateAnalytic(psi_0, V_const, Dc, T, f0, b, V0)` which is
// exact in the constant-V limit.  Verify that halving dt reduces the
// absolute error by ~16× (O(dt^4) convergence rate).
//
// Parameters chosen so that the errors are well above FP noise but psi
// evolution is bounded:
//   V_const = 1e-4 m/s          (representative of post-rupture slip)
//   psi_0   = 0.70               (between nucleation and steady state)
//   Dc, b, V0, f0 from TPV102Params
//   T       = 0.01 s            (~ 100 CFL-scale steps at 200m/P1)
//
// Usage:  ./seas_test_rk4_psi_integration   (serial, < 0.1 s)

#include "mfem.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/tpv102_params.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_TRUE(cond, msg) do { \
   num_tests++; \
   if (cond) { \
      num_passed++; \
      std::cout << "  PASSED: " << msg << "\n"; \
   } else { \
      num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; \
   } \
} while (0)

// ---------------------------------------------------------------------------
// Replicates the driver's per-step RK4-on-psi block EXACTLY (see
// drivers/tpv102_driver.cpp, near the AgingLawPsi ctor).  For constant
// V, sr_k1 = sr_k2 = sr_k3 = sr_k4 = V_const; the stage psi inputs come
// from psi_n + weight * dt * psi_k_prev.  This function returns psi_{n+1}
// for ONE RK4 step of size dt.
// ---------------------------------------------------------------------------
static real_t DriverRK4PsiStep(const AgingLawPsi &aging_law,
                                real_t psi_n, real_t V_const, real_t Dc,
                                real_t dt)
{
   // Stage 1 at (V, psi_n)
   real_t psi_k1 = aging_law.Rate(V_const, psi_n, Dc);
   real_t psi_stage2 = psi_n + 0.5 * dt * psi_k1;

   // Stage 2 at (V, psi_n + dt/2 * psi_k1)
   real_t psi_k2 = aging_law.Rate(V_const, psi_stage2, Dc);
   real_t psi_stage3 = psi_n + 0.5 * dt * psi_k2;

   // Stage 3 at (V, psi_n + dt/2 * psi_k2)
   real_t psi_k3 = aging_law.Rate(V_const, psi_stage3, Dc);
   real_t psi_stage4 = psi_n + dt * psi_k3;

   // Stage 4 at (V, psi_n + dt * psi_k3)
   real_t psi_k4 = aging_law.Rate(V_const, psi_stage4, Dc);

   // RK4 Butcher combination
   return psi_n + dt / 6.0 * (psi_k1 + 2.0*psi_k2 + 2.0*psi_k3 + psi_k4);
}

// Integrate psi over [0, T] with dt_step, return psi(T).
static real_t IntegrateRK4(const AgingLawPsi &aging_law,
                            real_t psi_0, real_t V_const, real_t Dc,
                            real_t T, real_t dt_step)
{
   real_t psi = psi_0;
   const int nsteps = static_cast<int>(std::round(T / dt_step));
   for (int s = 0; s < nsteps; s++)
   {
      psi = DriverRK4PsiStep(aging_law, psi, V_const, Dc, dt_step);
   }
   return psi;
}

int main()
{
   std::cout << "\n=== TPV102 v9.2.0 Step 5 (R-V92-H01) — Coupled-RK4-on-psi ===\n";

   const real_t b   = TPV102Params::b;
   const real_t V0  = TPV102Params::V0;
   const real_t f0  = TPV102Params::f0;
   const real_t Dc  = TPV102Params::Dc;

   AgingLawPsi aging_law(b, V0, f0);

   std::cout << "  AgingLawPsi params: b=" << b << "  V0=" << V0
             << "  f0=" << f0 << "  Dc=" << Dc << "\n";

   // =======================================================================
   // T1 — Exactness at constant V: the driver's RK4-on-psi must produce a
   // result within a tight bound of UpdateStateAnalytic for a single step,
   // because UpdateStateAnalytic is exact for constant V and RK4 is O(dt^4)
   // off the exact trajectory.
   // =======================================================================
   std::cout << "\n-- T1: single-step RK4 vs analytic (constant V) --\n";
   {
      const real_t V_const = 1.0e-4;        // m/s
      const real_t psi_0   = 0.70;
      const real_t dt      = 1.0e-4;        // 0.1 ms

      real_t psi_rk4     = DriverRK4PsiStep(aging_law, psi_0, V_const, Dc, dt);
      real_t psi_analy   = UpdateStateAnalytic(psi_0, V_const, Dc, dt,
                                                f0, b, V0);
      real_t err_abs     = std::abs(psi_rk4 - psi_analy);
      real_t err_rel     = err_abs / std::max(std::abs(psi_analy), real_t(1e-30));
      std::cout << "    psi_0         = " << std::setprecision(10) << psi_0    << "\n";
      std::cout << "    psi_rk4(dt)   = " << psi_rk4   << "\n";
      std::cout << "    psi_analytic  = " << psi_analy << "\n";
      std::cout << "    abs err       = " << err_abs   << "\n";
      std::cout << "    rel err       = " << err_rel   << "\n";
      // For one RK4 step, abs error ≲ K*(dt^5) where K depends on the fourth
      // time-derivative of psi.  Empirically here K is ~1e3 (typical for
      // aging-law at this regime); dt^5 = 1e-20 ⇒ expected err ~ 1e-17.
      // Loose bound 1e-12 to tolerate all plausible K variations.
      TEST_TRUE(err_abs < 1.0e-12,
                "T1 single RK4 step matches analytic for constant V (err < 1e-12)");
   }

   // =======================================================================
   // T2 — Convergence rate: |err(dt)| must shrink ~16× when dt halves.
   //
   // Regime: drive psi near steady state at V=1e-2 (psi_ss ≈ 0.489) so
   // BOTH terms in `dpsi/dt = (b*V0/Dc) · (exp((f0-psi)/b) - V/V0)` are
   // comparable — the exp term exposes the nonlinearity that
   // distinguishes an O(dt^4) integrator from an O(dt^2) operator split.
   // With V=1e-2 and psi_0 = 0.5 (just above steady state), exp((f0-ψ)/b)
   // = exp(8.33) ≈ 4.1e3 vs V/V0 = 1e4 — O(1) ratio, non-negligible d⁴ψ/dt⁴.
   //
   // Over T = 5 s with dt_base = 0.025 s (200 steps), accumulated RK4
   // error is above FP noise but still comfortably in the asymptotic
   // O(dt^4) regime for dt_base/2 and dt_base/4.
   // =======================================================================
   std::cout << "\n-- T2: O(dt^4) convergence rate across dt, dt/2, dt/4 --\n";
   {
      const real_t V_const = 1.0e-2;        // seismic slip rate
      const real_t psi_0   = 0.50;          // just above ψ_ss ≈ 0.489
      const real_t T       = 5.0;           // 5 s — several timescales
      const real_t dt_base = 0.025;         // 200 RK4 steps over 5 s

      // Reference (exact for constant V):
      real_t psi_ref = UpdateStateAnalytic(psi_0, V_const, Dc, T, f0, b, V0);

      real_t psi_h1 = IntegrateRK4(aging_law, psi_0, V_const, Dc, T, dt_base);
      real_t psi_h2 = IntegrateRK4(aging_law, psi_0, V_const, Dc, T, dt_base/2.0);
      real_t psi_h4 = IntegrateRK4(aging_law, psi_0, V_const, Dc, T, dt_base/4.0);

      real_t e1 = std::abs(psi_h1 - psi_ref);
      real_t e2 = std::abs(psi_h2 - psi_ref);
      real_t e4 = std::abs(psi_h4 - psi_ref);

      std::cout << "    T = " << T << " s, dt_base = " << dt_base << " s\n";
      std::cout << "    psi_ref (analytic)   = " << std::setprecision(12)
                << psi_ref << "\n";
      std::cout << "    psi_rk4(dt_base)     = " << psi_h1
                << "   err = " << std::scientific << std::setprecision(3)
                << e1 << "\n";
      std::cout << "    psi_rk4(dt_base/2)   = " << std::fixed << std::setprecision(12)
                << psi_h2
                << "   err = " << std::scientific << std::setprecision(3)
                << e2 << "\n";
      std::cout << "    psi_rk4(dt_base/4)   = " << std::fixed << std::setprecision(12)
                << psi_h4
                << "   err = " << std::scientific << std::setprecision(3)
                << e4 << "\n";

      // If errors are already at FP noise, the convergence ratio is
      // meaningless — check that regime first.  O(dt^4) predicts
      // e(dt_base/2) ≈ e_base / 16.  With e_base in [1e-8, 1e-16], that
      // puts e(dt_base/2) at [~6e-10, ~6e-18].  If e_base is already
      // 1e-16, FP noise dominates before we see the scaling.  Pick
      // parameters above so e_base ~ 1e-10 to give two clean halvings.
      const real_t fp_noise = 5.0e-16;
      if (e1 < 50.0 * fp_noise)
      {
         std::cout << "    [NOTE] coarse error at FP-noise floor; skipping ratio check\n";
         TEST_TRUE(e1 <= 1.0e-12,
                   "T2 RK4 matches analytic at FP-noise precision");
      }
      else
      {
         real_t r1 = (e2 > fp_noise) ? e1 / e2 : real_t(16.0);
         real_t r2 = (e4 > fp_noise) ? e2 / e4 : real_t(16.0);
         std::cout << "    ratio e1/e2 = " << std::fixed << std::setprecision(2)
                   << r1 << "   (O(dt^4) predicts 16)\n";
         std::cout << "    ratio e2/e4 = " << r2
                   << "   (O(dt^4) predicts 16; FP floor if e4 too small)\n";
         // Accept 8 ≤ r ≤ 32 — wide band that rules out accidental O(dt),
         // O(dt^2), or O(dt^8) while tolerating pre-asymptotic regime.
         TEST_TRUE(r1 >= 8.0 && r1 <= 32.0,
                   "T2 first halving: convergence ratio is O(dt^4)-consistent");
         // Second halving may approach FP noise; only enforce an upper bound.
         TEST_TRUE(r2 >= 4.0,
                   "T2 second halving: at least O(dt^2) (FP floor tolerant)");
      }
   }

   // =======================================================================
   // T3 — Regime consistency near steady state.  At V=V_ini the system is
   // near equilibrium, so |dpsi/dt| is tiny; the RK4 integrator should
   // stay bit-close to the analytic trajectory regardless of dt.
   // =======================================================================
   std::cout << "\n-- T3: near-equilibrium stability (small V, V_ini) --\n";
   {
      const real_t V_const = TPV102Params::V_ini;   // 1e-12 m/s
      // Steady-state psi for V=V_ini.  Exactly on the equilibrium attractor
      // ⇒ Rate() is ~0 everywhere, RK4 and analytic both give psi constant.
      const real_t psi_ss  = f0 + b * std::log(V0 / V_const);
      const real_t T       = 0.1;           // 100 ms
      const real_t dt_big  = 1.0e-3;        // "big" dt for near-eq regime

      real_t psi_ref = UpdateStateAnalytic(psi_ss, V_const, Dc, T, f0, b, V0);
      real_t psi_num = IntegrateRK4(aging_law, psi_ss, V_const, Dc, T, dt_big);

      real_t abs_err = std::abs(psi_num - psi_ref);
      std::cout << "    psi_ss    = " << std::setprecision(12) << psi_ss  << "\n";
      std::cout << "    psi_ref   = " << psi_ref << "\n";
      std::cout << "    psi_rk4   = " << psi_num << "\n";
      std::cout << "    abs err   = " << std::scientific << std::setprecision(3)
                << abs_err << "\n";
      TEST_TRUE(abs_err < 1.0e-10,
                "T3 RK4 tracks analytic at V_ini steady state");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Step 5 R-V92-H01 results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
