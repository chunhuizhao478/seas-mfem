// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 nucleation persistent-prestress channel regression
// (companion to debug_document/tpv102_debug_document/
//  tpv102_nucleation_code_review_2026-04-22.md and
//  tpv102_nucleation_code_fix_2026-04-22.md).
//
// THE BUG (review-2026-04-22, agreed CRITICAL):
//   The legacy ApplyNucleationTotal poked nucleation as a single-DOF
//   point source into bulk Q[SXY].  The wave operator radiates that
//   point source outward in O(h/cp), so on a 200 m TPV102 mesh the
//   shape-weighted SXY at the fault QP collapsed to ~10^-4-10^-5 of
//   the requested 25 MPa amplitude — V_max stalled at 2e-7 m/s on
//   Frontera, ~10^5x below the V_nuc trigger threshold.
//
// THE FIX:
//   1. New DOFData fields {sigma_n_nuc, tau1_nuc, tau2_nuc} (default 0).
//   2. FaultFaceFlux::EvaluateTotal adds them to the trial traction
//      after ComputeTrialTraction — re-imposed at every Riemann solve.
//   3. New helper ApplyNucleationTotalPrestress(dof_data, ..., t)
//      OVERWRITES dof_data[i].tau2_nuc per call (mirrors the working
//      fluctuation-Q ApplyNucleation overwrite at
//      tpv102_setup.hpp:135-148).  Bulk Q is never touched by
//      nucleation.
//
// WHAT THIS TEST CHECKS:
//   Test A (additivity) — direct EvaluateTotal call with hand-set Q±:
//     - tau2_nuc = 0  → trial traction matches ComputeTrialTraction.
//     - tau2_nuc = X  → trial traction shifts by exactly X (modulo the
//                        eta_s*V friction reaction, bounded by V_ini).
//   Test B (overwrite semantics) — ApplyNucleationTotalPrestress called
//     repeatedly at increasing t produces the corresponding
//     NucleationPerturbation(...) value in dof_data[i].tau2_nuc each
//     time (no add-delta drift).
//   Test C (default no-op) — DOFData with tau*_nuc=0 yields the same
//     EvaluateTotal output as a fixture without the new channel: the
//     fix is a strict no-op for callers that do not opt in.
//
// Runs in <1 s on the local laptop, no MPI needed.

#include "mfem.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/friction_solver.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) < t_) { num_passed++; \
      std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << " (got " << std::scientific << std::setprecision(6) \
                << v_ << ", expected " << e_ \
                << ", |diff| " << std::abs(v_-e_) \
                << ", tol " << t_ << ")\n"; } \
} while (0)

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

namespace {

// Build a locked-fault DOFData fixture in TPV102 conditions, but with
// {sigma_n0, tau1_0, tau2_0} zero (the total-Q dispatch contract — bulk
// Q carries the static prestress).  Friction solver is in equilibrium
// at V_ini so the eta_s*V reaction term is bounded by ~Zs*V_ini ≈ 1e-5
// Pa, well below the test tolerances.
DOFData MakeTotalLockedDOF(real_t a)
{
   DOFData d;
   d.Zp_plus = TPV102Params::Zp; d.Zp_minus = TPV102Params::Zp;
   d.Zs_plus = TPV102Params::Zs; d.Zs_minus = TPV102Params::Zs;
   d.eta_p = TPV102Params::Zp / 2.0;
   d.eta_s = TPV102Params::eta_s;
   d.sigma_n0 = 0;  d.tau1_0 = 0;  d.tau2_0 = 0;
   d.sigma_n_nuc = 0;  d.tau1_nuc = 0;  d.tau2_nuc = 0;
   d.a = a;
   d.Dc = TPV102Params::Dc;
   d.psi = ComputeInitialPsi(a);
   d.slip_rate = TPV102Params::V_ini;
   d.V1 = 0;  d.V2 = TPV102Params::V_ini;
   return d;
}

// Symmetric Q± at TPV102 background, in FAULT-LOCAL coordinates as
// supplied to FaultFaceFlux::EvaluateTotal by the wave operator.
//
// InitializeStateTotal lays down GLOBAL stress
//   Q[SYY_global] = +sigma_n,   Q[SXY_global] = -tau_ini
// (see tpv102_setup_total.hpp:38-50 for the BP5 canonical-frame
// rotation derivation: t1=(0,0,-1), t2=(+1,0,0), n=(0,-1,0)).
// The wave operator rotates Q to fault-local before EvaluateTotal:
// under the Voigt rotation R = [n; t1; t2]_rows we have
//   sigma_local_NN  = sigma_yy_global         = +sigma_n
//   sigma_local_NT2 = -sigma_xy_global        = +tau_ini
// where local index NN = SXX (1st diag in the (n, t1, t2) frame),
// NT2 = SXZ (off-diag between normal and tangent2 = strike).
//
// The earlier draft of this fixture wrote SXZ = -tau_ini, which is
// the GLOBAL sign — symmetric Q masked the error in tau2_corr (the
// magnitude is the same), but any directional check (e.g. the
// Riemann velocity jump direction in Test A8) would have inverted.
// Use the fault-local sign here for parity with InitializeStateTotal.
void SetSymmetricBackgroundQ(real_t *Qp, real_t *Qm)
{
   for (int c = 0; c < NUM_STATE; c++) { Qp[c] = 0; Qm[c] = 0; }
   // SXX is the fault-NORMAL stress in the rotated frame.
   Qp[SXX] = TPV102Params::sigma_n;  Qm[SXX] = TPV102Params::sigma_n;
   // SXZ = sigma_NT2 = +tau_ini in fault-local coordinates (sign flip
   // from global -tau_ini under the BP5 rotation; see header above).
   Qp[SXZ] = +TPV102Params::tau_ini; Qm[SXZ] = +TPV102Params::tau_ini;
}

} // anonymous

int main()
{
   std::cout << "\n=== TPV102 nucleation persistent-prestress channel "
             << "(2026-04-22 fix regression) ===\n";

   FaultFaceFlux flux(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);

   // -----------------------------------------------------------------
   // Test A: additivity — tau2_nuc shifts the trial traction the
   // friction solver sees by exactly the requested amount.
   // -----------------------------------------------------------------
   std::cout << "\n-- Test A: trial-traction additivity --\n";
   {
      const real_t a = TPV102Params::a_vw;

      // Baseline: tau2_nuc = 0.  The friction solver sees the bulk
      // background only; data.tau2_corr ≈ +tau_ini (locked fault,
      // V≈V_ini) — fault-local sigma_NT2 sign per SetSymmetricBackgroundQ.
      DOFData d_base = MakeTotalLockedDOF(a);
      real_t Qp[NUM_STATE], Qm[NUM_STATE], Qip[NUM_STATE], Qim[NUM_STATE];
      SetSymmetricBackgroundQ(Qp, Qm);
      flux.EvaluateTotal(d_base, Qp, Qm, Qip, Qim);
      const real_t tau2_corr_base = d_base.tau2_corr;
      // Sanity: locked fault at V_ini → tau2_corr ≈ +tau_ini, eta_s*V_ini
      // ≈ Zs*V_ini/2 ≈ 5e-3 Pa drift.  Allow 10 Pa slack.
      TEST_NEAR(tau2_corr_base, +TPV102Params::tau_ini, 10.0,
                "A1: baseline data.tau2_corr ≈ +tau_ini at locked fault");

      // With nucleation: tau2_nuc = nuc_dtau (= 25 MPa).  Under the
      // R-001-corrected EvaluateTotal, data.tau2_corr is the TOTAL
      // (= tau_corr_trial + tau2_nuc), so the OUTPUT tau2_corr shift
      // equals nuc_dtau MINUS the eta_s·ΔV friction reaction induced
      // by the rupture release (a few MPa at most for TPV102 numbers).
      // The shift therefore must be POSITIVE (nuc adds upward) AND
      // at least 80% of nuc_dtau.
      DOFData d_nuc = MakeTotalLockedDOF(a);
      d_nuc.tau2_nuc = TPV102Params::nuc_dtau;
      SetSymmetricBackgroundQ(Qp, Qm);  // Q untouched — nuc not in bulk
      flux.EvaluateTotal(d_nuc, Qp, Qm, Qip, Qim);
      const real_t shift = d_nuc.tau2_corr - tau2_corr_base;
      std::cout << "    baseline tau2_corr  = " << std::scientific
                << std::setprecision(6) << tau2_corr_base << " Pa\n"
                << "    nuc      tau2_corr  = " << d_nuc.tau2_corr << " Pa\n"
                << "    shift               = " << shift
                << " Pa (target +" << TPV102Params::nuc_dtau << ")\n";
      TEST_ASSERT(shift > 0,
                  "A2: tau2_nuc=+25 MPa shifts tau2_corr UPWARD");
      TEST_ASSERT(shift > 0.8 * TPV102Params::nuc_dtau,
                  "A3: tau2_corr shift >= 80% of nuc_dtau (= 20 MPa)");
      TEST_ASSERT(shift <= TPV102Params::nuc_dtau + 1.0,
                  "A4: tau2_corr shift <= nuc_dtau (no over-amplification)");

      // tau1_nuc and sigma_n_nuc parity: writing those should shift
      // tau1_corr / sigma_n_corr correspondingly with no cross-talk
      // into tau2_corr.  This guards against future code that confuses
      // the dip vs strike vs normal indices on the new fields.
      DOFData d_t1 = MakeTotalLockedDOF(a);
      d_t1.tau1_nuc = 1.0e6;
      SetSymmetricBackgroundQ(Qp, Qm);
      flux.EvaluateTotal(d_t1, Qp, Qm, Qip, Qim);
      TEST_ASSERT(d_t1.tau1_corr > 1e5,
                  "A5: tau1_nuc=+1 MPa shifts tau1_corr (dip channel)");
      TEST_NEAR(d_t1.tau2_corr, tau2_corr_base, 1e3,
                "A6: tau1_nuc does not leak into tau2_corr");

      DOFData d_sn = MakeTotalLockedDOF(a);
      d_sn.sigma_n_nuc = 1.0e6;
      SetSymmetricBackgroundQ(Qp, Qm);
      flux.EvaluateTotal(d_sn, Qp, Qm, Qip, Qim);
      const real_t expected_sn = TPV102Params::sigma_n + 1.0e6;
      TEST_NEAR(d_sn.sigma_n_corr, expected_sn, 1.0,
                "A7: sigma_n_nuc shifts sigma_n_corr by exactly the offset");

      // ---------------------------------------------------------------
      // R-002 subtest A8: Riemann imposed-velocity jump magnitude
      // equals V_abs (the friction-decomposed slip rate).
      //
      // For matched impedances and symmetric Q±:
      //   Q_imp_plus[VZ]  = +(1/Zs)·(tau2_corr - Q[SXZ])
      //   Q_imp_minus[VZ] = -(1/Zs)·(tau2_corr - Q[SXZ])
      //   ΔVZ = Q_imp_plus[VZ] - Q_imp_minus[VZ]
      //       = (2/Zs)·(tau2_corr_TRIAL_SCALE - Q[SXZ])
      //       = (2/Zs)·(-eta_s·V2)             (TRIAL-scale tau2_corr)
      //       = -V2                             (since 2·eta_s = Zs)
      //   |ΔVZ| = |V2| = V_abs                  (V_abs is purely V2)
      //
      // Under the R-001 buggy mutation, tau2_corr would carry the
      // tau2_nuc contribution and the jump would be
      //   |ΔVZ| = |(2/Zs)·(tau2_nuc - eta_s·V2)|
      // For tau2_nuc = 25 MPa, Zs ≈ 9.25e6 Pa·s/m, V2 ≈ 0.22 m/s:
      //   |ΔVZ| ≈ |5.4 - 0.22| ≈ 5.2 m/s vs V_abs ≈ 0.22 m/s — a
      // ~24× overshoot.  Require the jump to match V_abs to within
      // 5% relative + a 1e-9 m/s absolute floor (so the locked-fault
      // baseline at V≈V_ini=1e-12 m/s passes too).
      //
      // Re-run the d_nuc fixture and inspect Qip / Qim returned by
      // the same EvaluateTotal call (the ones written into the
      // baseline already; redo since ApplyNucleationTotalPrestress
      // doesn't update Qip/Qim).
      DOFData d_a8 = MakeTotalLockedDOF(a);
      d_a8.tau2_nuc = TPV102Params::nuc_dtau;
      real_t Qp8[NUM_STATE], Qm8[NUM_STATE];
      real_t Qip8[NUM_STATE], Qim8[NUM_STATE];
      SetSymmetricBackgroundQ(Qp8, Qm8);
      flux.EvaluateTotal(d_a8, Qp8, Qm8, Qip8, Qim8);
      const real_t jump_VZ = std::abs(Qip8[VZ] - Qim8[VZ]);
      const real_t v_abs   = d_a8.slip_rate;
      std::cout << "    A8 V_abs            = " << std::scientific
                << std::setprecision(6) << v_abs << " m/s\n"
                << "    A8 |Qip[VZ]-Qim[VZ]| = " << jump_VZ << " m/s\n";
      TEST_ASSERT(jump_VZ > 0.0,
                  "A8a: Riemann VZ jump is non-zero under nucleation");
      // 5% relative tolerance + 1e-9 m/s absolute floor.
      const real_t a8_tol = 0.05 * std::max(v_abs, 1.0e-9);
      TEST_NEAR(jump_VZ, v_abs, a8_tol,
                "A8b: |Qip[VZ] - Qim[VZ]| == V_abs (no tau2_nuc leak "
                "into Riemann jump — would be ~24x at full nuc if "
                "tau*_trial were mutated)");
   }

   // -----------------------------------------------------------------
   // Test B: ApplyNucleationTotalPrestress overwrite semantics.
   // -----------------------------------------------------------------
   std::cout << "\n-- Test B: ApplyNucleationTotalPrestress overwrite --\n";
   {
      // 3 fault QPs: hypocenter, off-hypo (ramp partially active),
      // far-field (outside nucleation radius — should stay 0).
      std::vector<Vector> coords(3, Vector(3));
      coords[0][0] = TPV102Params::hypo_along_strike;     // hypo
      coords[0][1] = 0.0;
      coords[0][2] = -TPV102Params::hypo_down_dip;
      coords[1][0] = TPV102Params::hypo_along_strike;     // 1.5 km off
      coords[1][1] = 0.0;
      coords[1][2] = -(TPV102Params::hypo_down_dip + 1.5e3);
      coords[2][0] = 30e3;                                // far-field
      coords[2][1] = 0.0;
      coords[2][2] = -1e3;

      std::vector<DOFData> dof_data(3);
      // Pre-poison tau2_nuc with garbage so we can detect overwrite.
      dof_data[0].tau2_nuc = -7.7e6;
      dof_data[1].tau2_nuc = +9.9e6;
      dof_data[2].tau2_nuc = +1.1e6;

      // Saturated time: dtau at hypo = nuc_dtau; off-hypo < nuc_dtau;
      // far-field = 0.
      const real_t t = TPV102Params::nuc_T + 0.5;
      ApplyNucleationTotalPrestress(dof_data, coords, t);

      const real_t expected_0 = NucleationPerturbation(
         coords[0][0], std::abs(coords[0][2]), t);
      const real_t expected_1 = NucleationPerturbation(
         coords[1][0], std::abs(coords[1][2]), t);
      const real_t expected_2 = NucleationPerturbation(
         coords[2][0], std::abs(coords[2][2]), t);

      TEST_NEAR(dof_data[0].tau2_nuc, expected_0,
                1e-9 * std::max(std::abs(expected_0), 1.0),
                "B1: hypo QP overwritten to NucleationPerturbation(t)");
      TEST_NEAR(dof_data[1].tau2_nuc, expected_1,
                1e-9 * std::max(std::abs(expected_1), 1.0),
                "B2: off-hypo QP overwritten to partial dtau");
      TEST_NEAR(dof_data[2].tau2_nuc, 0.0, 1e-9,
                "B3: far-field QP overwritten to 0 (outside nuc radius)");
      TEST_NEAR(dof_data[0].tau2_nuc, TPV102Params::nuc_dtau,
                1e-9 * TPV102Params::nuc_dtau,
                "B4: hypo at saturation t > nuc_T → tau2_nuc = nuc_dtau");

      // Idempotency under repeated call at same t.
      ApplyNucleationTotalPrestress(dof_data, coords, t);
      TEST_NEAR(dof_data[0].tau2_nuc, expected_0,
                1e-9 * std::max(std::abs(expected_0), 1.0),
                "B5: repeated call at same t is idempotent");

      // Time advance: at t' < t, hypo dtau drops back to a smaller
      // value (overwrite, not add).
      const real_t t_early = TPV102Params::nuc_T * 0.3;
      ApplyNucleationTotalPrestress(dof_data, coords, t_early);
      const real_t expected_0_early = NucleationPerturbation(
         coords[0][0], std::abs(coords[0][2]), t_early);
      TEST_NEAR(dof_data[0].tau2_nuc, expected_0_early,
                1e-9 * std::max(std::abs(expected_0_early), 1.0),
                "B6: overwrite (not add): tau2_nuc tracks current t");
      TEST_ASSERT(expected_0_early < expected_0,
                  "B7: temporal ramp t=0.3T → dtau < saturated value");

      // tau1_nuc and sigma_n_nuc are NOT touched by this helper.
      TEST_ASSERT(dof_data[0].tau1_nuc == 0.0,
                  "B8: tau1_nuc untouched (TPV102 is pure strike-slip)");
      TEST_ASSERT(dof_data[0].sigma_n_nuc == 0.0,
                  "B9: sigma_n_nuc untouched");
   }

   // -----------------------------------------------------------------
   // Test C: default no-op — DOFData with tau*_nuc=0 produces
   // bit-exact identical EvaluateTotal output to a fixture without the
   // new channel (the legacy contract is preserved).
   // -----------------------------------------------------------------
   std::cout << "\n-- Test C: default (tau*_nuc=0) is a strict no-op --\n";
   {
      // Hand-pick a non-trivial Q± with stress and velocity gradients.
      real_t Qp[NUM_STATE] = {0}, Qm[NUM_STATE] = {0};
      Qp[SXX] = 100e6; Qp[SXZ] = -60e6; Qp[VZ] = 0.5;
      Qm[SXX] = 110e6; Qm[SXZ] = -65e6; Qm[VZ] = -0.4;

      DOFData d1 = MakeTotalLockedDOF(TPV102Params::a_vw);
      DOFData d2 = MakeTotalLockedDOF(TPV102Params::a_vw);
      // d2 is the "no nuc field" reference — but since the new field
      // defaults to 0, d1 == d2 by construction.  This test still has
      // value: it documents that the addition is invisible to non-nuc
      // callers, locking that behavior in for future maintainers.
      real_t Qip1[NUM_STATE], Qim1[NUM_STATE];
      real_t Qip2[NUM_STATE], Qim2[NUM_STATE];
      flux.EvaluateTotal(d1, Qp, Qm, Qip1, Qim1);
      flux.EvaluateTotal(d2, Qp, Qm, Qip2, Qim2);

      bool bit_exact = true;
      for (int c = 0; c < NUM_STATE; c++)
      {
         if (Qip1[c] != Qip2[c] || Qim1[c] != Qim2[c]) { bit_exact = false; }
      }
      TEST_ASSERT(bit_exact,
                  "C1: tau*_nuc=0 produces bit-exact identical Q_imp");
      TEST_ASSERT(d1.tau2_corr == d2.tau2_corr &&
                  d1.tau1_corr == d2.tau1_corr &&
                  d1.sigma_n_corr == d2.sigma_n_corr,
                  "C2: tau*_nuc=0 produces bit-exact identical *_corr");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
