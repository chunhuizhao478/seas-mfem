// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "fault_face_flux.hpp"
#include "seas_diag_rank.hpp"
#include "tpv205_friction.hpp"
#include "../spatial/code/spatial_friction.hpp"  // Phase H.6: LSWFrictionCoefficient_ForcedRupture
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace mfem
{
namespace seas
{

#ifdef SEAS_TEST_INTERNAL
// TEST-ONLY cross-rank seed knob [Pa] (see header).  Zero => no-op / byte-exact.
real_t FaultFaceFlux::s_seas_test_tau2_trial_perturb_pa = 0.0;
// TEST-ONLY: disable the Phase-2 reconcile (negative-test leg).  False => runs.
bool FaultFaceFlux::s_seas_test_disable_reconcile = false;
#endif

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
FaultFaceFlux::FaultFaceFlux(real_t rho, real_t cp, real_t cs)
   : rho_(rho), cp_(cp), cs_(cs)
{
   Zp_ = rho * cp;
   Zs_ = rho * cs;
}

// ---------------------------------------------------------------------------
// Trial traction: Eq. (7a-c)
// In fault-local coordinates: x = normal, y = tangent1, z = tangent2.
// So SXX = sigma_nn, SXY = tau_1 (nt1), SXZ = tau_2 (nt2).
// VX = v_n, VY = v_t1, VZ = v_t2.
// ---------------------------------------------------------------------------
void FaultFaceFlux::ComputeTrialTraction(const DOFData &data,
                                         const real_t *Q_plus,
                                         const real_t *Q_minus,
                                         real_t &sigma_n_trial,
                                         real_t &tau1_trial,
                                         real_t &tau2_trial)
{
   real_t invZp_plus = 1.0 / data.Zp_plus;
   real_t invZp_minus = 1.0 / data.Zp_minus;
   real_t invZs_plus = 1.0 / data.Zs_plus;
   real_t invZs_minus = 1.0 / data.Zs_minus;

   // Eq. (7a): sigma_n^trial = eta_p * (v_n^- - v_n^+ + sigma_n^+/Zp^+ + sigma_n^-/Zp^-)
   sigma_n_trial = data.eta_p * (Q_minus[VX] - Q_plus[VX]
                                 + Q_plus[SXX] * invZp_plus
                                 + Q_minus[SXX] * invZp_minus);

   // Eq. (7b): tau_1^trial = eta_s * (v_t1^- - v_t1^+ + tau_1^+/Zs^+ + tau_1^-/Zs^-)
   tau1_trial = data.eta_s * (Q_minus[VY] - Q_plus[VY]
                              + Q_plus[SXY] * invZs_plus
                              + Q_minus[SXY] * invZs_minus);

   // Eq. (7c): tau_2^trial = eta_s * (v_t2^- - v_t2^+ + tau_2^+/Zs^+ + tau_2^-/Zs^-)
   tau2_trial = data.eta_s * (Q_minus[VZ] - Q_plus[VZ]
                              + Q_plus[SXZ] * invZs_plus
                              + Q_minus[SXZ] * invZs_minus);

#ifdef SEAS_TEST_INTERNAL
   // TEST-ONLY cross-rank seed (Phase 1 of the reconcile plan): when one rank
   // sets the knob, add a small absolute [Pa] nudge to its strike-channel
   // trial traction — the deterministic stand-in for the ~1e-14 (relative)
   // shared-face interpolation gap that tips two ranks onto opposite slip/lock
   // branches at the LSW kink.  Zero default => byte-exact; whole block
   // compiles out in production.
   tau2_trial += s_seas_test_tau2_trial_perturb_pa;
#endif
}

// ---------------------------------------------------------------------------
// Round-12 Patch 1: stage helpers.  `Evaluate` below composes
// ComputeStageState + BuildImposedState + WriteBackState; the split
// exists so wave_operator.inl can face-average one stage field without
// duplicating any physics.  The arithmetic ordering inside
// ComputeStageState → CompleteFromTrial → CompleteFromTheta →
// CompleteFromVabs is identical to the pre-refactor `Evaluate`, so the
// composition is byte-identical on baseline inputs.
// ---------------------------------------------------------------------------
void FaultFaceFlux::ComputeStageState(const DOFData &data,
                                      const real_t *Q_plus,
                                      const real_t *Q_minus,
                                      EvalStageState &s,
                                      FrictionSolver::Method method) const
{
   // Step 1: Trial traction (Eq. 7)
   ComputeTrialTraction(data, Q_plus, Q_minus,
                        s.sigma_n_trial, s.tau1_trial, s.tau2_trial);

   // Diagnostic (TPV104 normal-traction freeze).  When env var
   // SEAS_TPV104_FREEZE_SIGMA_N is set (and not "0"), pin the
   // resulting `sigma_n_total = data.sigma_n0 + data.sigma_n_nuc +
   // s.sigma_n_trial` to the override value by setting
   //     s.sigma_n_trial = override - data.sigma_n0 - data.sigma_n_nuc
   // which works in BOTH modes:
   //   - fluctuation-Q (TPV104, data.sigma_n0 = 120 MPa, nuc = 0):
   //       s.sigma_n_trial = 0 → sigma_n_total = 120 MPa, Q_imp[SXX] = 0
   //   - total-Q (TPV102 v9.3.0, data.sigma_n0 = 0, nuc = 0):
   //       s.sigma_n_trial = 120e6 → sigma_n_total = 120 MPa,
   //       Q_imp[SXX] = 120e6 (full physical stress)
   //   SEAS_TPV104_FREEZE_SIGMA_N=120e6  -> sigma_n_total = 120 MPa
   //   SEAS_TPV104_FREEZE_SIGMA_N=1      -> 120 MPa (default)
   //   unset / =0                         -> baseline (no override)
   // Sign: positive = compression (TPV104Params::sigma_n convention).
   //
   // Earlier revision pinned s.sigma_n_trial directly to the override,
   // which silently doubled sigma_n_total in fluctuation mode (TPV104
   // job 7677661 produced sigma_n_corr = 240 MPa).
   {
      const char *freeze = std::getenv("SEAS_TPV104_FREEZE_SIGMA_N");
      if (freeze && freeze[0] != '\0' &&
          !(freeze[0] == '0' && freeze[1] == '\0'))
      {
         char *endp = nullptr;
         const real_t parsed = std::strtod(freeze, &endp);
         const real_t target =
            (endp != freeze && parsed > 0.0) ? parsed : 120.0e6;
         s.sigma_n_trial = target - data.sigma_n0 - data.sigma_n_nuc;
      }
   }

   CompleteFromTrial(data, s, method);

#ifdef SEAS_DIAG_FAULT_FLUX
   // R-501 probe: full per-call diagnostic at the hypocenter QP.
   // Gated by SEAS_DIAG_FAULT_FLUX (build) AND SEAS_DIAG_V1_DRIFT
   // (runtime) AND data.diag_print (single DOF).  Prints AFTER
   // CompleteFromTrial so all derived quantities (V_abs, V1, V2,
   // tau1_corr, tau2_corr) are populated.  One line per Evaluate
   // call on the hypo QP — no MPI, single rank.  Sequence number
   // `seq` is a process-local counter that lets post-processing
   // sort calls in temporal order even when ADER substeps interleave.
   if (data.diag_print)
   {
      static const bool probe_enabled = []() {
         const char *e = std::getenv("SEAS_DIAG_V1_DRIFT");
         return e && e[0] != '\0' &&
                !(e[0] == '0' && e[1] == '\0');
      }();
      if (probe_enabled)
      {
         static long long seq = 0;
         ++seq;
         std::fprintf(stderr,
            "[R-501] seq=%lld  tau1_trial=%+.17e  tau2_trial=%+.17e  "
            "V_abs=%+.17e  Theta=%+.17e  V1=%+.17e  V2=%+.17e  "
            "tau1_corr=%+.17e  tau2_corr=%+.17e  "
            "Qp_VY=%+.17e  Qm_VY=%+.17e  Qp_SXY=%+.17e  Qm_SXY=%+.17e\n",
            seq, s.tau1_trial, s.tau2_trial,
            s.V_abs, s.Theta, s.V1, s.V2,
            s.tau1_corr, s.tau2_corr,
            Q_plus[VY], Q_minus[VY], Q_plus[SXY], Q_minus[SXY]);
      }
   }
#endif
}

void FaultFaceFlux::CompleteFromTrial(const DOFData &data,
                                      EvalStageState &s,
                                      FrictionSolver::Method method) const
{
   // v9.4.0 Commit 1: Total traction = static pre-stress
   // (`tau*_0` / `sigma_n0`) + persistent nucleation channel
   // (`tau*_nuc` / `sigma_n_nuc`) + dynamic trial from bulk Q
   // (`tau*_trial`).
   s.sigma_n_total = data.sigma_n0 + data.sigma_n_nuc + s.sigma_n_trial;
   s.tau1_total    = data.tau1_0   + data.tau1_nuc    + s.tau1_trial;
   s.tau2_total    = data.tau2_0   + data.tau2_nuc    + s.tau2_trial;

   // Traction magnitude Θ = sqrt(tau1_total² + tau2_total²)
   s.Theta = std::sqrt(s.tau1_total * s.tau1_total
                       + s.tau2_total * s.tau2_total);

   CompleteFromTheta(data, s, method);
}

void FaultFaceFlux::CompleteFromTheta(const DOFData &data,
                                      EvalStageState &s,
                                      FrictionSolver::Method method) const
{
   // Step 2: Solve friction equation for |V̂| (Eq. 8)
   s.V_abs = 0.0;
   if (s.Theta > 0.0)
   {
      s.V_abs = solver_.Solve(s.Theta, data.psi, std::abs(s.sigma_n_total),
                              data.eta_s, data.a, method);
   }
   CompleteFromVabs(data, s);
}

void FaultFaceFlux::CompleteFromVabs(const DOFData &data,
                                     EvalStageState &s) const
{
   // Step 3: Decompose slip rate into components (Eq. 9).  Defaults
   // mirror the pre-refactor code: V1 = V2 = 0 and tau*_corr =
   // tau*_trial when the friction solver returns V_abs = 0 (or
   // Theta = 0 at entry).
   s.V1 = 0.0; s.V2 = 0.0;
   s.tau1_corr = s.tau1_trial;
   s.tau2_corr = s.tau2_trial;
   // Normal traction is unchanged by slip (Eq. 10), regardless of
   // whether the friction solver ran.
   s.sigma_n_corr = s.sigma_n_trial;

   if (s.Theta > 0.0 && s.V_abs > 0.0)
   {
      // Friction strength.  σ_n strength floor (sliver-blowup plan
      // 2026-05-26): when enabled (floor >= 0) the strength's σ_n is
      // max(σ_n, floor) so a tensile excursion saturates at the floor
      // instead of feeding back; disabled (sentinel < 0) keeps the
      // historical |σ_n| ⇒ byte-exact for the TPV102/104 regressions.
      // Only the explicit strength is floored; the friction-solver σ_n
      // argument (CompleteFromTheta:201) is intentionally NOT floored in
      // v1 (plan §Phase 3 decision 2).
      real_t C = std::exp(data.psi / data.a) / (2.0 * FrictionSolver::V0);
      real_t f_V = data.a * std::asinh(s.V_abs * C);
      const real_t sn_str = (sigma_n_strength_floor_ >= 0.0)
                            ? std::max(s.sigma_n_total, sigma_n_strength_floor_)
                            : std::abs(s.sigma_n_total);
      real_t strength = sn_str * f_V;

      // H2 experiment (env-var SEAS_FORCE_V1_ZERO).  For pure-strike-slip
      // configurations (TPV104: tau1_0 = 0, tau1_nuc = 0, ideal V1 ≡ 0)
      // the closed-loop V1 amplification is structurally unstable: a 1-ULP
      // perturbation in tau1_total grows by ~1.07 per step during rupture,
      // reaching cm-scale slip_dip in ~1 s.  Locking V1 = 0 routes all
      // V_abs into V2 and zeros the dip-channel friction reaction,
      // breaking the loop at the friction-decomposition layer.  This is
      // the antiplane formulation; correct only for strike-slip.
      static const bool v1_zero_lock = []() {
         const char *e = std::getenv("SEAS_FORCE_V1_ZERO");
         return e && e[0] != '\0' &&
                !(e[0] == '0' && e[1] == '\0');
      }();
      if (v1_zero_lock)
      {
         s.V1 = 0.0;
         s.V2 = (s.tau2_total >= 0.0 ? +1.0 : -1.0) * s.V_abs;
         s.tau1_corr = 0.0;
         s.tau2_corr = s.tau2_trial - data.eta_s * s.V2;
      }
      else
      {
         // Slip rate decomposition (Eq. 9)
         s.V1 = s.V_abs * (s.tau1_total) / (strength + data.eta_s * s.V_abs);
         s.V2 = s.V_abs * (s.tau2_total) / (strength + data.eta_s * s.V_abs);

         // Step 4: Corrected traction (Eq. 10)
         s.tau1_corr = s.tau1_trial - data.eta_s * s.V1;
         s.tau2_corr = s.tau2_trial - data.eta_s * s.V2;
      }
   }
}

void FaultFaceFlux::BuildImposedState(const DOFData &data,
                                      const EvalStageState &s,
                                      const real_t *Q_plus,
                                      const real_t *Q_minus,
                                      real_t *Q_imp_plus,
                                      real_t *Q_imp_minus) const
{
   // Step 5: Construct imposed states (Eq. 11-12).  Initialize with
   // current state so non-normal stresses (SYY/SZZ/SYZ) and unused
   // velocity components carry through unchanged.
   std::memcpy(Q_imp_minus, Q_minus, NUM_STATE * sizeof(real_t));
   std::memcpy(Q_imp_plus,  Q_plus,  NUM_STATE * sizeof(real_t));

   // Minus side (Eq. 11a-d): v^{-,imp} = v^- - (1/Z)(sigma_corr - sigma^-)
   const real_t invZp_m = 1.0 / data.Zp_minus;
   const real_t invZs_m = 1.0 / data.Zs_minus;

   Q_imp_minus[VX] = Q_minus[VX] - invZp_m * (s.sigma_n_corr - Q_minus[SXX]);
   Q_imp_minus[VY] = Q_minus[VY] - invZs_m * (s.tau1_corr    - Q_minus[SXY]);
   Q_imp_minus[VZ] = Q_minus[VZ] - invZs_m * (s.tau2_corr    - Q_minus[SXZ]);

   // Plus side (Eq. 12a-d): v^{+,imp} = v^+ + (1/Z)(sigma_corr - sigma^+)
   const real_t invZp_p = 1.0 / data.Zp_plus;
   const real_t invZs_p = 1.0 / data.Zs_plus;

   Q_imp_plus[VX] = Q_plus[VX] + invZp_p * (s.sigma_n_corr - Q_plus[SXX]);
   Q_imp_plus[VY] = Q_plus[VY] + invZs_p * (s.tau1_corr    - Q_plus[SXY]);
   Q_imp_plus[VZ] = Q_plus[VZ] + invZs_p * (s.tau2_corr    - Q_plus[SXZ]);

   // Both sides: imposed stress = corrected traction (Eq. 11d/12d)
   Q_imp_minus[SXX] = s.sigma_n_corr;
   Q_imp_minus[SXY] = s.tau1_corr;
   Q_imp_minus[SXZ] = s.tau2_corr;

   Q_imp_plus[SXX]  = s.sigma_n_corr;
   Q_imp_plus[SXY]  = s.tau1_corr;
   Q_imp_plus[SXZ]  = s.tau2_corr;
}

void FaultFaceFlux::WriteBackState(DOFData &data,
                                   const EvalStageState &s) const
{
   // v9.4.0 Commit 1: store TOTAL physical traction (static +
   // nucleation + trial-scale corrected).  Riemann imposed-state
   // construction above uses TRIAL-scale `tau*_corr` to match
   // SeisSol's `tractionResults.traction*`.
   data.slip_rate = s.V_abs;
   data.V1 = s.V1;
   data.V2 = s.V2;
   data.tau1_corr    = data.tau1_0   + data.tau1_nuc    + s.tau1_corr;
   data.tau2_corr    = data.tau2_0   + data.tau2_nuc    + s.tau2_corr;
   data.sigma_n_corr = data.sigma_n0 + data.sigma_n_nuc + s.sigma_n_corr;
}

// ---------------------------------------------------------------------------
// Full Evaluate pipeline: Eq. (7) → (8) → (9) → (10) → (11)-(12)
// ---------------------------------------------------------------------------
void FaultFaceFlux::Evaluate(DOFData &data,
                             const real_t *Q_plus, const real_t *Q_minus,
                             real_t *Q_imp_plus, real_t *Q_imp_minus,
                             FrictionSolver::Method method) const
{
   // v9.2.0 F01+F02 invariant (REVIEW R-V92-H07): the driver's coupled
   // RK4 on (Q, psi) in `drivers/tpv102_driver.cpp` reads `data.psi` at
   // entry to this call as the stage-local psi and relies on it being
   // UNCHANGED on return.  A future refactor that adds a psi write here
   // would silently invalidate the RK4 stage arithmetic (the driver
   // would double-integrate psi).  Keep this guard whenever psi is not
   // a formal input to the Riemann solver.
#ifndef NDEBUG
   const real_t psi_at_entry = data.psi;
#endif

   // v9.0.0 Pelties-9 per-side flux (see
   // debug_document/tpv102_debug_document/tpv102_debug_v9.0.0_seissol_flux_comparison.md
   // §10.1 and §18 R-F08).  The call site in wave_operator.inl applies a
   // SINGLE GodunovFlux instance to both sides of the fault, which is
   // correct only when A_plus == A_minus (homogeneous material).  Guard
   // against silent use on a bimaterial face.
   // R-008: allow ~1-ULP drift — reject only when the two sides are
   // observably different materials, not when they're the same material
   // reached via two different expression chains (e.g. a future driver
   // computing Zp_plus, Zp_minus from separate rho*cp expressions).
   auto homog_ok = [](real_t a, real_t b)
   {
      return std::abs(a - b) <=
             1e-12 * std::max(std::abs(a), std::abs(b));
   };
   MFEM_VERIFY(homog_ok(data.Zp_plus, data.Zp_minus) &&
               homog_ok(data.Zs_plus, data.Zs_minus),
               "Bimaterial fault face detected (Zp_plus=" << data.Zp_plus
               << " Zp_minus=" << data.Zp_minus
               << " Zs_plus=" << data.Zs_plus
               << " Zs_minus=" << data.Zs_minus
               << ").  v9.0.0 Pelties-9 per-side flux assumes "
               "homogeneous material.  Extend GodunovFlux to per-side A "
               "before running this configuration.");

   // Round-12 Patch 1: compose the stage helpers.  Byte-identical to
   // the pre-refactor inline implementation.
   EvalStageState s;
   ComputeStageState(data, Q_plus, Q_minus, s, method);

#ifdef SEAS_DIAG_FAULT_FLUX
   // C-1 EVAL: v9.0.0 §0.5 checkpoint — printf only, no MPI.  Prints only
   // on DOFs the driver flagged as diagnostic (typically the hypocenter QP).
   if (data.diag_print)
   {
      std::fprintf(stderr,
         "[C-1 EVAL] rank=%d  tau1_trial=%+.3e Pa  tau2_trial=%+.3e Pa  "
         "|Q_plus[VY]|=%.3e  |Q_plus[SXY]|=%.3e  |Q_plus[VZ]|=%.3e  "
         "|Q_plus[SXZ]|=%.3e  psi=%.3e\n",
         g_seas_my_rank,
         s.tau1_trial, s.tau2_trial,
         std::abs(Q_plus[VY]),  std::abs(Q_plus[SXY]),
         std::abs(Q_plus[VZ]),  std::abs(Q_plus[SXZ]),
         data.psi);

      // C-1n NORMAL: trace the σ_n_trial decomposition to identify what
      // drives the trial-normal-stress perturbation we observe at the
      // hypocenter (~0.5 MPa during rupture transit).  σ_n_trial =
      // etaP * ((Q_minus[VX] − Q_plus[VX]) + Q_plus[SXX]/Zp+
      //         Q_minus[SXX]/Zp_neig).  All Q components are in
      // FAULT-LOCAL coords.  The two terms are physically:
      //   v_jump_term   = etaP * (Q_minus[VX] − Q_plus[VX])
      //                 = etaP × the fault-normal velocity jump.  Should
      //                   be ≡0 for pure strike-slip with mirror-
      //                   symmetric mesh.
      //   stress_term   = etaP * (Q_plus[SXX]/Zp + Q_minus[SXX]/Zp_neig)
      //                 = etaP × twice the fault-normal stress average.
      //                   For fluctuation-Q, both Q[SXX] should be ≡0
      //                   at fault QPs (pre-stress lives in DOFData,
      //                   not in Q).  Non-zero ⇒ wave operator's flux
      //                   back-feed has populated SXX_fluctuation, OR
      //                   the imposed-state back-flux is leaking
      //                   normal-stress into bulk Q.
      const real_t inv_Zp_p = 1.0 / data.Zp_plus;
      const real_t inv_Zp_m = 1.0 / data.Zp_minus;
      const real_t v_jump_term =
         data.eta_p * (Q_minus[VX] - Q_plus[VX]);
      const real_t stress_term =
         data.eta_p * (Q_plus[SXX] * inv_Zp_p
                     + Q_minus[SXX] * inv_Zp_m);
      std::fprintf(stderr,
         "[C-1n NORMAL] rank=%d  sigma_n_trial=%+.4e Pa  "
         "v_jump_term=%+.4e Pa  stress_term=%+.4e Pa  | "
         "Q_plus[VX]=%+.4e  Q_minus[VX]=%+.4e  "
         "Q_plus[SXX]=%+.4e  Q_minus[SXX]=%+.4e\n",
         g_seas_my_rank,
         s.sigma_n_trial, v_jump_term, stress_term,
         Q_plus[VX], Q_minus[VX], Q_plus[SXX], Q_minus[SXX]);
   }
#endif

   BuildImposedState(data, s, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus);
   WriteBackState(data, s);

#ifndef NDEBUG
   // R-V92-H07: psi must be pristine for the driver's coupled RK4 to be
   // correct.  Bit-exact equality is the right check — no arithmetic on
   // psi has happened between entry and exit.
   MFEM_ASSERT(data.psi == psi_at_entry,
               "FaultFaceFlux::Evaluate mutated data.psi "
               "(before = " << psi_at_entry << ", after = " << data.psi
               << ").  The driver's coupled-RK4-on-psi integrator assumes "
               "this function is psi-pure.  See plan Step 5 + REVIEW "
               "R-V92-H07.");
#endif
}

// ---------------------------------------------------------------------------
// Total-stress pipeline: Evaluate operating on Q that already carries the
// pre-stress tensor in every bulk DOF (v9.3.0 Phase 3, I-06 part A).
//
// Differences from Evaluate:
//   - Inputs Q_plus / Q_minus are TOTAL (pre-stress + fluctuation).
//     ComputeTrialTraction on TOTAL Q therefore produces a TOTAL trial
//     traction directly — no additional +sigma_n0/+tau_i_0 shift needed.
//   - DOFData.sigma_n0 / tau1_0 / tau2_0 are assumed ZEROED by the driver
//     under the Phase 4 migration so pre-stress is not double-counted.
//   - DOFData.sigma_n_corr / tau1_corr / tau2_corr are stored as TOTAL
//     directly (no sigma_n0 add).  This matches SeisSol's
//     RateAndState.h:222-240 convention.
//
// The imposed-state construction (Step 4) is byte-identical to the
// fluctuation path: the differences (sigma_n_corr - Q_{plus,minus}[SXX])
// are invariant under a constant pre-stress shift of both sides.
//
// R-003 guards: #ifndef NDEBUG psi-invariance entry/exit MFEM_ASSERT +
// SEAS_DIAG_FAULT_FLUX diag print block after ComputeTrialTraction.
// ---------------------------------------------------------------------------
void FaultFaceFlux::EvaluateTotal(DOFData &data,
                                  const real_t *Q_plus, const real_t *Q_minus,
                                  real_t *Q_imp_plus, real_t *Q_imp_minus,
                                  FrictionSolver::Method method) const
{
#ifndef NDEBUG
   // R-003 (review-incorporation plan): the driver's coupled-RK4-on-psi
   // integrator requires this function to be psi-pure.  Enforce the same
   // bit-exact guard as in Evaluate.
   const real_t psi_at_entry = data.psi;

   // R-I06-007 (review-round-2): EvaluateTotal's correctness contract is
   // that bulk Q carries the pre-stress, so DOFData's pre-stress fields
   // MUST be zero — the fluctuation-mode seeding done by
   // InitializeFaultDOFs must be followed by ZeroDOFDataPreStressTotal
   // before any EvaluateTotal call.  A future refactor that removes the
   // zeroing (or a "unification" that folds pre-stress back into
   // ComputeTrialTraction) would silently double-count the pre-stress
   // in the trial traction.  Assert the contract at runtime.
   MFEM_ASSERT(data.sigma_n0 == 0.0 &&
               data.tau1_0   == 0.0 &&
               data.tau2_0   == 0.0,
               "EvaluateTotal contract violated: DOFData pre-stress "
               "fields must be zero (sigma_n0=" << data.sigma_n0 <<
               ", tau1_0=" << data.tau1_0 <<
               ", tau2_0=" << data.tau2_0 <<
               ").  Total-stress drivers must call "
               "ZeroDOFDataPreStressTotal after InitializeFaultDOFs.");
#endif

   // Homogeneous check — identical to Evaluate.  v9.0.0 Pelties-9 per-side
   // flux assumes A_plus == A_minus on a fault face.
   auto homog_ok = [](real_t a, real_t b)
   {
      return std::abs(a - b) <=
             1e-12 * std::max(std::abs(a), std::abs(b));
   };
   MFEM_VERIFY(homog_ok(data.Zp_plus, data.Zp_minus) &&
               homog_ok(data.Zs_plus, data.Zs_minus),
               "Bimaterial fault face detected (Zp_plus=" << data.Zp_plus
               << " Zp_minus=" << data.Zp_minus
               << " Zs_plus=" << data.Zs_plus
               << " Zs_minus=" << data.Zs_minus
               << ").  EvaluateTotal assumes homogeneous material.  "
               "Extend GodunovFlux / FaultFaceFlux to per-side A before "
               "running this configuration.");

   // Step 1: Trial traction from bulk Q.  Under total-Q, bulk Q
   //         carries the static background prestress, so this trial is
   //         "background + wave" — but does NOT yet include the
   //         persistent nucleation channel.
   real_t sigma_n_trial, tau1_trial, tau2_trial;
   ComputeTrialTraction(data, Q_plus, Q_minus,
                        sigma_n_trial, tau1_trial, tau2_trial);

   // Step 1b: friction inputs = trial + persistent-nucleation channel.
   // Mirrors SeisSol's RateAndState.h:222-240 pattern, where
   // `totalTraction* = initialStressInFaultCS* + faultStresses.traction*`
   // is a temporary used ONLY for the friction solve and the
   // slip-rate divisor — never fed back into `tractionResults`
   // (which is the analog of our `tau*_corr` consumed by the Riemann
   // imposed state below).  Mutating `tau*_trial` here would leak
   // `tau*_nuc` into the Riemann velocity-jump
   // `(2/Zs)·(tau_corr - Q_bulk[SXY])`, producing a ~24× spurious
   // shear pulse into bulk Q every step at full nucleation
   // (tau2_nuc=25 MPa, V_abs≈0.22 m/s, Zs≈9.25 MPa·s/m).  Keep the
   // two scales separate.
   const real_t sigma_n_fric = sigma_n_trial + data.sigma_n_nuc;
   const real_t tau1_fric    = tau1_trial    + data.tau1_nuc;
   const real_t tau2_fric    = tau2_trial    + data.tau2_nuc;

#ifdef SEAS_DIAG_FAULT_FLUX
   // R-003: mirror Evaluate's C-1 checkpoint for total mode.  Same
   // diag_print gate, same fprintf semantics (no MPI).  Prints only
   // on DOFs the driver flagged as diagnostic.  Print the friction-
   // input total (what the solver actually sees); a separate trial
   // value can be reconstructed as tau*_fric - data.tau*_nuc if needed.
   if (data.diag_print)
   {
      std::fprintf(stderr,
         "[C-1 EVAL-TOTAL] rank=%d  tau1_fric=%+.3e Pa  "
         "tau2_fric=%+.3e Pa  sigma_n_fric=%+.3e Pa  psi=%.3e  "
         "a=%.3e  eta_s=%.3e  (tau2_nuc=%+.3e, sigma_n_nuc=%+.3e)\n",
         g_seas_my_rank, tau1_fric, tau2_fric, sigma_n_fric, data.psi,
         data.a, data.eta_s, data.tau2_nuc, data.sigma_n_nuc);
   }
#endif

   // Step 2: Solve friction equation for |V| on the FRICTION-INPUT
   //         total (trial + nuc).  Theta is the magnitude of the
   //         total tangential traction the friction law must balance.
   const real_t Theta = std::sqrt(tau1_fric * tau1_fric
                                + tau2_fric * tau2_fric);
   real_t V_abs = 0.0;
   if (Theta > 0.0)
   {
      V_abs = solver_.Solve(Theta, data.psi, std::abs(sigma_n_fric),
                            data.eta_s, data.a, method);
   }

   // Step 3: Slip-rate decomposition uses the TOTAL traction in the
   //         numerator (SeisSol numerator is `totalTraction*`).
   //         Riemann-side corrected traction tau*_corr stays on the
   //         TRIAL scale (no nuc).
   real_t V1 = 0.0, V2 = 0.0;
   real_t tau1_corr = tau1_trial, tau2_corr = tau2_trial;
   if (Theta > 0.0 && V_abs > 0.0)
   {
      real_t C = std::exp(data.psi / data.a) / (2.0 * FrictionSolver::V0);
      real_t f_V = data.a * std::asinh(V_abs * C);
      // σ_n strength floor (sliver-blowup plan 2026-05-26): enabled
      // (floor >= 0) ⇒ max(σ_n_fric, floor); disabled (sentinel < 0) ⇒
      // historical |σ_n_fric| (byte-exact TPV102/104).  Only the explicit
      // strength is floored; the friction-solver σ_n argument (:548) is
      // intentionally NOT floored in v1 (plan §Phase 3 decision 2).
      const real_t sn_str = (sigma_n_strength_floor_ >= 0.0)
                            ? std::max(sigma_n_fric, sigma_n_strength_floor_)
                            : std::abs(sigma_n_fric);
      real_t strength = sn_str * f_V;
      V1 = V_abs * tau1_fric / (strength + data.eta_s * V_abs);
      V2 = V_abs * tau2_fric / (strength + data.eta_s * V_abs);
      tau1_corr = tau1_trial - data.eta_s * V1;
      tau2_corr = tau2_trial - data.eta_s * V2;
   }

   // Step 4: Imposed states.  `*_corr` here is TRIAL-scale (matches
   //         SeisSol `tractionResults.traction*`), so the Riemann
   //         velocity jump `(2/Zs)·(tau_corr - Q_bulk[SXY])` carries
   //         only the friction reaction `eta_s·V` — the persistent
   //         nucleation channel does NOT radiate through bulk Q.
   const real_t sigma_n_corr = sigma_n_trial;
   std::memcpy(Q_imp_minus, Q_minus, NUM_STATE * sizeof(real_t));
   std::memcpy(Q_imp_plus,  Q_plus,  NUM_STATE * sizeof(real_t));

   const real_t invZp_m = 1.0 / data.Zp_minus;
   const real_t invZs_m = 1.0 / data.Zs_minus;
   const real_t invZp_p = 1.0 / data.Zp_plus;
   const real_t invZs_p = 1.0 / data.Zs_plus;

   Q_imp_minus[VX] = Q_minus[VX] - invZp_m * (sigma_n_corr - Q_minus[SXX]);
   Q_imp_minus[VY] = Q_minus[VY] - invZs_m * (tau1_corr    - Q_minus[SXY]);
   Q_imp_minus[VZ] = Q_minus[VZ] - invZs_m * (tau2_corr    - Q_minus[SXZ]);
   Q_imp_plus[VX]  = Q_plus[VX]  + invZp_p * (sigma_n_corr - Q_plus[SXX]);
   Q_imp_plus[VY]  = Q_plus[VY]  + invZs_p * (tau1_corr    - Q_plus[SXY]);
   Q_imp_plus[VZ]  = Q_plus[VZ]  + invZs_p * (tau2_corr    - Q_plus[SXZ]);

   Q_imp_minus[SXX] = sigma_n_corr;  Q_imp_plus[SXX] = sigma_n_corr;
   Q_imp_minus[SXY] = tau1_corr;     Q_imp_plus[SXY] = tau1_corr;
   Q_imp_minus[SXZ] = tau2_corr;     Q_imp_plus[SXZ] = tau2_corr;

   // Step 5: Update DOFData — store TOTAL (trial + nuc) for station
   //         output and post-step diagnostics.  Mirrors the
   //         fluctuation-path `data.tau*_corr = data.tau*_0 + tau*_corr`
   //         pattern: post-call DOFData carries the full physical
   //         traction the fault is bearing, while the Riemann path
   //         above operated on the trial scale.
   data.slip_rate    = V_abs;
   data.V1           = V1;
   data.V2           = V2;
   data.tau1_corr    = tau1_corr   + data.tau1_nuc;
   data.tau2_corr    = tau2_corr   + data.tau2_nuc;
   data.sigma_n_corr = sigma_n_corr + data.sigma_n_nuc;

#ifndef NDEBUG
   MFEM_ASSERT(data.psi == psi_at_entry,
               "FaultFaceFlux::EvaluateTotal mutated data.psi "
               "(before = " << psi_at_entry << ", after = " << data.psi
               << ").  R-V92-H07: the driver's coupled-RK4-on-psi "
               "integrator assumes this function is psi-pure.");
#endif
}

// ---------------------------------------------------------------------------
// ADER I-05 Phase 5: time-integrated Riemann solve (fluctuation variant).
// ---------------------------------------------------------------------------
// I± = ∫_0^{dt} Q±(τ) dτ  ⇒  Q̄± = I±/dt.  Call the standard Evaluate on
// the time-averaged state (which now represents the mean Q over [t_n,
// t_n+dt]); the friction solver returns tractions / slip-rate that are
// themselves time-averaged in the same sense.  Finally rescale the
// imposed Q-states back to time-integrated form: I_imp± = dt · Q_imp±.
//
// Proof of O(dt²) consistency with a midpoint-RK4 call: let
// Q(τ) = Q_m + O(dt),  Q_m = Q(t_n + dt/2).  Then Q̄ = Q_m + O(dt²) by
// Simpson's rule on smooth Q, so feeding Q̄ into the nonlinear friction
// equation gives the same root as feeding Q_m up to O(dt²) (Brent's
// root of a Lipschitz-smooth residual commutes with O(dt²) input
// perturbations).  Imposed states scale linearly in the inputs modulo
// the friction-law nonlinearity, again O(dt²).  Plan §Phase 5 §4.
// ---------------------------------------------------------------------------
void FaultFaceFlux::EvaluateADER(DOFData &data,
                                 const real_t *I_plus, const real_t *I_minus,
                                 real_t dt,
                                 real_t *I_imp_plus, real_t *I_imp_minus,
                                 FrictionSolver::Method method) const
{
   MFEM_VERIFY(dt > 0.0,
               "FaultFaceFlux::EvaluateADER: dt must be > 0, got " << dt);

   real_t Q_avg_plus[NUM_STATE], Q_avg_minus[NUM_STATE];
   const real_t inv_dt = 1.0 / dt;
   for (int c = 0; c < NUM_STATE; c++)
   {
      Q_avg_plus[c]  = I_plus[c]  * inv_dt;
      Q_avg_minus[c] = I_minus[c] * inv_dt;
   }

   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];

   // Round-11 FREEZE-A hook (test-only).  When SEAS_TEST_FREEZE_A=1,
   // bypass the friction Evaluate entirely: leave `data` unchanged
   // and pass the time-averaged bulk state through as the imposed
   // state.  This makes the fault face apply no flux correction
   // ("fault locked") so the pepper-guard test can distinguish
   // "asymmetry driven by fault friction / feedback" (closes under
   // FREEZE-A) from "asymmetry generated by the MFEM wave update
   // even with the fault locked" (persists under FREEZE-A).
   // Minimal hook — production builds with SEAS_TEST_FREEZE_A unset
   // are byte-identical to the pre-hook path.
   const char *freeze_a = std::getenv("SEAS_TEST_FREEZE_A");
   if (freeze_a && freeze_a[0] == '1')
   {
      for (int c = 0; c < NUM_STATE; c++)
      {
         Q_imp_plus[c]  = Q_avg_plus[c];
         Q_imp_minus[c] = Q_avg_minus[c];
      }
   }
   else
   {
      Evaluate(data, Q_avg_plus, Q_avg_minus, Q_imp_plus, Q_imp_minus, method);
   }

   for (int c = 0; c < NUM_STATE; c++)
   {
      I_imp_plus[c]  = Q_imp_plus[c]  * dt;
      I_imp_minus[c] = Q_imp_minus[c] * dt;
   }
}

// ---------------------------------------------------------------------------
// ADER I-05 Phase 5 + v9.3.0 §Phase 7: time-integrated Riemann solve
// (total-stress variant).  Same 1/dt ↔ dt scaling wrapping pattern as
// EvaluateADER, but dispatches to EvaluateTotal so TPV102's post-I-06
// total-Q path has a direct ADER entry point without the workaround the
// v9.3.0 plan describes.
// ---------------------------------------------------------------------------
void FaultFaceFlux::EvaluateADERTotal(DOFData &data,
                                      const real_t *I_plus_tot,
                                      const real_t *I_minus_tot,
                                      real_t dt,
                                      real_t *I_imp_plus_tot,
                                      real_t *I_imp_minus_tot,
                                      FrictionSolver::Method method) const
{
   MFEM_VERIFY(dt > 0.0,
               "FaultFaceFlux::EvaluateADERTotal: dt must be > 0, got " << dt);

   real_t Q_avg_plus[NUM_STATE], Q_avg_minus[NUM_STATE];
   const real_t inv_dt = 1.0 / dt;
   for (int c = 0; c < NUM_STATE; c++)
   {
      Q_avg_plus[c]  = I_plus_tot[c]  * inv_dt;
      Q_avg_minus[c] = I_minus_tot[c] * inv_dt;
   }

   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
   EvaluateTotal(data, Q_avg_plus, Q_avg_minus, Q_imp_plus, Q_imp_minus,
                 method);

   for (int c = 0; c < NUM_STATE; c++)
   {
      I_imp_plus_tot[c]  = Q_imp_plus[c]  * dt;
      I_imp_minus_tot[c] = Q_imp_minus[c] * dt;
   }
}

// ---------------------------------------------------------------------------
// REVIEW R-016: LSW counterpart to EvaluateADER.  Wraps the I/dt → I·dt
// pattern around the closed-form LSW solve so the wave operator's fault
// dispatch (interior + R-1600 shared-fault fallback) can route TPV205
// through correct physics instead of Brent on the rate-and-state law.
// Reads ONLY the LSW-native fields data.lsw_mu_s / lsw_mu_d / lsw_d_c;
// data.a / data.psi / data.Dc are not consumed.
// ---------------------------------------------------------------------------
void FaultFaceFlux::EvaluateADER_LSW(DOFData &data,
                                     const real_t *I_plus,
                                     const real_t *I_minus,
                                     real_t dt,
                                     real_t *I_imp_plus,
                                     real_t *I_imp_minus) const
{
   MFEM_VERIFY(dt > 0.0,
               "FaultFaceFlux::EvaluateADER_LSW: dt must be > 0, got " << dt);

   // Misuse guard: every TPV205 QP populates at least one LSW field
   // strictly > 0 (μ_s_barrier ≥ 10000 outside rupture area, μ_s = 0.677
   // inside; μ_d = 0.525; d_c = 0.4).  All-zero fields means the caller
   // routed a non-LSW DOFData here — abort loudly rather than silently
   // computing strength = 0 ⇒ unconstrained sliding.
   MFEM_VERIFY(data.lsw_mu_s > 0.0 || data.lsw_mu_d > 0.0
               || data.lsw_d_c > 0.0,
               "FaultFaceFlux::EvaluateADER_LSW: all LSW-native fields "
               "are zero (lsw_mu_s=" << data.lsw_mu_s
               << " lsw_mu_d=" << data.lsw_mu_d
               << " lsw_d_c=" << data.lsw_d_c
               << ").  This DOFData was not initialized by "
               "InitializeFaultDOFs_TPV205; the wave operator dispatched "
               "the LSW path on rate-and-state data.");

   // Homogeneous-material check (mirrors Evaluate / EvaluateTotal).
   // v9.0.0 Pelties-9 per-side flux assumes A_plus == A_minus.
   auto homog_ok = [](real_t a, real_t b)
   {
      return std::abs(a - b) <=
             static_cast<real_t>(1e-12) * std::max(std::abs(a), std::abs(b));
   };
   MFEM_VERIFY(homog_ok(data.Zp_plus, data.Zp_minus) &&
               homog_ok(data.Zs_plus, data.Zs_minus),
               "FaultFaceFlux::EvaluateADER_LSW: bimaterial fault face "
               "(Zp_plus=" << data.Zp_plus << " Zp_minus=" << data.Zp_minus
               << " Zs_plus=" << data.Zs_plus << " Zs_minus=" << data.Zs_minus
               << ").  Extend per-side handling before running this "
               "configuration.");

   // Step 0: Q̄± = I±/dt.
   real_t Q_avg_plus[NUM_STATE], Q_avg_minus[NUM_STATE];
   const real_t inv_dt = static_cast<real_t>(1.0) / dt;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      Q_avg_plus[c]  = I_plus[c]  * inv_dt;
      Q_avg_minus[c] = I_minus[c] * inv_dt;
   }

   // Step 1: trial traction (pure helper).
   EvalStageState s;
   ComputeTrialTraction(data, Q_avg_plus, Q_avg_minus,
                        s.sigma_n_trial, s.tau1_trial, s.tau2_trial);

   // Step 2: total traction.  TPV205 has zero nucleation channels
   // (sigma_n_nuc / tau*_nuc kept at init defaults of 0); pre-stress
   // (sigma_n0 / tau*_0) lives in DOFData.
   s.sigma_n_total = data.sigma_n0 + data.sigma_n_nuc + s.sigma_n_trial;
   s.tau1_total    = data.tau1_0   + data.tau1_nuc    + s.tau1_trial;
   s.tau2_total    = data.tau2_0   + data.tau2_nuc    + s.tau2_trial;
   s.Theta         = std::sqrt(s.tau1_total * s.tau1_total
                              + s.tau2_total * s.tau2_total);

   // Step 3: μ_eff(δ) at the slip magnitude carried in DOFData.  Reads
   // ONLY the LSW-native fields; the strength-barrier short-circuit
   // (R-002) lives inside the helper.
   const real_t delta = std::sqrt(data.slip1 * data.slip1
                                  + data.slip2 * data.slip2);
   const real_t mu_eff = LSWFrictionCoefficient_TPV205(delta,
                                                       data.lsw_mu_s,
                                                       data.lsw_mu_d,
                                                       data.lsw_d_c);

   // Step 4: closed-form LSW solve.  Sets s.V_abs, s.V{1,2}, s.tau{1,2}_corr.
   // R-003 barrier short-circuit (V = 0 in the barrier zone regardless of
   // σ_n sign) lives inside the helper.  Trailing arg = σ_n strength floor
   // (sliver-blowup plan 2026-05-26); disabled (sentinel < 0) maps to 0.0
   // ⇒ byte-exact `max(σ_n,0)`.
   SolveLSW_TPV205(s.tau1_trial, s.tau2_trial,
                   s.tau1_total, s.tau2_total,
                   s.sigma_n_total, data.eta_s,
                   mu_eff,
                   s.V_abs, s.V1, s.V2,
                   s.tau1_corr, s.tau2_corr,
                   SigmaNStrengthFloorForLSW(),
                   data.lsw_cohesion);   // Phase 10 (TPV31) additive C0 (0 ⇒ byte-exact)

   // σ_n is unaffected by friction — TRIAL-scale value matches the
   // rate-and-state path's CompleteFromVabs convention.
   s.sigma_n_corr = s.sigma_n_trial;

   // R-001 (final review): slip accumulation is the iterator's
   // responsibility (Tpv205SubStepIterator::StepOneQP_ updates
   // data.slip{1,2} once per sub-step for every dof_data entry,
   // including shared QPs).  The wave-operator dispatch (interior +
   // R-1601 shared-fault fallback) calls EvaluateADER_LSW AFTER the
   // iterator on shared QPs, so adding `data.slip{1,2} += V*dt` here
   // would double-count slip on shared faces and break the rupture
   // physics at np > 1.  EvaluateADER_LSW now only computes I_imp /
   // V / τ_corr / σ_n_corr; slip evolution is owned exclusively by
   // the iterator.

   // Step 6: imposed Q-state (Eq. 11-12 of the FaultFaceFlux pipeline).
   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
   BuildImposedState(data, s, Q_avg_plus, Q_avg_minus,
                     Q_imp_plus, Q_imp_minus);

   // Step 7: write back V/slip_rate/τ*_corr/σ_n_corr to DOFData (TOTAL
   // physical traction — pre + nuc + trial-scale corrected).
   WriteBackState(data, s);

   // Step 8: rescale Q_imp back to time-integrated form.
   for (int c = 0; c < NUM_STATE; ++c)
   {
      I_imp_plus[c]  = Q_imp_plus[c]  * dt;
      I_imp_minus[c] = Q_imp_minus[c] * dt;
   }
}

// ---------------------------------------------------------------------------
// Phase H.6 of spatial_dynamic_rupture_plan.md (rev-3): LSW Riemann solve
// with the TPV26/27 time-dependent forced-rupture friction coefficient.
//
// Implementation is a verbatim clone of EvaluateADER_LSW above with
// exactly one substitution at Step 3: the `LSWFrictionCoefficient_TPV205`
// call becomes `LSWFrictionCoefficient_ForcedRupture(delta, ..., t_now,
// data.T_forced_rupture, data.t0_decay_forced)`.  All other steps
// (trial-traction, total-traction, SolveLSW_TPV205, BuildImposedState,
// WriteBackState, I_imp rescaling) are byte-equivalent to the LSW
// path, so the BYTE-EXACT contract on TPV205-init DOFData
// (T_forced_rupture = 1.0e9, t0_decay_forced = 0 ⇒ f_2 = 0 ⇒ μ reduces
// to plain LSW) holds.
// ---------------------------------------------------------------------------
void FaultFaceFlux::EvaluateADER_LSW_ForcedRupture(
   DOFData &data,
   const real_t *I_plus,
   const real_t *I_minus,
   real_t dt,
   real_t t_now,
   real_t *I_imp_plus,
   real_t *I_imp_minus) const
{
   MFEM_VERIFY(dt > 0.0,
               "FaultFaceFlux::EvaluateADER_LSW_ForcedRupture: dt must be > 0, "
               "got " << dt);

   MFEM_VERIFY(data.lsw_mu_s > 0.0 || data.lsw_mu_d > 0.0
               || data.lsw_d_c > 0.0,
               "FaultFaceFlux::EvaluateADER_LSW_ForcedRupture: all LSW-native "
               "fields are zero (lsw_mu_s=" << data.lsw_mu_s
               << " lsw_mu_d=" << data.lsw_mu_d
               << " lsw_d_c=" << data.lsw_d_c
               << ").  This DOFData was not initialized for LSW; the wave "
               "operator dispatched the LSW_ForcedRupture path on rate-and-"
               "state data.");

   auto homog_ok = [](real_t a, real_t b)
   {
      return std::abs(a - b) <=
             static_cast<real_t>(1e-12) * std::max(std::abs(a), std::abs(b));
   };
   MFEM_VERIFY(homog_ok(data.Zp_plus, data.Zp_minus) &&
               homog_ok(data.Zs_plus, data.Zs_minus),
               "FaultFaceFlux::EvaluateADER_LSW_ForcedRupture: bimaterial "
               "fault face (Zp_plus=" << data.Zp_plus
               << " Zp_minus=" << data.Zp_minus
               << " Zs_plus=" << data.Zs_plus
               << " Zs_minus=" << data.Zs_minus
               << ").  Extend per-side handling before running this "
               "configuration.");

   real_t Q_avg_plus[NUM_STATE], Q_avg_minus[NUM_STATE];
   const real_t inv_dt = static_cast<real_t>(1.0) / dt;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      Q_avg_plus[c]  = I_plus[c]  * inv_dt;
      Q_avg_minus[c] = I_minus[c] * inv_dt;
   }

   EvalStageState s;
   ComputeTrialTraction(data, Q_avg_plus, Q_avg_minus,
                        s.sigma_n_trial, s.tau1_trial, s.tau2_trial);

   s.sigma_n_total = data.sigma_n0 + data.sigma_n_nuc + s.sigma_n_trial;
   s.tau1_total    = data.tau1_0   + data.tau1_nuc    + s.tau1_trial;
   s.tau2_total    = data.tau2_0   + data.tau2_nuc    + s.tau2_trial;
   s.Theta         = std::sqrt(s.tau1_total * s.tau1_total
                              + s.tau2_total * s.tau2_total);

   // Phase H.6 substitution: time-dependent friction coefficient that
   // reads the new D-4 DOFData fields T_forced_rupture / t0_decay_forced.
   const real_t delta = std::sqrt(data.slip1 * data.slip1
                                  + data.slip2 * data.slip2);
   const real_t mu_eff = spatial::LSWFrictionCoefficient_ForcedRupture(
      delta, data.lsw_mu_s, data.lsw_mu_d, data.lsw_d_c,
      t_now, data.T_forced_rupture, data.t0_decay_forced);

   // Trailing arg = σ_n strength floor (sliver-blowup plan 2026-05-26);
   // disabled (sentinel < 0) maps to 0.0 ⇒ byte-exact `max(σ_n,0)`.
   SolveLSW_TPV205(s.tau1_trial, s.tau2_trial,
                   s.tau1_total, s.tau2_total,
                   s.sigma_n_total, data.eta_s,
                   mu_eff,
                   s.V_abs, s.V1, s.V2,
                   s.tau1_corr, s.tau2_corr,
                   SigmaNStrengthFloorForLSW(),
                   data.lsw_cohesion);   // Phase 10 (TPV31) additive C0 (0 ⇒ byte-exact)

   s.sigma_n_corr = s.sigma_n_trial;

   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
   BuildImposedState(data, s, Q_avg_plus, Q_avg_minus,
                     Q_imp_plus, Q_imp_minus);

   WriteBackState(data, s);

   for (int c = 0; c < NUM_STATE; ++c)
   {
      I_imp_plus[c]  = Q_imp_plus[c]  * dt;
      I_imp_minus[c] = Q_imp_minus[c] * dt;
   }
}

} // namespace seas
} // namespace mfem
