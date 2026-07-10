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

#ifndef MFEM_SEAS_FAULT_FACE_FLUX_HPP
#define MFEM_SEAS_FAULT_FACE_FLUX_HPP

#include "mfem.hpp"
#include "wave_state.hpp"
#include "friction_solver.hpp"

#include <limits>

namespace mfem
{
namespace seas
{

/// Per-DOF fault data: impedances, initial stress, state variable, slip rate.
struct DOFData
{
   real_t Zp_plus = 0, Zp_minus = 0;    ///< P-impedance on ± sides
   real_t Zs_plus = 0, Zs_minus = 0;    ///< S-impedance on ± sides
   real_t eta_p = 0, eta_s = 0;          ///< Harmonic mean impedances
   real_t sigma_n0 = 0;                  ///< Background normal stress (>0 compression)
   real_t tau1_0 = 0, tau2_0 = 0;        ///< Background shear pre-stress
   /// Nucleation prestress (TPV102 total-Q persistent-driver channel).
   /// Read by FaultFaceFlux::EvaluateTotal each call and added to the
   /// trial traction so the persistent driver is re-imposed at every
   /// friction solve (SeisSol initialStressInFaultCS analog).  Distinct
   /// from {sigma_n0, tau1_0, tau2_0}: those carry STATIC background
   /// prestress and are zeroed under the total-Q dispatch contract
   /// (background prestress lives in bulk Q for total-Q).  tau*_nuc
   /// carries the time-varying nucleation perturbation that must NOT
   /// be poked into bulk Q (the wave operator radiates point sources
   /// away in O(h/cp), so a nucleation amplitude written to bulk Q
   /// dilutes ~10^4-10^5x before the friction solver sees it).
   /// For TPV102 only tau2_nuc is written (pure strike-slip);
   /// sigma_n_nuc and tau1_nuc are kept zero by default.
   real_t sigma_n_nuc = 0;               ///< Nucleation normal-stress driver
   real_t tau1_nuc = 0, tau2_nuc = 0;    ///< Nucleation shear-traction driver
   real_t a = 0.004;                      ///< Direct effect parameter
   real_t b = std::numeric_limits<real_t>::quiet_NaN();
                                          ///< RS state-evolution parameter (per-DOF; set by
                                          ///< every RS init path: TPV102 = TPV102Params::b,
                                          ///< SAFS-RS = rs.b(i)).  R-028: NaN default (NOT
                                          ///< 0.0).  LSW DOFData never reads it; an RS path
                                          ///< that forgot to set it propagates NaN through
                                          ///< the aging-law UpdateStateAnalytic (psi -> NaN),
                                          ///< caught by the blow-up / equilibrium checks.  A
                                          ///< 0.0 default would be SILENT-WRONG: for psi<f0
                                          ///< the analytic update snaps psi to f0 with no
                                          ///< NaN.  SeedEquilibriumPsi_RS also guards d.b>0.
   real_t Dc = 0.14;                      ///< Critical slip distance [m]
   real_t psi = 0;                        ///< State variable (logarithmic)
   real_t slip_rate = 0;                  ///< Current slip rate |V| [m/s]
   /// Max |V| over the CURRENT macro step's sub-step writes (iterator
   /// sub-steps + shared-fault macro solve).  slip_rate records only the
   /// last sub-step, so it is blind to an intermediate-node spike that
   /// still integrates into slip; this is the honest peak the blow-up
   /// monitor should sample.  TRANSIENT: reset to 0 each macro step by the
   /// driver, recomputed every step — intentionally NOT serialized in the
   /// checkpoint (PLAN_speckle_slip_runaway Phase 1, R-004 substrate).
   real_t slip_rate_substep_max = 0;      ///< Per-macro-step max |V| [m/s]
   /// MOST-TENSILE (minimum) sub-step normal traction sigma_n_total [Pa]
   /// (sliver_blowup speckle diag 2026-05-26).  The end-of-step sigma_n_corr
   /// can RECOVER to compressive after a sub-step tensile transient triggers
   /// the speckle, so this captures the worst sub-step value any QP saw.
   /// TRANSIENT and NOT serialized.  Accumulates the min over the whole DIAG
   /// INTERVAL (reset to numeric_limits::max() by the driver's [DIAG-SIGN]
   /// block AFTER it reports, NOT every macro step -- so a tensile transient
   /// on a step the diag does not sample is still captured at the next print).
   real_t sigma_n_substep_min = std::numeric_limits<real_t>::max();
   real_t V1 = 0, V2 = 0;               ///< Slip rate components [m/s] (from Eq. 9)
   real_t slip1 = 0, slip2 = 0;          ///< Accumulated slip components

   // Corrected traction from Riemann solver (populated by FaultFaceFlux::Evaluate)
   real_t tau1_corr = 0, tau2_corr = 0;  ///< Corrected tangential traction [Pa]
   real_t sigma_n_corr = 0;              ///< Corrected normal traction [Pa]

   // (Part C, TPV6/7) Per-side IMPOSED (split-node) particle velocity in the
   // fault-local frame: [0] = normal (VX), [1] = dip (VY), [2] = strike (VZ).
   // The two sides differ (bi-material): v_imp_plus - v_imp_minus = slip rate.
   // Written by Evaluate / EvaluateLSW / EvaluateADER_LSW (the imposed Godunov
   // state — the per-side velocity an on-fault SCEC TPV6 station reports); read
   // by the TPV6 per-side station writer (dynamic/tpv6_stations.hpp).  Default 0
   // => byte-exact for every non-TPV6 problem (other writers never read them).
   real_t v_imp_plus[3]  = {0.0, 0.0, 0.0};
   real_t v_imp_minus[3] = {0.0, 0.0, 0.0};

   // -----------------------------------------------------------------------
   // Linear slip-weakening (LSW) parameters — separate slot from the
   // rate-and-state `a / psi / Dc` so that:
   //   1. A code path that reads `data.a / data.psi / data.Dc` for
   //      rate-and-state friction (Evaluate, EvaluateTotal, the TPV104
   //      iterator, every Brent / Newton call) cannot accidentally
   //      consume LSW values and produce off-by-orders-of-magnitude
   //      strengths (REVIEW R-016).
   //   2. A code path that reads these LSW fields cannot accidentally
   //      consume rate-and-state values: zero defaults make
   //      `LSWFrictionCoefficient_TPV205(δ, 0, 0, 0)` deterministic
   //      and `EvaluateADER_LSW` aborts loudly via its own guard.
   //
   // Populated only by InitializeFaultDOFs_TPV205; pre-stress / nucleation /
   // impedance fields above are shared with rate-and-state callers.
   real_t lsw_mu_s = 0.0;   ///< LSW static friction μ_s (≥ mu_s_barrier ⇒ barrier QP)
   real_t lsw_mu_d = 0.0;   ///< LSW dynamic friction μ_d
   real_t lsw_d_c  = 0.0;   ///< LSW slip-weakening critical distance d_c [m]
   real_t lsw_cohesion = 0.0;  ///< Phase 10 (TPV31): LSW cohesion C0 [Pa],
                               ///< ADDITIVE to the strength (mu_eff·σ_n + C0).
                               ///< Default 0 ⇒ byte-exact for TPV205 / other LSW.

   // Phase H.7 of spatial_dynamic_rupture_plan.md (rev-3): TPV26/27
   // gradual forced-rupture fields.  In-class defaults are designed so
   // that any DOFData not explicitly initialised for forced rupture
   // produces f_2(t) = 0 for every physically reachable simulation
   // time, which makes
   // `EvaluateADER_LSW_ForcedRupture` (Phase H.6) bit-equivalent to
   // the existing `EvaluateADER_LSW` on TPV205-init DOFData (the
   // TPV205 byte-exact contract for the LSW path).
   //
   // - T_forced_rupture = 1.0e9 s is the "never forced" sentinel — for
   //   any physically reachable t the f_2(t) branch always returns 0.
   // - t0_decay_forced  = 0.0   is safe in combination with the
   //   sentinel above because the `t < T_forced` branch in
   //   LSWFrictionCoefficient_ForcedRupture fires unconditionally for
   //   the sentinel; the t0_decay-as-divisor branch is never entered.
   //
   // SAFS dynamic-rupture drivers OVERWRITE both fields per-DOF via
   // `InitializeFaultDOFs_Spatial` using values from
   // `SpatialFrictionResolver::ResolveForcedRupture`.
   real_t T_forced_rupture = 1.0e9;  ///< time-of-forced-rupture [s]
   real_t t0_decay_forced  = 0.0;    ///< forced-rupture decay time [s]

#ifdef SEAS_DIAG_FAULT_FLUX
   // v9.0.0 §0.5 DIAG gate for the C-1 / C-2 / C-3 bisection checkpoints.
   // Set true on exactly one hypocenter DOF (and optionally one off-hypo
   // witness) by the driver at init; all other DOFs keep diag_print=false
   // and the fprintf blocks are skipped.  Entire field is compiled out in
   // production builds so struct layout matches pre-change byte-for-byte.
   bool diag_print = false;
#endif
};

/// Per-QP working state for the `FaultFaceFlux::Evaluate` pipeline
/// (Round-12 Patch 1).  Exposes every intermediate of the friction
/// solve so that callers can (a) compute the full chain via
/// `ComputeStageState`, (b) overwrite one stage (e.g. face-averaged
/// trial traction in `wave_operator.inl`), and (c) drive the
/// downstream recompute via one of the four `CompleteFrom*` helpers.
/// Fields are populated in stage order; early stages are valid after
/// `ComputeTrialTraction` alone, later stages after their respective
/// helper runs.
struct EvalStageState
{
   real_t sigma_n_trial = 0, tau1_trial = 0, tau2_trial = 0;
   real_t sigma_n_total = 0, tau1_total = 0, tau2_total = 0;
   real_t Theta = 0;
   real_t V_abs = 0;
   real_t V1 = 0, V2 = 0;
   real_t sigma_n_corr = 0, tau1_corr = 0, tau2_corr = 0;
};

/// @brief Fault-face Riemann solver for dynamic rupture.
///
/// Implements the trial-and-correction approach (plan Section 2.4):
/// 1. Compute trial traction from Godunov state (Eq. 7)
/// 2. Solve friction equation for slip rate (Eq. 8-9)
/// 3. Compute corrected traction (Eq. 10)
/// 4. Construct imposed states (Eq. 11-12)
///
/// Reuses FaultBasis for coordinate transforms and FrictionSolver for Eq. 8.
/// Reference: SeisSol FrictionSolverCommon.h, de la Puente et al. (2009).
class FaultFaceFlux
{
public:
   /// @brief Construct with material parameters.
   ///
   /// @param[in] rho  Density [kg/m³].
   /// @param[in] cp  P-wave speed [m/s].
   /// @param[in] cs  S-wave speed [m/s].
   FaultFaceFlux(real_t rho, real_t cp, real_t cs);

   /// @brief Compute trial traction from Q± states in fault-local coordinates.
   ///
   /// Implements Eq. (7a-c) from the plan:
   ///   σ_n^trial = η_p * (v_n⁻ - v_n⁺ + σ_n⁺/Z_p⁺ + σ_n⁻/Z_p⁻)
   ///   τ_1^trial = η_s * (v_t1⁻ - v_t1⁺ + τ_1⁺/Z_s⁺ + τ_1⁻/Z_s⁻)
   ///   τ_2^trial = η_s * (v_t2⁻ - v_t2⁺ + τ_2⁺/Z_s⁺ + τ_2⁻/Z_s⁻)
   ///
   /// @param[in] data  Per-DOF impedance data.
   /// @param[in] Q_plus  Rotated state on + side (9 components, fault-local).
   /// @param[in] Q_minus  Rotated state on − side (9 components, fault-local).
   /// @param[out] sigma_n_trial  Trial normal stress.
   /// @param[out] tau1_trial  Trial tangential traction (component 1).
   /// @param[out] tau2_trial  Trial tangential traction (component 2).
   static void ComputeTrialTraction(const DOFData &data,
                                    const real_t *Q_plus,
                                    const real_t *Q_minus,
                                    real_t &sigma_n_trial,
                                    real_t &tau1_trial,
                                    real_t &tau2_trial);

#ifdef SEAS_TEST_INTERNAL
   /// TEST-ONLY (compiled out unless `SEAS_TEST_INTERNAL`): inject a controlled
   /// absolute [Pa] cross-rank seed into the strike-channel trial traction,
   /// mimicking the ~1e-14 (relative) DG shared-face interpolation gap that
   /// desyncs two ranks at the LSW slip-onset kink (Phase-0 result, job
   /// 7747304).  Set NONZERO on exactly ONE rank; `ComputeTrialTraction` then
   /// adds this value to `tau2_trial` every call.  An ABSOLUTE Pa (not a ULP):
   /// at the deterministic prestress knife's edge the unperturbed trial is
   /// exactly 0, where a ULP nudge is a useless denormal — a small absolute
   /// nudge (~1e-6 Pa = ~1.6e-14 of a 60 MPa traction, the real seed scale)
   /// is the faithful, controllable stand-in.  Consumed only by
   /// `test_shared_fault_reconcile_cross_rank`.  Zero default => no-op.
   static real_t s_seas_test_tau2_trial_perturb_pa;

   /// TEST-ONLY (compiled out unless `SEAS_TEST_INTERNAL`): when true,
   /// `ComputeADERSharedFaceFluxRHS` SKIPS the Phase-2 cross-rank reconcile
   /// (the boss-broadcast) and assembles from each rank's own un-reconciled
   /// friction state.  Drives the NEGATIVE leg of
   /// `test_shared_fault_reconcile_cross_rank`, proving the R-101 guard still
   /// trips on a real desync (and that the reconcile is what fixes it).  Must
   /// be set IDENTICALLY on all ranks — the reconcile's MPI exchange is
   /// collective, so a per-rank disable would deadlock.  False default =>
   /// reconcile runs (production behaviour).
   static bool s_seas_test_disable_reconcile;
#endif

   /// @brief Full fault-face Riemann solver pipeline.
   ///
   /// Given Q± in fault-local coordinates, computes the imposed states
   /// Q^{±,imp} that satisfy friction and the characteristic compatibility
   /// relations. Returns the Godunov flux for Elem1 and Elem2.
   ///
   /// Pipeline: Eq. (7) → (8) → (9) → (10) → (11)-(12).
   ///
   /// v9.4.0 Commit 1: the friction-input total traction is
   ///   tau*_total = data.tau*_0 + data.tau*_nuc + tau*_trial
   /// where `tau*_trial` comes from bulk Q via ComputeTrialTraction,
   /// `data.tau*_0` carries the static background pre-stress (split-
   /// prestress fluctuation dispatch: bulk Q carries only the dynamic
   /// fluctuation), and `data.tau*_nuc` carries the time-varying
   /// nucleation driver overwritten per step by the driver.  The
   /// Riemann imposed state uses `tau*_corr = tau*_trial - eta_s·V*`
   /// (TRIAL scale) so the velocity jump only carries the friction
   /// reaction — the persistent channel does NOT radiate through
   /// bulk Q every step.  On user output `data.tau*_corr` is stored
   /// as TOTAL = tau*_0 + tau*_nuc + tau*_corr.
   ///
   /// For callers that leave `tau*_nuc` / `sigma_n_nuc` at their
   /// default 0 (all non-TPV102 callers, including BP5), behavior is
   /// byte-identical to the pre-v9.4.0 `Evaluate`.
   ///
   /// @param[in,out] data  Per-DOF data (psi, slip_rate updated).
   /// @param[in] Q_plus  State on + side (fault-local, 9 components).
   /// @param[in] Q_minus  State on − side (fault-local, 9 components).
   /// @param[out] Q_imp_plus  Imposed state on + side (9 components).
   /// @param[out] Q_imp_minus  Imposed state on − side (9 components).
   /// @param[in] method  Friction solver method.
   void Evaluate(DOFData &data,
                 const real_t *Q_plus, const real_t *Q_minus,
                 real_t *Q_imp_plus, real_t *Q_imp_minus,
                 FrictionSolver::Method method = FrictionSolver::Method::Brent) const;

   /// @brief Instantaneous LSW counterpart to `Evaluate` — the `dt→0` limit of
   /// `EvaluateADER_LSW` (no I-form, no `dt`): operates directly on the bulk
   /// states `Q_plus`/`Q_minus` (fault-local, NUM_STATE components).
   ///
   /// Reads ONLY the LSW-native fields `data.lsw_mu_s / lsw_mu_d / lsw_d_c /
   /// lsw_cohesion`; `data.a / psi / Dc` are NOT consumed.  SLIP-STATELESS: reads
   /// `data.slip1 / slip2` to form δ = √(slip1²+slip2²) for μ(δ), but NEVER writes
   /// them — the coupled-RK-on-slip stepper (`AdvanceRKCoupledLSW_Spatial`)
   /// integrates slip.  This mirrors `Evaluate` being ψ-stateless (the R-V92-H07
   /// invariant) and `EvaluateADER_LSW`'s R-001 slip-invariant.
   ///
   /// Writes `data.{slip_rate,V1,V2,tau1_corr,tau2_corr,sigma_n_corr}` (the
   /// tau*_corr / sigma_n_corr fields carry TOTAL physical traction per
   /// `WriteBackState`).  Closed-form physics is the SCEC TPV5 §7-11 solve via
   /// `LSWFrictionCoefficient_TPV205` + `SolveLSW_TPV205` (no Brent iteration).
   ///
   /// @param[in,out] data         Per-DOF state (slip_rate/V/tau*_corr updated).
   /// @param[in]  Q_plus          State on + side (fault-local, NUM_STATE).
   /// @param[in]  Q_minus         State on − side (fault-local, NUM_STATE).
   /// @param[out] Q_imp_plus      Imposed state on + side (NUM_STATE).
   /// @param[out] Q_imp_minus     Imposed state on − side (NUM_STATE).
   void EvaluateLSW(DOFData &data,
                    const real_t *Q_plus, const real_t *Q_minus,
                    real_t *Q_imp_plus, real_t *Q_imp_minus) const;

   /// Round-12 Patch 1 stage helpers — split the body of `Evaluate` into
   /// reusable pieces so face-averaging experiments in
   /// `wave_operator.inl` can substitute per-QP stage values with face
   /// averages without duplicating any physics.  Composing
   /// `ComputeStageState` + `BuildImposedState` + `WriteBackState` is
   /// byte-identical to calling `Evaluate` directly on the same inputs.
   ///
   /// Full baseline chain: trial → total → Θ → V_abs → V1/V2 → tau*_corr.
   /// Leaves `data` unchanged.
   void ComputeStageState(const DOFData &data,
                          const real_t *Q_plus, const real_t *Q_minus,
                          EvalStageState &s,
                          FrictionSolver::Method method =
                             FrictionSolver::Method::Brent) const;

   /// Completion helpers: assume a specific stage field in `s` is
   /// already overwritten (e.g. by face-averaging), then recompute
   /// every downstream stage.  Do not touch earlier stages.
   ///
   /// `CompleteFromTrial` — trial traction was modified; recompute
   /// total → Θ → V_abs → V1/V2 → corrected.
   void CompleteFromTrial(const DOFData &data,
                          EvalStageState &s,
                          FrictionSolver::Method method =
                             FrictionSolver::Method::Brent) const;

   /// `CompleteFromTheta` — Θ was modified (e.g. face-averaged);
   /// recompute V_abs → V1/V2 → corrected.  Assumes total traction
   /// (sigma_n_total, tau*_total) is already valid.
   void CompleteFromTheta(const DOFData &data,
                          EvalStageState &s,
                          FrictionSolver::Method method =
                             FrictionSolver::Method::Brent) const;

   /// `CompleteFromVabs` — V_abs was modified; recompute V1/V2 →
   /// corrected.  Assumes total traction and Θ are already valid.
   void CompleteFromVabs(const DOFData &data, EvalStageState &s) const;

   /// Construct the imposed states (Eq. 11-12) from a completed
   /// `EvalStageState`.  Pure function on `s` + per-side bulk Q.
   void BuildImposedState(const DOFData &data,
                          const EvalStageState &s,
                          const real_t *Q_plus, const real_t *Q_minus,
                          real_t *Q_imp_plus, real_t *Q_imp_minus) const;

   /// Write `slip_rate`, `V1`, `V2`, `tau1_corr`, `tau2_corr`, and
   /// `sigma_n_corr` onto `data` (same convention as `Evaluate`: the
   /// tau*_corr / sigma_n_corr fields carry TOTAL physical traction,
   /// i.e. pre-stress + nucleation + trial-scale corrected).
   void WriteBackState(DOFData &data, const EvalStageState &s) const;

   /// Total-stress variant of Evaluate (I-06 Phase 3).  Expects Q_plus,
   /// Q_minus to carry TOTAL stresses (pre-stress + fluctuation); outputs
   /// Q_imp_± also carry TOTAL stresses.  On a symmetric-pre-stress fault
   /// with homogeneous material, mathematically equivalent to
   /// Evaluate(data, Q_fluc_±, ...) + pre-stress shift (plan §3).
   ///
   /// DOFData.sigma_n_corr, tau1_corr, tau2_corr are stored as TOTAL.
   /// Under the v9.3.0 migration, DOFData.sigma_n0, tau1_0, tau2_0 are
   /// ZEROED by the driver (Phase 4) so that EvaluateTotal does not
   /// double-count pre-stress through ComputeTrialTraction on TOTAL Q.
   ///
   /// R-003 guards preserved: psi-invariance is asserted at entry/exit
   /// (the driver's coupled RK4-on-psi integrator relies on psi-purity),
   /// and the SEAS_DIAG_FAULT_FLUX diagnostic print block mirrors the
   /// fluctuation path so C-1 / C-2 / C-3 bisection still functions.
   ///
   /// @param[in,out] data  Per-DOF state.  sigma_n_corr, tau1_corr,
   ///                      tau2_corr stored as TOTAL values.
   /// @param[in]  Q_plus  Total-stress state on + side (fault-local).
   /// @param[in]  Q_minus Total-stress state on - side (fault-local).
   /// @param[out] Q_imp_plus  Total-stress imposed state on + side.
   /// @param[out] Q_imp_minus Total-stress imposed state on - side.
   /// @param[in]  method Friction solver choice.
   void EvaluateTotal(DOFData &data,
                      const real_t *Q_plus, const real_t *Q_minus,
                      real_t *Q_imp_plus, real_t *Q_imp_minus,
                      FrictionSolver::Method method = FrictionSolver::Method::Brent) const;

   /// ADER Phase 5: time-integrated friction solve (fluctuation-Q variant).
   ///
   /// Inputs I± = ∫_0^{dt} Q±(τ) dτ in fault-local coordinates.  Converts
   /// to time-averaged Q̄± = I±/dt, calls `Evaluate` once on the averaged
   /// state (the ADER one-shot replacement for the 4 per-stage RK4 calls),
   /// and rescales the imposed outputs back to time-integrated form:
   ///   I_imp± = dt · Q_imp±.
   ///
   /// After the call, `data.{slip_rate,V1,V2,tau1_corr,tau2_corr,sigma_n_corr}`
   /// carry the TIME-AVERAGED values over [t_n, t_n+dt], as specified by
   /// plan §Phase 5 §4.  `data.psi` is NOT updated — the caller must call
   /// `UpdateStateAnalytic` with the returned time-averaged slip rate.
   ///
   /// In the dt→0 limit, `EvaluateADER(data, I±, dt)` ≈
   /// `Evaluate(data, Q±(t_n + dt/2), dt) · dt` to O(dt²)
   /// (averaging-then-solving vs solving-at-midpoint for the nonlinear
   /// friction law).
   ///
   /// @param[in,out] data  Per-DOF state.  Slip-rate & traction fields
   ///                      updated to time-averaged values.
   /// @param[in]  I_plus   Time-integrated + side state (9 components,
   ///                      fault-local).
   /// @param[in]  I_minus  Time-integrated − side state.
   /// @param[in]  dt       Time step.  Must be > 0.
   /// @param[out] I_imp_plus   Time-integrated imposed + state.
   /// @param[out] I_imp_minus  Time-integrated imposed − state.
   /// @param[in]  method   Friction solver choice.
   void EvaluateADER(DOFData &data,
                     const real_t *I_plus, const real_t *I_minus,
                     real_t dt,
                     real_t *I_imp_plus, real_t *I_imp_minus,
                     FrictionSolver::Method method = FrictionSolver::Method::Brent) const;

   /// ADER Phase 5 + v9.3.0 §Phase 7 follow-up: time-integrated friction
   /// solve with TOTAL-stress inputs.  Mirrors `EvaluateADER` but wraps
   /// `EvaluateTotal` internally so the ADER corrector can drive the
   /// post-I-06 TPV102 dispatch path without the explicit 1/dt-then-dt
   /// workaround mentioned in the v9.3.0 plan.
   ///
   /// @param[in,out] data  Per-DOF state; slip/traction fields updated
   ///                      to time-averaged TOTAL values.
   /// @param[in]  I_plus_tot  Total-stress time-integrated + state.
   /// @param[in]  I_minus_tot Total-stress time-integrated − state.
   /// @param[in]  dt          Time step (> 0).
   /// @param[out] I_imp_plus_tot   Total-stress imposed + state.
   /// @param[out] I_imp_minus_tot  Total-stress imposed − state.
   /// @param[in]  method      Friction solver choice.
   void EvaluateADERTotal(DOFData &data,
                          const real_t *I_plus_tot, const real_t *I_minus_tot,
                          real_t dt,
                          real_t *I_imp_plus_tot, real_t *I_imp_minus_tot,
                          FrictionSolver::Method method = FrictionSolver::Method::Brent) const;

   /// @brief LSW counterpart to `EvaluateADER` — time-integrated Riemann
   /// solve for the linear slip-weakening law (SCEC TPV5 / TPV205).
   ///
   /// I-form wrap pattern identical to `EvaluateADER`:
   ///   Q̄± = I±/dt, run the LSW closed form on Q̄±, scale I_imp = dt·Q_imp.
   /// The friction physics is the SCEC TPV5 §7-11 closed form:
   ///   μ(δ) via `LSWFrictionCoefficient_TPV205`,
   ///   V_abs / V1 / V2 / τ*_corr via `SolveLSW_TPV205`,
   /// reading the LSW-native fields `data.lsw_mu_s / lsw_mu_d / lsw_d_c`.
   /// `data.a / data.psi / data.Dc` are NOT read.
   ///
   /// On return `data.{slip_rate, V1, V2, tau1_corr, tau2_corr,
   /// sigma_n_corr}` carry the values from the closed-form solve at the
   /// macro-step's time-averaged Q, matching `EvaluateADER`'s station-
   /// output convention; `data.tau*_corr` and `data.sigma_n_corr` are
   /// TOTAL (= pre-stress + nuc + trial-scale corrected) per
   /// `WriteBackState`. `data.psi` is NOT touched.
   ///
   /// R-001 (final review): `data.slip{1,2}` is NOT touched by this
   /// call — slip evolution is owned exclusively by
   /// `Tpv205SubStepIterator::StepOneQP_`, which integrates `slip{1,2}
   /// += V{1,2} * dt_sub` once per sub-step over every fault QP.
   /// HISTORY: under the retired R-1601 shared-fault fallback
   /// (superseded by unify-plan Phase 2; see debug_document/
   /// tpv104_debug_document/R1601_root_cause_2026-07-10.md) this
   /// function ran AFTER the iterator on shared QPs, so accumulating
   /// slip here double-counted at np > 1.  The rule REMAINS
   /// load-bearing on the one-shot (no-buffer) path and preserves the
   /// iterator's exclusive ownership of slip.  See REVIEW R-001.
   ///
   /// @param[in,out] data   Per-DOF state.
   /// @param[in]  I_plus    Time-integrated + side state (NUM_STATE).
   /// @param[in]  I_minus   Time-integrated − side state (NUM_STATE).
   /// @param[in]  dt        Time step (> 0).
   /// @param[out] I_imp_plus   Time-integrated imposed + state (NUM_STATE).
   /// @param[out] I_imp_minus  Time-integrated imposed − state (NUM_STATE).
   void EvaluateADER_LSW(DOFData &data,
                         const real_t *I_plus, const real_t *I_minus,
                         real_t dt,
                         real_t *I_imp_plus, real_t *I_imp_minus) const;

   /// @brief Phase H.6 of spatial_dynamic_rupture_plan.md (rev-3):
   /// LSW Riemann solve with the TPV26/27 time-dependent forced-
   /// rupture friction coefficient.
   ///
   /// Same as `EvaluateADER_LSW` but the friction coefficient is
   ///   μ(t, δ) = μ_s + (μ_d − μ_s) · max(f_1(δ), f_2(t, T, t_0))
   /// where `f_2` ramps from 0 → 1 over the interval
   /// `[T_forced_rupture, T_forced_rupture + t0_decay_forced]`.  Reads
   /// the two new `DOFData` fields (`T_forced_rupture`,
   /// `t0_decay_forced`).
   ///
   /// BYTE-EXACT CONTRACT: with the in-class `DOFData` defaults
   /// (`T_forced_rupture = 1.0e9`, `t0_decay_forced = 0`),
   /// `f_2(t) = 0` for any physically reachable simulation time
   /// (< 1.0e8 s), so the result reduces to plain LSW and equals
   /// `EvaluateADER_LSW` to the last bit.  This preserves the TPV205
   /// byte-exact gate: TPV205 stays on `FaultFrictionLaw::LSW` (the
   /// existing path), and a hypothetical mistaken routing of TPV205
   /// data through this method still produces TPV205-identical output.
   ///
   /// Slip-accumulation invariant from `EvaluateADER_LSW` applies
   /// verbatim: this method does NOT update `data.slip{1,2}`.
   ///
   /// @param[in,out] data    Per-DOF state.  Reads
   ///                        `T_forced_rupture` / `t0_decay_forced`;
   ///                        writes V{1,2}, tau{1,2}_corr, sigma_n_corr.
   /// @param[in]  I_plus     Time-integrated + side state (NUM_STATE).
   /// @param[in]  I_minus    Time-integrated − side state (NUM_STATE).
   /// @param[in]  dt         Time step (> 0).
   /// @param[in]  t_now      Macro-step simulation time at which to
   ///                        evaluate the f_2(t) factor.  Pass
   ///                        `wave.GetTime()` from the dispatch site.
   /// @param[out] I_imp_plus  Time-integrated imposed + state.
   /// @param[out] I_imp_minus Time-integrated imposed − state.
   void EvaluateADER_LSW_ForcedRupture(DOFData &data,
                                       const real_t *I_plus,
                                       const real_t *I_minus,
                                       real_t dt,
                                       real_t t_now,
                                       real_t *I_imp_plus,
                                       real_t *I_imp_minus) const;

   /// Access the friction solver.
   const FrictionSolver &GetSolver() const { return solver_; }

   /// Set the compressive normal-stress strength floor [Pa] (sliver-
   /// blowup plan 2026-05-26).  `v < 0` is the disabled sentinel: each
   /// friction law keeps its exact current strength expression (LSW
   /// `max(σ_n,0)`, RS `|σ_n|`) ⇒ byte-exact for the TPV/BP5
   /// regressions.  `v >= 0` floors the σ_n that enters the SHEAR
   /// STRENGTH (only) at `v`, so below `v` compression the strength is
   /// the constant `μ·v` rather than the spurious tensile `0` (LSW) or
   /// `|σ_n|` (RS).  The written-back `sigma_n_corr` channel is
   /// untouched.  Set once from config by the driver; default disabled.
   void SetSigmaNStrengthFloor(real_t v) { sigma_n_strength_floor_ = v; }

   /// Raw floor value (the disabled sentinel `< 0` is preserved).
   real_t GetSigmaNStrengthFloor() const { return sigma_n_strength_floor_; }

   /// Floor value for the LSW path, with the disabled sentinel `< 0`
   /// mapped to `0.0`.  Passing `0.0` to `SolveLSW_TPV205`'s
   /// `sigma_n_floor` argument reproduces the current LSW `max(σ_n,0)`
   /// behavior byte-exactly, so this is the value the LSW dispatch sites
   /// forward when the floor is disabled.
   real_t SigmaNStrengthFloorForLSW() const
   { return sigma_n_strength_floor_ >= 0.0 ? sigma_n_strength_floor_ : 0.0; }

   /// (Unified bi-material plan, Part B / B2) Affirm that the CALLER converts the
   /// per-side imposed state to the bulk flux with PER-SIDE A (i.e. applies A_plus
   /// to Q_imp_plus and A_minus to Q_imp_minus — the matrix `BimaterialWaveOperator`
   /// does this via `FluxForElem_`).  The fault-flux MATH (ComputeTrialTraction /
   /// BuildImposedState) is already per-side-correct for unequal impedance (verified:
   /// `seas_test_bimaterial_fault_riemann`, B0 verdict
   /// debug_document/tpv6_debug_document/bimaterial_fault_verification_2026-06-06.md);
   /// the ONLY thing the homogeneity guards protected was a single-A conversion site.
   /// Default FALSE => the guards still abort on a bimaterial face (the scalar
   /// `WaveOperator` never sets this), so byte-exact behavior is preserved.  Set TRUE
   /// ONLY by the matrix operator.
   void SetPerSideFluxApplied(bool v) { per_side_flux_applied_ = v; }
   bool GetPerSideFluxApplied() const { return per_side_flux_applied_; }

   /// (Part C, TPV6) Capture the per-side imposed (split-node) velocity into
   /// DOFData for the TPV6 per-side station writer.  Maps the fault-local
   /// velocity components of the imposed Godunov states into
   /// DOFData::v_imp_{plus,minus} = [normal(VX), dip(VY), strike(VZ)].  Called
   /// by the Evaluate variants after they build Q_imp_{plus,minus}; pure store,
   /// no effect on the flux (other writers never read these fields).
   static void StoreImposedVelocity_(DOFData &d,
                                     const real_t *Q_imp_plus,
                                     const real_t *Q_imp_minus)
   {
      d.v_imp_plus[0]  = Q_imp_plus[VX];  d.v_imp_plus[1]  = Q_imp_plus[VY];
      d.v_imp_plus[2]  = Q_imp_plus[VZ];
      d.v_imp_minus[0] = Q_imp_minus[VX]; d.v_imp_minus[1] = Q_imp_minus[VY];
      d.v_imp_minus[2] = Q_imp_minus[VZ];
   }

private:
   real_t rho_, cp_, cs_;
   real_t Zp_, Zs_;  ///< Impedances (homogeneous)
   FrictionSolver solver_;
   /// Compressive σ_n strength floor [Pa].  `< 0` ⇒ disabled (current
   /// behavior).  See SetSigmaNStrengthFloor.
   real_t sigma_n_strength_floor_ = -1.0;
   /// (Part B / B2) caller-applies-per-side-A affirmation; default false ⇒ the
   /// bimaterial-fault homogeneity guards abort.  See SetPerSideFluxApplied.
   bool per_side_flux_applied_ = false;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_FACE_FLUX_HPP
