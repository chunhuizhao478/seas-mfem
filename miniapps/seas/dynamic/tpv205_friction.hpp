// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV5 / TPV205 linear slip-weakening (LSW) friction — closed-form
// slip-rate solve.
//
// LSW friction law (SCEC TPV5 §7-11):
//
//   μ(δ) = μ_s − (μ_s − μ_d) · min(δ / d_c, 1)
//   τ_strength = μ(δ) · σ_n
//
// Once the trial traction has been built from the bulk Q via
// FaultFaceFlux::ComputeTrialTraction, the radiation-damping balance
//
//   τ_corr = τ_total − η_s · V
//   |τ_corr| ≤ τ_strength    (with equality when the fault is sliding)
//
// gives the closed-form magnitude
//
//   V_abs = max(0, (|τ_total| − τ_strength) / η_s)
//
// and the direction of V is parallel to (τ1_total, τ2_total).
//
// This is fundamentally simpler than rate-and-state: there is no Newton
// iteration, no state variable evolution, no asinh/exp.  The Riemann
// imposed-state construction (FaultFaceFlux::BuildImposedState) is reused
// unchanged because it consumes (sigma_n_corr, tau1_corr, tau2_corr) in
// TRIAL scale, which we set per the same convention as the rate-and-state
// path.

#ifndef MFEM_SEAS_TPV205_FRICTION_HPP
#define MFEM_SEAS_TPV205_FRICTION_HPP

#include "mfem.hpp"
#include "../config/tpv205_params.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>  // std::getenv for the SEAS_NOOPENING diagnostic cap

namespace mfem
{
namespace seas
{

/// Effective LSW friction coefficient at slip magnitude `delta`.
///
/// μ(δ) = μ_s − (μ_s − μ_d) · min(δ/d_c, 1)
///   - δ ≤ 0       → μ_s
///   - 0 < δ < d_c → linear interpolation
///   - δ ≥ d_c     → μ_d
///
/// `mu_s` may be a very large finite value (e.g. 10000) to encode the
/// strength barrier outside the 30 km × 15 km rupture area; the function
/// is monotonically non-increasing in δ.
inline real_t LSWFrictionCoefficient_TPV205(real_t delta,
                                            real_t mu_s, real_t mu_d,
                                            real_t d_c)
{
   // R-002: SCEC TPV5 §11 strength barrier — `mu_s = mu_s_barrier`
   // (10000) is a sentinel meaning "rupture forbidden here for the
   // entire run".  Do NOT slip-weaken from this value; otherwise a QP
   // that ever accumulates δ ≥ d_c would silently collapse to μ_d
   // and let the rupture front escape the 30 km × 15 km rupture area.
   // Use a generous threshold so any future bump to the sentinel
   // (e.g. 1e5) still trips this branch.
   if (mu_s >= 0.5 * TPV205Params::mu_s_barrier) { return mu_s; }
   if (delta <= 0.0) { return mu_s; }
   if (delta >= d_c) { return mu_d; }
   return mu_s - (mu_s - mu_d) * (delta / d_c);
}

/// Solve the linear slip-weakening radiation-damping balance.
///
/// Inputs:
///   tau1_total, tau2_total : total tangential traction components [Pa]
///   sigma_n_total          : total normal stress [Pa] (positive = compression)
///   eta_s                  : radiation damping = ρ·c_s / 2 [Pa·s/m]
///   mu_eff                 : LSW μ(δ) at this QP (already evaluated)
///
/// Outputs:
///   V_abs : non-negative slip-rate magnitude [m/s]
///   V1, V2: slip-rate components (parallel to (τ1_total, τ2_total))
///   tau1_corr, tau2_corr : trial-scale corrected traction
///                          (`τ_corr = τ_trial − η_s·V`; the caller
///                           supplies tau*_trial and adds back the
///                           static pre-stress on writeback per the
///                           FaultFaceFlux convention)
///
/// All inputs and outputs are in fault-local coordinates (component 1 =
/// dip, component 2 = strike per the BP5 / Tandem canonical frame).
///
/// For LSW, σ_n is unaffected by friction: σ_n_corr = σ_n_trial (caller
/// sets that directly when filling EvalStageState).  This function does
/// NOT touch the normal-stress channel.
inline void SolveLSW_TPV205(real_t tau1_trial, real_t tau2_trial,
                            real_t tau1_total, real_t tau2_total,
                            real_t sigma_n_total, real_t eta_s,
                            real_t mu_eff,
                            real_t &V_abs, real_t &V1, real_t &V2,
                            real_t &tau1_corr, real_t &tau2_corr)
{
   MFEM_ASSERT(eta_s > 0.0,
               "SolveLSW_TPV205: eta_s must be positive; got " << eta_s);

   // Compute |τ_total|.  The slip-rate direction is parallel to the
   // total tangential traction.  TPV205 is pure strike-slip in steady
   // state, but during rupture both components can be non-zero.
   const real_t tau_abs = std::sqrt(tau1_total * tau1_total
                                    + tau2_total * tau2_total);

   // R-003: SCEC TPV5 §11 strength barrier is an ABSOLUTE LOCK — the
   // rock outside the 30 km × 15 km rupture area cannot rupture there
   // for the duration of the simulation, regardless of σ_n sign or
   // magnitude (including σ_n = 0 exactly, where both `max(σ_n, 0)`
   // and `|σ_n|` evaluate to 0 and would otherwise let τ_strength
   // collapse to 0).  Short-circuit V = 0 here so no path through the
   // closed-form solve can seed slip in the barrier.
   if (mu_eff >= 0.5 * TPV205Params::mu_s_barrier)
   {
      V_abs = 0.0;
      V1 = 0.0;
      V2 = 0.0;
      // No friction reaction — corrected traction equals trial traction
      // (the imposed-state Riemann then carries the trial values
      // unchanged, consistent with the velocity-discontinuity-zero
      // boundary condition on a locked face).
      tau1_corr = tau1_trial;
      tau2_corr = tau2_trial;
      return;
   }

   // DIAGNOSTIC no-opening cap (env SEAS_NOOPENING, default OFF).  Under
   // tension (sigma_n_total <= 0) the standard LSW law below free-slides
   // (sigma_n_pos clamps to 0 => tau_strength = 0 => V = |tau|/eta_s), which
   // — when the trial normal channel collapses tensile — drives the
   // unbounded slip runaway seen in the SAFS [SLIP] trace.  This env-gated
   // branch zeros the frictional slip under tension (the physical "fault
   // opens and decouples" limit) to TEST whether that free-slide path is the
   // operative amplifier (REVIEW_speckle_tension_analysis test 1).  Mirrors
   // the strength-barrier short-circuit above.  Default OFF => the branch is
   // never taken => byte-exact for the TPV205 regression.  Parsed once.
   static const bool s_no_opening = [] {
      const char *e = std::getenv("SEAS_NOOPENING");
      return e && e[0] && e[0] != '0'; }();
   if (s_no_opening && sigma_n_total <= 0.0)
   {
      V_abs = 0.0;
      V1 = 0.0;
      V2 = 0.0;
      tau1_corr = tau1_trial;
      tau2_corr = tau2_trial;
      return;
   }

   // Inside the rupture area the standard "fault opens under tension"
   // semantic applies — σ_n_pos clamps to 0 so a tensile transient
   // drops τ_strength to 0 and the fault freely slides.
   const real_t sigma_n_pos = std::max<real_t>(sigma_n_total, 0.0);
   const real_t tau_strength = mu_eff * sigma_n_pos;

   // Closed-form V_abs from radiation damping balance.
   if (tau_abs > tau_strength)
   {
      V_abs = (tau_abs - tau_strength) / eta_s;
   }
   else
   {
      V_abs = 0.0;
   }

   // Decompose V into (V1, V2) along (τ1_total, τ2_total).
   if (V_abs > 0.0 && tau_abs > 0.0)
   {
      const real_t inv_tau = 1.0 / tau_abs;
      V1 = V_abs * tau1_total * inv_tau;
      V2 = V_abs * tau2_total * inv_tau;
   }
   else
   {
      V1 = 0.0;
      V2 = 0.0;
   }

   // Trial-scale corrected traction (FaultFaceFlux::BuildImposedState
   // convention): τ_corr = τ_trial − η_s · V.  This drops the static
   // pre-stress and the (zero in TPV205) nucleation channel; on
   // writeback the driver/iterator adds them back to produce the
   // TOTAL physical traction for the station file.
   tau1_corr = tau1_trial - eta_s * V1;
   tau2_corr = tau2_trial - eta_s * V2;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV205_FRICTION_HPP
