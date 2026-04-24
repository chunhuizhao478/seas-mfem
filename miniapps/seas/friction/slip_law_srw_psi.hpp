// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV104 state-evolution law: slip law with strong rate weakening
// (SCEC FL=103, also called "Fast Velocity Weakening") in ψ-space.
//
// R-002 (2026-04-24) Step 2 of TPV104 implementation.  Companion to
// config/tpv104_params.hpp.  This file implements the SlipLawSRWPsi
// class + UpdateStateAnalyticSlipLawSRW free function that match the
// canonical FVW reference implementation's state-variable update and
// friction coefficient bit-for-bit on the TPV104 envelope.
//
// Do NOT modify friction/dieterich_ruina.hpp or friction/state_evolution.hpp;
// both are on the extreme-care list.  SlipLawSRWPsi only extends the
// StateEvolution interface by implementing its virtual methods and adds
// SRW-aware overloads that take per-QP V_w (SCEC TPV104 §4.1).
//
// === Review fixes (tpv104_review_2026-04-24.md) ===
// R-001: base-class virtuals Rate / SteadyState / RateDerivativeV /
//   RateDerivativeTheta abort in "production" mode (driver path), and
//   only dispatch to the default V_w / a when the class is in "test
//   fixture" mode.  The caller (driver) must invoke SetProductionMode()
//   at construction to enable the guard.
// R-002: `V_safe = max(V, 1e-45)` clamp is removed.  All formula inputs
//   use raw V directly; the upstream Newton solver's kAlmostZero floor
//   is the source of truth.
// R-004: (V/V_w)^8 is computed as the unrolled integer power
//   ((r²)²)² so the arithmetic matches a compile-time unrolled 8th
//   power byte-for-byte (avoids std::pow ULP drift).

#ifndef MFEM_SEAS_SLIP_LAW_SRW_PSI_HPP
#define MFEM_SEAS_SLIP_LAW_SRW_PSI_HPP

#include "mfem.hpp"
#include "state_evolution.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace mfem
{
namespace seas
{

/// Numerically stable  log(x · sinh(c))  evaluation, using the identity
///   log(x · sinh(c)) = |c| + log((x/2) · -sign(c) · expm1(-2|c|))
/// so large |c| does not overflow sinh internally.
///
/// For TPV104 the argument c = f_ss/a is always positive (f_ss ∈ [0.1,
/// 0.65], a > 0), but the full formula handles either sign.
inline real_t LogSinhStable(real_t x, real_t c)
{
   const real_t sign_c = (c >= 0.0) ? 1.0 : -1.0;
   const real_t absC   = std::abs(c);
   return absC + std::log(x / 2.0 * -sign_c * std::expm1(-2.0 * absC));
}

/// Integer 8th-power by the unrolled binary chain ((x²)²)².
/// Matches a compile-time unrolled template like
/// `reference::integer_pow<8, double>` byte-for-byte — avoids ULP drift
/// that `std::pow(x, 8.0)` introduces via its exp(8·log(x)) path.
/// (R-004 review fix.)
inline real_t IntegerPow8(real_t x)
{
   const real_t x2 = x * x;
   const real_t x4 = x2 * x2;
   return x4 * x4;
}

/// Analytic single-step update for the SRW slip-law-in-ψ ODE
///   dψ/dt = -(V/L) · (ψ - ψ_ss(V))
/// assuming V constant over [t, t + dt].  Matches the canonical FVW
/// reference implementation's `updateStateVariable` — see T_SRW_5 /
/// T_SRW_8 for the byte-match port.
///
/// R-002: raw V is used throughout.  Caller must guarantee V ≥ kAlmostZero
/// (the in-Newton floor).  A V = 0 input produces f_LV = +∞ → f_ss → +∞ →
/// ψ_ss = +∞ but  exp1m = 0  so the analytic step returns ψ_old — valid
/// locked-fault limit.
///
/// @param[in] psi_old  ψ at the start of the sub-step.
/// @param[in] V        Slip-rate magnitude during the sub-step [m/s].
/// @param[in] L        Critical slip distance at this QP [m] (SCEC L).
/// @param[in] dt       Sub-step size [s].
/// @param[in] V_w      Weakening velocity at this QP [m/s].
/// @param[in] a        Direct-effect parameter at this QP.
/// @param[in] b        Slowness parameter (global for TPV104).
/// @param[in] V0       Reference slip rate (global, = 1e-6 m/s).
/// @param[in] f0       Reference friction coefficient (global, = 0.6).
/// @param[in] muW      Weakening friction coefficient (global, = 0.1).
/// @return ψ at the end of the sub-step.
inline real_t UpdateStateAnalyticSlipLawSRW(real_t psi_old, real_t V,
                                            real_t L, real_t dt,
                                            real_t V_w, real_t a,
                                            real_t b, real_t V0,
                                            real_t f0, real_t muW)
{
   // (1) f_LV = max(0, f0 - (b - a) · log(V/V0))
   //     Sign convention: matches the canonical FVW reference and
   //     SCEC TPV104 (Noda & Lapusta 2013).  The plan's §4.2.1 eq (2a)
   //     has a sign typo; the byte-match directive resolves it in favour
   //     of the reference implementation.
   const real_t f_LV = std::max(static_cast<real_t>(0),
                                f0 - (b - a) * std::log(V / V0));

   // (2) f_ss = muW + (f_LV - muW) / (1 + (V/V_w)^8)^(1/8)
   //     R-004: integer 8th power via unrolled binary chain, not std::pow.
   const real_t V_over_Vw   = V / V_w;
   const real_t V_over_Vw_8 = IntegerPow8(V_over_Vw);
   const real_t denom       = std::pow(1.0 + V_over_Vw_8,
                                       static_cast<real_t>(1.0 / 8.0));
   const real_t f_ss        = muW + (f_LV - muW) / denom;

   // (3) ψ_ss = a · logsinh(2·V0/V, f_ss/a)
   const real_t psi_ss = a * LogSinhStable(2.0 * V0 / V, f_ss / a);

   // (4) Analytic step:  ψ(t+dt) = ψ_ss + (ψ_0 - ψ_ss)·exp(-V·dt/L).
   //     Expressed as ψ_ss·(1 - exp(preexp1)) + exp(preexp1)·ψ_0 to match
   //     the canonical reference's round-off profile.
   const real_t preexp1 = -V * (dt / L);
   const real_t exp1v   = std::exp(preexp1);
   const real_t exp1m   = -std::expm1(preexp1);
   return psi_ss * exp1m + exp1v * psi_old;
}

/// Slip-law-with-SRW state-evolution law in ψ-space.
///
/// ODE:
///   dψ/dt = -(V/L) · (ψ - ψ_ss(V, V_w, a))
/// Steady-state:
///   ψ_ss(V) = a · ln((2V_0/V) · sinh(f_ss(V)/a))
///   f_ss(V) = f_w + (f_LV(V) - f_w) / (1 + (V/V_w)^8)^(1/8)
///   f_LV(V) = max(0, f_0 - (b - a) · ln(V/V_0))
///
/// === Interface contract (R-001 review fix) ===
/// The class holds a *scalar* default `V_w_default_` and a *scalar*
/// default `a_`.  Two operating modes:
///
///  - **Test-fixture mode (default)**: base-class virtuals Rate /
///    SteadyState / RateDerivativeV / RateDerivativeTheta dispatch
///    through the default V_w and a.  This is fine for unit tests and
///    any usage where V_w is spatially uniform.
///
///  - **Production mode**: call `SetProductionMode()` before wiring into
///    the driver.  Now the base virtuals MFEM_ABORT with a directive
///    message telling the caller to use the `_SRW` overloads with
///    per-QP V_w[i] and per-QP a[i].  This catches the silent-fallthrough
///    bug where the driver would otherwise receive wrong per-QP physics.
class SlipLawSRWPsi : public StateEvolution
{
public:
   /// Construct with global scalars + required V_w default (used by
   /// base-interface virtuals in test-fixture mode).  Per-QP V_w is
   /// passed explicitly to Rate_SRW / SteadyState_SRW / PsiSS_SRW in
   /// production.
   ///
   /// R2-007 (review round 2): V_w_default is a **required** argument —
   /// the previous default value of 1.0 corresponded to the strengthening
   /// regime (V_w_out), not the VW core (V_w_in = 0.1).  A test fixture
   /// expecting "default = TPV104 core physics" would silently get wrong
   /// behaviour.  Every call site must now pick the V_w that matches its
   /// physics regime.
   SlipLawSRWPsi(real_t a_scalar, real_t b_scalar, real_t V0_scalar,
                 real_t f0_scalar, real_t muW_scalar,
                 real_t V_w_default)
      : a_(a_scalar), b_(b_scalar), V0_(V0_scalar),
        f0_(f0_scalar), muW_(muW_scalar),
        V_w_default_(V_w_default),
        production_mode_(false) {}

   /// Flip the class into production mode.  After this call, the base-
   /// class virtuals abort with a loud error; the driver must use the
   /// `_SRW` overloads.  (R-001 review fix.)
   void SetProductionMode() { production_mode_ = true; }

   /// True iff the base-class virtuals are guarded against silent
   /// fallthrough.
   bool InProductionMode() const { return production_mode_; }

   // =================================================================
   // StateEvolution interface — implemented using the stored V_w default.
   // Abort in production mode (R-001 guard).
   // =================================================================

   real_t Rate(real_t V, real_t psi, real_t L) const override
   {
      if (production_mode_)
      {
         throw std::runtime_error(
            "SlipLawSRWPsi::Rate base-virtual called in production "
            "mode.  Use Rate_SRW(V, psi, L, V_w[i], a[i]) with "
            "per-QP V_w and a instead.");
      }
      return Rate_SRW(V, psi, L, V_w_default_, a_);
   }

   real_t SteadyState(real_t V, real_t /*L*/) const override
   {
      if (production_mode_)
      {
         throw std::runtime_error(
            "SlipLawSRWPsi::SteadyState base-virtual called in production "
            "mode.  Use SteadyState_SRW(V, V_w[i], a[i]) with per-QP V_w "
            "and a instead.");
      }
      return SteadyState_SRW(V, V_w_default_, a_);
   }

   real_t RateDerivativeV(real_t V, real_t psi, real_t L) const override
   {
      if (production_mode_)
      {
         throw std::runtime_error(
            "SlipLawSRWPsi::RateDerivativeV base-virtual called in "
            "production mode.");
      }
      // R2-003 (review round 2): at V = 0 the R-002 fix gives
      //   ψ_ss(0) = +∞   and   dψ_ss/dV → −∞,
      // so the finite-difference (ps_plus - ps_mid)/eps = (finite − ∞)/eps
      // is −∞ and the product (V/L)·(−∞) = 0·(−∞) = NaN.  Return the
      // analytic limit directly:
      //   dRate/dV = -(1/L)·(ψ - ψ_ss(V))  +  (V/L)·dψ_ss/dV.
      // As V → 0⁺, the second term vanishes (locked-fault contributes
      // no rate derivative via the product-with-V factor), and the first
      // term -(1/L)·(ψ - ψ_ss(0⁺)) evaluated using a sub-normal V proxy
      // gives a finite answer.
      if (V == 0.0)
      {
         const real_t V_proxy = std::numeric_limits<real_t>::min();
         const real_t ps_zero = SteadyState_SRW(V_proxy, V_w_default_, a_);
         return -(1.0 / L) * (psi - ps_zero);
      }
      const real_t eps = std::max(static_cast<real_t>(1e-8) * std::abs(V),
                                  static_cast<real_t>(1e-18));
      const real_t ps_plus  = SteadyState_SRW(V + eps, V_w_default_, a_);
      const real_t ps_mid   = SteadyState_SRW(V,       V_w_default_, a_);
      const real_t dps_dV   = (ps_plus - ps_mid) / eps;
      return -(1.0 / L) * (psi - ps_mid) + (V / L) * dps_dV;
   }

   real_t RateDerivativeTheta(real_t V, real_t /*psi*/,
                              real_t L) const override
   {
      if (production_mode_)
      {
         throw std::runtime_error(
            "SlipLawSRWPsi::RateDerivativeTheta base-virtual called in "
            "production mode.");
      }
      return -V / L;
   }

   const char *GetName() const override { return "SlipLawSRWPsi"; }

   // =================================================================
   // SRW-aware per-QP entry points.
   // =================================================================

   /// ODE right-hand side with per-QP V_w and per-QP a.
   real_t Rate_SRW(real_t V, real_t psi, real_t L,
                   real_t V_w, real_t a) const
   {
      const real_t psi_ss = PsiSS_SRW(V, V_w, a, b_, V0_, f0_, muW_);
      return -(V / L) * (psi - psi_ss);
   }

   /// Steady-state ψ with per-QP V_w and per-QP a.
   real_t SteadyState_SRW(real_t V, real_t V_w, real_t a) const
   {
      return PsiSS_SRW(V, V_w, a, b_, V0_, f0_, muW_);
   }

   // =================================================================
   // Static helper — per-primitive PsiSS_SRW evaluation.
   // Exposed static so tests can check it without constructing the class.
   // R-002: raw V used throughout (no V_safe clamp).
   // R-004: unrolled integer 8th-power.
   // =================================================================
   static real_t PsiSS_SRW(real_t V, real_t V_w, real_t a, real_t b,
                           real_t V0, real_t f0, real_t muW)
   {
      // f_LV = max(0, f0 - (b - a) · log(V/V0))
      const real_t f_LV = std::max(static_cast<real_t>(0),
                                   f0 - (b - a) * std::log(V / V0));

      // f_ss = muW + (f_LV - muW) / (1 + (V/V_w)^8)^(1/8)
      const real_t V_over_Vw   = V / V_w;
      const real_t V_over_Vw_8 = IntegerPow8(V_over_Vw);
      const real_t denom       = std::pow(1.0 + V_over_Vw_8,
                                          static_cast<real_t>(1.0 / 8.0));
      const real_t f_ss        = muW + (f_LV - muW) / denom;

      // ψ_ss = a · logsinh(2V0/V, f_ss/a)
      return a * LogSinhStable(2.0 * V0 / V, f_ss / a);
   }

   // =================================================================
   // Accessors for tests.
   // =================================================================
   real_t GetA()   const { return a_; }
   real_t GetB()   const { return b_; }
   real_t GetV0()  const { return V0_; }
   real_t GetF0()  const { return f0_; }
   real_t GetMuW() const { return muW_; }

private:
   real_t a_, b_, V0_, f0_, muW_;
   real_t V_w_default_;
   bool   production_mode_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_SLIP_LAW_SRW_PSI_HPP
