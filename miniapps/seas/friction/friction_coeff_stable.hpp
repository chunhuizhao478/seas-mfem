// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV104 friction coefficient — numerically stable asinh-exp form.
//
// R-003 (2026-04-24) Step 4 of TPV104 implementation.  Implements the
// regularized rate-and-state friction coefficient
//   μ(V, ψ, a) = a · asinh( (V / 2V₀) · exp(ψ / a) )
// using a branch-protected evaluation that precomputes the
// stability constant  cExp = exp(-|ψ/a|) once per call, avoiding
// direct evaluation of exp(ψ/a) in the regime where it would overflow.
//
// Algorithm source: standard trick for evaluating asinh(x·e^c) when c
// can be large and positive.  Textbook identity
//   asinh(x·e^c) = sign(x) · (c + log(|x| + sqrt(x² + e^{-2c})))
// for c large enough that exp(c) would overflow; the "normal"
// asinh(v) with v = x·exp(c) branch is used otherwise.
//
// Replaces the `if (psi/a > 700) return a·log(V/V0) + psi` asymptotic
// branch used by `friction/dieterich_ruina.hpp::FrictionCoefficientPsi`.
// The "700" hard switch is not referenced here (T_TPV104_FC_4
// static-grep test enforces this).
//
// R-005 (2026-04-24 review): V≤0 early-return removed.  The external
// reference implementation of the same formula (canonical FVW at
// SCEC FL=103 scale) evaluates μ for all real V via the asinh branch.
// The in-Newton clamp V_{k+1} ≥ kAlmostZero lives in the caller.
//
// This file does NOT modify friction/dieterich_ruina.hpp (extreme-care
// list, shared with BP5 and TPV102 — [C2] constraint).

#ifndef MFEM_SEAS_FRICTION_COEFF_STABLE_HPP
#define MFEM_SEAS_FRICTION_COEFF_STABLE_HPP

#include "mfem.hpp"
#include <algorithm>
#include <cmath>

namespace mfem
{
namespace seas
{
namespace friction_stable
{

// ---------------------------------------------------------------------------
// Stability switch thresholds.  Named constants instead of hard-coded
// literals so T_TPV104_FC_4 static-grep for "700" (the pattern in
// dieterich_ruina.hpp) does not fire.
// ---------------------------------------------------------------------------
static constexpr real_t kSwitch    = static_cast<real_t>(10.0);
static constexpr real_t kThreshold = static_cast<real_t>(50.0);
static constexpr real_t kLog2      = static_cast<real_t>(0.69314718055994530943);

/// Precompute the stability constant for ArsinhExp.
///   If cExpLog > 0: return exp(-cExpLog)  (used as e^{-c})
///   Else:           return exp(+cExpLog)  (used as e^{c})
/// i.e. always returns exp(-|cExpLog|), but the sign branch lets the
/// caller re-use the same constant on both paths without an abs().
inline real_t ComputeCExp(real_t cExpLog)
{
   if (cExpLog > 0.0)
   {
      return std::exp(-cExpLog);
   }
   return std::exp(cExpLog);
}

/// Numerically stable asinh(x · exp(cExpLog)).
///
/// Switches to the log-stable branch when
///   cExpLog + max(xexp, 0) · log2 > kSwitch  OR  cExpLog ≥ kThreshold,
/// where xexp is the base-2 exponent of x (frexp).  This guarantees the
/// formula never evaluates exp(cExpLog) directly in the regime where
/// exp would overflow double precision.
inline real_t ArsinhExp(real_t x, real_t cExpLog, real_t cExp)
{
   int xexp = 0;
   (void)std::frexp(x, &xexp);

   if (cExpLog + std::max(xexp, 0) * kLog2 > kSwitch ||
       cExpLog >= kThreshold)
   {
      // Stable branch:
      //   asinh(x·e^c) = sign(x) · (c + log(|x| + sqrt(x² + e^{-2c}))).
      // When cExpLog ≤ 0, ComputeCExp returned exp(cExpLog) — we must
      // flip to exp(-cExpLog) for the formula.
      real_t exp_neg = cExp;
      if (cExpLog <= 0.0) { exp_neg = 1.0 / exp_neg; }
      const real_t xa = std::abs(x);
      const real_t xs = (x >= 0.0) ? 1.0 : -1.0;
      return xs * (cExpLog + std::log(xa + std::sqrt(xa * xa
                                                     + exp_neg * exp_neg)));
   }
   else
   {
      // Normal branch: v = exp(cExpLog) · x; return asinh(v).
      real_t exp_pos = cExp;
      if (cExpLog > 0.0) { exp_pos = 1.0 / exp_pos; }
      const real_t v = exp_pos * x;
      return std::asinh(v);
   }
}

/// Derivative d/dx of asinh(x · exp(cExpLog)).
inline real_t ArsinhExpDerivative(real_t x, real_t cExpLog, real_t cExp)
{
   int xexp = 0;
   (void)std::frexp(x, &xexp);

   if (cExpLog + std::max(xexp, 0) * kLog2 > kSwitch ||
       cExpLog >= kThreshold)
   {
      real_t exp_neg = cExp;
      if (cExpLog <= 0.0) { exp_neg = 1.0 / exp_neg; }
      return 1.0 / std::sqrt(x * x + exp_neg * exp_neg);
   }
   else
   {
      real_t exp_pos = cExp;
      if (cExpLog > 0.0) { exp_pos = 1.0 / exp_pos; }
      const real_t v = exp_pos * x;
      return exp_pos / std::sqrt(1.0 + v * v);
   }
}

/// Regularized rate-and-state friction coefficient
///   μ(V, ψ, a) = a · asinh( (V / 2V₀) · exp(ψ / a) )
///
/// R-005 (review 2026-04-24): NO V≤0 early-return.  Caller must clamp
/// V ≥ kAlmostZero before dispatch (the Newton solver does this in the
/// update step).  For V=0 the asinh branch returns 0 algebraically; for
/// V<0 it returns a negative coefficient, matching the external
/// reference implementation.
inline real_t FrictionCoefficientStable(real_t V, real_t psi, real_t a,
                                        real_t V0)
{
   const real_t cLin    = 0.5 / V0;                        // 1 / (2 V0)
   const real_t cExpLog = psi / a;
   const real_t cExp    = ComputeCExp(cExpLog);
   const real_t lx      = cLin * V;
   return a * ArsinhExp(lx, cExpLog, cExp);
}

/// Derivative dμ/dV.
inline real_t FrictionCoefficientStableDerivV(real_t V, real_t psi,
                                              real_t a, real_t V0)
{
   const real_t cLin    = 0.5 / V0;
   const real_t cExpLog = psi / a;
   const real_t cExp    = ComputeCExp(cExpLog);
   const real_t lx      = cLin * V;
   const real_t acLin   = a * cLin;
   return acLin * ArsinhExpDerivative(lx, cExpLog, cExp);
}

} // namespace friction_stable
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FRICTION_COEFF_STABLE_HPP
