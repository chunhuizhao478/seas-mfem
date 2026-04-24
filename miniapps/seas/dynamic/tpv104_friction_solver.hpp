// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV104 friction solver — Newton-Raphson on the regularized
// rate-and-state friction balance.
//
// R-004 (2026-04-24) Step 5 of TPV104 implementation.  Solves
//   g(V) = -(1/η_s) · ( |σ_n| · μ(V, ψ, a) − |τ| ) − V = 0
// via standard Newton-Raphson, using the Step-4
// `friction_stable::FrictionCoefficientStable` +
// `FrictionCoefficientStableDerivV` as μ(V, ψ) and dμ/dV, so the
// residual and gradient are bit-identical with the canonical FVW
// solver at every iterate.
//
// The existing `FrictionSolver::SolveNR` in `dynamic/friction_solver.cpp`
// uses the MFEM-native μ formula (from dieterich_ruina.hpp) and remains
// the CLI fallback (`--friction-solver=legacy-newton`).  This TPV104-
// specific variant is the production path; it keeps the BP5/TPV102
// Extreme-Care files untouched [C2].
//
// R-006 (review 2026-04-24): The initial-guess V = V_prev is NOT clamped
// on entry.  The in-loop update  V_{k+1} = max(kAlmostZero, V_k - g/g')
// handles the kAlmostZero floor per iterate.  Clamping the first guess
// would drop V_prev < kAlmostZero samples onto a different iterate
// trajectory than the canonical reference Newton.
//
// **WARNING**: Do NOT use this solver for BP5.  BP5 requires Brent
// (CLAUDE.md §"Friction Solver").  This file is TPV104-only.

#ifndef MFEM_SEAS_TPV104_FRICTION_SOLVER_HPP
#define MFEM_SEAS_TPV104_FRICTION_SOLVER_HPP

#include "mfem.hpp"
#include "../friction/friction_coeff_stable.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace mfem
{
namespace seas
{

// Sub-normal floor used by the external reference implementation; chosen
// empirically at 1e-45 for double precision.
inline constexpr real_t kTpv104AlmostZero = static_cast<real_t>(1e-45);

/// Newton-Raphson slip-rate solver for TPV104.
///
/// Solves  g(V) = -(1/η_s) · ( |σ_n| · μ(V, ψ, a) − |τ| ) − V = 0.
/// Converged when |g| < tol.  Default tol = 1e-8 matches the reference
/// implementation's newtonTolerance; plan §4.10 Step 5 allows tightening
/// to 1e-10 for tighter root accuracy.
///
/// Initial guess: V_prev (unclamped — see R-006).
/// Per-iterate update: V_{k+1} = max(kTpv104AlmostZero, V_k - g/g').
///
/// @param[in] tau_abs  |τ| (positive).
/// @param[in] psi      ψ (logarithmic state variable).
/// @param[in] sigma_n  |σ_n| (positive compression — TPV104 convention).
/// @param[in] eta_s    η_s = ρ·c_s / 2.  Must be > 0.
/// @param[in] a        Direct-effect parameter at this QP.
/// @param[in] V0       Reference slip rate (= 1e-6 m/s).
/// @param[in] V_prev   Previous-step V (the first-guess seed).
/// @param[in] max_iter Maximum Newton iterations.
/// @param[in] tol      Absolute tolerance on g.
/// @param[out] iterations   Number of iterations used (if non-null).
/// @param[out] has_converged True if |g| < tol at return (if non-null).
/// @return Slip-rate magnitude V.
inline real_t SolveSlipRateNewtonStable(
   real_t tau_abs, real_t psi, real_t sigma_n, real_t eta_s,
   real_t a, real_t V0,
   real_t V_prev,
   int max_iter = 60,
   real_t tol   = 1e-8,
   int *iterations = nullptr,
   bool *has_converged = nullptr)
{
   // R2-004 + R2-006 (review round 2): validate all inputs.  Silent
   // propagation of NaN or negative V_prev through the Newton loop would
   // corrupt μ (which is odd in V after R-005) and yield NaN without a
   // loud failure.  Throw std::runtime_error — portable across MFEM
   // builds without MFEM_USE_EXCEPTIONS.
   if (!std::isfinite(tau_abs) || tau_abs < 0.0)
   {
      throw std::runtime_error(
         "SolveSlipRateNewtonStable: tau_abs must be finite and "
         "non-negative; got " + std::to_string(tau_abs));
   }
   if (!std::isfinite(sigma_n))
   {
      throw std::runtime_error(
         "SolveSlipRateNewtonStable: sigma_n must be finite; got "
         + std::to_string(sigma_n));
   }
   if (!std::isfinite(eta_s) || eta_s <= 0.0)
   {
      throw std::runtime_error(
         "SolveSlipRateNewtonStable: eta_s must be finite and positive; "
         "got " + std::to_string(eta_s));
   }
   if (!std::isfinite(a) || a <= 0.0)
   {
      throw std::runtime_error(
         "SolveSlipRateNewtonStable: a must be finite and positive; "
         "got " + std::to_string(a));
   }
   if (!std::isfinite(V0) || V0 <= 0.0)
   {
      throw std::runtime_error(
         "SolveSlipRateNewtonStable: V0 must be finite and positive; "
         "got " + std::to_string(V0));
   }
   if (!std::isfinite(psi))
   {
      throw std::runtime_error(
         "SolveSlipRateNewtonStable: psi must be finite; got "
         + std::to_string(psi));
   }
   if (!std::isfinite(V_prev) || V_prev < 0.0)
   {
      throw std::runtime_error(
         "SolveSlipRateNewtonStable: V_prev must be finite and "
         "non-negative; got " + std::to_string(V_prev));
   }

   if (tau_abs <= 0.0)
   {
      if (iterations)    { *iterations = 0; }
      if (has_converged) { *has_converged = true; }
      return 0.0;
   }
   if (sigma_n <= 0.0)
   {
      // Fault in tension — no friction; V = tau / eta_s.
      if (iterations)    { *iterations = 0; }
      if (has_converged) { *has_converged = true; }
      return (eta_s > 0.0) ? (tau_abs / eta_s) : 0.0;
   }
   const real_t inv_eta_s = 1.0 / eta_s;
   const real_t sigma_n_abs = std::abs(sigma_n);

   // R-006: first guess unclamped — match reference.
   real_t V = V_prev;

   int it = 0;
   real_t g = 0.0;
   for (it = 0; it < max_iter; ++it)
   {
      const real_t mu = friction_stable::FrictionCoefficientStable(
                           V, psi, a, V0);
      g = -inv_eta_s * (sigma_n_abs * mu - tau_abs) - V;

      if (std::abs(g) < tol)
      {
         if (iterations)    { *iterations = it + 1; }
         if (has_converged) { *has_converged = true; }
         return V;
      }

      const real_t dmu = friction_stable::FrictionCoefficientStableDerivV(
                            V, psi, a, V0);
      const real_t dg  = -inv_eta_s * (sigma_n_abs * dmu) - 1.0;

      const real_t step = g / dg;
      V = std::max(kTpv104AlmostZero, V - step);
   }

   if (iterations)    { *iterations = max_iter; }
   if (has_converged) { *has_converged = std::abs(g) < tol; }
   return V;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV104_FRICTION_SOLVER_HPP
