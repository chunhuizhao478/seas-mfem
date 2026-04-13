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

#ifndef MFEM_SEAS_FRICTION_SOLVER_HPP
#define MFEM_SEAS_FRICTION_SOLVER_HPP

#include "mfem.hpp"
#include <cmath>
#include <functional>

namespace mfem
{
namespace seas
{

/// @brief Dual friction solver: Brent (CPU) and Newton-Raphson (GPU).
///
/// Solves the friction balance equation (plan Eq. 8):
///   Θ = |σ_n| * f(V̂, ψ) + η_s * V̂
/// where f(V, ψ) = a * asinh[(V / 2V₀) * exp(ψ/a)] is the regularized
/// rate-and-state friction coefficient.
///
/// Reference: Uphoff (2020) Eq. 4.57, SeisSol SlowVelocityWeakeningLaw.h.
class FrictionSolver
{
public:
   enum class Method { Brent, NewtonRaphson, HybridNRBisection };

   /// Reference slip rate V₀ [m/s].
   static constexpr real_t V0 = 1e-6;

   /// Solve for slip rate V̂ given trial traction Θ.
   ///
   /// @param[in] tau  Magnitude of total traction Θ = |τ_total|.
   /// @param[in] psi  State variable (logarithmic).
   /// @param[in] sigma_n  Effective normal stress |σ_n| (positive compression).
   /// @param[in] eta  S-wave impedance η_s = Z_s^+ Z_s^- / (Z_s^+ + Z_s^-).
   /// @param[in] a  Direct effect parameter.
   /// @param[in] method  Solver method (default: Brent for CPU).
   /// @return Slip rate V̂ [m/s].
   real_t Solve(real_t tau, real_t psi, real_t sigma_n,
                real_t eta, real_t a, Method method = Method::Brent) const;

   /// Brent's method (robust, guaranteed convergence). Plan Eq. 8.
   /// Bracket: [V_lo=0, V_hi=Θ/η_s]. R-002 fix: sign polarity.
   real_t SolveBrent(real_t tau, real_t psi, real_t sigma_n,
                     real_t eta, real_t a) const;

   /// Newton-Raphson (fast, for GPU). Plan Eq. 15.
   /// Uses analytical derivative df/dV̂ = a*C / sqrt(1 + (V̂*C)²).
   real_t SolveNR(real_t tau, real_t psi, real_t sigma_n,
                  real_t eta, real_t a) const;

   /// Hybrid NR+Bisection: try NR for 5 iterations, fall back to bisection.
   real_t SolveHybrid(real_t tau, real_t psi, real_t sigma_n,
                      real_t eta, real_t a) const;

   /// Evaluate the residual g(V) = |σ_n| * f(V, ψ) + η * V - Θ.
   /// (R-002 fix) g(0) = -Θ < 0, g(Θ/η) > 0 for valid bracket.
   static real_t Residual(real_t V, real_t tau, real_t psi,
                          real_t sigma_n, real_t eta, real_t a);

   /// Analytical derivative dg/dV = |σ_n| * df/dV + η.
   /// df/dV = a * C / sqrt(1 + (V*C)²), C = exp(ψ/a) / (2*V₀).
   static real_t ResidualDerivative(real_t V, real_t psi,
                                    real_t sigma_n, real_t eta, real_t a);

private:
   /// Brent's method root finder. Finds root of F in [a, b].
   /// Requires F(a) * F(b) < 0. Returns root.
   static real_t ZeroIn(real_t a, real_t b,
                        const std::function<real_t(real_t)> &F,
                        real_t tol = 0.0);
};

/// Aging law state variable update (plan Eq. 13).
/// Uses expm1() for numerical stability when V*dt/L is small.
///
/// ψ(t+Δt) = ψ(t) * exp(-V*Δt/L) + (L/V) * (1 - exp(-V*Δt/L))
///
/// @param[in] psi_old  Current state variable.
/// @param[in] V  Slip rate [m/s].
/// @param[in] Dc  Critical slip distance [m].
/// @param[in] dt  Time step [s].
/// @return Updated state variable.
inline real_t UpdateStateAnalytic(real_t psi_old, real_t V, real_t Dc, real_t dt)
{
   // For V→0: aging law dψ/dt = 1 - Vψ/L → dψ/dt ≈ 1, so ψ(t+dt) = ψ + dt.
   real_t x = V * dt / Dc;
   if (x < 1e-15) { return psi_old + dt; }
   real_t exp_neg_x = std::exp(-x);
   // Use expm1 for stability: 1 - exp(-x) = -expm1(-x)
   real_t one_minus_exp = -std::expm1(-x);
   return psi_old * exp_neg_x + (Dc / V) * one_minus_exp;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FRICTION_SOLVER_HPP
