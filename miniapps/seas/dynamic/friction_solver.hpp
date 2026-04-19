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
#include "../friction/dieterich_ruina.hpp"
#include "../friction/state_evolution.hpp"
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
/// Brent solve delegates to the proven QD DieterichRuinaFriction solver
/// (log10-V space, verified against Tandem across 62 debug iterations).
/// NR and Hybrid solvers are new (for GPU performance).
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
   /// QD friction solver (Tandem-verified log10-V Brent).
   /// Only V0 is used by SolveSlipRatePsi (a is passed as argument).
   /// Other constants (f0, b, Dc) are unused by the Brent solve.
   DieterichRuinaFriction qd_friction_{
      DieterichRuinaFriction::Constants{V0, 0.6, 0.012, 0.02}};
};

// UpdateStateAnalytic now lives in friction/state_evolution.hpp (R-005 fix).
// It's available here via the #include "../friction/state_evolution.hpp" above.

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FRICTION_SOLVER_HPP
