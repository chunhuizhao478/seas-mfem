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

#include "friction_solver.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

namespace mfem
{
namespace seas
{

// ---------------------------------------------------------------------------
// Residual: g(V) = |sigma_n| * f(V, psi) + eta * V - Theta
// (R-002 fix): g(0) = 0 + 0 - Theta = -Theta < 0
//              g(Theta/eta) = |sigma_n| * f(...) + Theta - Theta > 0
// This is OPPOSITE polarity from the QD Brent solver.
// ---------------------------------------------------------------------------
real_t FrictionSolver::Residual(real_t V, real_t tau, real_t psi,
                                real_t sigma_n, real_t eta, real_t a)
{
   // f(V, psi) = a * asinh[(V / 2V0) * exp(psi / a)]
   real_t C = std::exp(psi / a) / (2.0 * V0);
   real_t f = a * std::asinh(V * C);
   return sigma_n * f + eta * V - tau;
}

// ---------------------------------------------------------------------------
// Residual derivative: dg/dV = |sigma_n| * df/dV + eta
// df/dV = a * C / sqrt(1 + (V*C)^2)   (plan Eq. 15)
// ---------------------------------------------------------------------------
real_t FrictionSolver::ResidualDerivative(real_t V, real_t psi,
                                          real_t sigma_n, real_t eta, real_t a)
{
   real_t C = std::exp(psi / a) / (2.0 * V0);
   real_t VC = V * C;
   real_t df_dV = a * C / std::sqrt(1.0 + VC * VC);
   return sigma_n * df_dV + eta;
}

// ---------------------------------------------------------------------------
// Dispatch
// ---------------------------------------------------------------------------
real_t FrictionSolver::Solve(real_t tau, real_t psi, real_t sigma_n,
                             real_t eta, real_t a, Method method) const
{
   if (tau <= 0.0) { return 0.0; }

   switch (method)
   {
      case Method::Brent:    return SolveBrent(tau, psi, sigma_n, eta, a);
      case Method::NewtonRaphson: return SolveNR(tau, psi, sigma_n, eta, a);
      case Method::HybridNRBisection: return SolveHybrid(tau, psi, sigma_n, eta, a);
      default: return SolveBrent(tau, psi, sigma_n, eta, a);
   }
}

// ---------------------------------------------------------------------------
// Brent's method: delegates to proven QD solver (log10-V space, Tandem-verified)
// ---------------------------------------------------------------------------
real_t FrictionSolver::SolveBrent(real_t tau, real_t psi, real_t sigma_n,
                                  real_t eta, real_t a) const
{
   if (tau <= 0.0) { return 0.0; }
   return qd_friction_.SolveSlipRatePsi(tau, psi, sigma_n, eta, a);
}

// ---------------------------------------------------------------------------
// Newton-Raphson: fast convergence, for GPU
// ---------------------------------------------------------------------------
real_t FrictionSolver::SolveNR(real_t tau, real_t psi, real_t sigma_n,
                               real_t eta, real_t a) const
{
   if (tau <= 0.0) { return 0.0; }

   // Initial guess using f estimate from psi: f ≈ a * psi/a = psi (crude)
   real_t f_est = std::max(0.1, std::min(std::abs(psi), 1.0));
   real_t V = tau / (sigma_n * f_est + eta);

   for (int iter = 0; iter < 60; iter++)
   {
      real_t g = Residual(V, tau, psi, sigma_n, eta, a);
      if (std::abs(g) < 1e-8 * tau) { return V; }  // relative tolerance

      real_t dg = ResidualDerivative(V, psi, sigma_n, eta, a);
      real_t dV = g / dg;
      V = std::max(1e-45, V - dV);  // floor to prevent negative V
   }

   return V;  // may not be converged for extreme psi/a
}

// ---------------------------------------------------------------------------
// Hybrid NR+Bisection: try NR 5 iterations, fall back to bisection
// ---------------------------------------------------------------------------
real_t FrictionSolver::SolveHybrid(real_t tau, real_t psi, real_t sigma_n,
                                   real_t eta, real_t a) const
{
   if (tau <= 0.0) { return 0.0; }

   // Phase 1: Newton-Raphson (fast, 5 iterations)
   real_t f_est = std::max(0.1, std::min(std::abs(psi), 1.0));
   real_t V = tau / (sigma_n * f_est + eta);
   bool converged = false;

   for (int iter = 0; iter < 5; iter++)
   {
      real_t g = Residual(V, tau, psi, sigma_n, eta, a);
      if (std::abs(g) < 1e-8 * tau) { converged = true; break; }

      real_t dg = ResidualDerivative(V, psi, sigma_n, eta, a);
      V = std::max(1e-45, V - g / dg);
   }

   if (converged) { return V; }

   // Phase 2: Brent fallback (guaranteed convergence for all parameter ranges)
   return SolveBrent(tau, psi, sigma_n, eta, a);
}

// ZeroIn removed — Brent now delegates to DieterichRuinaFriction::SolveSlipRatePsi
// which uses the proven QD log10-V Brent solver (Tandem-verified).

} // namespace seas
} // namespace mfem
