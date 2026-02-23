// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-443211.
// All rights reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_DIETERICH_RUINA_HPP
#define MFEM_SEAS_DIETERICH_RUINA_HPP

#include "mfem.hpp"
#include "friction_law.hpp"
#include <cmath>
#include <limits>

namespace mfem
{
namespace seas
{

/// Dieterich-Ruina rate-and-state friction law (regularized form).
///
/// The regularized friction coefficient is:
///   f(V, theta) = a * asinh[(V / 2V0) * exp((f0 + b * ln(V0 * theta / Dc)) / a)]
///
/// This formulation:
/// - Regularizes the logarithmic form to handle V -> 0
/// - Provides a smooth friction coefficient for all positive V
/// - Reduces to the classical form for V >> V0
///
/// References:
/// - Rice, J.R., Lapusta, N., Ranjith, K. (2001). J. Mech. Phys. Solids.
/// - SCEC SEAS Benchmarks: https://strike.scec.org/cvws/seas/
class DieterichRuinaFriction : public FrictionLaw
{
public:
   /// Global constants for the friction law.
   struct Constants
   {
      real_t V0 = 1.0e-6;   ///< Reference slip rate [m/s]
      real_t f0 = 0.6;      ///< Reference friction coefficient [-]
      real_t b = 0.015;     ///< State evolution parameter [-]
      real_t Dc = 0.004;    ///< Critical slip distance [m] (BP2 value)
   };

   /// Constructor with constants.
   explicit DieterichRuinaFriction(const Constants &c) : cp_(c) {}

   /// Default constructor with BP2 default values.
   DieterichRuinaFriction()
      : DieterichRuinaFriction(Constants{1.0e-6, 0.6, 0.015, 0.004}) {}

   // =========================================================================
   // FrictionLaw interface implementation
   // =========================================================================

   /// Compute the friction coefficient f(V, theta).
   ///
   /// f(V, theta) = a * asinh[(V / 2V0) * exp((f0 + b * ln(V0 * theta / Dc)) / a)]
   ///
   /// @param[in] V Slip rate [m/s]
   /// @param[in] theta State variable [s]
   /// @param[in] a Rate-and-state direct effect parameter
   /// @return Friction coefficient [-]
   real_t FrictionCoefficient(real_t V, real_t theta, real_t a) const override
   {
      // Ensure positive values to avoid numerical issues
      V = std::max(V, V_min_);
      theta = std::max(theta, theta_min_);

      // f = a * asinh[(V / 2V0) * exp((f0 + b * ln(V0 * theta / Dc)) / a)]
      real_t log_arg = cp_.V0 * theta / cp_.Dc;
      real_t exp_arg = (cp_.f0 + cp_.b * std::log(log_arg)) / a;
      real_t sinh_arg = (V / (2.0 * cp_.V0)) * std::exp(exp_arg);

      return a * std::asinh(sinh_arg);
   }

   /// Compute df/dV for Newton solver.
   ///
   /// Using the chain rule:
   /// df/dV = a * (1 / sqrt(1 + arg^2)) * d(arg)/dV
   ///       = a * (1 / sqrt(1 + arg^2)) * (1 / 2V0) * exp(...)
   ///       = (1 / sqrt(1 + arg^2)) * (a / 2V0) * exp(...)
   ///
   /// Simplified using: sqrt(1 + sinh^2(x)) = cosh(x)
   /// For y = asinh(x), we have x = sinh(y) and dy/dx = 1/cosh(y)
   /// Since cosh(asinh(x)) = sqrt(1 + x^2):
   /// df/dV = a / (V * sqrt(1 + (2V0/V)^2 * exp(-2*(...)/a)))
   ///
   /// More directly: df/dV = (f/V) * (1 / (1 + exp(-2*f/a)))
   /// But let's use the straightforward form:
   real_t FrictionDerivativeV(real_t V, real_t theta, real_t a) const override
   {
      V = std::max(V, V_min_);
      theta = std::max(theta, theta_min_);

      real_t log_arg = cp_.V0 * theta / cp_.Dc;
      real_t exp_arg = (cp_.f0 + cp_.b * std::log(log_arg)) / a;
      real_t exp_val = std::exp(exp_arg);
      real_t sinh_arg = (V / (2.0 * cp_.V0)) * exp_val;

      // d(asinh(x))/dx = 1/sqrt(1 + x^2)
      // d(arg)/dV = exp_val / (2 * V0)
      real_t d_arg_dV = exp_val / (2.0 * cp_.V0);
      real_t d_asinh = 1.0 / std::sqrt(1.0 + sinh_arg * sinh_arg);

      return a * d_asinh * d_arg_dV;
   }

   /// Compute df/dtheta for implicit time stepping.
   ///
   /// df/dtheta = a * (1/sqrt(1 + arg^2)) * d(arg)/dtheta
   /// where d(arg)/dtheta = arg * (b / a) * (1 / theta)
   real_t FrictionDerivativeTheta(real_t V, real_t theta, real_t a) const override
   {
      V = std::max(V, V_min_);
      theta = std::max(theta, theta_min_);

      real_t log_arg = cp_.V0 * theta / cp_.Dc;
      real_t exp_arg = (cp_.f0 + cp_.b * std::log(log_arg)) / a;
      real_t exp_val = std::exp(exp_arg);
      real_t sinh_arg = (V / (2.0 * cp_.V0)) * exp_val;

      // d(arg)/dtheta = arg * (b / a) / theta
      real_t d_arg_dtheta = sinh_arg * (cp_.b / a) / theta;
      real_t d_asinh = 1.0 / std::sqrt(1.0 + sinh_arg * sinh_arg);

      return a * d_asinh * d_arg_dtheta;
   }

   /// Solve for slip rate V given stress tau and state theta.
   ///
   /// Solves: tau = sigma_n * f(V, theta) + eta * V
   ///
   /// Uses Newton-Raphson iteration with Brent's method fallback.
   /// The function F(V) = sigma_n * f(V, theta) + eta * V - tau = 0
   /// is monotonically increasing in V, guaranteeing a unique solution.
   ///
   /// @param[in] tau Total shear stress [Pa]
   /// @param[in] theta State variable [s]
   /// @param[in] sigma_n Normal stress [Pa]
   /// @param[in] eta Radiation damping coefficient [Pa*s/m]
   /// @param[in] a Rate-and-state direct effect parameter
   /// @param[out] iterations Optional: number of iterations used
   /// @return Slip rate V [m/s]
   real_t SolveSlipRate(real_t tau, real_t theta, real_t sigma_n,
                        real_t eta, real_t a,
                        int *iterations = nullptr) const override
   {
      // Handle fault in tension (sigma_n <= 0)
      // Per Tandem implementation: return viscous sliding if eta > 0
      if (sigma_n <= 0.0)
      {
         if (iterations) { *iterations = 0; }
         if (eta > 0.0)
         {
            return tau / eta;  // Viscous sliding
         }
         else
         {
            return 0.0;  // Cannot determine slip rate without friction or damping
         }
      }

      // Initial guess: start near the reference velocity
      real_t V = cp_.V0;

      // Bounds for bracketing (slip rate should be positive)
      real_t V_lo = V_min_;
      // Physical upper bound: when friction = 0, tau = eta * V
      // V_max = tau / eta (with fallback if eta is very small)
      // Following Tandem's approach for tighter bracketing
      real_t V_hi = (eta > 1e-6) ? (tau / eta) : 100.0;
      V_hi = std::min(V_hi, 100.0);  // Cap at 100 m/s for safety

      // Newton-Raphson iteration
      const int max_iter = 100;
      const real_t tol = 1.0e-12;

      int iter = 0;
      for (; iter < max_iter; ++iter)
      {
         real_t f = FrictionCoefficient(V, theta, a);
         real_t df_dV = FrictionDerivativeV(V, theta, a);

         real_t F = sigma_n * f + eta * V - tau;
         real_t dF_dV = sigma_n * df_dV + eta;

         // Newton update
         real_t dV = -F / dF_dV;

         // Limit step size to stay in reasonable bounds
         real_t V_new = V + dV;

         // Ensure positivity and reasonable bounds
         if (V_new < V_lo)
         {
            V_new = 0.5 * (V + V_lo);
         }
         else if (V_new > V_hi)
         {
            V_new = 0.5 * (V + V_hi);
         }

         // Check convergence
         real_t rel_change = std::abs(V_new - V) / std::max(V, V_min_);
         V = V_new;

         if (rel_change < tol || std::abs(F) < tol * tau)
         {
            break;
         }
      }

      if (iterations != nullptr)
      {
         *iterations = iter;
      }

      // Verify solution
      MFEM_ASSERT(iter < max_iter,
                  "Newton solver failed to converge for slip rate");

      return V;
   }

   /// Compute initial state theta(0) for stress equilibrium.
   ///
   /// Given tau0 and V_init, solve for theta such that:
   ///   tau0 = sigma_n * f(V_init, theta) + eta * V_init
   ///
   /// Rearranging the friction law:
   ///   f = tau0' / sigma_n, where tau0' = tau0 - eta * V_init
   ///   a * asinh(arg) = f
   ///   arg = sinh(f / a)
   ///   (V_init / 2V0) * exp((f0 + b * ln(V0 * theta / Dc)) / a) = sinh(f / a)
   ///   exp((f0 + b * ln(V0 * theta / Dc)) / a) = (2V0 / V_init) * sinh(f / a)
   ///   (f0 + b * ln(V0 * theta / Dc)) / a = ln[(2V0 / V_init) * sinh(f / a)]
   ///   b * ln(V0 * theta / Dc) = a * ln[(2V0 / V_init) * sinh(f / a)] - f0
   ///   ln(V0 * theta / Dc) = (a / b) * ln[(2V0 / V_init) * sinh(f / a)] - f0 / b
   ///   theta = (Dc / V0) * exp{(a / b) * ln[(2V0 / V_init) * sinh(f / a)] - f0 / b}
   ///
   /// @param[in] tau0 Pre-stress [Pa]
   /// @param[in] V_init Initial slip rate [m/s]
   /// @param[in] sigma_n Normal stress [Pa]
   /// @param[in] eta Radiation damping coefficient [Pa*s/m]
   /// @param[in] a Rate-and-state direct effect parameter
   /// @return Initial state variable theta [s]
   real_t InitialState(real_t tau0, real_t V_init, real_t sigma_n,
                       real_t eta, real_t a) const override
   {
      // Effective stress after removing radiation damping
      real_t tau_eff = tau0 - eta * V_init;

      // Friction coefficient at initial state
      real_t f = tau_eff / sigma_n;

      // Solve for theta using the inverted friction law
      // theta = (Dc / V0) * exp{(a / b) * ln[(2V0 / V_init) * sinh(f / a)] - f0 / b}

      // Protect against sinh overflow for large f/a
      // sinh(x) ≈ exp(x)/2 for large x; sinh(700) ≈ 5e303, near double max
      real_t f_over_a = f / a;
      real_t sinh_arg;
      if (f_over_a > 700.0)
      {
         // Use asymptotic approximation: sinh(x) ≈ exp(x)/2 for large x
         sinh_arg = 0.5 * std::exp(f_over_a);
      }
      else
      {
         sinh_arg = std::sinh(f_over_a);
      }
      real_t log_arg = (2.0 * cp_.V0 / V_init) * sinh_arg;

      // Handle potential numerical issues
      MFEM_ASSERT(log_arg > 0.0, "Invalid argument for logarithm in InitialState");

      real_t exp_arg = (a / cp_.b) * std::log(log_arg) - cp_.f0 / cp_.b;
      real_t theta = (cp_.Dc / cp_.V0) * std::exp(exp_arg);

      return theta;
   }

   // =========================================================================
   // Psi-space methods (logarithmic state variable)
   // =========================================================================
   // psi = f0 + b*ln(V0*theta/Dc)
   // These methods work directly with psi, avoiding exp/log round-trips.

   /// Convert theta to psi: psi = f0 + b*ln(V0*theta/Dc).
   real_t ThetaToPsi(real_t theta) const
   {
      theta = std::max(theta, theta_min_);
      return cp_.f0 + cp_.b * std::log(cp_.V0 * theta / cp_.Dc);
   }

   /// Convert psi to theta: theta = (Dc/V0)*exp((psi - f0)/b).
   real_t PsiToTheta(real_t psi) const
   {
      return (cp_.Dc / cp_.V0) * std::exp((psi - cp_.f0) / cp_.b);
   }

   /// Friction coefficient in psi-space:
   /// f(V, psi) = a * asinh[(V / 2V0) * exp(psi / a)]
   real_t FrictionCoefficientPsi(real_t V, real_t psi, real_t a) const
   {
      V = std::max(V, V_min_);
      real_t sinh_arg = (V / (2.0 * cp_.V0)) * std::exp(psi / a);
      return a * std::asinh(sinh_arg);
   }

   /// Derivative df/dV in psi-space:
   /// df/dV = a / (2V0) * exp(psi/a) / sqrt(1 + [(V/2V0)*exp(psi/a)]^2)
   real_t FrictionDerivativeVPsi(real_t V, real_t psi, real_t a) const
   {
      V = std::max(V, V_min_);
      real_t exp_val = std::exp(psi / a);
      real_t sinh_arg = (V / (2.0 * cp_.V0)) * exp_val;
      real_t d_arg_dV = exp_val / (2.0 * cp_.V0);
      real_t d_asinh = 1.0 / std::sqrt(1.0 + sinh_arg * sinh_arg);
      return a * d_asinh * d_arg_dV;
   }

   /// Solve for slip rate V given stress tau and state psi.
   ///
   /// Solves: tau = sigma_n * f(V, psi) + eta * V
   /// where f(V, psi) = a * asinh[(V / 2V0) * exp(psi / a)]
   real_t SolveSlipRatePsi(real_t tau, real_t psi, real_t sigma_n,
                           real_t eta, real_t a,
                           int *iterations = nullptr) const
   {
      if (sigma_n <= 0.0)
      {
         if (iterations) { *iterations = 0; }
         if (eta > 0.0) { return tau / eta; }
         else { return 0.0; }
      }

      real_t V = cp_.V0;
      real_t V_lo = V_min_;
      real_t V_hi = (eta > 1e-6) ? (tau / eta) : 100.0;
      V_hi = std::min(V_hi, 100.0);

      const int max_iter = 100;
      const real_t tol = 1.0e-12;

      int iter = 0;
      for (; iter < max_iter; ++iter)
      {
         real_t f = FrictionCoefficientPsi(V, psi, a);
         real_t df_dV = FrictionDerivativeVPsi(V, psi, a);

         real_t F = sigma_n * f + eta * V - tau;
         real_t dF_dV = sigma_n * df_dV + eta;

         real_t V_new = V - F / dF_dV;

         if (V_new < V_lo) { V_new = 0.5 * (V + V_lo); }
         else if (V_new > V_hi) { V_new = 0.5 * (V + V_hi); }

         real_t rel_change = std::abs(V_new - V) / std::max(V, V_min_);
         V = V_new;

         if (rel_change < tol || std::abs(F) < tol * tau) { break; }
      }

      if (iterations) { *iterations = iter; }
      MFEM_ASSERT(iter < max_iter,
                  "Newton solver failed to converge for slip rate (psi)");
      return V;
   }

   /// Compute initial psi from stress equilibrium.
   ///
   /// Given tau0 and V_init, solve for psi such that:
   ///   tau0 = sigma_n * f(V_init, psi) + eta * V_init
   ///
   /// f = (tau0 - eta*V_init) / sigma_n
   /// a * asinh(arg) = f  =>  arg = sinh(f/a)
   /// (V_init/2V0) * exp(psi/a) = sinh(f/a)
   /// psi = a * ln[2V0/V_init * sinh(f/a)]
   real_t InitialStatePsi(real_t tau0, real_t V_init, real_t sigma_n,
                          real_t eta, real_t a) const
   {
      real_t tau_eff = tau0 - eta * V_init;
      real_t f = tau_eff / sigma_n;

      real_t f_over_a = f / a;
      real_t sinh_val;
      if (f_over_a > 700.0)
      {
         sinh_val = 0.5 * std::exp(f_over_a);
      }
      else
      {
         sinh_val = std::sinh(f_over_a);
      }

      real_t log_arg = (2.0 * cp_.V0 / V_init) * sinh_val;
      MFEM_ASSERT(log_arg > 0.0,
                  "Invalid argument for logarithm in InitialStatePsi");

      return a * std::log(log_arg);
   }

   // =========================================================================
   // Accessors for constants
   // =========================================================================

   real_t GetV0() const override { return cp_.V0; }
   real_t GetF0() const override { return cp_.f0; }
   real_t GetB() const override { return cp_.b; }
   real_t GetDc() const override { return cp_.Dc; }

   /// Get reference to the constants structure.
   const Constants &GetConstants() const { return cp_; }

   /// Set the constants.
   void SetConstants(const Constants &c) { cp_ = c; }

private:
   Constants cp_;

   /// Minimum slip rate to avoid numerical issues
   static constexpr real_t V_min_ = 1.0e-30;

   /// Minimum state variable to avoid numerical issues
   static constexpr real_t theta_min_ = 1.0e-30;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DIETERICH_RUINA_HPP
