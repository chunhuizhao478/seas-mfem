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
#include <functional>
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
      if (V <= 0.0) { return 0.0; }
      theta = std::max(theta, theta_min_);

      // f = a * asinh[(V / 2V0) * exp((f0 + b * ln(V0 * theta / Dc)) / a)]
      real_t log_arg = cp_.V0 * theta / cp_.Dc;
      real_t exp_arg = (cp_.f0 + cp_.b * std::log(log_arg)) / a;

      // Safe evaluation for large exp_arg (equivalent to large psi/a)
      if (exp_arg > 700.0)
      {
         // asinh(x) ≈ log(2x) for large x
         // sinh_arg = (V/2V0)*exp(exp_arg), log(2*sinh_arg) = log(V/V0) + exp_arg
         return std::max(0.0, a * (std::log(V / cp_.V0) + exp_arg));
      }

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
      if (sigma_n <= 0.0)
      {
         if (iterations) { *iterations = 0; }
         if (eta > 0.0) { return tau / eta; }
         else { return 0.0; }
      }

      // Use Brent's method with bracket [0, V_hi]
      // Upper bound: when friction ≥ 0, V ≤ tau/eta. Cap for safety.
      real_t V_lo = 0.0;
      real_t V_hi = (eta > 0.0) ? (tau / eta) : 100.0;
      V_hi = std::min(V_hi, 100.0);

      auto residual = [&](real_t V) -> real_t
      {
         // FrictionCoefficient clamps V >= V_min_ internally
         return tau - sigma_n * FrictionCoefficient(V, theta, a) - eta * V;
      };

      real_t V = zeroIn(V_lo, V_hi, residual);

      if (iterations) { *iterations = 0; }
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

   /// Convert psi to theta with per-DOF Dc: theta = (Dc/V0)*exp((psi - f0)/b).
   real_t PsiToTheta(real_t psi, real_t Dc) const
   {
      return (Dc / cp_.V0) * std::exp((psi - cp_.f0) / cp_.b);
   }

   /// Friction coefficient in psi-space:
   /// f(V, psi) = a * asinh[(V / 2V0) * exp(psi / a)]
   ///
   /// Safe evaluation: handles V=0 and large psi/a without overflow.
   /// At V=0, f=0 by convention (no sliding, no friction force).
   /// For large psi/a (>=700), uses asinh(x) ≈ log(2x) for large x.
   real_t FrictionCoefficientPsi(real_t V, real_t psi, real_t a) const
   {
      if (V <= 0.0) { return 0.0; }

      real_t psi_over_a = psi / a;
      if (psi_over_a > 700.0)
      {
         // For large psi/a, exp(psi/a) overflows. Use:
         // asinh(x) = log(x + sqrt(x^2+1)) ≈ log(2x) for large x
         // sinh_arg = (V/2V0)*exp(psi/a), so log(2*sinh_arg) = log(V/V0) + psi/a
         // f = a * [log(V/V0) + psi/a] = a*log(V/V0) + psi
         real_t log_term = std::log(V / cp_.V0);
         return std::max(0.0, a * log_term + psi);
      }

      real_t sinh_arg = (V / (2.0 * cp_.V0)) * std::exp(psi_over_a);
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
   ///
   /// Uses Brent's method with bracket [0, tau/eta].
   /// - At V=0: R(0) = tau > 0
   /// - At V=tau/eta: R = -sigma_n*f < 0
   ///
   /// When psi is very negative (collapsed state variable during RK45
   /// intermediate stages), f → 0 and R(V_hi) ≈ 0 may become slightly
   /// positive due to floating-point error. In this degenerate case,
   /// friction is negligible and V ≈ tau/eta is returned.
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

      if (tau <= 0.0)
      {
         if (iterations) { *iterations = 0; }
         return 0.0;
      }

      // Brent's method with bracket [0, tau/eta]
      real_t V_lo = 0.0;
      real_t V_hi = tau / eta;

      auto residual = [&](real_t V) -> real_t
      {
         return tau - sigma_n * FrictionCoefficientPsi(V, psi, a) - eta * V;
      };

      real_t Fb = residual(V_hi);
      if (Fb >= 0.0)
      {
         // Friction is negligible (collapsed state variable, psi << 0).
         static int degen_count = 0;
         if (++degen_count <= 5)
         {
            std::cerr << "[WARNING] SolveSlipRatePsi degenerate case #"
                      << degen_count << ": psi=" << psi
                      << " a=" << a << " tau=" << tau
                      << " V=tau/eta=" << V_hi << "\n";
         }
         if (iterations) { *iterations = 0; }
         return V_hi;
      }

      real_t V = zeroIn(V_lo, V_hi, residual);

      if (iterations) { *iterations = 0; }
      return V;
   }

   /// Solve for 2-component slip rate given 2-component traction and scalar psi.
   ///
   /// Algorithm (following Tandem's DieterichRuinaAgeing::slip_rate):
   /// 1. tau_abs = ||tau_vec||
   /// 2. V_abs = SolveSlipRatePsi(tau_abs, psi, sigma_n, eta, a)
   /// 3. V_vec = (V_abs / tau_abs) * tau_vec
   ///
   /// Slip velocity is parallel to traction
   /// (slip occurs in the direction of driving stress).
   ///
   /// @param[in] tau_vec Traction vector (2 components) [Pa]
   /// @param[in] psi Logarithmic state variable [-]
   /// @param[in] sigma_n Normal stress [Pa]
   /// @param[in] eta Radiation damping coefficient [Pa*s/m]
   /// @param[in] a Rate-and-state direct effect parameter
   /// @param[out] V_vec Slip rate vector (2 components) [m/s]
   /// @param[out] iterations Optional: number of iterations used
   void SolveSlipRateVectorPsi(const real_t tau_vec[2], real_t psi,
                                real_t sigma_n, real_t eta, real_t a,
                                real_t V_vec[2],
                                int *iterations = nullptr) const
   {
      real_t tau_abs = std::sqrt(tau_vec[0] * tau_vec[0] +
                                 tau_vec[1] * tau_vec[1]);
      if (tau_abs <= 0.0)
      {
         V_vec[0] = 0.0;
         V_vec[1] = 0.0;
         if (iterations) { *iterations = 0; }
         return;
      }
      real_t V_abs = SolveSlipRatePsi(tau_abs, psi, sigma_n, eta, a, iterations);
      // v55 D8: anti-parallel to tau, matching Tandem DieterichRuinaBase.h:174.
      // The negation only affects the internal state variable S (cumulative slip).
      // GetSlip() negates S back to physical slip before the domain solver sees it.
      V_vec[0] = -(V_abs / tau_abs) * tau_vec[0];
      V_vec[1] = -(V_abs / tau_abs) * tau_vec[1];
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

   /// Minimum slip rate to avoid numerical issues (used in friction evaluation)
   static constexpr real_t V_min_ = 1.0e-30;

   /// Minimum state variable to avoid numerical issues
   static constexpr real_t theta_min_ = 1.0e-30;

   /// Brent's method for finding zeros of a function in [a, b].
   ///
   /// Port of Tandem's zeroIn (from Forsythe/Malcolm/Moler ZEROIN).
   /// F(a) and F(b) must have opposite signs (or F(a)==0).
   /// Returns x in [a,b] with |bracket| <= tol + 4*eps*|x|.
   static real_t zeroIn(real_t a, real_t b,
                         std::function<real_t(real_t)> F,
                         real_t tol = 0.0)
   {
      real_t eps = std::numeric_limits<real_t>::epsilon();
      real_t Fa = F(a);
      if (Fa == 0.0) { return a; }
      real_t Fb = F(b);

      MFEM_VERIFY(!std::isnan(Fa) && !std::isinf(Fa),
                  "zeroIn: NaN/Inf at a=" << a << ", F(a)=" << Fa);
      MFEM_VERIFY(!std::isnan(Fb) && !std::isinf(Fb),
                  "zeroIn: NaN/Inf at b=" << b << ", F(b)=" << Fb);
      MFEM_VERIFY(std::copysign(Fa, Fb) != Fa || Fb == 0.0,
                  "zeroIn: F(a) and F(b) must have different signs. "
                  << "a=" << a << " F(a)=" << Fa
                  << " b=" << b << " F(b)=" << Fb);

      real_t c = a, Fc = Fa;
      real_t d = b - a, e = d;

      while (Fb != 0.0)
      {
         if (std::copysign(Fb, Fc) == Fb)
         {
            c = a; Fc = Fa; d = b - a; e = d;
         }
         if (std::fabs(Fc) < std::fabs(Fb))
         {
            a = b; b = c; c = a;
            Fa = Fb; Fb = Fc; Fc = Fa;
         }
         real_t xm = 0.5 * (c - b);
         real_t tol1 = 2.0 * eps * std::fabs(b) + 0.5 * tol;
         if (std::fabs(xm) <= tol1 || Fb == 0.0) { break; }
         if (std::fabs(e) < tol1 || std::fabs(Fa) <= std::fabs(Fb))
         {
            d = xm; e = d;  // bisection
         }
         else
         {
            real_t s = Fb / Fa;
            real_t p, q;
            if (a != c)
            {
               // inverse quadratic interpolation
               real_t qq = Fa / Fc;
               real_t r = Fb / Fc;
               p = s * (2.0 * xm * qq * (qq - r) - (b - a) * (r - 1.0));
               q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
            }
            else
            {
               // linear interpolation
               p = 2.0 * xm * s;
               q = 1.0 - s;
            }
            if (p > 0) { q = -q; } else { p = -p; }
            if (2.0 * p < 3.0 * xm * q - std::fabs(tol1 * q) &&
                p < std::fabs(0.5 * e * q))
            {
               e = d; d = p / q;
            }
            else
            {
               d = xm; e = d;  // bisection
            }
         }
         a = b; Fa = Fb;
         if (std::fabs(d) > tol1) { b += d; }
         else { b += std::copysign(tol1, xm); }
         Fb = F(b);
         MFEM_VERIFY(!std::isnan(Fb) && !std::isinf(Fb),
                     "zeroIn: NaN/Inf at b=" << b << ", F(b)=" << Fb);
      }
      return b;
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DIETERICH_RUINA_HPP
