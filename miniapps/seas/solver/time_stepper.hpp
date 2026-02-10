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

#ifndef MFEM_SEAS_TIME_STEPPER_HPP
#define MFEM_SEAS_TIME_STEPPER_HPP

#include "mfem.hpp"
#include <algorithm>
#include <cmath>

namespace mfem
{
namespace seas
{

/// @brief Adaptive time stepper for SEAS simulations.
///
/// Adjusts the time step based on the maximum slip rate on the fault.
/// During interseismic periods (slow slip), large time steps are used.
/// During nucleation and coseismic phases (fast slip), small time steps
/// are used to resolve the dynamics.
///
/// Time step control strategy:
///   V < 1e-12 m/s  -> dt_max (very slow interseismic)
///   V > V_target   -> dt * V_target / V (rapid slip, reduce dt)
///   otherwise      -> min(dt * growth_factor, dt_max) (gradual increase)
class AdaptiveTimeStepper
{
public:
   /// @brief Default constructor with standard SEAS parameters.
   AdaptiveTimeStepper()
      : dt_min_(1e-6),
        dt_max_(3.15e7),     // ~1 year in seconds
        V_target_(1e-6),     // Target max slip rate for dt control
        growth_factor_(1.5), // Factor for increasing dt
        dt_(1e3)             // Initial time step (1000 seconds)
   {}

   // =========================================================================
   // Configuration
   // =========================================================================

   /// Set minimum time step [seconds].
   void SetDtMin(real_t dt_min) { dt_min_ = dt_min; }

   /// Set maximum time step [seconds].
   void SetDtMax(real_t dt_max) { dt_max_ = dt_max; }

   /// Set target slip rate for dt control [m/s].
   void SetTargetSlipRateMax(real_t V_target) { V_target_ = V_target; }

   /// Set growth factor for increasing dt during slow slip.
   void SetGrowthFactor(real_t factor) { growth_factor_ = factor; }

   /// Set initial time step [seconds].
   void SetInitialDt(real_t dt0) { dt_ = dt0; }

   // =========================================================================
   // Time step computation
   // =========================================================================

   /// @brief Compute new time step based on maximum slip rate.
   ///
   /// @param V_max Maximum slip rate on the fault [m/s]
   /// @param dt_current Current time step [seconds]
   /// @return New time step [seconds], clamped to [dt_min, dt_max]
   real_t ComputeNewDt(real_t V_max, real_t dt_current) const
   {
      real_t dt_new;

      if (V_max < 1e-12)
      {
         // Very slow slip: use maximum dt
         dt_new = dt_max_;
      }
      else if (V_max > V_target_)
      {
         // Rapid slip: reduce time step proportionally
         dt_new = dt_current * V_target_ / V_max;
      }
      else
      {
         // Slow slip: gradually increase time step
         dt_new = std::min(dt_current * growth_factor_, dt_max_);
      }

      // Clamp to bounds
      return std::max(dt_min_, std::min(dt_max_, dt_new));
   }

   // =========================================================================
   // Accessors
   // =========================================================================

   /// Get current time step.
   real_t GetDt() const { return dt_; }

   /// Update current time step.
   void SetDt(real_t dt) { dt_ = dt; }

   /// Get minimum time step.
   real_t GetDtMin() const { return dt_min_; }

   /// Get maximum time step.
   real_t GetDtMax() const { return dt_max_; }

   /// Get target slip rate.
   real_t GetTargetSlipRateMax() const { return V_target_; }

private:
   real_t dt_min_;          ///< Minimum time step [seconds]
   real_t dt_max_;          ///< Maximum time step [seconds]
   real_t V_target_;        ///< Target max slip rate for dt control [m/s]
   real_t growth_factor_;   ///< Factor for increasing dt during slow slip
   real_t dt_;              ///< Current time step [seconds]
};

/// @brief Dormand-Prince 5(4) embedded Runge-Kutta integrator with adaptive
///        step size control based on local truncation error estimation.
///
/// Uses the standard DOPRI5(4) pair (7 stages, FSAL). The 5th-order solution
/// advances the state; the difference from the 4th-order solution provides
/// a local truncation error estimate. Step size is adjusted so the weighted
/// L-infinity error norm stays below 1.
///
/// Designed following Tandem's approach (PETSc TS with `-ts_rk_type 5dp`,
/// `-ts_atol 1e-7`, `-ts_adapt_wnormtype infinity`).
class DormandPrinceRK45
{
public:
   DormandPrinceRK45()
      : atol_(1e-7),
        rtol_(1e-7),
        safety_(0.9),
        growth_max_(5.0),
        shrink_min_(0.2),
        dt_min_(1e-6),
        dt_max_(0.1 * 3.15576e7),  // 0.1 year
        dt_(1e3),
        initialized_(false),
        total_rejections_(0),
        diag_count_(0)
   {}

   // =========================================================================
   // Configuration
   // =========================================================================

   void SetAbsTol(real_t atol) { atol_ = atol; }
   void SetRelTol(real_t rtol) { rtol_ = rtol; }
   void SetSafety(real_t safety) { safety_ = safety; }
   void SetGrowthMax(real_t gmax) { growth_max_ = gmax; }
   void SetShrinkMin(real_t smin) { shrink_min_ = smin; }
   void SetDtMin(real_t dt_min) { dt_min_ = dt_min; }
   void SetDtMax(real_t dt_max) { dt_max_ = dt_max; }
   void SetDt(real_t dt) { dt_ = dt; }

   // =========================================================================
   // Initialization
   // =========================================================================

   /// @brief Initialize the integrator for a given operator.
   ///
   /// Must be called before the first Step(). Allocates stage vectors.
   void Init(TimeDependentOperator &op)
   {
      int n = op.Width();
      for (int i = 0; i < 7; i++) { k_[i].SetSize(n); }
      y_tmp_.SetSize(n);
      err_.SetSize(n);
      initialized_ = false;
   }

   /// @brief Take one adaptive step.
   ///
   /// Attempts a step of size dt_. On success, advances state and t,
   /// updates dt_ for the next step, and returns true. On failure
   /// (error too large), shrinks dt_ and returns false — the caller
   /// should retry without incrementing the step counter.
   ///
   /// @param op The ODE operator (dy/dt = op.Mult(y))
   /// @param state [in/out] State vector
   /// @param t [in/out] Current time (advanced on success)
   /// @param dt [out] The step size that was used (or attempted)
   /// @return true if step accepted, false if rejected
   bool Step(TimeDependentOperator &op, Vector &state, real_t &t, real_t &dt)
   {
      const int n = state.Size();
      dt = dt_;

      // Stage 1: use FSAL from previous step, or compute fresh
      if (!initialized_)
      {
         op.SetTime(t);
         op.Mult(state, k_[0]);
         initialized_ = true;
      }
      // else k_[0] = k_[6] from previous accepted step (FSAL)

      // Stage 2
      for (int i = 0; i < n; i++)
      {
         y_tmp_(i) = state(i) + dt * (a21 * k_[0](i));
      }
      op.SetTime(t + c2 * dt);
      op.Mult(y_tmp_, k_[1]);

      // Stage 3
      for (int i = 0; i < n; i++)
      {
         y_tmp_(i) = state(i) + dt * (a31 * k_[0](i) + a32 * k_[1](i));
      }
      op.SetTime(t + c3 * dt);
      op.Mult(y_tmp_, k_[2]);

      // Stage 4
      for (int i = 0; i < n; i++)
      {
         y_tmp_(i) = state(i) + dt * (a41 * k_[0](i) + a42 * k_[1](i)
                                       + a43 * k_[2](i));
      }
      op.SetTime(t + c4 * dt);
      op.Mult(y_tmp_, k_[3]);

      // Stage 5
      for (int i = 0; i < n; i++)
      {
         y_tmp_(i) = state(i) + dt * (a51 * k_[0](i) + a52 * k_[1](i)
                                       + a53 * k_[2](i) + a54 * k_[3](i));
      }
      op.SetTime(t + c5 * dt);
      op.Mult(y_tmp_, k_[4]);

      // Stage 6
      for (int i = 0; i < n; i++)
      {
         y_tmp_(i) = state(i) + dt * (a61 * k_[0](i) + a62 * k_[1](i)
                                       + a63 * k_[2](i) + a64 * k_[3](i)
                                       + a65 * k_[4](i));
      }
      op.SetTime(t + dt);
      op.Mult(y_tmp_, k_[5]);

      // Stage 7 (= 5th-order solution, also next k_[0] via FSAL)
      // y5 = state + dt * (b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
      // (b2 = 0, b7 = 0)
      for (int i = 0; i < n; i++)
      {
         y_tmp_(i) = state(i) + dt * (b1 * k_[0](i) + b3 * k_[2](i)
                                       + b4 * k_[3](i) + b5 * k_[4](i)
                                       + b6 * k_[5](i));
      }
      op.SetTime(t + dt);
      op.Mult(y_tmp_, k_[6]);

      // Error estimate: e = dt * (b - b*) . k
      // e_i = dt * (e1*k1_i + e3*k3_i + e4*k4_i + e5*k5_i + e6*k6_i + e7*k7_i)
      for (int i = 0; i < n; i++)
      {
         err_(i) = dt * (e1 * k_[0](i) + e3 * k_[2](i) + e4 * k_[3](i)
                         + e5 * k_[4](i) + e6 * k_[5](i) + e7 * k_[6](i));
      }

      // Compute weighted L-infinity error norm
      real_t err_norm = 0.0;
      int worst_idx = 0;
      for (int i = 0; i < n; i++)
      {
         real_t scale = atol_ + rtol_ * std::abs(y_tmp_(i));
         real_t ei = std::abs(err_(i)) / scale;
         if (ei > err_norm)
         {
            err_norm = ei;
            worst_idx = i;
         }
      }

      // Compute new dt using standard PI controller formula
      //   dt_new = safety * dt * err_norm^(-1/q), q = min(p, p*) = 5
      real_t dt_new;
      if (err_norm <= 0.0)
      {
         dt_new = dt * growth_max_;
      }
      else
      {
         dt_new = safety_ * dt * std::pow(err_norm, -0.2);
      }

      // Clamp growth/shrink factors
      dt_new = std::max(dt_new, shrink_min_ * dt);
      dt_new = std::min(dt_new, growth_max_ * dt);
      dt_new = std::max(dt_min_, std::min(dt_max_, dt_new));

      if (err_norm <= 1.0)
      {
         // Accept step
         state = y_tmp_;
         t += dt;
         dt_ = dt_new;

         // FSAL: k_[6] becomes k_[0] for next step
         k_[0] = k_[6];

         // Diagnostic: log when accepted but dt stuck near dt_min
         if (dt <= dt_min_ * 1.5 && dt_new <= dt_min_ * 1.5)
         {
            diag_count_++;
            if (diag_count_ <= 5 || diag_count_ % 10000 == 0)
            {
               int dof = worst_idx / 2;
               bool is_theta = (worst_idx % 2 == 1);
               std::cout << "[RK45 dt_min] step accepted, dt=" << dt << " s"
                         << ", err_norm=" << err_norm
                         << ", dt_new=" << dt_new
                         << ", worst: DOF " << dof
                         << (is_theta ? " (theta)" : " (slip)")
                         << ", |err|=" << std::abs(err_(worst_idx))
                         << ", |y|=" << std::abs(y_tmp_(worst_idx))
                         << ", diag_count=" << diag_count_
                         << "\n";
            }
         }
         else
         {
            diag_count_ = 0;
         }

         return true;
      }
      else
      {
         // Reject step
         dt_ = dt_new;
         initialized_ = true;  // keep k_[0] from this step start
         total_rejections_++;

         // Diagnostic: log when stuck at dt_min
         if (dt_ <= dt_min_ * 1.01)
         {
            int dof = worst_idx / 2;
            bool is_theta = (worst_idx % 2 == 1);
            real_t scale = atol_ + rtol_ * std::abs(y_tmp_(worst_idx));
            std::cout << "[RK45 STUCK] dt=" << dt << " s"
                      << ", err_norm=" << err_norm
                      << ", worst: DOF " << dof
                      << (is_theta ? " (theta)" : " (slip)")
                      << ", |err|=" << std::abs(err_(worst_idx))
                      << ", scale=" << scale
                      << ", |y|=" << std::abs(y_tmp_(worst_idx))
                      << ", |k0|=" << std::abs(k_[0](worst_idx))
                      << "\n";
         }

         return false;
      }
   }

   // =========================================================================
   // Accessors
   // =========================================================================

   real_t GetDt() const { return dt_; }
   real_t GetDtMin() const { return dt_min_; }
   real_t GetDtMax() const { return dt_max_; }
   int GetTotalRejections() const { return total_rejections_; }

private:
   real_t atol_;           ///< Absolute tolerance
   real_t rtol_;           ///< Relative tolerance
   real_t safety_;         ///< Safety factor for dt adjustment
   real_t growth_max_;     ///< Maximum dt growth factor per step
   real_t shrink_min_;     ///< Minimum dt shrink factor per step
   real_t dt_min_;         ///< Minimum allowed dt
   real_t dt_max_;         ///< Maximum allowed dt
   real_t dt_;             ///< Current time step
   bool initialized_;      ///< Whether k_[0] is valid from a previous step
   int total_rejections_;  ///< Total number of rejected steps
   int diag_count_;        ///< Counter for dt_min diagnostic messages

   Vector k_[7];  ///< Stage vectors
   Vector y_tmp_; ///< Temporary solution
   Vector err_;   ///< Error estimate

   // Dormand-Prince 5(4) coefficients
   // Nodes
   static constexpr real_t c2 = 1.0 / 5.0;
   static constexpr real_t c3 = 3.0 / 10.0;
   static constexpr real_t c4 = 4.0 / 5.0;
   static constexpr real_t c5 = 8.0 / 9.0;

   // A matrix (lower triangular)
   static constexpr real_t a21 = 1.0 / 5.0;
   static constexpr real_t a31 = 3.0 / 40.0;
   static constexpr real_t a32 = 9.0 / 40.0;
   static constexpr real_t a41 = 44.0 / 45.0;
   static constexpr real_t a42 = -56.0 / 15.0;
   static constexpr real_t a43 = 32.0 / 9.0;
   static constexpr real_t a51 = 19372.0 / 6561.0;
   static constexpr real_t a52 = -25360.0 / 2187.0;
   static constexpr real_t a53 = 64448.0 / 6561.0;
   static constexpr real_t a54 = -212.0 / 729.0;
   static constexpr real_t a61 = 9017.0 / 3168.0;
   static constexpr real_t a62 = -355.0 / 33.0;
   static constexpr real_t a63 = 46732.0 / 5247.0;
   static constexpr real_t a64 = 49.0 / 176.0;
   static constexpr real_t a65 = -5103.0 / 18656.0;

   // 5th order weights (b)
   static constexpr real_t b1 = 35.0 / 384.0;
   // b2 = 0
   static constexpr real_t b3 = 500.0 / 1113.0;
   static constexpr real_t b4 = 125.0 / 192.0;
   static constexpr real_t b5 = -2187.0 / 6784.0;
   static constexpr real_t b6 = 11.0 / 84.0;
   // b7 = 0

   // Error coefficients: e_i = b_i - b*_i
   static constexpr real_t e1 = 35.0/384.0 - 5179.0/57600.0;
   // e2 = 0
   static constexpr real_t e3 = 500.0/1113.0 - 7571.0/16695.0;
   static constexpr real_t e4 = 125.0/192.0 - 393.0/640.0;
   static constexpr real_t e5 = -2187.0/6784.0 + 92097.0/339200.0;
   static constexpr real_t e6 = 11.0/84.0 - 187.0/2100.0;
   static constexpr real_t e7 = -1.0/40.0;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TIME_STEPPER_HPP
