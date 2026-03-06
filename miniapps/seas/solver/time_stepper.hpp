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
#include "../common/mpi_context.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

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
        rtol_(1e-50),           // Match Tandem/PETSc: pure absolute tolerance
        safety_(0.9),
        reject_safety_(0.5),    // PETSc default: extra shrink after rejection
        growth_max_(10.0),      // PETSc default clip[1]
        shrink_min_(0.1),       // PETSc default clip[0]
        dt_min_(1e-6),
        dt_max_(0.5 * 3.15576e7),  // 0.5 year
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
   void SetRejectSafety(real_t rs) { reject_safety_ = rs; }
   void SetGrowthMax(real_t gmax) { growth_max_ = gmax; }
   void SetShrinkMin(real_t smin) { shrink_min_ = smin; }
   void SetDtMin(real_t dt_min) { dt_min_ = dt_min; }
   void SetDtMax(real_t dt_max) { dt_max_ = dt_max; }
   void SetDt(real_t dt) { dt_ = dt; }

   /// Enable verbose per-step diagnostics (error norm, worst DOF, etc.)
   void SetVerbose(bool v) { diag_verbose_ = v; }

   /// Set number of state components per fault node for diagnostic output.
   /// BP2: 2 (slip, theta), BP5: 3 (slip_dip, slip_strike, psi).
   void SetStatePerNode(int spn) { state_per_node_ = spn; }

   /// Use weighted RMS (2-norm) instead of L-infinity for error norm.
   /// More robust to outlier DOFs at MPI partition boundaries.
   void SetUse2Norm(bool v) { use_2norm_ = v; }

   /// @brief Set MPI context for parallel error norm reduction.
   ///
   /// In parallel, the error norm must be reduced across all ranks
   /// so that accept/reject decisions are consistent. Without this,
   /// ranks can diverge, causing deadlock in subsequent MPI collectives.
   void SetMPIContext(MPIContext *ctx) { mpi_ctx_ = ctx; }

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

      // Check for NaN in stage 2 (skip remaining stages on failure)
      if (!std::isfinite(NormL2(k_[1])))
      {
         dt_ = std::max(dt_min_, dt * shrink_min_);
         initialized_ = false;
         total_rejections_++;
         return false;
      }

      // Stage 3
      for (int i = 0; i < n; i++)
      {
         y_tmp_(i) = state(i) + dt * (a31 * k_[0](i) + a32 * k_[1](i));
      }
      op.SetTime(t + c3 * dt);
      op.Mult(y_tmp_, k_[2]);

      if (!std::isfinite(NormL2(k_[2])))
      {
         dt_ = std::max(dt_min_, dt * shrink_min_);
         initialized_ = false;
         total_rejections_++;
         return false;
      }

      // Stage 4
      for (int i = 0; i < n; i++)
      {
         y_tmp_(i) = state(i) + dt * (a41 * k_[0](i) + a42 * k_[1](i)
                                       + a43 * k_[2](i));
      }
      op.SetTime(t + c4 * dt);
      op.Mult(y_tmp_, k_[3]);

      if (!std::isfinite(NormL2(k_[3])))
      {
         dt_ = std::max(dt_min_, dt * shrink_min_);
         initialized_ = false;
         total_rejections_++;
         return false;
      }

      // Stage 5
      for (int i = 0; i < n; i++)
      {
         y_tmp_(i) = state(i) + dt * (a51 * k_[0](i) + a52 * k_[1](i)
                                       + a53 * k_[2](i) + a54 * k_[3](i));
      }
      op.SetTime(t + c5 * dt);
      op.Mult(y_tmp_, k_[4]);

      if (!std::isfinite(NormL2(k_[4])))
      {
         dt_ = std::max(dt_min_, dt * shrink_min_);
         initialized_ = false;
         total_rejections_++;
         return false;
      }

      // Stage 6
      for (int i = 0; i < n; i++)
      {
         y_tmp_(i) = state(i) + dt * (a61 * k_[0](i) + a62 * k_[1](i)
                                       + a63 * k_[2](i) + a64 * k_[3](i)
                                       + a65 * k_[4](i));
      }
      op.SetTime(t + dt);
      op.Mult(y_tmp_, k_[5]);

      if (!std::isfinite(NormL2(k_[5])))
      {
         dt_ = std::max(dt_min_, dt * shrink_min_);
         initialized_ = false;
         total_rejections_++;
         return false;
      }

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

      // Compute weighted error norm
      real_t err_norm = 0.0;
      int worst_idx = 0;
      if (use_2norm_)
      {
         // Weighted RMS (2-norm): sqrt(1/N * sum((err_i/scale_i)^2))
         // More robust to outlier DOFs than L-infinity.
         real_t sum_sq = 0.0;
         real_t max_ei = 0.0;
         for (int i = 0; i < n; i++)
         {
            real_t scale = atol_ + rtol_ * std::abs(y_tmp_(i));
            real_t ei = std::abs(err_(i)) / scale;
            sum_sq += ei * ei;
            if (ei > max_ei) { max_ei = ei; worst_idx = i; }
         }
         if (mpi_ctx_)
         {
            // Global sum for 2-norm, global max for worst_idx diagnostic
            int global_n = mpi_ctx_->GlobalSumInt(n);
            real_t global_sum = mpi_ctx_->GlobalSum(sum_sq);
            max_ei = mpi_ctx_->GlobalMax(max_ei);
            err_norm = std::sqrt(global_sum / global_n);
         }
         else
         {
            err_norm = (n > 0) ? std::sqrt(sum_sq / n) : 0.0;
         }
      }
      else
      {
         // Weighted L-infinity (matches Tandem/PETSc default)
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

         // In parallel, reduce err_norm across all ranks so that
         // accept/reject decisions are consistent.
         if (mpi_ctx_)
         {
            err_norm = mpi_ctx_->GlobalMax(err_norm);
         }
      }

      // Guard against NaN in error estimate.
      // NaN can arise from solver failure (e.g., CG divergence) in any RK stage.
      // Force rejection with aggressive dt reduction -- do NOT use the PI
      // controller formula since err_norm is meaningless.
      if (!std::isfinite(err_norm))
      {
         dt_ = std::max(dt_min_, dt * shrink_min_);
         initialized_ = false;  // Discard FSAL k_[0] -- it may contain NaN
         total_rejections_++;
         return false;
      }

      // Compute new dt using standard PI controller formula (PETSc TSAdaptBasic)
      //   dt_new = safety * dt * err_norm^(-1/order), order = 5 for DOPRI5(4)
      //   Matches PETSc: hfac_lte = safety * enorm^(-1/order)
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

      // Verbose diagnostic for every step (enabled by diag_verbose_)
      if (diag_verbose_ && (!mpi_ctx_ || mpi_ctx_->IsRoot()))
      {
         int dof = worst_idx / state_per_node_;
         int comp = worst_idx % state_per_node_;
         const char *comp_name = (comp == state_per_node_ - 1)
                                    ? "(theta/psi)" : "(slip)";
         real_t scale = atol_ + rtol_ * std::abs(y_tmp_(worst_idx));
         std::cout << (err_norm <= 1.0 ? "[ACCEPT]" : "[REJECT]")
                   << " dt=" << std::scientific << std::setprecision(3) << dt
                   << " err=" << err_norm
                   << " dt_new=" << dt_new
                   << " worst=DOF" << dof << comp_name
                   << " |err|=" << std::abs(err_(worst_idx))
                   << " scale=" << scale
                   << " |y|=" << std::abs(y_tmp_(worst_idx))
                   << "\n";
      }

      if (err_norm <= 1.0)
      {
         // Accept step
         state = y_tmp_;
         t += dt;
         dt_ = dt_new;

         // FSAL: k_[6] becomes k_[0] for next step
         k_[0] = k_[6];

         diag_count_ = 0;
         return true;
      }
      else
      {
         // Reject step: apply extra safety factor (PETSc reject_safety)
         dt_ = dt_new * reject_safety_;
         dt_ = std::max(dt_min_, dt_);
         initialized_ = true;  // keep k_[0] from this step start
         total_rejections_++;

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

   /// @brief Check if FSAL stage k_[0] is valid from a previous step.
   bool IsInitialized() const { return initialized_; }

   /// @brief Get the first stage vector (for checkpoint).
   const Vector &GetK0() const { return k_[0]; }

   /// @brief Restore FSAL state from checkpoint.
   ///
   /// Sets k_[0] and marks the integrator as initialized so the next
   /// Step() reuses k_[0] instead of recomputing it.
   void RestoreFSAL(const Vector &k0)
   {
      k_[0] = k0;
      initialized_ = true;
   }

private:
   MPIContext *mpi_ctx_ = nullptr;  ///< MPI context for parallel error reduction

   real_t atol_;           ///< Absolute tolerance
   real_t rtol_;           ///< Relative tolerance
   real_t safety_;         ///< Safety factor for dt adjustment
   real_t reject_safety_;  ///< Extra shrink factor after rejection (PETSc default 0.5)
   real_t growth_max_;     ///< Maximum dt growth factor per step
   real_t shrink_min_;     ///< Minimum dt shrink factor per step
   real_t dt_min_;         ///< Minimum allowed dt
   real_t dt_max_;         ///< Maximum allowed dt
   real_t dt_;             ///< Current time step
   bool initialized_;      ///< Whether k_[0] is valid from a previous step
   bool diag_verbose_ = false; ///< Verbose per-step diagnostics
   bool use_2norm_ = false;    ///< Use RMS (2-norm) instead of L-inf for error
   int total_rejections_;  ///< Total number of rejected steps
   int diag_count_;        ///< Counter for dt_min diagnostic messages
   int state_per_node_ = 2; ///< State components per node (2=BP2, 3=BP5)

   Vector k_[7];  ///< Stage vectors
   Vector y_tmp_; ///< Temporary solution
   Vector err_;   ///< Error estimate

   /// @brief Compute L2 norm, checking for NaN.
   ///
   /// In parallel, a local NaN may not show up in the local Norml2()
   /// because only a subset of DOFs are on each rank. We check locally
   /// and broadcast via GlobalMax if available.
   real_t NormL2(const Vector &v) const
   {
      real_t local_norm = v.Norml2();
      if (mpi_ctx_)
      {
         // If any rank has NaN, propagate it: max(NaN, x) behavior varies,
         // so explicitly check and broadcast.
         int local_nan = std::isfinite(local_norm) ? 0 : 1;
         int global_nan = mpi_ctx_->GlobalSumInt(local_nan);
         if (global_nan > 0) { return std::numeric_limits<real_t>::quiet_NaN(); }
         return mpi_ctx_->GlobalMax(local_norm);
      }
      return local_norm;
   }

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
