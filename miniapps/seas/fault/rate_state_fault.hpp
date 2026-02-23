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

#ifndef MFEM_SEAS_RATE_STATE_FAULT_HPP
#define MFEM_SEAS_RATE_STATE_FAULT_HPP

#include "mfem.hpp"
#include "fault_geometry.hpp"
#include "../friction/friction_law.hpp"
#include "../friction/dieterich_ruina.hpp"
#include "../friction/state_evolution.hpp"
#include "../config/bp2_params.hpp"
#include "../common/mpi_context.hpp"

#include <algorithm>
#include <cmath>

namespace mfem
{
namespace seas
{

/// @brief Rate-state fault operator for quasi-dynamic earthquake cycle simulation.
///
/// This class manages the fault state evolution and interfaces with the domain
/// solver. It follows Tandem's RateAndState.h architecture with key methods:
/// - PreInit(): Set initial slip values (before domain solve)
/// - Init(): Compute initial θ from stress equilibrium (after first domain solve)
/// - ComputeRHS(): Compute time derivatives of state variables
///
/// State layout: [slip_0, theta_0, slip_1, theta_1, ...]
/// - Each fault node has 2 state variables: slip and theta (state variable)
/// - slip: accumulated fault slip [m]
/// - theta: state variable from rate-and-state friction law [s]
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
template <typename MeshType = Mesh>
class RateStateFaultOperator
{
public:
   /// State layout constants
   static constexpr int StatePerNode = 2;  ///< slip + theta per node
   static constexpr int SlipIndex = 0;     ///< Index of slip in per-node state
   static constexpr int ThetaIndex = 1;    ///< Index of theta in per-node state

   /// @brief Constructor.
   ///
   /// @param geom Fault geometry (provides depths and depth-dependent params)
   /// @param friction Friction law (DieterichRuinaFriction)
   /// @param evolution State evolution law (AgingLaw, SlipLaw, AgingLawPsi, etc.)
   /// @param params BP2 benchmark parameters
   /// @param mpi_ctx MPI context (nullptr for serial)
   /// @param use_psi If true, integrate in psi-space (logarithmic state variable)
   RateStateFaultOperator(FaultGeometry<MeshType> *geom,
                          FrictionLaw *friction,
                          StateEvolution *evolution,
                          const BP2Params &params,
                          MPIContext *mpi_ctx = nullptr,
                          bool use_psi = false)
      : geom_(geom),
        friction_(friction),
        evolution_(evolution),
        params_(params),
        mpi_ctx_(mpi_ctx),
        num_nodes_(geom ? geom->NumFaultDOFs() : 0),
        tau0_(0.0),
        V_max_(0.0),
        use_psi_(use_psi),
        dr_friction_(dynamic_cast<DieterichRuinaFriction*>(friction))
   {
      if (use_psi_)
      {
         MFEM_ASSERT(dr_friction_ != nullptr,
                     "Psi-space integration requires DieterichRuinaFriction");
      }
      if (num_nodes_ > 0)
      {
         slip_rate_.SetSize(num_nodes_);
         slip_rate_ = params_.V_init;
      }
   }

   // =========================================================================
   // State size information
   // =========================================================================

   /// Total state vector size (2 values per node: slip + theta)
   int StateSize() const { return num_nodes_ * StatePerNode; }

   /// Number of fault nodes
   int NumNodes() const { return num_nodes_; }

   // =========================================================================
   // Methods following Tandem RateAndState.h interface
   // =========================================================================

   /// @brief Pre-initialize: set initial slip values.
   ///
   /// Called before the first domain solve. Sets initial slip to zero
   /// and theta to a placeholder value.
   ///
   /// Following Tandem RateAndState.h pre_init()
   ///
   /// @param[out] state State vector to initialize [StateSize()]
   void PreInit(Vector &state)
   {
      MFEM_ASSERT(state.Size() == StateSize(),
                  "State vector has wrong size");

      for (int i = 0; i < num_nodes_; i++)
      {
         // Initial slip = 0
         state(i * StatePerNode + SlipIndex) = 0.0;

         // Placeholder state (will be computed in Init after first domain solve)
         if (use_psi_)
         {
            // psi_ss = f0 + b*ln(V0/V_init)
            state(i * StatePerNode + ThetaIndex) =
               evolution_->SteadyState(params_.V_init, params_.Dc);
         }
         else
         {
            // theta_ss = Dc/V_init
            state(i * StatePerNode + ThetaIndex) = params_.Dc / params_.V_init;
         }
      }
   }

   /// @brief Initialize: compute initial θ from stress equilibrium.
   ///
   /// Called after the first domain solve with zero slip. Computes the
   /// initial state variable θ such that stress equilibrium is satisfied:
   ///   τ₀ + τ_qs = σ_n · f(V_init, θ) + η · V_init
   ///
   /// Following Tandem RateAndState.h init()
   ///
   /// @param[in] traction Traction from domain solve (τ_qs) [NumNodes()]
   /// @param[in,out] state State vector to update [StateSize()]
   /// @return Maximum initial slip rate
   real_t Init(const Vector &traction, Vector &state)
   {
      MFEM_ASSERT(traction.Size() == num_nodes_,
                  "Traction vector has wrong size");
      MFEM_ASSERT(state.Size() == StateSize(),
                  "State vector has wrong size");

      // Compute pre-stress τ₀ from BP2 parameters
      tau0_ = params_.tau0();

      V_max_ = 0.0;
      const Vector &a_values = geom_->GetAValues();
      const Vector &eta_values = geom_->GetEtaValues();
      const Vector &depths = geom_->GetDepths();

      for (int i = 0; i < num_nodes_; i++)
      {
         if (depths(i) < -params_.Wf)
         {
            // Below Wf: prescribed plate rate, keep placeholder θ
            slip_rate_(i) = params_.Vp;
            continue;
         }

         // Total stress = pre-stress + quasi-static traction
         real_t tau = tau0_ + traction(i);
         real_t a = a_values(i);
         real_t eta = eta_values(i);

         if (use_psi_)
         {
            // Compute initial psi from stress equilibrium
            real_t psi0 = dr_friction_->InitialStatePsi(tau, params_.V_init,
                                                         params_.sigma_n, eta, a);
            state(i * StatePerNode + ThetaIndex) = psi0;

            real_t V = dr_friction_->SolveSlipRatePsi(tau, psi0,
                                                       params_.sigma_n, eta, a);
            slip_rate_(i) = V;
            V_max_ = std::max(V_max_, V);
         }
         else
         {
            // Compute initial θ from stress equilibrium
            real_t theta0 = friction_->InitialState(tau, params_.V_init,
                                                     params_.sigma_n, eta, a);
            state(i * StatePerNode + ThetaIndex) = theta0;

            real_t V = friction_->SolveSlipRate(tau, theta0,
                                                 params_.sigma_n, eta, a);
            slip_rate_(i) = V;
            V_max_ = std::max(V_max_, V);
         }
      }

      return V_max_;
   }

   /// @brief Compute RHS: time derivatives of state variables.
   ///
   /// Called during time stepping. Computes:
   ///   dslip/dt = V (slip rate from stress balance)
   ///   dtheta/dt = G(V, θ) (from state evolution law)
   ///
   /// Following Tandem RateAndState.h rhs()
   ///
   /// @param[in] traction Traction from domain solve (τ_qs) [NumNodes()]
   /// @param[in] state Current state vector [StateSize()]
   /// @param[out] rate Time derivatives of state [StateSize()]
   /// @return Maximum slip rate
   real_t ComputeRHS(const Vector &traction, const Vector &state, Vector &rate)
   {
      MFEM_ASSERT(traction.Size() == num_nodes_,
                  "Traction vector has wrong size");
      MFEM_ASSERT(state.Size() == StateSize(),
                  "State vector has wrong size");
      MFEM_ASSERT(rate.Size() == StateSize(),
                  "Rate vector has wrong size");

      V_max_ = 0.0;
      const Vector &a_values = geom_->GetAValues();
      const Vector &eta_values = geom_->GetEtaValues();
      const Vector &depths = geom_->GetDepths();

      for (int i = 0; i < num_nodes_; i++)
      {
         if (depths(i) < -params_.Wf)
         {
            // Below Wf: prescribed plate rate, no state evolution
            rate(i * StatePerNode + SlipIndex) = params_.Vp;
            rate(i * StatePerNode + ThetaIndex) = 0.0;
            slip_rate_(i) = params_.Vp;
            continue;
         }

         // Get current state variable (theta or psi depending on mode)
         real_t state_var = state(i * StatePerNode + ThetaIndex);

         // Total stress = pre-stress + quasi-static traction
         real_t tau = tau0_ + traction(i);
         real_t a = a_values(i);
         real_t eta = eta_values(i);

         real_t V;
         if (use_psi_)
         {
            V = dr_friction_->SolveSlipRatePsi(tau, state_var,
                                                params_.sigma_n, eta, a);
         }
         else
         {
            V = friction_->SolveSlipRate(tau, state_var,
                                          params_.sigma_n, eta, a);
         }
         slip_rate_(i) = V;
         V_max_ = std::max(V_max_, V);

         // dslip/dt = V
         rate(i * StatePerNode + SlipIndex) = V;

         // d(state_var)/dt from evolution law (works for both theta and psi)
         rate(i * StatePerNode + ThetaIndex) =
            evolution_->Rate(V, state_var, params_.Dc);
      }

      return V_max_;
   }

   // =========================================================================
   // State access methods
   // =========================================================================

   /// @brief Extract slip from state vector.
   ///
   /// @param[in] state Full state vector [StateSize()]
   /// @param[out] slip Slip at each node [NumNodes()]
   void GetSlip(const Vector &state, Vector &slip) const
   {
      MFEM_ASSERT(state.Size() == StateSize(), "State vector has wrong size");
      slip.SetSize(num_nodes_);
      for (int i = 0; i < num_nodes_; i++)
      {
         slip(i) = state(i * StatePerNode + SlipIndex);
      }
   }

   /// @brief Extract theta from state vector.
   ///
   /// In psi-space mode, converts psi -> theta so that I/O output
   /// always produces physical theta values (SCEC-format log10(theta)).
   ///
   /// @param[in] state Full state vector [StateSize()]
   /// @param[out] theta State variable (physical theta) at each node [NumNodes()]
   void GetTheta(const Vector &state, Vector &theta) const
   {
      MFEM_ASSERT(state.Size() == StateSize(), "State vector has wrong size");
      theta.SetSize(num_nodes_);
      if (use_psi_)
      {
         for (int i = 0; i < num_nodes_; i++)
         {
            real_t psi = state(i * StatePerNode + ThetaIndex);
            theta(i) = dr_friction_->PsiToTheta(psi);
         }
      }
      else
      {
         for (int i = 0; i < num_nodes_; i++)
         {
            theta(i) = state(i * StatePerNode + ThetaIndex);
         }
      }
   }

   /// @brief Set slip in state vector.
   ///
   /// @param[in] slip Slip values to set [NumNodes()]
   /// @param[out] state State vector to modify [StateSize()]
   void SetSlip(const Vector &slip, Vector &state)
   {
      MFEM_ASSERT(slip.Size() == num_nodes_, "Slip vector has wrong size");
      MFEM_ASSERT(state.Size() == StateSize(), "State vector has wrong size");
      for (int i = 0; i < num_nodes_; i++)
      {
         state(i * StatePerNode + SlipIndex) = slip(i);
      }
   }

   /// @brief Set theta in state vector.
   ///
   /// @param[in] theta State variable values to set [NumNodes()]
   /// @param[out] state State vector to modify [StateSize()]
   void SetTheta(const Vector &theta, Vector &state)
   {
      MFEM_ASSERT(theta.Size() == num_nodes_, "Theta vector has wrong size");
      MFEM_ASSERT(state.Size() == StateSize(), "State vector has wrong size");
      for (int i = 0; i < num_nodes_; i++)
      {
         state(i * StatePerNode + ThetaIndex) = theta(i);
      }
   }

   // =========================================================================
   // Accessors
   // =========================================================================

   /// Get slip rate from last RHS evaluation.
   const Vector &GetSlipRate() const { return slip_rate_; }

   /// Get maximum slip rate from last RHS evaluation (local).
   real_t GetMaxSlipRate() const { return V_max_; }

   /// Get global maximum slip rate (reduced across all MPI ranks).
   /// In serial mode, returns the same as GetMaxSlipRate().
   real_t GetGlobalMaxSlipRate() const
   {
      if (mpi_ctx_)
      {
         return mpi_ctx_->GlobalMax(V_max_);
      }
      return V_max_;
   }

   /// Get pre-stress τ₀ [Pa].
   real_t GetTau0() const { return tau0_; }

   /// @brief Initialize pre-stress from BP2 parameters (for restart).
   ///
   /// Normally tau0_ is set by Init() during SetInitialCondition().
   /// On restart, we skip SetInitialCondition() but still need tau0_.
   void InitPreStress() { tau0_ = params_.tau0(); }

   /// Get fault geometry.
   const FaultGeometry<MeshType> *GetGeometry() const { return geom_; }

   /// Get friction law.
   const FrictionLaw *GetFrictionLaw() const { return friction_; }

   /// Get state evolution law.
   const StateEvolution *GetEvolution() const { return evolution_; }

   /// Get BP2 parameters.
   const BP2Params &GetParams() const { return params_; }

   // =========================================================================
   // Verification and output
   // =========================================================================

   /// @brief Verify stress equilibrium at each node.
   ///
   /// Checks that: τ = σ_n · f(V, θ) + η · V for all nodes.
   ///
   /// @param[in] traction Traction from domain solve [NumNodes()]
   /// @param[in] state Current state vector [StateSize()]
   /// @param[in] tol Relative tolerance for equilibrium check
   /// @return Maximum relative error in stress balance
   real_t VerifyStressEquilibrium(const Vector &traction, const Vector &state,
                                  real_t tol = 1e-10) const
   {
      const Vector &a_values = geom_->GetAValues();
      const Vector &eta_values = geom_->GetEtaValues();
      const Vector &depths = geom_->GetDepths();

      real_t max_rel_error = 0.0;

      for (int i = 0; i < num_nodes_; i++)
      {
         // Skip below-Wf DOFs (prescribed loading, no friction law)
         if (depths(i) < -params_.Wf) { continue; }

         real_t state_var = state(i * StatePerNode + ThetaIndex);
         real_t tau = tau0_ + traction(i);
         real_t a = a_values(i);
         real_t eta = eta_values(i);

         real_t V, f;
         if (use_psi_)
         {
            V = dr_friction_->SolveSlipRatePsi(tau, state_var,
                                                params_.sigma_n, eta, a);
            f = dr_friction_->FrictionCoefficientPsi(V, state_var, a);
         }
         else
         {
            V = friction_->SolveSlipRate(tau, state_var,
                                          params_.sigma_n, eta, a);
            f = friction_->FrictionCoefficient(V, state_var, a);
         }
         real_t tau_computed = params_.sigma_n * f + eta * V;

         real_t rel_error = std::abs(tau - tau_computed) /
                            std::max(std::abs(tau), 1.0);
         max_rel_error = std::max(max_rel_error, rel_error);
      }

      return max_rel_error;
   }

   /// @brief Print fault state information.
   void PrintState(const Vector &state, std::ostream &os = mfem::out) const
   {
      os << "Fault State Summary:\n";
      os << "  tau0 = " << tau0_ / 1e6 << " MPa\n";
      os << "  V_max = " << V_max_ << " m/s\n";

      if (num_nodes_ > 0)
      {
         // Extract slip and theta
         Vector slip, theta;
         GetSlip(state, slip);
         GetTheta(state, theta);

         os << "  Slip range: [" << slip.Min() << ", " << slip.Max() << "] m\n";
         os << "  Theta range: [" << theta.Min() << ", " << theta.Max() << "] s\n";
         os << "  Slip rate range: [" << slip_rate_.Min() << ", "
            << slip_rate_.Max() << "] m/s\n";
      }
   }

   /// @brief Set cached slip rate from checkpoint data.
   void SetSlipRate(const Vector &V)
   {
      slip_rate_ = V;
      V_max_ = V.Normlinf();
   }

   /// Whether psi-space integration is active.
   bool UsePsi() const { return use_psi_; }

private:
   FaultGeometry<MeshType> *geom_;
   FrictionLaw *friction_;
   StateEvolution *evolution_;
   BP2Params params_;
   MPIContext *mpi_ctx_ = nullptr;

   int num_nodes_;      ///< Number of fault DOFs
   real_t tau0_;        ///< Pre-stress [Pa]
   real_t V_max_;       ///< Maximum slip rate from last evaluation

   bool use_psi_ = false;  ///< If true, state variable is psi instead of theta
   DieterichRuinaFriction *dr_friction_ = nullptr;  ///< Downcast for psi methods

   Vector slip_rate_;   ///< Cached slip rate from last RHS [NumNodes()]
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_RATE_STATE_FAULT_HPP
