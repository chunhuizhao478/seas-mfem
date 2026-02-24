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

#ifndef MFEM_SEAS_OPERATOR_HPP
#define MFEM_SEAS_OPERATOR_HPP

#include "mfem.hpp"
#include "../domain/antiplane_operator.hpp"
#include "../fault/rate_state_fault.hpp"
#include "../common/seas_types.hpp"
#include "../common/mpi_context.hpp"

#include <memory>
#include <cmath>

namespace mfem
{
namespace seas
{

/// @brief SEAS quasi-dynamic time integration operator.
///
/// Couples the domain solver (Phase 2) with the fault operator (Phase 3)
/// for quasi-dynamic earthquake cycle simulation.
///
/// Inherits from TimeDependentOperator for use with MFEM ODE solvers
/// (RK4Solver, etc.).
///
/// State vector layout: [slip_0, theta_0, slip_1, theta_1, ..., slip_{n-1}, theta_{n-1}]
/// Rate vector layout:  [V_0, dtheta_0/dt, V_1, dtheta_1/dt, ..., V_{n-1}, dtheta_{n-1}/dt]
///
/// The coupling flow in each Mult() call:
/// 1. Extract slip from state vector
/// 2. Solve domain problem with slip BC -> displacement u
/// 3. Compute traction on fault from u -> tau_qs
/// 4. Compute fault RHS (slip rate V and dtheta/dt) from stress balance
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
template <typename MeshType = Mesh>
class SEASQuasiDynamicOperator : public TimeDependentOperator
{
public:
   using DomainOpType = AntiplaneDomainOperator<MeshType>;
   using GridFuncType = typename DomainOpType::GridFuncType;

   /// @brief Construct the SEAS quasi-dynamic operator.
   ///
   /// @param domain Domain operator (Phase 2) - owned externally
   /// @param fault Fault operator (Phase 3) - owned externally
   /// @param mpi_ctx MPI context for parallel reductions (optional)
   SEASQuasiDynamicOperator(DomainOpType *domain,
                             RateStateFaultOperator<MeshType> *fault,
                             MPIContext *mpi_ctx = nullptr);

   /// Destructor
   ~SEASQuasiDynamicOperator() override = default;

   /// @brief Initialize state vector at t=0.
   ///
   /// Follows Tandem's two-phase initialization:
   /// 1. PreInit: Set slip=0, theta=placeholder
   /// 2. Solve domain with zero slip to get initial traction
   /// 3. Init: Compute theta from stress equilibrium
   /// 4. Verify initial slip rate matches V_init
   ///
   /// @param[out] state State vector to initialize [2 * num_fault_dofs]
   void SetInitialCondition(Vector &state);

   /// @brief Compute d(state)/dt = RHS(t, state).
   ///
   /// This is the main ODE function called by MFEM ODE solvers (e.g., RK4Solver).
   ///
   /// @param[in] state Current state [slip_0, theta_0, slip_1, theta_1, ...]
   /// @param[out] rate Time derivatives [V_0, dtheta_0/dt, V_1, dtheta_1/dt, ...]
   void Mult(const Vector &state, Vector &rate) const override;

   /// @brief Get current displacement solution.
   const GridFuncType &GetDisplacement() const { return *u_gf_; }

   /// @brief Set displacement from checkpoint data.
   void SetDisplacement(const Vector &u) { *u_gf_ = u; }

   /// @brief Get traction at fault from last evaluation.
   const Vector &GetTraction() const { return traction_; }

   /// @brief Get maximum slip rate from last evaluation (global in parallel).
   real_t GetMaxSlipRate() const
   {
      if (mpi_ctx_)
      {
         return fault_->GetGlobalMaxSlipRate();
      }
      return fault_->GetMaxSlipRate();
   }

   /// @brief Get the domain operator.
   const DomainOpType *GetDomain() const { return domain_; }

   /// @brief Get the fault operator.
   const RateStateFaultOperator<MeshType> *GetFault() const { return fault_; }

private:
   DomainOpType *domain_;
   RateStateFaultOperator<MeshType> *fault_;
   MPIContext *mpi_ctx_ = nullptr;

   /// Displacement grid function (solution of domain problem)
   std::unique_ptr<GridFuncType> u_gf_;

   /// Work vectors (mutable for use in const Mult)
   mutable Vector slip_;
   mutable Vector traction_;
};

// ============================================================================
// Implementation
// ============================================================================

template <typename MeshType>
SEASQuasiDynamicOperator<MeshType>::SEASQuasiDynamicOperator(
   DomainOpType *domain,
   RateStateFaultOperator<MeshType> *fault,
   MPIContext *mpi_ctx)
   : TimeDependentOperator(fault->StateSize()),
     domain_(domain), fault_(fault), mpi_ctx_(mpi_ctx)
{
   MFEM_VERIFY(domain_ != nullptr, "Domain operator must not be null");
   MFEM_VERIFY(fault_ != nullptr, "Fault operator must not be null");

   // Allocate displacement grid function on the domain FE space
   u_gf_ = std::make_unique<GridFuncType>(&domain_->GetFESpace());
   *u_gf_ = 0.0;

   // Allocate work vectors
   slip_.SetSize(fault_->NumNodes());
   traction_.SetSize(fault_->NumNodes());
}

template <typename MeshType>
void SEASQuasiDynamicOperator<MeshType>::SetInitialCondition(Vector &state)
{
   MFEM_VERIFY(state.Size() == fault_->StateSize(),
               "State vector size mismatch: got " << state.Size()
               << ", expected " << fault_->StateSize());

   // Phase 1: Pre-initialize (slip=0, theta=placeholder)
   fault_->PreInit(state);

   // Phase 2: Solve domain with zero slip to get initial traction
   fault_->GetSlip(state, slip_);
   domain_->Solve(0.0, slip_, *u_gf_);
   domain_->ComputeTraction(*u_gf_, slip_, traction_);

   // Phase 3: Initialize theta from stress equilibrium
   //   tau0 + traction = sigma_n * f(V_init, theta) + eta * V_init
   real_t V_max = fault_->Init(traction_, state);

   // Phase 4: Verify initial slip rate
   // Re-solve domain and recompute traction with updated state
   fault_->GetSlip(state, slip_);
   domain_->Solve(0.0, slip_, *u_gf_);
   domain_->ComputeTraction(*u_gf_, slip_, traction_);

   // Compute RHS to populate slip rates
   Vector rate_temp(fault_->StateSize());
   fault_->ComputeRHS(traction_, state, rate_temp);

   // Use global V_max in parallel, local in serial
   V_max = GetMaxSlipRate();

   const BP2Params &params = fault_->GetParams();

   // Verify stress equilibrium (local check, each rank verifies its own DOFs)
   real_t eq_error = fault_->VerifyStressEquilibrium(traction_, state);
   if (mpi_ctx_)
   {
      eq_error = mpi_ctx_->GlobalMax(eq_error);
   }
   MFEM_VERIFY(eq_error < 1e-6,
               "Initial stress equilibrium error too large: " << eq_error);

   // Verify initial slip rate is close to V_init
   // Allow generous tolerance since the domain solve with zero slip
   // gives near-zero traction, so V should be close to V_init
   real_t V_rel_err = std::abs(V_max - params.V_init) /
                      std::max(params.V_init, 1e-30);
   MFEM_VERIFY(V_rel_err < 0.1,
               "Initial V_max = " << V_max
               << " differs from V_init = " << params.V_init
               << " by " << V_rel_err * 100 << "%");
}

template <typename MeshType>
void SEASQuasiDynamicOperator<MeshType>::Mult(
   const Vector &state, Vector &rate) const
{
   // 1. Extract slip from state vector
   fault_->GetSlip(state, slip_);

   // Debug: check for NaN in slip before domain solve
#ifndef NDEBUG
   for (int i = 0; i < slip_.Size(); i++)
   {
      MFEM_VERIFY(std::isfinite(slip_(i)),
         "NaN/Inf in slip at DOF " << i << ": " << slip_(i));
   }
#endif

   // 2. Solve domain problem with slip BC
   //    ∇²u = 0 with [[u]] = slip on fault (all x=0 faces)
   domain_->Solve(t, slip_, *u_gf_);

   // 3. Compute traction at fault from displacement
   //    τ_qs = μ * {{∂u/∂x}} + μ * κ * h⁻¹ * ([[u]] - δ)
   domain_->ComputeTraction(*u_gf_, slip_, traction_);

   // 4. Compute fault RHS (slip rate and state rate)
   //    Stress balance: τ₀ + τ_qs = σ_n * f(V, θ) + η * V
   //    State evolution: dθ/dt = G(V, θ)
   fault_->ComputeRHS(traction_, state, rate);
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_OPERATOR_HPP
