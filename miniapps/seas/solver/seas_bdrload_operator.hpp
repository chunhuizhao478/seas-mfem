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

#ifndef MFEM_SEAS_BDRLOAD_OPERATOR_HPP
#define MFEM_SEAS_BDRLOAD_OPERATOR_HPP

/// @file seas_bdrload_operator.hpp
/// @brief SEAS quasi-dynamic operator with far-field Dirichlet loading.
///
/// This operator couples AntiplaneBdrLoadOperator (domain with far-field
/// Dirichlet BCs) with RateStateFaultOperator (fault) for BP1 simulation.
///
/// Matches Tandem's BP1 boundary condition approach:
/// - Far-field left/right boundaries: Dirichlet u = sign(x) * Vp/2 * t
/// - Free surface and bottom: Natural (zero traction)
/// - Below-Wf fault nodes: prescribed V=Vp (equivalent to Tandem's
///   Dirichlet on the deep extension below the fault)
/// - The domain solve uses the current time `t` for the far-field BC

#include "mfem.hpp"
#include "../domain/antiplane_bdrload_operator.hpp"
#include "../fault/rate_state_fault.hpp"
#include "../common/seas_types.hpp"
#include "../common/mpi_context.hpp"

#include <memory>
#include <cmath>

namespace mfem
{
namespace seas
{

/// @brief Adapter for RateStateFaultOperator used with far-field Dirichlet
/// loading. Matches Tandem's BP1 approach.
///
/// Below-Wf nodes are prescribed at V=Vp (same as base operator), matching
/// Tandem's Dirichlet BC on the deep extension (Physical Curve 5).
/// Above-Wf nodes use rate-and-state friction as usual.
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
template <typename MeshType = Mesh>
class BdrLoadFaultAdapter
{
public:
   static constexpr int StatePerNode = RateStateFaultOperator<MeshType>::StatePerNode;
   static constexpr int SlipIndex = RateStateFaultOperator<MeshType>::SlipIndex;
   static constexpr int ThetaIndex = RateStateFaultOperator<MeshType>::ThetaIndex;

   /// @brief Construct adapter wrapping an existing fault operator.
   ///
   /// @param fault The underlying fault operator (owned externally)
   BdrLoadFaultAdapter(RateStateFaultOperator<MeshType> *fault)
      : fault_(fault),
        num_nodes_(fault->NumNodes()),
        V_max_(0.0)
   {
      MFEM_VERIFY(fault_ != nullptr, "Fault operator must not be null");
      slip_rate_.SetSize(num_nodes_);
      slip_rate_ = fault_->GetParams().V_init;
   }

   int StateSize() const { return fault_->StateSize(); }
   int NumNodes() const { return num_nodes_; }

   /// @brief Pre-initialize: delegates to underlying fault operator.
   void PreInit(Vector &state) { fault_->PreInit(state); }

   /// @brief Initialize: compute initial theta from stress equilibrium.
   ///
   /// Above-Wf nodes: rate-and-state friction initialization.
   /// Below-Wf nodes: prescribed V=Vp (matching Tandem's deep extension
   /// Dirichlet BC).
   real_t Init(const Vector &traction, Vector &state)
   {
      MFEM_ASSERT(traction.Size() == num_nodes_,
                  "Traction vector has wrong size");
      MFEM_ASSERT(state.Size() == StateSize(),
                  "State vector has wrong size");

      const BP2Params &params = fault_->GetParams();
      const FaultGeometry<MeshType> *geom = fault_->GetGeometry();
      const FrictionLaw *friction = fault_->GetFrictionLaw();
      const auto *dr_friction =
         dynamic_cast<const DieterichRuinaFriction*>(friction);
      bool use_psi = fault_->UsePsi();

      // Compute pre-stress (same as base)
      real_t tau0 = params.tau0();

      const Vector &a_values = geom->GetAValues();
      const Vector &eta_values = geom->GetEtaValues();
      const Vector &depths = geom->GetDepths();

      V_max_ = 0.0;

      for (int i = 0; i < num_nodes_; i++)
      {
         if (depths(i) < -params.Wf)
         {
            // Below Wf: prescribed plate rate (matches Tandem's Dirichlet
            // on deep extension below fault)
            slip_rate_(i) = params.Vp;
            continue;
         }

         real_t tau = tau0 + traction(i);
         real_t a = a_values(i);
         real_t eta = eta_values(i);

         if (use_psi)
         {
            real_t psi0 = dr_friction->InitialStatePsi(tau, params.V_init,
                                                        params.sigma_n, eta, a);
            state(i * StatePerNode + ThetaIndex) = psi0;

            real_t V = dr_friction->SolveSlipRatePsi(tau, psi0,
                                                      params.sigma_n, eta, a);
            slip_rate_(i) = V;
            V_max_ = std::max(V_max_, V);
         }
         else
         {
            real_t theta0 = friction->InitialState(tau, params.V_init,
                                                    params.sigma_n, eta, a);
            state(i * StatePerNode + ThetaIndex) = theta0;

            real_t V = friction->SolveSlipRate(tau, theta0,
                                                params.sigma_n, eta, a);
            slip_rate_(i) = V;
            V_max_ = std::max(V_max_, V);
         }
      }

      return V_max_;
   }

   /// @brief Compute RHS: above-Wf nodes use rate-and-state, below-Wf
   /// nodes are prescribed at V=Vp (matching Tandem's deep extension).
   real_t ComputeRHS(const Vector &traction, const Vector &state, Vector &rate)
   {
      MFEM_ASSERT(traction.Size() == num_nodes_,
                  "Traction vector has wrong size");
      MFEM_ASSERT(state.Size() == StateSize(),
                  "State vector has wrong size");
      MFEM_ASSERT(rate.Size() == StateSize(),
                  "Rate vector has wrong size");

      const BP2Params &params = fault_->GetParams();
      const FaultGeometry<MeshType> *geom = fault_->GetGeometry();
      const FrictionLaw *friction = fault_->GetFrictionLaw();
      const StateEvolution *evolution = fault_->GetEvolution();
      const auto *dr_friction =
         dynamic_cast<const DieterichRuinaFriction*>(friction);
      bool use_psi = fault_->UsePsi();

      real_t tau0 = params.tau0();
      const Vector &a_values = geom->GetAValues();
      const Vector &eta_values = geom->GetEtaValues();
      const Vector &depths = geom->GetDepths();

      V_max_ = 0.0;

      for (int i = 0; i < num_nodes_; i++)
      {
         if (depths(i) < -params.Wf)
         {
            // Below Wf: prescribed plate rate, no state evolution
            // (matches Tandem's Dirichlet on deep extension)
            rate(i * StatePerNode + SlipIndex) = params.Vp;
            rate(i * StatePerNode + ThetaIndex) = 0.0;
            slip_rate_(i) = params.Vp;
            continue;
         }

         real_t state_var = state(i * StatePerNode + ThetaIndex);
         real_t tau = tau0 + traction(i);
         real_t a = a_values(i);
         real_t eta = eta_values(i);

         real_t V;
         if (use_psi)
         {
            V = dr_friction->SolveSlipRatePsi(tau, state_var,
                                               params.sigma_n, eta, a);
         }
         else
         {
            V = friction->SolveSlipRate(tau, state_var,
                                         params.sigma_n, eta, a);
         }
         slip_rate_(i) = V;
         V_max_ = std::max(V_max_, V);

         // dslip/dt = V
         rate(i * StatePerNode + SlipIndex) = V;

         // d(state_var)/dt from evolution law
         rate(i * StatePerNode + ThetaIndex) =
            evolution->Rate(V, state_var, params.Dc);
      }

      // Update the underlying fault operator's cached slip rate for
      // GetMaxSlipRate() and GetSlipRate() to work correctly.
      fault_->SetSlipRate(slip_rate_);

      return V_max_;
   }

   // Delegation methods
   void GetSlip(const Vector &state, Vector &slip) const
   {
      fault_->GetSlip(state, slip);
   }

   void GetTheta(const Vector &state, Vector &theta) const
   {
      fault_->GetTheta(state, theta);
   }

   const Vector &GetSlipRate() const { return slip_rate_; }
   real_t GetMaxSlipRate() const { return V_max_; }

   real_t GetGlobalMaxSlipRate() const
   {
      const FaultGeometry<MeshType> *geom = fault_->GetGeometry();
      MPIContext *mpi_ctx = geom->GetMPIContext();
      if (mpi_ctx)
      {
         return mpi_ctx->GlobalMax(V_max_);
      }
      return V_max_;
   }

   real_t GetTau0() const { return fault_->GetParams().tau0(); }
   void InitPreStress() { fault_->InitPreStress(); }

   const FaultGeometry<MeshType> *GetGeometry() const
   {
      return fault_->GetGeometry();
   }

   const BP2Params &GetParams() const { return fault_->GetParams(); }

   RateStateFaultOperator<MeshType> *GetUnderlying() { return fault_; }
   const RateStateFaultOperator<MeshType> *GetUnderlying() const { return fault_; }

   /// @brief Verify stress equilibrium at each node (above-Wf only).
   real_t VerifyStressEquilibrium(const Vector &traction, const Vector &state,
                                   real_t tol = 1e-10) const
   {
      const BP2Params &params = fault_->GetParams();
      const FaultGeometry<MeshType> *geom = fault_->GetGeometry();
      const FrictionLaw *friction = fault_->GetFrictionLaw();
      const auto *dr_friction =
         dynamic_cast<const DieterichRuinaFriction*>(friction);
      bool use_psi = fault_->UsePsi();

      real_t tau0 = params.tau0();
      const Vector &a_values = geom->GetAValues();
      const Vector &eta_values = geom->GetEtaValues();
      const Vector &depths = geom->GetDepths();

      real_t max_rel_error = 0.0;

      for (int i = 0; i < num_nodes_; i++)
      {
         // Skip below-Wf nodes (prescribed V=Vp, no stress equilibrium)
         if (depths(i) < -params.Wf) { continue; }

         real_t state_var = state(i * StatePerNode + ThetaIndex);
         real_t tau = tau0 + traction(i);
         real_t a = a_values(i);
         real_t eta = eta_values(i);

         real_t V, f;
         if (use_psi)
         {
            V = dr_friction->SolveSlipRatePsi(tau, state_var,
                                               params.sigma_n, eta, a);
            f = dr_friction->FrictionCoefficientPsi(V, state_var, a);
         }
         else
         {
            V = friction->SolveSlipRate(tau, state_var,
                                         params.sigma_n, eta, a);
            f = friction->FrictionCoefficient(V, state_var, a);
         }
         real_t tau_computed = params.sigma_n * f + eta * V;

         real_t rel_error = std::abs(tau - tau_computed) /
                            std::max(std::abs(tau), 1.0);
         max_rel_error = std::max(max_rel_error, rel_error);
      }

      return max_rel_error;
   }

   void SetSlipRate(const Vector &V)
   {
      slip_rate_ = V;
      V_max_ = V.Normlinf();
      fault_->SetSlipRate(V);
   }

private:
   RateStateFaultOperator<MeshType> *fault_;
   int num_nodes_;
   real_t V_max_;
   Vector slip_rate_;
};

// ============================================================================
// SEAS Boundary-Load Operator
// ============================================================================

/// @brief SEAS quasi-dynamic operator with far-field Dirichlet loading.
///
/// Same coupling pattern as SEASQuasiDynamicOperator but:
/// - Uses AntiplaneBdrLoadOperator (domain with far-field Dirichlet BC)
/// - Uses BdrLoadFaultAdapter (below-Wf nodes prescribed at V=Vp)
/// - Time `t` matters in Solve() for the far-field BC
///
/// State vector layout: [slip_0, theta_0, slip_1, theta_1, ...]
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
template <typename MeshType = Mesh>
class SEASBdrLoadOperator : public TimeDependentOperator
{
public:
   using DomainOpType = AntiplaneBdrLoadOperator<MeshType>;
   using FaultAdapterType = BdrLoadFaultAdapter<MeshType>;
   using GridFuncType = typename DomainOpType::GridFuncType;

   /// @brief Construct the SEAS boundary-load operator.
   ///
   /// @param domain Domain operator with bottom Dirichlet BC (owned externally)
   /// @param fault_adapter Fault adapter for all-RS-node behavior (owned externally)
   /// @param mpi_ctx MPI context for parallel reductions (optional)
   SEASBdrLoadOperator(DomainOpType *domain,
                       FaultAdapterType *fault_adapter,
                       MPIContext *mpi_ctx = nullptr);

   ~SEASBdrLoadOperator() override = default;

   /// @brief Initialize state vector at t=0.
   void SetInitialCondition(Vector &state);

   /// @brief Compute d(state)/dt = RHS(t, state).
   void Mult(const Vector &state, Vector &rate) const override;

   const GridFuncType &GetDisplacement() const { return *u_gf_; }
   void SetDisplacement(const Vector &u) { *u_gf_ = u; }
   const Vector &GetTraction() const { return traction_; }

   real_t GetMaxSlipRate() const
   {
      if (mpi_ctx_)
      {
         return fault_adapter_->GetGlobalMaxSlipRate();
      }
      return fault_adapter_->GetMaxSlipRate();
   }

   const DomainOpType *GetDomain() const { return domain_; }
   const FaultAdapterType *GetFaultAdapter() const { return fault_adapter_; }

private:
   DomainOpType *domain_;
   FaultAdapterType *fault_adapter_;
   MPIContext *mpi_ctx_ = nullptr;

   std::unique_ptr<GridFuncType> u_gf_;

   mutable Vector slip_;
   mutable Vector traction_;
};

// ============================================================================
// Implementation
// ============================================================================

template <typename MeshType>
SEASBdrLoadOperator<MeshType>::SEASBdrLoadOperator(
   DomainOpType *domain,
   FaultAdapterType *fault_adapter,
   MPIContext *mpi_ctx)
   : TimeDependentOperator(fault_adapter->StateSize()),
     domain_(domain), fault_adapter_(fault_adapter), mpi_ctx_(mpi_ctx)
{
   MFEM_VERIFY(domain_ != nullptr, "Domain operator must not be null");
   MFEM_VERIFY(fault_adapter_ != nullptr, "Fault adapter must not be null");

   u_gf_ = std::make_unique<GridFuncType>(&domain_->GetFESpace());
   *u_gf_ = 0.0;

   slip_.SetSize(fault_adapter_->NumNodes());
   traction_.SetSize(fault_adapter_->NumNodes());
}

template <typename MeshType>
void SEASBdrLoadOperator<MeshType>::SetInitialCondition(Vector &state)
{
   MFEM_VERIFY(state.Size() == fault_adapter_->StateSize(),
               "State vector size mismatch: got " << state.Size()
               << ", expected " << fault_adapter_->StateSize());

   // Phase 1: Pre-initialize (slip=0, theta=placeholder)
   fault_adapter_->PreInit(state);

   // Phase 2: Solve domain with zero slip at t=0
   // At t=0, the bottom Dirichlet BC is zero, so traction should be ~zero
   fault_adapter_->GetSlip(state, slip_);
   domain_->Solve(0.0, slip_, *u_gf_);
   domain_->ComputeTraction(*u_gf_, slip_, traction_);

   // Phase 3: Initialize theta from stress equilibrium
   real_t V_max = fault_adapter_->Init(traction_, state);

   // Phase 4: Verify initial slip rate
   fault_adapter_->GetSlip(state, slip_);
   domain_->Solve(0.0, slip_, *u_gf_);
   domain_->ComputeTraction(*u_gf_, slip_, traction_);

   // Compute RHS to populate slip rates
   Vector rate_temp(fault_adapter_->StateSize());
   fault_adapter_->ComputeRHS(traction_, state, rate_temp);

   V_max = GetMaxSlipRate();

   const BP2Params &params = fault_adapter_->GetParams();

   // Verify stress equilibrium
   real_t eq_error = fault_adapter_->VerifyStressEquilibrium(traction_, state);
   if (mpi_ctx_)
   {
      eq_error = mpi_ctx_->GlobalMax(eq_error);
   }
   MFEM_VERIFY(eq_error < 1e-6,
               "Initial stress equilibrium error too large: " << eq_error);

   // Verify initial slip rate is close to V_init
   real_t V_rel_err = std::abs(V_max - params.V_init) /
                      std::max(params.V_init, 1e-30);
   MFEM_VERIFY(V_rel_err < 0.1,
               "Initial V_max = " << V_max
               << " differs from V_init = " << params.V_init
               << " by " << V_rel_err * 100 << "%");
}

template <typename MeshType>
void SEASBdrLoadOperator<MeshType>::Mult(
   const Vector &state, Vector &rate) const
{
   // 1. Extract slip from state
   fault_adapter_->GetSlip(state, slip_);

#ifndef NDEBUG
   for (int i = 0; i < slip_.Size(); i++)
   {
      MFEM_VERIFY(std::isfinite(slip_(i)),
         "NaN/Inf in slip at DOF " << i << ": " << slip_(i));
   }
#endif

   // 2. Solve domain with slip BC and time-dependent bottom Dirichlet BC
   //    The `t` member (inherited from TimeDependentOperator) is set by ODE solver
   domain_->Solve(t, slip_, *u_gf_);

   // 3. Compute traction at fault
   domain_->ComputeTraction(*u_gf_, slip_, traction_);

   // 4. Compute fault RHS (ALL nodes use rate-and-state)
   fault_adapter_->ComputeRHS(traction_, state, rate);
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BDRLOAD_OPERATOR_HPP
