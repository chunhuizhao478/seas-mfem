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
#include "../domain/domain_operator.hpp"
#include "../domain/antiplane_operator.hpp"
#include "../domain/elasticity_operator.hpp"
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
/// The coupling flow in each Mult() call:
/// 1. Extract slip from state vector
/// 2. Solve domain problem with slip BC -> displacement u
/// 3. Compute traction on fault from u -> tau_qs
/// 4. Compute fault RHS (slip rate V and dtheta/dt) from stress balance
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
/// @tparam DomainOpType Domain operator type (default: AntiplaneDomainOperator)
/// @tparam FaultOpType Fault operator type (default: RateStateFaultOperator<MeshType>)
template <typename MeshType = Mesh,
          typename DomainOpType = AntiplaneDomainOperator<MeshType>,
          typename FaultOpType = RateStateFaultOperator<MeshType>>
class SEASQuasiDynamicOperator : public TimeDependentOperator
{
public:
   using GridFuncType = typename DomainOpType::GridFuncType;

   /// @brief Construct the SEAS quasi-dynamic operator.
   ///
   /// @param domain Domain operator (Phase 2) - owned externally
   /// @param fault Fault operator (Phase 3) - owned externally
   /// @param mpi_ctx MPI context for parallel reductions (optional)
   SEASQuasiDynamicOperator(DomainOpType *domain,
                             FaultOpType *fault,
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
   /// @param[out] state State vector to initialize [fault->StateSize()]
   void SetInitialCondition(Vector &state);

   /// @brief Compute d(state)/dt = RHS(t, state).
   ///
   /// This is the main ODE function called by MFEM ODE solvers (e.g., RK4Solver).
   ///
   /// @param[in] state Current state [StateSize()]
   /// @param[out] rate Time derivatives [StateSize()]
   void Mult(const Vector &state, Vector &rate) const override;

   /// @brief PETSc TS explicit RHS callback.
   ///
   /// PetscODESolver routes explicit RHS evaluations through ExplicitMult().
   /// This operator is explicit and already defines its RHS in Mult(), so
   /// forward both code paths to the same implementation.
   void ExplicitMult(const Vector &state, Vector &rate) const override
   {
      Mult(state, rate);
   }

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
   const FaultOpType *GetFault() const { return fault_; }

   /// v51: Zero dip traction component after ComputeTraction, before friction law.
   /// Tests whether the 21% cross-component contamination in {σ·n} causes the dip offset.
   void SetZeroDipTraction(bool v) { zero_dip_traction_ = v; }

   /// v51: Dump dip/strike traction ratio during coseismic (V_max > threshold).
   void SetDiagCoseismicDip(bool v, real_t v_threshold = 0.1)
   {
      diag_coseismic_dip_ = v;
      coseismic_v_threshold_ = v_threshold;
   }

   /// v51: Use elastic normal stress (sigma_n from displacement field)
   /// instead of constant sigma_n. Matches Tandem's DieterichRuinaAgeing.
   void SetElasticSigmaN(bool v) { elastic_sigma_n_ = v; }

private:
   DomainOpType *domain_;
   FaultOpType *fault_;
   MPIContext *mpi_ctx_ = nullptr;

   /// Displacement grid function (solution of domain problem)
   std::unique_ptr<GridFuncType> u_gf_;

   /// Work vectors (mutable for use in const Mult)
   mutable Vector slip_;
   mutable Vector traction_;
   mutable Vector normal_traction_;  // v51: elastic T_n for sigma_n feedback

   // v51 flags
   bool zero_dip_traction_ = false;
   // v54: default ON to match Tandem's DieterichRuinaBase.h:87
   // (sigma_n = -sn_elastic + SnPre). Previously off by default (v51).
   bool elastic_sigma_n_ = true;
   bool diag_coseismic_dip_ = false;
   mutable bool diag_coseismic_dip_done_ = false;
   real_t coseismic_v_threshold_ = 0.1;
};

// ============================================================================
// Implementation
// ============================================================================

template <typename MeshType, typename DomainOpType, typename FaultOpType>
SEASQuasiDynamicOperator<MeshType, DomainOpType, FaultOpType>::SEASQuasiDynamicOperator(
   DomainOpType *domain,
   FaultOpType *fault,
   MPIContext *mpi_ctx)
   : TimeDependentOperator(fault->StateSize()),
     domain_(domain), fault_(fault), mpi_ctx_(mpi_ctx)
{
   MFEM_VERIFY(domain_ != nullptr, "Domain operator must not be null");
   MFEM_VERIFY(fault_ != nullptr, "Fault operator must not be null");

   // Allocate displacement grid function on the domain FE space
   u_gf_ = std::make_unique<GridFuncType>(&domain_->GetFESpace());
   *u_gf_ = 0.0;

   // Allocate work vectors (sized for slip/traction components)
   slip_.SetSize(fault_->SlipSize());
   traction_.SetSize(fault_->TractionSize());
}

template <typename MeshType, typename DomainOpType, typename FaultOpType>
void SEASQuasiDynamicOperator<MeshType, DomainOpType, FaultOpType>::SetInitialCondition(Vector &state)
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
   // V_max is recomputed below after the verification re-solve
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

   // Verify stress equilibrium (local check, each rank verifies its own DOFs)
   real_t eq_error = fault_->VerifyStressEquilibrium(traction_, state);
   if (mpi_ctx_)
   {
      eq_error = mpi_ctx_->GlobalMax(eq_error);
   }
   MFEM_VERIFY(eq_error < 1e-6,
               "Initial stress equilibrium error too large: " << eq_error);

   // Log initial slip rate comparison (use global V_ref in parallel).
   // For BP5-QD, the nucleation zone has δτ overstress so V_max > V_nuc
   // is expected and correct — do not assert.
   {
      real_t V_ref = fault_->GetReferenceVInit();
      if (mpi_ctx_)
      {
         V_ref = mpi_ctx_->GlobalMax(V_ref);
      }
      if (mpi_ctx_ == nullptr || mpi_ctx_->IsRoot())
      {
         mfem::out << "  Initial V_max = " << V_max
                   << ", V_ref (max V_init) = " << V_ref << "\n";
      }
   }
}

template <typename MeshType, typename DomainOpType, typename FaultOpType>
void SEASQuasiDynamicOperator<MeshType, DomainOpType, FaultOpType>::Mult(
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
   domain_->Solve(t, slip_, *u_gf_);

   // 3. Compute traction at fault from displacement
   // v51: optionally compute elastic normal traction for sigma_n feedback
   domain_->ComputeTraction(*u_gf_, slip_, traction_,
                             elastic_sigma_n_ ? &normal_traction_ : nullptr);

   // v51: Zero dip traction component (index 0 of each DOF's [dip, strike] pair)
   if (zero_dip_traction_)
   {
      for (int i = 0; i < traction_.Size() / 2; i++)
      {
         traction_(2 * i) = 0.0;  // tau_dip = 0
      }
   }

   // v51: Dump dip/strike traction ratio during coseismic phase
   if (diag_coseismic_dip_ && !diag_coseismic_dip_done_)
   {
      // Check if V_max exceeds coseismic threshold
      real_t v_max = fault_->GetMaxSlipRate();
      if (mpi_ctx_) { v_max = fault_->GetGlobalMaxSlipRate(); }
      if (v_max > coseismic_v_threshold_)
      {
         int n_dofs = traction_.Size() / 2;
         real_t sum_ratio = 0.0;
         real_t max_ratio = 0.0;
         int count = 0;
         real_t max_tau_dip = 0.0, max_tau_strike = 0.0;
         int max_ratio_dof = -1;

         for (int i = 0; i < n_dofs; i++)
         {
            real_t td = std::abs(traction_(2*i));      // |tau_dip|
            real_t ts = std::abs(traction_(2*i + 1));   // |tau_strike|
            if (ts > 1e3)  // only count DOFs with significant strike traction (> 1 kPa)
            {
               real_t ratio = td / ts;
               sum_ratio += ratio;
               count++;
               if (ratio > max_ratio)
               {
                  max_ratio = ratio;
                  max_ratio_dof = i;
                  max_tau_dip = traction_(2*i);
                  max_tau_strike = traction_(2*i + 1);
               }
            }
         }

         if (count > 0)
         {
            bool is_root = !mpi_ctx_ || mpi_ctx_->IsRoot();
            if (is_root)
            {
               mfem::out << "[COSEISMIC-DIP] V_max=" << v_max
                  << " n_active=" << count
                  << " mean_|td/ts|=" << sum_ratio / count
                  << " max_|td/ts|=" << max_ratio
                  << " @DOF=" << max_ratio_dof
                  << " tau_dip=" << max_tau_dip
                  << " tau_strike=" << max_tau_strike
                  << std::endl;
            }
            diag_coseismic_dip_done_ = true;
         }
      }
   }

   // 4. Compute fault RHS (slip rate and state rate)
   // v51: pass elastic normal traction for sigma_n feedback (nullptr = use constant)
   fault_->ComputeRHS(traction_, state, rate,
                       elastic_sigma_n_ ? &normal_traction_ : nullptr);
}

// BP2 type alias (uses default template arguments)
using BP2SEASOp = SEASQuasiDynamicOperator<Mesh>;

// BP5 type aliases
using BP5DomainOp = ElasticityDomainOperator<Mesh>;
using BP5SEASOp   = SEASQuasiDynamicOperator<Mesh, BP5DomainOp, BP5FaultOp>;

#ifdef MFEM_USE_MPI
// Parallel BP5 type alias
using PBP5SEASOp = SEASQuasiDynamicOperator<ParMesh,
                      ElasticityDomainOperator<ParMesh>,
                      RateStateFaultOperator<ParMesh, 2>>;
#endif

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_OPERATOR_HPP
