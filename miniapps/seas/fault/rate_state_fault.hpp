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
#include "../config/bp5_params.hpp"
#include "../common/mpi_context.hpp"

#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>

namespace mfem
{
namespace seas
{

/// @brief Rate-state fault operator for quasi-dynamic earthquake cycle simulation.
///
/// This class manages the fault state evolution and interfaces with the domain
/// solver. It follows Tandem's RateAndState.h architecture with key methods:
/// - PreInit(): Set initial slip values (before domain solve)
/// - Init(): Compute initial theta/psi from stress equilibrium (after first domain solve)
/// - ComputeRHS(): Compute time derivatives of state variables
///
/// State layout depends on SlipComponents:
///   SlipComponents=1 (BP1/BP2): [slip_0, theta_0, slip_1, theta_1, ...]
///   SlipComponents=2 (BP5):     [s_dip_0, s_strike_0, psi_0, s_dip_1, s_strike_1, psi_1, ...]
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
/// @tparam SlipComponents Number of slip components (1 for antiplane, 2 for 3D)
template <typename MeshType = Mesh, int SlipComponents = 1>
class RateStateFaultOperator
{
public:
   /// State layout constants
   static constexpr int NumSlipComp = SlipComponents;
   static constexpr int StatePerNode = SlipComponents + 1;
   static constexpr int SlipIndex = 0;     ///< First slip component index
   static constexpr int ThetaIndex = SlipComponents; ///< Theta/psi index (last)
   static constexpr int PsiIndex = SlipComponents;   ///< Alias for ThetaIndex

   // =========================================================================
   // BP2 Constructor (SlipComponents=1)
   // =========================================================================

   /// @brief Constructor for antiplane/BP2 (scalar slip).
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
      static_assert(SlipComponents == 1,
                    "BP2Params constructor requires SlipComponents=1");
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
   // BP5 Constructor (SlipComponents=2)
   // =========================================================================

   /// @brief Constructor for 3D elasticity/BP5 (vector slip).
   ///
   /// BP5 always uses psi-space integration.
   ///
   /// @param geom Fault geometry (provides 2D coordinates and spatially varying params)
   /// @param friction DieterichRuina friction law (required for vector solver)
   /// @param evolution State evolution law
   /// @param params BP5 benchmark parameters
   /// @param mpi_ctx MPI context (nullptr for serial)
   RateStateFaultOperator(FaultGeometry<MeshType> *geom,
                          DieterichRuinaFriction *friction,
                          StateEvolution *evolution,
                          const BP5Params &params,
                          MPIContext *mpi_ctx = nullptr)
      : geom_(geom),
        friction_(friction),
        evolution_(evolution),
        mpi_ctx_(mpi_ctx),
        num_nodes_(geom ? geom->NumFaultDOFs() : 0),
        tau0_(0.0),
        V_max_(0.0),
        use_psi_(true),  // BP5 always uses psi
        dr_friction_(friction),
        bp5_params_(params),
        sigma_n_bp5_(params.sigma_n),
        Vp_bp5_(params.Vp),
        Wf_bp5_(params.Wf)
   {
      static_assert(SlipComponents == 2,
                    "BP5Params constructor requires SlipComponents=2");
      MFEM_ASSERT(dr_friction_ != nullptr,
                  "BP5 requires DieterichRuinaFriction");
      MFEM_ASSERT(geom_ != nullptr && geom_->IsBP5(),
                  "BP5 constructor requires BP5 FaultGeometry");

      if (num_nodes_ > 0)
      {
         // Slip rate: 2 components per node
         slip_rate_.SetSize(2 * num_nodes_);

         // Initialize with V_init from geometry
         const Vector &V_init = geom_->GetVInit();
         MFEM_ASSERT(V_init.Size() == 2 * num_nodes_,
                     "V_init size mismatch");
         for (int i = 0; i < 2 * num_nodes_; i++)
         {
            slip_rate_(i) = V_init(i);
         }

         // Cache per-DOF parameters from FaultGeometry
         Dc_values_ = geom_->GetDcValues();
         tau_pre_ = geom_->GetTauPre();
         V_init_values_ = geom_->GetVInit();
      }
   }

   // =========================================================================
   // State size information
   // =========================================================================

   /// Total state vector size (StatePerNode values per node)
   int StateSize() const { return num_nodes_ * StatePerNode; }

   /// Total slip vector size (SlipComponents per node)
   int SlipSize() const { return num_nodes_ * SlipComponents; }

   /// Total traction vector size (SlipComponents per node)
   int TractionSize() const { return num_nodes_ * SlipComponents; }

   /// Number of fault nodes
   int NumNodes() const { return num_nodes_; }

   // =========================================================================
   // Methods following Tandem RateAndState.h interface
   // =========================================================================

   /// @brief Pre-initialize: set initial slip values.
   ///
   /// Called before the first domain solve. Sets initial slip to zero
   /// and theta/psi to a placeholder value.
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
         // Initial slip = 0 (all components)
         for (int c = 0; c < SlipComponents; c++)
         {
            state(i * StatePerNode + c) = 0.0;
         }

         // Placeholder state (will be computed in Init after first domain solve)
         if constexpr (SlipComponents == 1)
         {
            if (use_psi_)
            {
               // psi_ss = f0 + b*ln(V0/V_init)
               state(i * StatePerNode + PsiIndex) =
                  evolution_->SteadyState(params_.V_init, params_.Dc);
            }
            else
            {
               // theta_ss = Dc/V_init
               state(i * StatePerNode + ThetaIndex) = params_.Dc / params_.V_init;
            }
         }
         else
         {
            // BP5: always psi-space
            real_t V_abs_init = std::sqrt(
               V_init_values_(2*i) * V_init_values_(2*i) +
               V_init_values_(2*i+1) * V_init_values_(2*i+1));
            real_t Dc = Dc_values_(i);
            state(i * StatePerNode + PsiIndex) =
               evolution_->SteadyState(std::max(V_abs_init, 1e-30), Dc);
         }
      }
   }

   /// @brief Initialize: compute initial theta/psi from stress equilibrium.
   ///
   /// Called after the first domain solve with zero slip. Computes the
   /// initial state variable such that stress equilibrium is satisfied.
   ///
   /// Following Tandem RateAndState.h init()
   ///
   /// @param[in] traction Traction from domain solve [TractionSize()]
   /// @param[in,out] state State vector to update [StateSize()]
   /// @return Maximum initial slip rate
   real_t Init(const Vector &traction, Vector &state)
   {
      MFEM_ASSERT(traction.Size() == TractionSize(),
                  "Traction vector has wrong size");
      MFEM_ASSERT(state.Size() == StateSize(),
                  "State vector has wrong size");

      V_max_ = 0.0;
      const Vector &a_values = geom_->GetAValues();
      const Vector &eta_values = geom_->GetEtaValues();
      const Vector &depths = geom_->GetDepths();

      if constexpr (SlipComponents == 1)
      {
         // ---- BP2 scalar path (unchanged) ----
         tau0_ = params_.tau0();

         for (int i = 0; i < num_nodes_; i++)
         {
            if (depths(i) < -params_.Wf)
            {
               // Below Wf: prescribed plate rate, keep placeholder theta
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
               // Compute initial theta from stress equilibrium
               real_t theta0 = friction_->InitialState(tau, params_.V_init,
                                                        params_.sigma_n, eta, a);
               state(i * StatePerNode + ThetaIndex) = theta0;

               real_t V = friction_->SolveSlipRate(tau, theta0,
                                                    params_.sigma_n, eta, a);
               slip_rate_(i) = V;
               V_max_ = std::max(V_max_, V);
            }
         }
      }
      else
      {
         // ---- BP5 vector path ----
         for (int i = 0; i < num_nodes_; i++)
         {
            // Below fault zone check (shouldn't happen with proper fault detection)
            if (depths(i) > Wf_bp5_ + 1.0)
            {
               slip_rate_(2*i) = 0.0;
               slip_rate_(2*i+1) = Vp_bp5_;
               continue;
            }

            // Vector stress: tau_pre + elastic traction
            real_t tau_vec[2] = {tau_pre_(2*i) + traction(2*i),
                                 tau_pre_(2*i+1) + traction(2*i+1)};
            real_t tau_abs = std::sqrt(tau_vec[0]*tau_vec[0] +
                                       tau_vec[1]*tau_vec[1]);
            real_t V_abs_init = std::sqrt(
               V_init_values_(2*i) * V_init_values_(2*i) +
               V_init_values_(2*i+1) * V_init_values_(2*i+1));

            real_t a = a_values(i);
            real_t eta = eta_values(i);

            real_t psi0;
            if (scec_psi_init_)
            {
               // SCEC Eq. 18: psi(0) = f0 + b*ln(V0/V_init) everywhere
               // delta_tau is genuine overstress, not absorbed into state
               psi0 = bp5_params_.psi_init();
            }
            else
            {
               // Tandem-style: absorb delta_tau into psi via InitialStatePsi
               // System starts in equilibrium at V = V_init (no immediate earthquake)
               psi0 = dr_friction_->InitialStatePsi(
                  tau_abs, V_abs_init, sigma_n_bp5_, eta, a);
            }
            state(i * StatePerNode + PsiIndex) = psi0;

            // Verify by solving vector equation
            real_t V_vec[2];
            dr_friction_->SolveSlipRateVectorPsi(
               tau_vec, psi0, sigma_n_bp5_, eta, a, V_vec);
            slip_rate_(2*i) = V_vec[0];
            slip_rate_(2*i+1) = V_vec[1];
            real_t V_abs = std::sqrt(V_vec[0]*V_vec[0] + V_vec[1]*V_vec[1]);
            V_max_ = std::max(V_max_, V_abs);
         }
      }

      return V_max_;
   }

   /// @brief Compute RHS: time derivatives of state variables.
   ///
   /// Called during time stepping. Computes slip rates and state evolution.
   ///
   /// Following Tandem RateAndState.h rhs()
   ///
   /// @param[in] traction Traction from domain solve [TractionSize()]
   /// @param[in] state Current state vector [StateSize()]
   /// @param[out] rate Time derivatives of state [StateSize()]
   /// @return Maximum slip rate
   real_t ComputeRHS(const Vector &traction, const Vector &state, Vector &rate,
                     const Vector *normal_traction = nullptr)
   {
      MFEM_ASSERT(traction.Size() == TractionSize(),
                  "Traction vector has wrong size");
      MFEM_ASSERT(state.Size() == StateSize(),
                  "State vector has wrong size");
      MFEM_ASSERT(rate.Size() == StateSize(),
                  "Rate vector has wrong size");
      if (normal_traction)
      {
         MFEM_ASSERT(normal_traction->Size() == num_nodes_,
                     "Normal traction vector has wrong size: "
                     << normal_traction->Size() << " vs " << num_nodes_);
      }

      V_max_ = 0.0;
      const Vector &a_values = geom_->GetAValues();
      const Vector &eta_values = geom_->GetEtaValues();
      const Vector &depths = geom_->GetDepths();

      for (int i = 0; i < num_nodes_; i++)
      {
         if constexpr (SlipComponents == 1)
         {
            // ---- BP2 scalar path (unchanged) ----
            if (depths(i) < -params_.Wf)
            {
               // Below Wf: prescribed plate rate, no state evolution
               rate(i * StatePerNode + SlipIndex) = params_.Vp;
               rate(i * StatePerNode + ThetaIndex) = 0.0;
               slip_rate_(i) = params_.Vp;
               continue;
            }

            real_t state_var = state(i * StatePerNode + ThetaIndex);
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

            rate(i * StatePerNode + SlipIndex) = V;
            rate(i * StatePerNode + ThetaIndex) =
               evolution_->Rate(V, state_var, params_.Dc);
         }
         else
         {
            // ---- BP5 vector path ----
            if (depths(i) > Wf_bp5_ + 1.0)
            {
               rate(i * StatePerNode + 0) = 0.0;
               rate(i * StatePerNode + 1) = Vp_bp5_;
               rate(i * StatePerNode + PsiIndex) = 0.0;
               slip_rate_(2*i) = 0.0;
               slip_rate_(2*i+1) = Vp_bp5_;
               continue;
            }

            real_t psi = state(i * StatePerNode + PsiIndex);
            real_t tau_vec[2] = {tau_pre_(2*i) + traction(2*i),
                                 tau_pre_(2*i+1) + traction(2*i+1)};
            real_t a = a_values(i);
            real_t eta = eta_values(i);
            real_t Dc = Dc_values_(i);

            // v51: elastic sigma_n feedback (matches Tandem DieterichRuinaAgeing.h:86)
            // sigma_n_eff = sigma_n_pre + T_n_elastic
            // NormalStress returns positive for compression → sigma_n stays near 25 MPa
            real_t sigma_n_eff = sigma_n_bp5_;
            if (normal_traction)
            {
               sigma_n_eff = sigma_n_bp5_ + (*normal_traction)(i);
               // Safety: ensure sigma_n stays positive (physical requirement)
               sigma_n_eff = std::max(sigma_n_eff, 0.1 * sigma_n_bp5_);
            }

            real_t V_vec[2];
            dr_friction_->SolveSlipRateVectorPsi(
               tau_vec, psi, sigma_n_eff, eta, a, V_vec);

            real_t V_abs = std::sqrt(V_vec[0]*V_vec[0] + V_vec[1]*V_vec[1]);

            rate(i * StatePerNode + 0) = V_vec[0];
            rate(i * StatePerNode + 1) = V_vec[1];
            rate(i * StatePerNode + PsiIndex) =
               evolution_->Rate(V_abs, psi, Dc);

            slip_rate_(2*i) = V_vec[0];
            slip_rate_(2*i+1) = V_vec[1];
            V_max_ = std::max(V_max_, V_abs);
         }
      }

      // Traction monitoring (BP5 only)
      if constexpr (SlipComponents == 2)
      {
         if (monitor_interval_ > 0)
         {
            monitor_call_count_++;
            if (monitor_call_count_ % monitor_interval_ == 0)
            {
               // Auto-pick stations if not specified
               std::vector<int> stations = monitor_stations_;
               if (stations.empty() && num_nodes_ > 0)
               {
                  stations.push_back(0);
                  if (num_nodes_ > 1) { stations.push_back(num_nodes_ / 2); }
                  if (num_nodes_ > 2) { stations.push_back(num_nodes_ - 1); }
               }
               for (int idx : stations)
               {
                  if (idx < 0 || idx >= num_nodes_) { continue; }
                  real_t tp_d = tau_pre_(2*idx);
                  real_t tp_s = tau_pre_(2*idx+1);
                  real_t tr_d = traction(2*idx);
                  real_t tr_s = traction(2*idx+1);
                  real_t tot_d = tp_d + tr_d;
                  real_t tot_s = tp_s + tr_s;
                  std::cout << "[TRACTION] call=" << monitor_call_count_
                            << " node=" << idx
                            << " tau_pre=(" << tp_d << "," << tp_s << ")"
                            << " traction=(" << tr_d << "," << tr_s << ")"
                            << " total=(" << tot_d << "," << tot_s << ")"
                            << " |total|=" << std::sqrt(tot_d*tot_d + tot_s*tot_s)
                            << " V=" << std::sqrt(
                                 slip_rate_(2*idx)*slip_rate_(2*idx) +
                                 slip_rate_(2*idx+1)*slip_rate_(2*idx+1))
                            << "\n";
               }
            }
         }
      }

      return V_max_;
   }

   // =========================================================================
   // State access methods
   // =========================================================================

   /// @brief Extract slip from state vector.
   ///
   /// @param[in] state Full state vector [StateSize()]
   /// @param[out] slip Slip at each node [SlipSize()]
   void GetSlip(const Vector &state, Vector &slip) const
   {
      MFEM_ASSERT(state.Size() == StateSize(), "State vector has wrong size");
      slip.SetSize(SlipSize());
      for (int i = 0; i < num_nodes_; i++)
      {
         for (int c = 0; c < SlipComponents; c++)
         {
            slip(i * SlipComponents + c) = state(i * StatePerNode + c);
         }
      }
   }

   /// @brief Extract theta from state vector.
   ///
   /// In psi-space mode, converts psi -> theta so that I/O output
   /// always produces physical theta values.
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
            real_t psi = state(i * StatePerNode + PsiIndex);
            if constexpr (SlipComponents == 2)
            {
               // BP5: use per-DOF Dc for correct theta conversion
               real_t Dc = Dc_values_(i);
               theta(i) = dr_friction_->PsiToTheta(psi, Dc);
            }
            else
            {
               theta(i) = dr_friction_->PsiToTheta(psi);
            }
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
   /// @param[in] slip Slip values to set [SlipSize()]
   /// @param[out] state State vector to modify [StateSize()]
   void SetSlip(const Vector &slip, Vector &state)
   {
      MFEM_ASSERT(slip.Size() == SlipSize(), "Slip vector has wrong size");
      MFEM_ASSERT(state.Size() == StateSize(), "State vector has wrong size");
      for (int i = 0; i < num_nodes_; i++)
      {
         for (int c = 0; c < SlipComponents; c++)
         {
            state(i * StatePerNode + c) = slip(i * SlipComponents + c);
         }
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
   /// Size: SlipComponents * NumNodes()
   const Vector &GetSlipRate() const { return slip_rate_; }

   /// Get maximum slip rate from last RHS evaluation (local).
   /// For vector case, this is max ||V_i||.
   real_t GetMaxSlipRate() const { return V_max_; }

   /// Get global maximum slip rate (reduced across all MPI ranks).
   real_t GetGlobalMaxSlipRate() const
   {
      if (mpi_ctx_)
      {
         return mpi_ctx_->GlobalMax(V_max_);
      }
      return V_max_;
   }

   /// Get pre-stress tau0 [Pa] (BP2 only).
   real_t GetTau0() const { return tau0_; }

   /// @brief Initialize pre-stress from BP2 parameters (for restart).
   void InitPreStress()
   {
      if constexpr (SlipComponents == 1)
      {
         tau0_ = params_.tau0();
      }
      // BP5: tau_pre_ is set in constructor, no action needed
   }

   /// Get fault geometry.
   const FaultGeometry<MeshType> *GetGeometry() const { return geom_; }

   /// Get friction law.
   const FrictionLaw *GetFrictionLaw() const { return friction_; }

   /// Get state evolution law.
   const StateEvolution *GetEvolution() const { return evolution_; }

   /// Get BP2 parameters (SlipComponents=1 only).
   const BP2Params &GetParams() const
   {
      static_assert(SlipComponents == 1,
                    "GetParams() is only valid for SlipComponents=1. "
                    "Use GetBP5Params() for SlipComponents=2.");
      return params_;
   }

   /// Get BP5 parameters (SlipComponents=2 only).
   const BP5Params &GetBP5Params() const { return bp5_params_; }

   /// Get a reference V_init value for verification (scalar).
   /// For BP5, returns max |V_init| across all DOFs (accounts for nucleation zone).
   real_t GetReferenceVInit() const
   {
      if constexpr (SlipComponents == 1)
      {
         return params_.V_init;
      }
      else
      {
         real_t V_ref = 0.0;
         for (int i = 0; i < num_nodes_; i++)
         {
            real_t V = std::sqrt(V_init_values_(2*i) * V_init_values_(2*i) +
                                 V_init_values_(2*i+1) * V_init_values_(2*i+1));
            V_ref = std::max(V_ref, V);
         }
         return V_ref;
      }
   }

   /// Get normal stress.
   real_t GetSigmaN() const
   {
      if constexpr (SlipComponents == 1)
      {
         return params_.sigma_n;
      }
      else
      {
         return sigma_n_bp5_;
      }
   }

   // =========================================================================
   // Verification and output
   // =========================================================================

   /// @brief Verify stress equilibrium at each node.
   ///
   /// @param[in] traction Traction from domain solve [TractionSize()]
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
         if constexpr (SlipComponents == 1)
         {
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
         else
         {
            if (depths(i) > Wf_bp5_ + 1.0) { continue; }

            real_t psi = state(i * StatePerNode + PsiIndex);
            real_t tau_vec[2] = {tau_pre_(2*i) + traction(2*i),
                                 tau_pre_(2*i+1) + traction(2*i+1)};
            real_t tau_abs = std::sqrt(tau_vec[0]*tau_vec[0] +
                                       tau_vec[1]*tau_vec[1]);
            real_t a = a_values(i);
            real_t eta = eta_values(i);

            real_t V_vec[2];
            dr_friction_->SolveSlipRateVectorPsi(
               tau_vec, psi, sigma_n_bp5_, eta, a, V_vec);
            real_t V_abs = std::sqrt(V_vec[0]*V_vec[0] + V_vec[1]*V_vec[1]);
            real_t f = dr_friction_->FrictionCoefficientPsi(V_abs, psi, a);
            real_t tau_computed = sigma_n_bp5_ * f + eta * V_abs;

            real_t rel_error = std::abs(tau_abs - tau_computed) /
                               std::max(tau_abs, 1.0);
            max_rel_error = std::max(max_rel_error, rel_error);
         }
      }

      return max_rel_error;
   }

   /// @brief Verify initial slip rate is reasonable.
   ///
   /// @param[in] V_max Maximum slip rate from initialization
   void VerifyInitialSlipRate(real_t V_max) const
   {
      real_t V_ref = GetReferenceVInit();
      real_t V_rel_err = std::abs(V_max - V_ref) /
                         std::max(V_ref, 1e-30);
      MFEM_VERIFY(V_rel_err < 0.5,
                  "Initial V_max = " << V_max
                  << " differs from V_init = " << V_ref
                  << " by " << V_rel_err * 100 << "%");
   }

   /// @brief Print fault state information.
   void PrintState(const Vector &state, std::ostream &os = mfem::out) const
   {
      os << "Fault State Summary:\n";
      if constexpr (SlipComponents == 1)
      {
         os << "  tau0 = " << tau0_ / 1e6 << " MPa\n";
      }
      os << "  V_max = " << V_max_ << " m/s\n";

      if (num_nodes_ > 0)
      {
         Vector slip, theta;
         GetSlip(state, slip);
         GetTheta(state, theta);

         os << "  Theta range: [" << theta.Min() << ", " << theta.Max() << "] s\n";

         if constexpr (SlipComponents == 1)
         {
            os << "  Slip range: [" << slip.Min() << ", " << slip.Max() << "] m\n";
            os << "  Slip rate range: [" << slip_rate_.Min() << ", "
               << slip_rate_.Max() << "] m/s\n";
         }
         else
         {
            real_t slip_min = 1e30, slip_max = 0.0;
            real_t Vr_min = 1e30, Vr_max = 0.0;
            for (int i = 0; i < num_nodes_; i++)
            {
               real_t s = std::sqrt(slip(2*i)*slip(2*i) +
                                    slip(2*i+1)*slip(2*i+1));
               real_t v = std::sqrt(slip_rate_(2*i)*slip_rate_(2*i) +
                                    slip_rate_(2*i+1)*slip_rate_(2*i+1));
               slip_min = std::min(slip_min, s);
               slip_max = std::max(slip_max, s);
               Vr_min = std::min(Vr_min, v);
               Vr_max = std::max(Vr_max, v);
            }
            os << "  |Slip| range: [" << slip_min << ", " << slip_max << "] m\n";
            os << "  |V| range: [" << Vr_min << ", " << Vr_max << "] m/s\n";
         }
      }
   }

   /// @brief Set cached slip rate from checkpoint data.
   void SetSlipRate(const Vector &V)
   {
      slip_rate_ = V;
      if constexpr (SlipComponents == 1)
      {
         V_max_ = V.Normlinf();
      }
      else
      {
         V_max_ = 0.0;
         for (int i = 0; i < num_nodes_; i++)
         {
            real_t V_abs = std::sqrt(V(2*i)*V(2*i) + V(2*i+1)*V(2*i+1));
            V_max_ = std::max(V_max_, V_abs);
         }
      }
   }

   /// Whether psi-space integration is active.
   bool UsePsi() const { return use_psi_; }

   /// Set psi initialization mode for BP5.
   /// If true (default), use SCEC-correct psi = f0 + b*ln(V0/V_init).
   /// If false, absorb delta_tau into psi via InitialStatePsi (matches Tandem).
   void SetScecPsiInit(bool scec) { scec_psi_init_ = scec; }

   /// Enable traction monitoring at a few stations every N steps.
   /// @param interval Log every N calls to ComputeRHS (0 = disabled)
   /// @param station_indices Fault DOF indices to monitor (empty = auto-pick 3)
   void SetTractionMonitoring(int interval,
                               const std::vector<int> &station_indices = {})
   {
      monitor_interval_ = interval;
      monitor_stations_ = station_indices;
      monitor_call_count_ = 0;
   }

private:
   FaultGeometry<MeshType> *geom_;
   FrictionLaw *friction_;
   StateEvolution *evolution_;
   MPIContext *mpi_ctx_ = nullptr;

   int num_nodes_;      ///< Number of fault DOFs
   real_t tau0_;        ///< Pre-stress [Pa] (BP2 scalar)
   real_t V_max_;       ///< Maximum slip rate from last evaluation

   bool use_psi_ = false;  ///< If true, state variable is psi instead of theta
   bool scec_psi_init_ = false;  ///< If false (default), Tandem InitialStatePsi; if true, SCEC fixed psi
   DieterichRuinaFriction *dr_friction_ = nullptr;  ///< Downcast for psi methods

   Vector slip_rate_;   ///< Cached slip rate [SlipComponents * NumNodes()]

   // BP2 parameters (only used when SlipComponents == 1)
   BP2Params params_;

   // BP5 parameters (only used when SlipComponents == 2)
   BP5Params bp5_params_;
   real_t sigma_n_bp5_ = 0.0;
   real_t Vp_bp5_ = 0.0;
   real_t Wf_bp5_ = 0.0;
   Vector Dc_values_;       ///< Per-DOF critical slip distance [NumNodes()]
   Vector tau_pre_;         ///< Per-DOF pre-stress [2 * NumNodes()] (BP5 only)
   Vector V_init_values_;   ///< Per-DOF initial velocity [2 * NumNodes()] (BP5 only)

   // Traction monitoring
   int monitor_interval_ = 0;        ///< Log every N ComputeRHS calls (0=off)
   mutable int monitor_call_count_ = 0;
   std::vector<int> monitor_stations_;  ///< Fault DOF indices to monitor
};

/// Convenience type aliases for common template instantiations.
using BP2FaultOp = RateStateFaultOperator<Mesh, 1>;
using BP5FaultOp = RateStateFaultOperator<Mesh, 2>;

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_RATE_STATE_FAULT_HPP
