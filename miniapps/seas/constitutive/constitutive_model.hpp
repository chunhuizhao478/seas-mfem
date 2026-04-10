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

#ifndef MFEM_SEAS_CONSTITUTIVE_MODEL_HPP
#define MFEM_SEAS_CONSTITUTIVE_MODEL_HPP

#include "mfem.hpp"

namespace mfem
{
namespace seas
{

/// @brief Abstract interface for constitutive models.
///
/// Handles three tiers of constitutive complexity:
///   Tier 1 (Algebraic):  Linear elastic — no state, no evolution.
///   Tier 2 (Evolving):   Damage-breakage — state ODEs integrated by RK45.
///   Tier 3 (Non-local):  Gradient-regularized damage — state + FE gradients.
///
/// Phase 3 implements Tier 1 only (LinearElastic). The interface is designed
/// so Tiers 2-3 work without breaking changes — all evolution/non-local
/// methods have default no-op implementations.
///
/// Voigt notation convention (engineering strain):
///   Index map: 0=xx, 1=yy, 2=zz, 3=xy, 4=yz, 5=xz
///   Shear strains are ENGINEERING: gamma_xy = 2*epsilon_xy
///   For isotropic: C(0,0) = lambda+2*mu, C(3,3) = mu (NOT 2*mu)
class ConstitutiveModel
{
public:
   virtual ~ConstitutiveModel() = default;

   // =========================================================================
   // Tier 1: Core stress-strain interface (required for ALL models)
   // =========================================================================

   /// Internal state variables per quadrature point.
   /// Linear elastic: 0. Damage-breakage: 8.
   virtual int NumInternalVars() const = 0;

   /// Whether Newton iteration is required (nonlinear stress-strain relation).
   virtual bool IsNonlinear() const = 0;

   /// Compute stress from strain.
   /// @param epsilon  Strain in Voigt notation [6]
   /// @param int_vars Internal state variables [NumInternalVars()] (nullptr if 0)
   /// @param sigma    Output stress in Voigt notation [6]
   virtual void ComputeStress(const real_t *epsilon, const real_t *int_vars,
                              real_t *sigma) const = 0;

   /// Compute 6x6 tangent stiffness d(sigma)/d(epsilon) in engineering Voigt.
   /// @param epsilon  Strain in Voigt notation [6] (nullptr for linear models)
   /// @param int_vars Internal state variables [NumInternalVars()] (nullptr if 0)
   /// @param C_tang   Output 6x6 tangent matrix
   virtual void ComputeTangent(const real_t *epsilon, const real_t *int_vars,
                               DenseMatrix &C_tang) const = 0;

   /// Update internal state after a converged Newton step.
   /// For algebraic plasticity: return mapping. For evolving models: no-op.
   virtual void UpdateState(const real_t *epsilon,
                            real_t *int_vars) const = 0;

   /// Max eigenvalue of C. Used for DG penalty and CFL.
   virtual real_t GetPenaltyModulus() const = 0;

   /// Max wave speed sqrt(max_eig(C)/rho). Used for CFL and radiation damping.
   virtual real_t GetMaxWaveSpeed(real_t rho) const = 0;

   // =========================================================================
   // Tier 2: Time-evolving state variables (default: no evolution)
   // =========================================================================

   /// Whether state variables evolve via ODEs that must be time-integrated.
   virtual bool HasStateEvolution() const { return false; }

   /// Compute d(int_vars)/dt for time integration.
   virtual void ComputeStateRates(const real_t *epsilon,
                                  const real_t *int_vars,
                                  real_t *rates) const {}

   // =========================================================================
   // Tier 3: Non-local models (default: local only)
   // =========================================================================

   /// Number of state variables requiring FE-space gradient computation.
   virtual int NumNonLocalVars() const { return 0; }

   /// Compute stress including non-local gradient terms.
   /// Default: delegates to ComputeStress (ignoring gradient).
   virtual void ComputeStressNonLocal(const real_t *epsilon,
                                      const real_t *int_vars,
                                      const real_t *grad_state,
                                      real_t *sigma) const
   {
      ComputeStress(epsilon, int_vars, sigma);
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CONSTITUTIVE_MODEL_HPP
