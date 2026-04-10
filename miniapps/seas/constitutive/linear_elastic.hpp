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

#ifndef MFEM_SEAS_LINEAR_ELASTIC_HPP
#define MFEM_SEAS_LINEAR_ELASTIC_HPP

#include "constitutive_model.hpp"
#include <cmath>

namespace mfem
{
namespace seas
{

/// @brief Isotropic linear elastic constitutive model.
///
/// Implements Tier 1 only. No internal state variables, no evolution.
/// Stress: sigma = lambda * tr(epsilon) * I + 2 * mu * epsilon
/// Tangent: 6x6 isotropic Voigt stiffness (constant).
class LinearElastic : public ConstitutiveModel
{
public:
   LinearElastic(real_t lambda, real_t mu)
      : lambda_(lambda), mu_(mu)
   {
      MFEM_VERIFY(mu > 0, "LinearElastic: mu must be positive, got " << mu);
   }

   int NumInternalVars() const override { return 0; }
   bool IsNonlinear() const override { return false; }

   void ComputeStress(const real_t *epsilon, const real_t *int_vars,
                      real_t *sigma) const override
   {
      // sigma = lambda * tr(eps) * I + 2 * mu * eps
      // Voigt: eps = [exx, eyy, ezz, gamma_xy, gamma_yz, gamma_xz]
      //   where gamma = 2*epsilon (engineering strain)
      // sigma = [sxx, syy, szz, tau_xy, tau_yz, tau_xz]
      real_t tr_eps = epsilon[0] + epsilon[1] + epsilon[2];
      sigma[0] = lambda_ * tr_eps + 2.0 * mu_ * epsilon[0];
      sigma[1] = lambda_ * tr_eps + 2.0 * mu_ * epsilon[1];
      sigma[2] = lambda_ * tr_eps + 2.0 * mu_ * epsilon[2];
      // Engineering shear: gamma_xy = 2*eps_xy, tau_xy = 2*mu*eps_xy = mu*gamma_xy
      sigma[3] = mu_ * epsilon[3];
      sigma[4] = mu_ * epsilon[4];
      sigma[5] = mu_ * epsilon[5];
   }

   void ComputeTangent(const real_t *epsilon, const real_t *int_vars,
                       DenseMatrix &C_tang) const override
   {
      C_tang.SetSize(6, 6);
      C_tang = 0.0;

      // Normal components
      real_t lam2mu = lambda_ + 2.0 * mu_;
      C_tang(0, 0) = lam2mu;
      C_tang(1, 1) = lam2mu;
      C_tang(2, 2) = lam2mu;

      C_tang(0, 1) = lambda_;
      C_tang(0, 2) = lambda_;
      C_tang(1, 0) = lambda_;
      C_tang(1, 2) = lambda_;
      C_tang(2, 0) = lambda_;
      C_tang(2, 1) = lambda_;

      // Shear components (engineering Voigt: C(3,3) = mu, NOT 2*mu)
      C_tang(3, 3) = mu_;
      C_tang(4, 4) = mu_;
      C_tang(5, 5) = mu_;
   }

   void UpdateState(const real_t *epsilon, real_t *int_vars) const override
   {
      // No-op: no internal state
   }

   real_t GetPenaltyModulus() const override
   {
      return lambda_ + 2.0 * mu_;
   }

   real_t GetMaxWaveSpeed(real_t rho) const override
   {
      MFEM_VERIFY(rho > 0, "GetMaxWaveSpeed: rho must be positive");
      return std::sqrt((lambda_ + 2.0 * mu_) / rho);
   }

   // Convenience accessors (isotropic-only, not on base interface)
   real_t GetLambda() const { return lambda_; }
   real_t GetMu() const { return mu_; }

private:
   real_t lambda_;
   real_t mu_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_LINEAR_ELASTIC_HPP
