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

#ifndef MFEM_SEAS_GODUNOV_FLUX_HPP
#define MFEM_SEAS_GODUNOV_FLUX_HPP

#include "mfem.hpp"
#include "wave_state.hpp"

namespace mfem
{
namespace seas
{

/// @brief Godunov upwind flux for the 3D isotropic elastic wave equation.
///
/// Precomputes split flux matrices A_x^+ and A_x^- for the x-direction
/// Jacobian from material parameters (lambda, mu, rho). For a general face
/// with unit normal n, the flux is computed via:
///   1. Rotate Q to face-normal frame (n becomes x-axis)
///   2. Apply A_x^+ to self-state, A_x^- to neighbor-state
///   3. Rotate result back to global frame
///
/// Reference: de la Puente et al. (2009), Section 2.2;
///            SeisSol ElasticSetup.h:146-165.
class GodunovFlux
{
public:
   /// Construct from material parameters (homogeneous, constant).
   /// @param[in] lambda  First Lame parameter [Pa]
   /// @param[in] mu  Shear modulus [Pa]
   /// @param[in] rho  Density [kg/m^3]
   GodunovFlux(real_t lambda, real_t mu, real_t rho);

   /// @name Flux computations
   ///@{

   /// Interior face flux: F_h = A_n^+ Q_self + A_n^- Q_nbr (Eq. 4).
   /// @param[in] nor  Unit outward normal from Elem1 (3 components)
   /// @param[in] Q_self  State on the self (Elem1) side (9 components)
   /// @param[in] Q_nbr  State on the neighbor (Elem2) side (9 components)
   /// @param[out] F_h  Numerical flux (9 components)
   void Interior(const real_t *nor, const real_t *Q_self,
                 const real_t *Q_nbr, real_t *F_h) const;

   /// First-order absorbing BC flux: F_abs = A_n^+ Q_self (Eq. 5).
   /// Sets incoming waves to zero (no reflection for normal incidence).
   void Absorbing(const real_t *nor, const real_t *Q_self,
                  real_t *F_h) const;

   /// Free-surface BC flux (Eq. 6).
   /// Mirrors stress components with odd number of normal indices via Gamma,
   /// enforcing sigma . n = 0 at the surface.
   void FreeSurface(const real_t *nor, const real_t *Q_self,
                    real_t *F_h) const;

   ///@}

   /// @name Accessors
   ///@{
   real_t GetCp() const { return cp_; }
   real_t GetCs() const { return cs_; }
   real_t GetLambda() const { return lambda_; }
   real_t GetMu() const { return mu_; }
   real_t GetRho() const { return rho_; }
   real_t GetZp() const { return Zp_; }
   real_t GetZs() const { return Zs_; }

   /// Get the x-direction Jacobian matrix A (for testing).
   const DenseMatrix &GetAx() const { return Ax_; }
   /// Get A_x^+ split flux matrix (for testing).
   const DenseMatrix &GetAxPlus() const { return Ax_plus_; }
   /// Get A_x^- split flux matrix (for testing).
   const DenseMatrix &GetAxMinus() const { return Ax_minus_; }
   /// Build the 9x9 Jacobian matrix for direction dir (0=x, 1=y, 2=z).
   ///
   /// NOTE: this method must remain public.
   /// tests/unit/test_godunov_interior_equal_sides_identity.cpp
   /// (see plan §15.4 T-E in tpv102_debug_v9.0.0_seissol_flux_comparison.md)
   /// reconstructs A externally to verify the identity
   ///    flux_.Interior(n, Q, Q) == T . A_x . T^{-1} . Q
   /// which underpins the Pelties 2012 eq. (9) per-side fix at
   /// wave_operator.inl:811-887 and :1232-1244.
   void BuildJacobian(int dir, DenseMatrix &A) const;

   /// Compute the 9x9 rotation matrix T^{-1} (global -> face-local).
   static void BuildRotationInverse(const real_t *nor, const real_t *t1,
                                    const real_t *t2, DenseMatrix &Tinv);

   /// Compute the 9x9 inverse rotation T (face-local -> global).
   static void BuildRotation(const real_t *nor, const real_t *t1,
                             const real_t *t2, DenseMatrix &T);

   /// Build an orthonormal frame (nor, t1, t2) from a unit normal.
   static void BuildFrame(const real_t *nor, real_t *t1, real_t *t2);
   ///@}

private:
   real_t lambda_, mu_, rho_;
   real_t cp_, cs_;
   real_t Zp_, Zs_;  ///< P-wave and S-wave impedance

   DenseMatrix Ax_;        ///< x-direction Jacobian (9x9)
   DenseMatrix Ax_plus_;   ///< Positive split flux A_x^+ (9x9)
   DenseMatrix Ax_minus_;  ///< Negative split flux A_x^- (9x9)

   /// Apply split flux in the rotated frame: F_rot = A_x^+ Q_self_rot + A_x^- Q_nbr_rot.
   void ApplySplitFlux(const real_t *Q_self_rot, const real_t *Q_nbr_rot,
                       real_t *F_rot) const;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_GODUNOV_FLUX_HPP
