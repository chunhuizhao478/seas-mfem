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

#include <array>

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
   /// Correct for FLUCTUATION-Q drivers (the ambient background is zero,
   /// so Q_ghost = 0 preserves equilibrium).  Under the v9.3.0 total-Q
   /// TPV102 path the ambient background is the pre-stress tensor, not
   /// zero — use `AbsorbingTotal` instead.
   void Absorbing(const real_t *nor, const real_t *Q_self,
                  real_t *F_h) const;

   /// Absorbing BC flux for total-Q drivers (I-06 migration).
   /// Replaces the zero ghost with a caller-supplied background state
   /// `Q_bg` (typically the uniform TPV102 pre-stress tensor in global
   /// coordinates).  When the local state equals the background, this
   /// flux reduces to the interior identity `F = A_n Q_bg`, preserving
   /// uniform pre-stress as a true equilibrium.  Perturbations
   /// `Q_self - Q_bg` are upwind-damped by the standard A^- projection.
   ///
   /// @param[in]  nor      Unit outward face normal (3 components).
   /// @param[in]  Q_self   Local state (9 components, global frame).
   /// @param[in]  Q_bg     Background state on the ghost side (9 components).
   /// @param[out] F_h      Numerical flux (9 components, global frame).
   void AbsorbingTotal(const real_t *nor, const real_t *Q_self,
                       const real_t *Q_bg, real_t *F_h) const;

   /// Free-surface BC flux (Eq. 6).
   /// Mirrors stress components with odd number of normal indices via Gamma,
   /// enforcing sigma . n = 0 at the surface.  Correct for FLUCTUATION-Q
   /// drivers (ambient background = 0).  Under total-Q with a tilted free
   /// surface, the pre-stress traction sigma_bg . n is non-zero and this
   /// flux radiates it as a spurious outgoing wave — use
   /// `FreeSurfaceTotal` instead.
   void FreeSurface(const real_t *nor, const real_t *Q_self,
                    real_t *F_h) const;

   /// Free-surface BC flux for total-Q drivers (I-06 R-I06-005).  Enforces
   /// `(sigma - sigma_bg) . n = 0`, i.e. zero FLUCTUATION traction at the
   /// surface.  Decomposes `Q_self = Q_bg + Q_pert` and applies the
   /// gamma-mirror to `Q_pert` only; the reconstructed ghost is
   /// `Q_bg + gamma * Q_pert`.  On a uniform pre-stress initial condition
   /// the flux reduces to the interior identity F = A Q_bg and equilibrium
   /// is preserved — no spurious radiation from tilted free-surface faces.
   /// TPV102's horizontal z=0 surface has `sigma_bg . n = 0` so this flux
   /// numerically agrees with `FreeSurface` there; the fix is inert in
   /// production but required for any future topography.
   void FreeSurfaceTotal(const real_t *nor, const real_t *Q_self,
                         const real_t *Q_bg, real_t *F_h) const;

   /// Free-surface flux via Godunov characteristic projection (I-04).
   /// Replaces the gamma-mirror ghost-cell state with a direct construction
   /// of the sigma.n=0 imposed state in the rotated frame.  Algebraically
   /// equivalent to `FreeSurface` at ANY tilt on a flat facet: let
   /// delta := Q_ghost_gamma - Q_god; then A^- * delta = 0 by construction
   /// (delta lies entirely in the right-going + zero-mode subspace of R),
   /// so the two paths produce bit-equal flux modulo FP noise from the
   /// rotate -> ApplySplitFlux -> rotate-back chain.
   ///
   /// Added for SeisSol API parity and as the integration point for the
   /// v9.3.1 ADER free-surface variant.  NOT a corner-pump remedy:
   /// both variants call the same `BuildFrame` for the tangent basis,
   /// so any Gram-Schmidt sensitivity in `BuildFrame` (e.g. the up-vector
   /// swap at |dot| > 0.9) affects both paths identically.
   ///
   /// R-004 note: the compliance projector S = -R_{21} R_{11}^{-1} is
   /// inlined using Zp_ / Zs_ already stored by the ctor — no cached
   /// member needed.  See the SIGN NOTE in `godunov_flux.cpp` for MFEM's
   /// eigenvector sign convention (S_MFEM = +diag(1/Zp, 1/Zs, 1/Zs),
   /// opposite of the v9.3.0 plan pseudocode's SeisSol convention).
   ///
   /// @param[in]  nor     Unit outward face normal (3 components)
   /// @param[in]  Q_self  State on the local side (9 components, global frame)
   /// @param[out] F_h     Numerical flux (9 components, global frame)
   void FreeSurfaceGodunov(const real_t *nor, const real_t *Q_self,
                           real_t *F_h) const;

   /// Godunov-projection free-surface flux for total-Q drivers
   /// (I-06 R-I06-005 companion to FreeSurfaceTotal).  Applies the
   /// characteristic projection to the FLUCTUATION `Q_pert = Q_self -
   /// Q_bg` rather than to `Q_self` itself, and adds the background
   /// traction back so uniform `Q = Q_bg` is preserved.
   void FreeSurfaceGodunovTotal(const real_t *nor, const real_t *Q_self,
                                const real_t *Q_bg, real_t *F_h) const;

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

   /// ADER I-05 Phase 2: "star" Jacobian matrix used by the CK recursion.
   ///
   /// In the usual ADER-DG formulation the CK recursion acts on
   /// REFERENCE-frame derivatives and would use reference-frame Jacobians
   /// `A^*_d = J · A_d` per element (J = element Jacobian).  This
   /// implementation sidesteps that: `WaveOperator::ApplySpatialDerivative`
   /// computes PHYSICAL-frame derivatives via
   /// `FiniteElement::CalcPhysDShape` (which already applies the `J^{-T}`
   /// transform to the reference-frame shape gradients), so the CK
   /// recursion can use the PHYSICAL-frame material Jacobians `A_d`
   /// directly.  For the homogeneous-isotropic material TPV102 uses, these
   /// physical-frame Jacobians are exactly the matrices returned by
   /// `BuildJacobian(d, ·)` and are the same on every element — so we
   /// precompute them once at ctor time and expose them here as a cheap
   /// const reference.
   ///
   /// The name "star" is retained for parity with SeisSol / ADER-DG
   /// literature, not because this implementation uses reference-frame
   /// Jacobians.  Future extension to heterogeneous material or bent
   /// elements will need either per-element star matrices or a
   /// reference-frame CK path.
   ///
   /// @param[in] dir  Spatial direction in {0, 1, 2}.
   /// @return          9x9 (physical-frame) Jacobian for direction `dir`.
   const DenseMatrix &GetReferenceStarMatrix(int dir) const;

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

   /// ADER I-05 Phase 2: cached reference star matrices [A_x, A_y, A_z].
   /// Precomputed in the ctor via BuildJacobian(d, .).  Homogeneous
   /// isotropic simplex assumption: ref_star_[d] == global-frame Jacobian
   /// (see GetReferenceStarMatrix doc comment).
   std::array<DenseMatrix, 3> ref_star_;

   /// Apply split flux in the rotated frame: F_rot = A_x^+ Q_self_rot + A_x^- Q_nbr_rot.
   void ApplySplitFlux(const real_t *Q_self_rot, const real_t *Q_nbr_rot,
                       real_t *F_rot) const;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_GODUNOV_FLUX_HPP
