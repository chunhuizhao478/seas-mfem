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

#include "godunov_flux.hpp"
#include <cmath>
#include <cstring>

namespace mfem
{
namespace seas
{

// ---------------------------------------------------------------------------
// Construction: precompute A_x, A_x^+, A_x^-
// ---------------------------------------------------------------------------
GodunovFlux::GodunovFlux(real_t lambda, real_t mu, real_t rho)
   : lambda_(lambda), mu_(mu), rho_(rho),
     Ax_(NUM_STATE, NUM_STATE),
     Ax_plus_(NUM_STATE, NUM_STATE),
     Ax_minus_(NUM_STATE, NUM_STATE)
{
   MFEM_VERIFY(mu > 0, "Shear modulus mu must be positive, got " << mu);
   MFEM_VERIFY(rho > 0, "Density rho must be positive, got " << rho);
   MFEM_VERIFY(lambda + 2.0*mu > 0,
               "lambda+2*mu must be positive, got " << lambda + 2.0*mu);

   cp_ = std::sqrt((lambda + 2.0*mu) / rho);
   cs_ = std::sqrt(mu / rho);
   Zp_ = rho * cp_;
   Zs_ = rho * cs_;

   // Build x-direction Jacobian
   BuildJacobian(0, Ax_);

   // Compute A_x^+ and A_x^- via eigenvector decomposition.
   //
   // For the x-direction Jacobian:
   //   Eigenvalues: {+cp, +cs, +cs, 0, 0, 0, -cs, -cs, -cp}
   //   Eigenvectors: known analytically (see plan Section 2.1.3)
   //
   // We build R (columns = right eigenvectors), invert to get R^{-1},
   // then: A_x^+ = R * diag(cp,cs,cs,0,0,0,0,0,0) * R^{-1}
   //        A_x^- = R * diag(0,0,0,0,0,0,-cs,-cs,-cp) * R^{-1}

   real_t lp = lambda + 2.0*mu;  // = rho * cp^2
   real_t inv_rho = 1.0 / rho;

   // Build eigenvector matrix R (9x9, columns = eigenvectors)
   DenseMatrix R(NUM_STATE, NUM_STATE);
   R = 0.0;

   // Column 0: lambda = +cp (rightgoing P-wave)
   R(SXX, 0) = -lp;   R(SYY, 0) = -lambda;  R(SZZ, 0) = -lambda;
   R(VX, 0) = cp_;

   // Column 1: lambda = +cs (rightgoing S-wave, y-polarized)
   R(SXY, 1) = -mu;    R(VY, 1) = cs_;

   // Column 2: lambda = +cs (rightgoing S-wave, z-polarized)
   R(SXZ, 2) = -mu;    R(VZ, 2) = cs_;

   // Columns 3-5: lambda = 0 (zero modes)
   R(SYY, 3) = 1.0;  // sigma_yy
   R(SZZ, 4) = 1.0;  // sigma_zz
   R(SYZ, 5) = 1.0;  // sigma_yz

   // Column 6: lambda = -cs (leftgoing S-wave, y-polarized)
   R(SXY, 6) = mu;     R(VY, 6) = cs_;

   // Column 7: lambda = -cs (leftgoing S-wave, z-polarized)
   R(SXZ, 7) = mu;     R(VZ, 7) = cs_;

   // Column 8: lambda = -cp (leftgoing P-wave)
   R(SXX, 8) = lp;    R(SYY, 8) = lambda;   R(SZZ, 8) = lambda;
   R(VX, 8) = cp_;

   // Invert R to get R^{-1}
   DenseMatrix Rinv(R);
   DenseMatrixInverse Rinv_solver(Rinv);
   DenseMatrix Rinv_mat(NUM_STATE, NUM_STATE);
   Rinv_solver.GetInverseMatrix(Rinv_mat);

   // Build diagonal eigenvalue matrices
   // Lambda_plus  = diag(cp, cs, cs, 0, 0, 0, 0, 0, 0)
   // Lambda_minus = diag(0, 0, 0, 0, 0, 0, -cs, -cs, -cp)
   DenseMatrix Lambda_plus(NUM_STATE, NUM_STATE);
   DenseMatrix Lambda_minus(NUM_STATE, NUM_STATE);
   Lambda_plus = 0.0;
   Lambda_minus = 0.0;

   Lambda_plus(0, 0) = cp_;
   Lambda_plus(1, 1) = cs_;
   Lambda_plus(2, 2) = cs_;

   Lambda_minus(6, 6) = -cs_;
   Lambda_minus(7, 7) = -cs_;
   Lambda_minus(8, 8) = -cp_;

   // A_x^+ = R * Lambda_plus * R^{-1}
   // A_x^- = R * Lambda_minus * R^{-1}
   DenseMatrix temp(NUM_STATE, NUM_STATE);

   // temp = Lambda_plus * R^{-1}
   Mult(Lambda_plus, Rinv_mat, temp);
   // Ax_plus_ = R * temp
   mfem::Mult(R, temp, Ax_plus_);

   // temp = Lambda_minus * R^{-1}
   Mult(Lambda_minus, Rinv_mat, temp);
   // Ax_minus_ = R * temp
   mfem::Mult(R, temp, Ax_minus_);

   // ADER I-05 Phase 2: precompute the physical-frame material Jacobians
   // used by the CK recursion.  Because `ApplySpatialDerivative` returns
   // PHYSICAL-frame derivatives (via `CalcPhysDShape`, which absorbs the
   // reference→physical J^{-T} transform), the CK recursion `L(Q) =
   // -Σ_d A_d ∂_d Q` can use `BuildJacobian(d, ·)` directly as `A_d`
   // without an additional reference-frame transform.  For the
   // homogeneous isotropic material TPV102 uses, the result is the same
   // 9x9 matrix on every element — cache it once.
   // See GetReferenceStarMatrix's doc-comment for the "star" naming.
   for (int d = 0; d < 3; d++)
   {
      ref_star_[d].SetSize(NUM_STATE, NUM_STATE);
      BuildJacobian(d, ref_star_[d]);
   }
}

// ---------------------------------------------------------------------------
// ADER I-05 Phase 2: accessor for cached reference star matrices.
// ---------------------------------------------------------------------------
const DenseMatrix &GodunovFlux::GetReferenceStarMatrix(int dir) const
{
   MFEM_VERIFY(dir >= 0 && dir < 3,
               "GetReferenceStarMatrix: dir must be in {0,1,2}, got " << dir);
   return ref_star_[dir];
}

// ---------------------------------------------------------------------------
// Build the 9x9 Jacobian matrix for direction dir
// ---------------------------------------------------------------------------
void GodunovFlux::BuildJacobian(int dir, DenseMatrix &A) const
{
   A.SetSize(NUM_STATE, NUM_STATE);
   A = 0.0;

   real_t lp = lambda_ + 2.0*mu_;
   real_t inv_rho = 1.0 / rho_;

   // Map direction to velocity and stress indices.
   // For dir=0 (x): normal stress = SXX, shear stresses = SXY, SXZ
   // For dir=1 (y): normal stress = SYY, shear stresses = SXY, SYZ
   // For dir=2 (z): normal stress = SZZ, shear stresses = SXZ, SYZ

   // Velocity component along this direction
   int v_dir = VX + dir;  // VX=6 for dir=0, VY=7 for dir=1, VZ=8 for dir=2

   // Normal stress component (sigma_{dir,dir})
   int s_nn = dir;  // SXX=0 for dir=0, SYY=1 for dir=1, SZZ=2 for dir=2

   // Map to shear stress indices
   // sigma_{dir,0}, sigma_{dir,1}, sigma_{dir,2} for each direction
   // Using Voigt: SXY=3 -> (0,1), SYZ=4 -> (1,2), SXZ=5 -> (0,2)
   auto voigt = [](int i, int j) -> int {
      if (i == j) return i;  // SXX=0, SYY=1, SZZ=2
      if ((i == 0 && j == 1) || (i == 1 && j == 0)) return SXY;
      if ((i == 1 && j == 2) || (i == 2 && j == 1)) return SYZ;
      if ((i == 0 && j == 2) || (i == 2 && j == 0)) return SXZ;
      return -1;
   };

   // Stress equations: d(sigma_ij)/dt = -C_{ijkl} d(v_k)/d(x_l)
   // For flux in direction dir: only d/d(x_dir) contributes
   // A[sigma_ij][v_k] = -C_{ij,k,dir}
   //
   // For isotropic: C_{ijkl} = lambda delta_ij delta_kl + mu(delta_ik delta_jl + delta_il delta_jk)
   // So C_{ij,k,dir} = lambda delta_ij delta_{k,dir} + mu(delta_{i,k} delta_{j,dir} + delta_{i,dir} delta_{j,k})

   for (int i = 0; i < 3; i++)
   {
      for (int j = i; j < 3; j++)
      {
         int s_ij = voigt(i, j);
         for (int k = 0; k < 3; k++)
         {
            real_t c_val = 0.0;
            if (i == j && k == dir) { c_val += lambda_; }
            if (i == k && j == dir) { c_val += mu_; }
            if (i == dir && j == k) { c_val += mu_; }
            if (std::abs(c_val) > 0.0)
            {
               A(s_ij, VX + k) = -c_val;
            }
         }
      }
   }

   // Momentum equations: rho d(v_k)/dt = d(sigma_{k,dir})/d(x_dir)
   // A[v_k][sigma_{k,dir}] = -1/rho
   for (int k = 0; k < 3; k++)
   {
      int s_kd = voigt(k, dir);
      A(VX + k, s_kd) = -inv_rho;
   }
}

// ---------------------------------------------------------------------------
// Rotation matrices for stress and velocity
// ---------------------------------------------------------------------------

// Build T^{-1}: global → face-local rotation (9x9)
// Face-local: x → normal, y → t1, z → t2
void GodunovFlux::BuildRotationInverse(const real_t *n, const real_t *t1,
                                       const real_t *t2, DenseMatrix &Tinv)
{
   Tinv.SetSize(NUM_STATE, NUM_STATE);
   Tinv = 0.0;

   // Basis vectors as rows of the 3x3 rotation matrix Q
   // Q maps global → local: v_local = Q * v_global
   // Row 0 = normal, Row 1 = t1, Row 2 = t2
   real_t Q[3][3] = {
      {n[0],  n[1],  n[2]},
      {t1[0], t1[1], t1[2]},
      {t2[0], t2[1], t2[2]}
   };

   // Velocity block (rows 6-8, cols 6-8): v_local = Q * v_global
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
         Tinv(VX + i, VX + j) = Q[i][j];

   // Stress block: sigma'_ab = Q_ai Q_bj sigma_ij
   // Using Voigt notation: sigma'_{Voigt(a,b)} = sum over Voigt(i,j) of
   //   (appropriate combination of Q_ai Q_bj) * sigma_{Voigt(i,j)}
   //
   // Map Voigt pairs: 0=(0,0), 1=(1,1), 2=(2,2), 3=(0,1), 4=(1,2), 5=(0,2)
   int vi[6] = {0, 1, 2, 0, 1, 0};
   int vj[6] = {0, 1, 2, 1, 2, 2};

   for (int ab = 0; ab < 6; ab++)
   {
      int a = vi[ab], b = vj[ab];
      for (int ij = 0; ij < 6; ij++)
      {
         int i = vi[ij], j = vj[ij];
         real_t val = Q[a][i] * Q[b][j];
         if (i != j) { val += Q[a][j] * Q[b][i]; }
         Tinv(ab, ij) = val;
      }
   }
}

// Build T: face-local → global rotation (9x9)
void GodunovFlux::BuildRotation(const real_t *n, const real_t *t1,
                                const real_t *t2, DenseMatrix &T)
{
   T.SetSize(NUM_STATE, NUM_STATE);
   T = 0.0;

   // Q^T maps local → global: v_global = Q^T * v_local
   real_t Q[3][3] = {
      {n[0],  n[1],  n[2]},
      {t1[0], t1[1], t1[2]},
      {t2[0], t2[1], t2[2]}
   };

   // Velocity block: v_global = Q^T * v_local
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
         T(VX + i, VX + j) = Q[j][i];

   // Stress block: inverse rotation
   // sigma_{Voigt(i,j)} = sum_{Voigt(a,b)} T_sigma[ij][ab] * sigma'_{Voigt(a,b)}
   // T_sigma is the inverse of Tinv_sigma. For orthogonal Q, T_sigma = Tinv_sigma^T
   // when using the same Voigt conventions... but it's cleaner to just use Q^T.
   //
   // sigma_ij = Q^T_ia Q^T_jb sigma'_ab = Q_ai Q_bj sigma'_ab (same formula, swap a↔i,b↔j)
   int vi[6] = {0, 1, 2, 0, 1, 0};
   int vj[6] = {0, 1, 2, 1, 2, 2};

   for (int ij = 0; ij < 6; ij++)
   {
      int i = vi[ij], j = vj[ij];
      for (int ab = 0; ab < 6; ab++)
      {
         int a = vi[ab], b = vj[ab];
         real_t val = Q[a][i] * Q[b][j];
         if (a != b) { val += Q[b][i] * Q[a][j]; }
         T(ij, ab) = val;
      }
   }
}

// Build orthonormal tangent frame from a unit normal
void GodunovFlux::BuildFrame(const real_t *nor, real_t *t1, real_t *t2)
{
   // Choose initial vector not parallel to nor
   real_t up[3] = {0.0, 0.0, 1.0};
   real_t dot = nor[0]*up[0] + nor[1]*up[1] + nor[2]*up[2];
   if (std::abs(dot) > 0.9)
   {
      up[0] = 1.0; up[1] = 0.0; up[2] = 0.0;
   }

   // t1 = up x nor (Gram-Schmidt)
   t1[0] = up[1]*nor[2] - up[2]*nor[1];
   t1[1] = up[2]*nor[0] - up[0]*nor[2];
   t1[2] = up[0]*nor[1] - up[1]*nor[0];
   real_t len = std::sqrt(t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2]);
   t1[0] /= len; t1[1] /= len; t1[2] /= len;

   // t2 = nor x t1
   t2[0] = nor[1]*t1[2] - nor[2]*t1[1];
   t2[1] = nor[2]*t1[0] - nor[0]*t1[2];
   t2[2] = nor[0]*t1[1] - nor[1]*t1[0];
}

// ---------------------------------------------------------------------------
// Apply precomputed split flux in the rotated (face-normal) frame
// ---------------------------------------------------------------------------
void GodunovFlux::ApplySplitFlux(const real_t *Q_self_rot,
                                 const real_t *Q_nbr_rot,
                                 real_t *F_rot) const
{
   // F_rot = A_x^+ * Q_self_rot + A_x^- * Q_nbr_rot
   for (int i = 0; i < NUM_STATE; i++)
   {
      real_t sum = 0.0;
      for (int j = 0; j < NUM_STATE; j++)
      {
         sum += Ax_plus_(i, j) * Q_self_rot[j]
              + Ax_minus_(i, j) * Q_nbr_rot[j];
      }
      F_rot[i] = sum;
   }
}

// ---------------------------------------------------------------------------
// Interior flux: F_h = A_n^+ Q_self + A_n^- Q_nbr
// ---------------------------------------------------------------------------
void GodunovFlux::Interior(const real_t *nor, const real_t *Q_self,
                           const real_t *Q_nbr, real_t *F_h) const
{
   // 1. Build orthonormal frame
   real_t t1[3], t2[3];
   BuildFrame(nor, t1, t2);

   // 2. Build rotation matrices
   DenseMatrix Tinv(NUM_STATE, NUM_STATE);
   DenseMatrix T(NUM_STATE, NUM_STATE);
   BuildRotationInverse(nor, t1, t2, Tinv);
   BuildRotation(nor, t1, t2, T);

   // 3. Rotate states to face-normal frame
   real_t Q_self_rot[NUM_STATE], Q_nbr_rot[NUM_STATE];
   Tinv.Mult(Q_self, Q_self_rot);
   Tinv.Mult(Q_nbr, Q_nbr_rot);

   // 4. Apply split flux in rotated frame
   real_t F_rot[NUM_STATE];
   ApplySplitFlux(Q_self_rot, Q_nbr_rot, F_rot);

   // 5. Rotate back to global frame
   T.Mult(F_rot, F_h);
}

// ---------------------------------------------------------------------------
// Central (non-dissipative) flux for Zhang et al. 2023 Mixed-Flux dispatch.
// F_h = 0.5 · A_n · (Q_self + Q_nbr) — same rotation as Interior, but the
// rotated-frame inner step uses (Ax_plus_ + Ax_minus_) (= full A_x_face_local)
// instead of the eigenvalue-split upwind combination.
// ---------------------------------------------------------------------------
void GodunovFlux::Central(const real_t *nor, const real_t *Q_self,
                          const real_t *Q_nbr, real_t *F_h) const
{
   // 1. Build orthonormal frame (same as Interior).
   real_t t1[3], t2[3];
   BuildFrame(nor, t1, t2);

   // 2. Build rotation matrices (same as Interior).
   DenseMatrix Tinv(NUM_STATE, NUM_STATE);
   DenseMatrix T(NUM_STATE, NUM_STATE);
   BuildRotationInverse(nor, t1, t2, Tinv);
   BuildRotation(nor, t1, t2, T);

   // 3. Rotate states to face-normal frame (same as Interior).
   real_t Q_self_rot[NUM_STATE], Q_nbr_rot[NUM_STATE];
   Tinv.Mult(Q_self, Q_self_rot);
   Tinv.Mult(Q_nbr, Q_nbr_rot);

   // 4. Central flux in rotated frame:
   //    F_rot = 0.5 · (Ax_plus_ + Ax_minus_) · (Q_self_rot + Q_nbr_rot)
   //    Note: (Ax_plus_ + Ax_minus_) is the FULL face-normal Jacobian A_x
   //    in the face-rotated frame, NOT the global GetReferenceStarMatrix(0)
   //    (which is A_x in global Cartesian coords).
   real_t Q_sum[NUM_STATE];
   for (int i = 0; i < NUM_STATE; i++)
   {
      Q_sum[i] = Q_self_rot[i] + Q_nbr_rot[i];
   }
   real_t F_rot[NUM_STATE];
   for (int i = 0; i < NUM_STATE; i++)
   {
      real_t s = 0.0;
      for (int j = 0; j < NUM_STATE; j++)
      {
         s += (Ax_plus_(i, j) + Ax_minus_(i, j)) * Q_sum[j];
      }
      F_rot[i] = 0.5 * s;
   }

   // 5. Rotate back to global frame (same as Interior).
   T.Mult(F_rot, F_h);
}

// ---------------------------------------------------------------------------
// Absorbing BC flux: F_abs = A_n^+ Q_self (Eq. 5)
// ---------------------------------------------------------------------------
void GodunovFlux::Absorbing(const real_t *nor, const real_t *Q_self,
                            real_t *F_h) const
{
   // Same as Interior with Q_nbr = 0
   real_t Q_zero[NUM_STATE];
   std::memset(Q_zero, 0, NUM_STATE * sizeof(real_t));
   Interior(nor, Q_self, Q_zero, F_h);
}

// ---------------------------------------------------------------------------
// AbsorbingTotal (I-06 migration): background-aware absorbing BC.
// ---------------------------------------------------------------------------
// Addresses REVIEW.md R-I06-001.  Under the Phase 4 total-stress migration,
// bulk Q carries the pre-stress tensor at every DOF, including DOFs
// adjacent to absorbing (lateral, bottom, far-Y) boundary faces.  The
// original Absorbing flux uses Q_ghost = 0, so at equilibrium (uniform
// Q = Q_pre with no waves) the flux reduces to A^+ Q_pre -- NOT equal
// to the interior identity A Q_pre, so a constant-in-time source of
// outgoing waves radiates the pre-stress at every absorbing face.
//
// AbsorbingTotal uses Q_ghost = Q_bg (the supplied background).  When
// Q_self = Q_bg, Interior(Q_bg, Q_bg) = A Q_bg, which gives the same
// value at every face of a closed element and sums to zero after the
// standard surface-integral -- equilibrium is preserved.  Perturbations
// Q_self - Q_bg are upwind-damped by A^-, so the BC still absorbs
// outgoing waves relative to the background.
// ---------------------------------------------------------------------------
void GodunovFlux::AbsorbingTotal(const real_t *nor, const real_t *Q_self,
                                 const real_t *Q_bg, real_t *F_h) const
{
   Interior(nor, Q_self, Q_bg, F_h);
}

// ---------------------------------------------------------------------------
// FreeSurfaceTotal (I-06 R-I06-005): background-aware gamma-mirror free
// surface.  Enforces (sigma - sigma_bg).n = 0, i.e. ZERO FLUCTUATION
// traction, rather than sigma.n = 0.  On a uniform Q = Q_bg field the
// flux reduces to the interior identity F = A Q_bg and equilibrium is
// preserved — no spurious radiation from tilted free-surface faces with
// non-zero pre-stress traction.
//
// Construction: Q_pert = Q_self - Q_bg carries the fluctuation; apply
// the fluctuation-path gamma-mirror to Q_pert alone to get the ghost
// fluctuation; reconstruct Q_ghost = Q_bg + gamma*Q_pert.  Then the
// standard Godunov flux F = A^+ Q_self + A^- Q_ghost.
//
// Equivalence with FreeSurface on a flat horizontal free surface
// (TPV102's z=0): sigma_bg.n = 0 there, so Q_bg traction rotated to
// face-local has zero sigma_nn / sigma_nt1 / sigma_nt2; the
// gamma-mirror of Q_pert is also the gamma-mirror of Q_self; the two
// fluxes agree.  This matches the "fix inert for TPV102 production"
// property the reviewer noted.
// ---------------------------------------------------------------------------
void GodunovFlux::FreeSurfaceTotal(const real_t *nor, const real_t *Q_self,
                                   const real_t *Q_bg, real_t *F_h) const
{
   // 1. Build orthonormal frame.
   real_t t1[3], t2[3];
   BuildFrame(nor, t1, t2);

   // 2. Build rotation matrices.
   DenseMatrix Tinv(NUM_STATE, NUM_STATE);
   DenseMatrix T(NUM_STATE, NUM_STATE);
   BuildRotationInverse(nor, t1, t2, Tinv);
   BuildRotation(nor, t1, t2, T);

   // 3. Rotate Q_self and Q_bg only (R-003: Q_pert_rot = Q_self_rot -
   //    Q_bg_rot algebraically, by linearity of Tinv; saves one 9x9
   //    mat-vec vs rotating a separately-formed Q_pert).
   real_t Q_self_rot[NUM_STATE];
   real_t Q_bg_rot  [NUM_STATE];
   Tinv.Mult(Q_self, Q_self_rot);
   Tinv.Mult(Q_bg,   Q_bg_rot);
   real_t Q_pert_rot[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
   {
      Q_pert_rot[c] = Q_self_rot[c] - Q_bg_rot[c];
   }

   // 4. Gamma-mirror of the fluctuation.  In rotated frame:
   //      SXX (sigma_nn)  → flip (odd-normal-index)
   //      SYY (sigma_t1t1) → keep
   //      SZZ (sigma_t2t2) → keep
   //      SXY (sigma_nt1) → flip
   //      SYZ (sigma_t1t2) → keep
   //      SXZ (sigma_nt2) → flip
   //      VX, VY, VZ       → keep
   static const real_t gamma[NUM_STATE] = {-1, 1, 1, -1, 1, -1, 1, 1, 1};

   // 5. Ghost = Q_bg + gamma * Q_pert in rotated frame.
   real_t Q_ghost_rot[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
   {
      Q_ghost_rot[c] = Q_bg_rot[c] + gamma[c] * Q_pert_rot[c];
   }

   // 6. Apply split flux in rotated frame.
   real_t F_rot[NUM_STATE];
   ApplySplitFlux(Q_self_rot, Q_ghost_rot, F_rot);

   // 7. Rotate back to global frame.
   T.Mult(F_rot, F_h);
}

// ---------------------------------------------------------------------------
// FreeSurfaceGodunovTotal (I-06 R-I06-005): background-aware Godunov-
// projection free surface.  Same idea as FreeSurfaceTotal but uses the
// characteristic projection on the fluctuation Q_pert.
// ---------------------------------------------------------------------------
void GodunovFlux::FreeSurfaceGodunovTotal(const real_t *nor,
                                          const real_t *Q_self,
                                          const real_t *Q_bg,
                                          real_t *F_h) const
{
   // 1. Build frame.
   real_t t1[3], t2[3];
   BuildFrame(nor, t1, t2);

   // 2. Rotation matrices.
   DenseMatrix Tinv(NUM_STATE, NUM_STATE);
   DenseMatrix T(NUM_STATE, NUM_STATE);
   BuildRotationInverse(nor, t1, t2, Tinv);
   BuildRotation(nor, t1, t2, T);

   // 3. Rotate Q_self and Q_bg into face-local frame; compute
   //    fluctuation Q_pert_rot.
   real_t Q_self_rot[NUM_STATE];
   real_t Q_bg_rot  [NUM_STATE];
   Tinv.Mult(Q_self, Q_self_rot);
   Tinv.Mult(Q_bg,   Q_bg_rot);
   real_t Q_pert_rot[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
   {
      Q_pert_rot[c] = Q_self_rot[c] - Q_bg_rot[c];
   }

   // 4. Godunov projection on the fluctuation: sigma_pert.n = 0 at the
   //    surface, velocity perturbed by Z^{-1} . sigma_pert.  MFEM sign
   //    convention (see FreeSurfaceGodunov): velocity increment has
   //    MINUS sign.
   const real_t invZp = 1.0 / Zp_;
   const real_t invZs = 1.0 / Zs_;
   real_t Q_god_pert_rot[NUM_STATE];
   std::memcpy(Q_god_pert_rot, Q_pert_rot, NUM_STATE * sizeof(real_t));
   Q_god_pert_rot[SXX] = 0.0;
   Q_god_pert_rot[SXY] = 0.0;
   Q_god_pert_rot[SXZ] = 0.0;
   Q_god_pert_rot[VX]  = Q_pert_rot[VX] - invZp * Q_pert_rot[SXX];
   Q_god_pert_rot[VY]  = Q_pert_rot[VY] - invZs * Q_pert_rot[SXY];
   Q_god_pert_rot[VZ]  = Q_pert_rot[VZ] - invZs * Q_pert_rot[SXZ];

   // 5. Ghost = Q_bg + Q_god_pert in rotated frame.
   real_t Q_god_rot[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
   {
      Q_god_rot[c] = Q_bg_rot[c] + Q_god_pert_rot[c];
   }

   // 6. Apply split flux.
   real_t F_rot[NUM_STATE];
   ApplySplitFlux(Q_self_rot, Q_god_rot, F_rot);

   // 7. Rotate back.
   T.Mult(F_rot, F_h);
}

// ---------------------------------------------------------------------------
// Free-surface BC flux (Eq. 6)
// ---------------------------------------------------------------------------
void GodunovFlux::FreeSurface(const real_t *nor, const real_t *Q_self,
                              real_t *F_h) const
{
   // 1. Build orthonormal frame
   real_t t1[3], t2[3];
   BuildFrame(nor, t1, t2);

   // 2. Build rotation matrices
   DenseMatrix Tinv(NUM_STATE, NUM_STATE);
   DenseMatrix T(NUM_STATE, NUM_STATE);
   BuildRotationInverse(nor, t1, t2, Tinv);
   BuildRotation(nor, t1, t2, T);

   // 3. Rotate to face-normal frame
   real_t Q_rot[NUM_STATE];
   Tinv.Mult(Q_self, Q_rot);

   // 4. Create ghost state: Gamma mirrors stress with odd normal indices
   // In rotated frame (x = normal direction):
   //   SXX (0): sigma_nn  → flip  (1 normal index in each of nn)
   //   SYY (1): sigma_t1t1 → keep
   //   SZZ (2): sigma_t2t2 → keep
   //   SXY (3): sigma_nt1  → flip  (1 normal index)
   //   SYZ (4): sigma_t1t2 → keep
   //   SXZ (5): sigma_nt2  → flip  (1 normal index)
   //   VX (6): v_n → keep
   //   VY (7): v_t1 → keep
   //   VZ (8): v_t2 → keep
   static const real_t gamma[NUM_STATE] = {-1, 1, 1, -1, 1, -1, 1, 1, 1};

   real_t Q_ghost_rot[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
   {
      Q_ghost_rot[c] = gamma[c] * Q_rot[c];
   }

   // 5. Apply split flux in rotated frame
   real_t F_rot[NUM_STATE];
   ApplySplitFlux(Q_rot, Q_ghost_rot, F_rot);

   // 6. Rotate back to global frame
   T.Mult(F_rot, F_h);
}

// ---------------------------------------------------------------------------
// Free-surface BC flux via Godunov characteristic projection (I-04).
// ---------------------------------------------------------------------------
// PURPOSE (addresses REVIEW.md R-I04-004): this variant exists for
// SeisSol API parity and as the integration point for the v9.3.1 ADER
// free-surface path.  It is NOT a corner-pump remedy — both this path
// and `FreeSurface` call the same `BuildFrame`, so any Gram-Schmidt
// sensitivity at near-axis-aligned normals affects them identically.
//
// EQUIVALENCE WITH gamma-MIRROR: algebraically exact at any tilt on a
// flat facet.  Let delta := Q_ghost_gamma - Q_god in the rotated frame.
// Then delta lies entirely in span{R(:,0..5)} (right-going + zero modes),
// so A_x^- * delta = 0 and the two flux expressions are bit-equal modulo
// FP noise from the rotate -> ApplySplitFlux -> rotate-back chain.
//
// The compliance projector S = -R_{21} R_{11}^{-1} in the rotated frame is
// diagonal for isotropic elasticity.  SIGN NOTE (deviation from v9.3.0 plan
// pseudocode): the plan Phase 1 section writes S = diag(-1/Zp, -1/Zs, -1/Zs)
// and "Q_god[vel] = Q_self[vel] + (1/Zp) * Q_self[sigma_nn]" (plus sign).
// That pseudocode assumes SeisSol's eigenvector sign convention where the
// +cp mode has +lp at SXX.  MFEM's R (see godunov_flux.cpp ctor: R(SXX,0) =
// -lp) uses the OPPOSITE sign, so the derived compliance block is
//      S_MFEM = +diag(1/Zp, 1/Zs, 1/Zs)
// and the velocity update is
//      Q_god[vel] = Q_self[vel] + S_MFEM * (0 - Q_self[trac])
//                 = Q_self[vel] - S_MFEM * Q_self[trac]        (MINUS sign)
//
// Both sign conventions give the SAME physical Godunov state — what
// matters is that the implementation's S sign matches its R sign.  The
// minus form below is what passes the equivalence-with-gamma-mirror test
// on a flat free surface (proof: the Riemann interface state
// Q_L + P^-(Q_R - Q_L) with Q_R = gamma*Q_L has Q_god[VX] = Q_L[VX] -
// (1/Zp) * Q_L[SXX] under MFEM's eigenvector scaling).
//
// Rotated-frame components (x = normal, y = t1, z = t2):
//    Q_god[SXX] = 0            (sigma_nn)
//    Q_god[SXY] = 0            (sigma_nt1)
//    Q_god[SXZ] = 0            (sigma_nt2)
//    Q_god[VX]  = Q_self[VX] - (1/Zp) * Q_self[SXX]
//    Q_god[VY]  = Q_self[VY] - (1/Zs) * Q_self[SXY]
//    Q_god[VZ]  = Q_self[VZ] - (1/Zs) * Q_self[SXZ]
// Passive modes (SYY, SZZ, SYZ) pass through unchanged: they are zero-
// mode eigenvectors of A_x and do not contribute to A^-.
//
// R-004: inlined invZp/invZs; no cached compliance-block member.
// ---------------------------------------------------------------------------
void GodunovFlux::FreeSurfaceGodunov(const real_t *nor, const real_t *Q_self,
                                     real_t *F_h) const
{
   // 1. Build orthonormal frame (matches Interior / FreeSurface).
   real_t t1[3], t2[3];
   BuildFrame(nor, t1, t2);

   // 2. Build rotation matrices.
   DenseMatrix Tinv(NUM_STATE, NUM_STATE);
   DenseMatrix T(NUM_STATE, NUM_STATE);
   BuildRotationInverse(nor, t1, t2, Tinv);
   BuildRotation(nor, t1, t2, T);

   // 3. Rotate Q_self into face-local frame.
   real_t Q_rot[NUM_STATE];
   Tinv.Mult(Q_self, Q_rot);

   // 4. Godunov projection in MFEM's R convention (see SIGN NOTE above).
   const real_t invZp = 1.0 / Zp_;
   const real_t invZs = 1.0 / Zs_;

   real_t Q_god_rot[NUM_STATE];
   std::memcpy(Q_god_rot, Q_rot, NUM_STATE * sizeof(real_t));
   Q_god_rot[SXX] = 0.0;                       // sigma_nn  = 0
   Q_god_rot[SXY] = 0.0;                       // sigma_nt1 = 0
   Q_god_rot[SXZ] = 0.0;                       // sigma_nt2 = 0
   Q_god_rot[VX]  = Q_rot[VX] - invZp * Q_rot[SXX];
   Q_god_rot[VY]  = Q_rot[VY] - invZs * Q_rot[SXY];
   Q_god_rot[VZ]  = Q_rot[VZ] - invZs * Q_rot[SXZ];

   // 5. Apply split flux in rotated frame with the Godunov-projected ghost.
   real_t F_rot[NUM_STATE];
   ApplySplitFlux(Q_rot, Q_god_rot, F_rot);

   // 6. Rotate back to global frame.
   T.Mult(F_rot, F_h);
}

} // namespace seas
} // namespace mfem
