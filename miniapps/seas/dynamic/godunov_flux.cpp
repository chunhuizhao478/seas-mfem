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

} // namespace seas
} // namespace mfem
