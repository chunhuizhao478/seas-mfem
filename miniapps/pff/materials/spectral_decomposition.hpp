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

#ifndef MFEM_PFF_SPECTRAL_DECOMPOSITION_HPP
#define MFEM_PFF_SPECTRAL_DECOMPOSITION_HPP

#include "mfem.hpp"
#include "linalg/kernels.hpp"
#include <cmath>
#include <algorithm>

namespace mfem
{

namespace pff
{

/** @brief Spectral decomposition of strain tensor for tension-compression split.
 *
 * Implements the Miehe et al. (2010) spectral decomposition of the strain tensor
 * into positive (tensile) and negative (compressive) parts:
 *
 *   epsilon = epsilon_+ + epsilon_-
 *
 * where:
 *   epsilon_+ = sum_{i: lambda_i > 0} lambda_i * n_i (x) n_i
 *   epsilon_- = sum_{i: lambda_i < 0} lambda_i * n_i (x) n_i
 *
 * with lambda_i and n_i being eigenvalues and eigenvectors of epsilon.
 *
 * The strain energy is then split as:
 *   psi = psi_+ + psi_-
 *   psi_+ = lambda/2 * <tr(epsilon)>_+^2 + mu * tr(epsilon_+^2)
 *   psi_- = lambda/2 * <tr(epsilon)>_-^2 + mu * tr(epsilon_-^2)
 *
 * where <x>_+ = max(x, 0) and <x>_- = min(x, 0).
 *
 * @tparam dim Spatial dimension (2 or 3)
 */
template <int dim>
class SpectralDecomposition
{
public:
   /** @brief Decompose strain tensor into positive and negative parts.
    *
    * @param[in] strain Input strain tensor (dim x dim symmetric matrix)
    * @param[out] strain_pos Positive (tensile) part of strain
    * @param[out] strain_neg Negative (compressive) part of strain
    */
   static void Decompose(const DenseMatrix &strain,
                         DenseMatrix &strain_pos,
                         DenseMatrix &strain_neg);

   /** @brief Compute positive (tensile) strain energy density.
    *
    * psi_+ = lambda/2 * <tr(epsilon)>_+^2 + mu * tr(epsilon_+^2)
    *
    * @param[in] strain Input strain tensor
    * @param[in] lambda First Lame parameter
    * @param[in] mu Second Lame parameter (shear modulus)
    * @return Positive strain energy density
    */
   static real_t PositiveEnergy(const DenseMatrix &strain,
                                real_t lambda, real_t mu);

   /** @brief Compute negative (compressive) strain energy density.
    *
    * psi_- = lambda/2 * <tr(epsilon)>_-^2 + mu * tr(epsilon_-^2)
    *
    * @param[in] strain Input strain tensor
    * @param[in] lambda First Lame parameter
    * @param[in] mu Second Lame parameter (shear modulus)
    * @return Negative strain energy density
    */
   static real_t NegativeEnergy(const DenseMatrix &strain,
                                real_t lambda, real_t mu);

   /** @brief Compute stress from positive strain part (for degradation).
    *
    * sigma_+ = lambda * <tr(epsilon)>_+ * I + 2 * mu * epsilon_+
    *
    * @param[in] strain Input strain tensor
    * @param[in] lambda First Lame parameter
    * @param[in] mu Second Lame parameter
    * @param[out] stress_pos Positive stress tensor
    */
   static void PositiveStress(const DenseMatrix &strain,
                              real_t lambda, real_t mu,
                              DenseMatrix &stress_pos);

   /** @brief Compute stress from negative strain part (undegraded).
    *
    * sigma_- = lambda * <tr(epsilon)>_- * I + 2 * mu * epsilon_-
    *
    * @param[in] strain Input strain tensor
    * @param[in] lambda First Lame parameter
    * @param[in] mu Second Lame parameter
    * @param[out] stress_neg Negative stress tensor
    */
   static void NegativeStress(const DenseMatrix &strain,
                              real_t lambda, real_t mu,
                              DenseMatrix &stress_neg);

   /** @brief Compute the derivative of positive stress w.r.t. strain.
    *
    * This is needed for the Jacobian of the elasticity problem.
    * Returns the 4th order tensor C_+ such that:
    *   dsigma_+ = C_+ : depsilon
    *
    * @param[in] strain Input strain tensor
    * @param[in] lambda First Lame parameter
    * @param[in] mu Second Lame parameter
    * @param[out] tangent The 4th order tangent tensor (dim^2 x dim^2 matrix)
    */
   static void PositiveStressTangent(const DenseMatrix &strain,
                                     real_t lambda, real_t mu,
                                     DenseMatrix &tangent);

   /** @brief Compute the derivative of negative stress w.r.t. strain.
    *
    * @param[in] strain Input strain tensor
    * @param[in] lambda First Lame parameter
    * @param[in] mu Second Lame parameter
    * @param[out] tangent The 4th order tangent tensor (dim^2 x dim^2 matrix)
    */
   static void NegativeStressTangent(const DenseMatrix &strain,
                                     real_t lambda, real_t mu,
                                     DenseMatrix &tangent);

private:
   /// Macaulay bracket: <x>_+ = max(x, 0)
   MFEM_HOST_DEVICE static real_t MacaulayPos(real_t x)
   {
      return (x > 0.0) ? x : 0.0;
   }

   /// Macaulay bracket: <x>_- = min(x, 0)
   MFEM_HOST_DEVICE static real_t MacaulayNeg(real_t x)
   {
      return (x < 0.0) ? x : 0.0;
   }

   /// Heaviside function: H(x) = 1 if x > 0, else 0
   MFEM_HOST_DEVICE static real_t Heaviside(real_t x)
   {
      return (x > 0.0) ? 1.0 : 0.0;
   }

   /// Compute eigenvalues and eigenvectors
   static void ComputeEigen(const DenseMatrix &A,
                            Vector &eigenvalues,
                            DenseMatrix &eigenvectors);
};

// Implementation for 2D
template <>
inline void SpectralDecomposition<2>::ComputeEigen(const DenseMatrix &A,
                                                   Vector &eigenvalues,
                                                   DenseMatrix &eigenvectors)
{
   MFEM_ASSERT(A.Height() == 2 && A.Width() == 2, "Matrix must be 2x2");

   eigenvalues.SetSize(2);
   eigenvectors.SetSize(2, 2);

   real_t lambda[2];
   real_t vec[4];

   // Use MFEM's optimized 2x2 eigenvalue solver
   kernels::CalcEigenvalues<2>(A.Data(), lambda, vec);

   eigenvalues[0] = lambda[0];
   eigenvalues[1] = lambda[1];

   // Eigenvectors are stored column-wise
   eigenvectors(0, 0) = vec[0];
   eigenvectors(1, 0) = vec[1];
   eigenvectors(0, 1) = vec[2];
   eigenvectors(1, 1) = vec[3];
}

// Implementation for 3D
template <>
inline void SpectralDecomposition<3>::ComputeEigen(const DenseMatrix &A,
                                                   Vector &eigenvalues,
                                                   DenseMatrix &eigenvectors)
{
   MFEM_ASSERT(A.Height() == 3 && A.Width() == 3, "Matrix must be 3x3");

   eigenvalues.SetSize(3);
   eigenvectors.SetSize(3, 3);

   real_t lambda[3];
   real_t vec[9];

   // Use MFEM's optimized 3x3 eigenvalue solver
   kernels::CalcEigenvalues<3>(A.Data(), lambda, vec);

   eigenvalues[0] = lambda[0];
   eigenvalues[1] = lambda[1];
   eigenvalues[2] = lambda[2];

   // Eigenvectors are stored column-wise
   for (int i = 0; i < 3; i++)
   {
      eigenvectors(0, i) = vec[3*i + 0];
      eigenvectors(1, i) = vec[3*i + 1];
      eigenvectors(2, i) = vec[3*i + 2];
   }
}

template <int dim>
void SpectralDecomposition<dim>::Decompose(const DenseMatrix &strain,
                                           DenseMatrix &strain_pos,
                                           DenseMatrix &strain_neg)
{
   MFEM_ASSERT(strain.Height() == dim && strain.Width() == dim,
               "Strain tensor must be dim x dim");

   strain_pos.SetSize(dim, dim);
   strain_neg.SetSize(dim, dim);
   strain_pos = 0.0;
   strain_neg = 0.0;

   Vector eigenvalues;
   DenseMatrix eigenvectors;
   ComputeEigen(strain, eigenvalues, eigenvectors);

   // Reconstruct positive and negative parts
   // epsilon_+ = sum_{i: lambda_i > 0} lambda_i * n_i (x) n_i
   // epsilon_- = sum_{i: lambda_i < 0} lambda_i * n_i (x) n_i
   for (int k = 0; k < dim; k++)
   {
      real_t lambda_k = eigenvalues[k];
      real_t lambda_pos = MacaulayPos(lambda_k);
      real_t lambda_neg = MacaulayNeg(lambda_k);

      for (int i = 0; i < dim; i++)
      {
         for (int j = 0; j < dim; j++)
         {
            real_t n_i = eigenvectors(i, k);
            real_t n_j = eigenvectors(j, k);
            strain_pos(i, j) += lambda_pos * n_i * n_j;
            strain_neg(i, j) += lambda_neg * n_i * n_j;
         }
      }
   }
}

template <int dim>
real_t SpectralDecomposition<dim>::PositiveEnergy(const DenseMatrix &strain,
                                                  real_t lambda, real_t mu)
{
   // psi_+ = lambda/2 * <tr(epsilon)>_+^2 + mu * tr(epsilon_+^2)

   Vector eigenvalues;
   DenseMatrix eigenvectors;
   ComputeEigen(strain, eigenvalues, eigenvectors);

   // Compute trace
   real_t tr_strain = 0.0;
   for (int i = 0; i < dim; i++)
   {
      tr_strain += eigenvalues[i];
   }

   real_t tr_pos = MacaulayPos(tr_strain);

   // Compute tr(epsilon_+^2) = sum of (max(lambda_i, 0))^2
   real_t tr_eps_pos_sq = 0.0;
   for (int i = 0; i < dim; i++)
   {
      real_t lambda_pos = MacaulayPos(eigenvalues[i]);
      tr_eps_pos_sq += lambda_pos * lambda_pos;
   }

   return 0.5 * lambda * tr_pos * tr_pos + mu * tr_eps_pos_sq;
}

template <int dim>
real_t SpectralDecomposition<dim>::NegativeEnergy(const DenseMatrix &strain,
                                                  real_t lambda, real_t mu)
{
   // psi_- = lambda/2 * <tr(epsilon)>_-^2 + mu * tr(epsilon_-^2)

   Vector eigenvalues;
   DenseMatrix eigenvectors;
   ComputeEigen(strain, eigenvalues, eigenvectors);

   // Compute trace
   real_t tr_strain = 0.0;
   for (int i = 0; i < dim; i++)
   {
      tr_strain += eigenvalues[i];
   }

   real_t tr_neg = MacaulayNeg(tr_strain);

   // Compute tr(epsilon_-^2) = sum of (min(lambda_i, 0))^2
   real_t tr_eps_neg_sq = 0.0;
   for (int i = 0; i < dim; i++)
   {
      real_t lambda_neg = MacaulayNeg(eigenvalues[i]);
      tr_eps_neg_sq += lambda_neg * lambda_neg;
   }

   return 0.5 * lambda * tr_neg * tr_neg + mu * tr_eps_neg_sq;
}

template <int dim>
void SpectralDecomposition<dim>::PositiveStress(const DenseMatrix &strain,
                                                real_t lambda, real_t mu,
                                                DenseMatrix &stress_pos)
{
   // sigma_+ = lambda * <tr(epsilon)>_+ * I + 2 * mu * epsilon_+

   stress_pos.SetSize(dim, dim);
   stress_pos = 0.0;

   Vector eigenvalues;
   DenseMatrix eigenvectors;
   ComputeEigen(strain, eigenvalues, eigenvectors);

   // Compute trace
   real_t tr_strain = 0.0;
   for (int i = 0; i < dim; i++)
   {
      tr_strain += eigenvalues[i];
   }

   real_t tr_pos = MacaulayPos(tr_strain);

   // Add lambda * <tr(epsilon)>_+ * I
   for (int i = 0; i < dim; i++)
   {
      stress_pos(i, i) += lambda * tr_pos;
   }

   // Add 2 * mu * epsilon_+
   DenseMatrix strain_pos(dim, dim), strain_neg(dim, dim);
   Decompose(strain, strain_pos, strain_neg);

   for (int i = 0; i < dim; i++)
   {
      for (int j = 0; j < dim; j++)
      {
         stress_pos(i, j) += 2.0 * mu * strain_pos(i, j);
      }
   }
}

template <int dim>
void SpectralDecomposition<dim>::NegativeStress(const DenseMatrix &strain,
                                                real_t lambda, real_t mu,
                                                DenseMatrix &stress_neg)
{
   // sigma_- = lambda * <tr(epsilon)>_- * I + 2 * mu * epsilon_-

   stress_neg.SetSize(dim, dim);
   stress_neg = 0.0;

   Vector eigenvalues;
   DenseMatrix eigenvectors;
   ComputeEigen(strain, eigenvalues, eigenvectors);

   // Compute trace
   real_t tr_strain = 0.0;
   for (int i = 0; i < dim; i++)
   {
      tr_strain += eigenvalues[i];
   }

   real_t tr_neg = MacaulayNeg(tr_strain);

   // Add lambda * <tr(epsilon)>_- * I
   for (int i = 0; i < dim; i++)
   {
      stress_neg(i, i) += lambda * tr_neg;
   }

   // Add 2 * mu * epsilon_-
   DenseMatrix strain_pos(dim, dim), strain_neg_part(dim, dim);
   Decompose(strain, strain_pos, strain_neg_part);

   for (int i = 0; i < dim; i++)
   {
      for (int j = 0; j < dim; j++)
      {
         stress_neg(i, j) += 2.0 * mu * strain_neg_part(i, j);
      }
   }
}

template <int dim>
void SpectralDecomposition<dim>::PositiveStressTangent(const DenseMatrix &strain,
                                                       real_t lambda, real_t mu,
                                                       DenseMatrix &tangent)
{
   // The tangent tensor C_+ is computed as:
   // C_+_{ijkl} = d(sigma_+)_{ij} / d(epsilon)_{kl}
   //
   // For spectral decomposition, this involves derivatives of the
   // projection operators, which can be complex when eigenvalues are close.
   //
   // Simplified form (assuming distinct eigenvalues):
   // C_+ = lambda * H(tr(eps)) * I (x) I
   //     + 2 * mu * sum_{i: lambda_i > 0} P_i (x) P_i
   //     + correction terms for eigenvalue derivatives

   tangent.SetSize(dim * dim, dim * dim);
   tangent = 0.0;

   Vector eigenvalues;
   DenseMatrix eigenvectors;
   ComputeEigen(strain, eigenvalues, eigenvectors);

   // Compute trace
   real_t tr_strain = 0.0;
   for (int i = 0; i < dim; i++)
   {
      tr_strain += eigenvalues[i];
   }

   real_t H_tr = Heaviside(tr_strain);

   // Lambda term: lambda * H(tr(eps)) * I (x) I
   // (I (x) I)_{ijkl} = delta_ij * delta_kl
   for (int i = 0; i < dim; i++)
   {
      for (int k = 0; k < dim; k++)
      {
         // Index mapping: (i,j) -> i*dim + j
         int ij = i * dim + i;  // j = i for I
         int kl = k * dim + k;  // l = k for I
         tangent(ij, kl) += lambda * H_tr;
      }
   }

   // Mu term: 2 * mu * sum_{i: lambda_i > 0} P_i (x) P_i
   // where P_i = n_i (x) n_i is the projection onto eigenvector i
   for (int m = 0; m < dim; m++)
   {
      if (eigenvalues[m] > 0.0)
      {
         for (int i = 0; i < dim; i++)
         {
            for (int j = 0; j < dim; j++)
            {
               for (int k = 0; k < dim; k++)
               {
                  for (int l = 0; l < dim; l++)
                  {
                     int ij = i * dim + j;
                     int kl = k * dim + l;
                     real_t n_i = eigenvectors(i, m);
                     real_t n_j = eigenvectors(j, m);
                     real_t n_k = eigenvectors(k, m);
                     real_t n_l = eigenvectors(l, m);
                     tangent(ij, kl) += 2.0 * mu * n_i * n_j * n_k * n_l;
                  }
               }
            }
         }
      }
   }

   // Add cross-derivative terms for distinct eigenvalues
   // These terms arise from d(n_i)/d(epsilon) contributions
   const real_t tol = 1e-10;
   for (int a = 0; a < dim; a++)
   {
      for (int b = a + 1; b < dim; b++)
      {
         real_t diff = eigenvalues[a] - eigenvalues[b];
         if (std::abs(diff) > tol)
         {
            real_t lambda_a_pos = MacaulayPos(eigenvalues[a]);
            real_t lambda_b_pos = MacaulayPos(eigenvalues[b]);
            real_t coeff = 2.0 * mu * (lambda_a_pos - lambda_b_pos) / diff;

            for (int i = 0; i < dim; i++)
            {
               for (int j = 0; j < dim; j++)
               {
                  for (int k = 0; k < dim; k++)
                  {
                     for (int l = 0; l < dim; l++)
                     {
                        int ij = i * dim + j;
                        int kl = k * dim + l;
                        real_t n_ai = eigenvectors(i, a);
                        real_t n_aj = eigenvectors(j, a);
                        real_t n_bi = eigenvectors(i, b);
                        real_t n_bj = eigenvectors(j, b);
                        real_t n_ak = eigenvectors(k, a);
                        real_t n_al = eigenvectors(l, a);
                        real_t n_bk = eigenvectors(k, b);
                        real_t n_bl = eigenvectors(l, b);

                        // Symmetric part of (n_a (x) n_b + n_b (x) n_a) (x) (n_a (x) n_b + n_b (x) n_a)
                        real_t term = 0.5 * (n_ai * n_bj + n_bi * n_aj) *
                                      0.5 * (n_ak * n_bl + n_bk * n_al);
                        tangent(ij, kl) += coeff * term;
                     }
                  }
               }
            }
         }
      }
   }
}

template <int dim>
void SpectralDecomposition<dim>::NegativeStressTangent(const DenseMatrix &strain,
                                                       real_t lambda, real_t mu,
                                                       DenseMatrix &tangent)
{
   // Similar to PositiveStressTangent but for negative eigenvalues

   tangent.SetSize(dim * dim, dim * dim);
   tangent = 0.0;

   Vector eigenvalues;
   DenseMatrix eigenvectors;
   ComputeEigen(strain, eigenvalues, eigenvectors);

   // Compute trace
   real_t tr_strain = 0.0;
   for (int i = 0; i < dim; i++)
   {
      tr_strain += eigenvalues[i];
   }

   real_t H_tr_neg = (tr_strain < 0.0) ? 1.0 : 0.0;

   // Lambda term: lambda * H(-tr(eps)) * I (x) I
   for (int i = 0; i < dim; i++)
   {
      for (int k = 0; k < dim; k++)
      {
         int ij = i * dim + i;
         int kl = k * dim + k;
         tangent(ij, kl) += lambda * H_tr_neg;
      }
   }

   // Mu term: 2 * mu * sum_{i: lambda_i < 0} P_i (x) P_i
   for (int m = 0; m < dim; m++)
   {
      if (eigenvalues[m] < 0.0)
      {
         for (int i = 0; i < dim; i++)
         {
            for (int j = 0; j < dim; j++)
            {
               for (int k = 0; k < dim; k++)
               {
                  for (int l = 0; l < dim; l++)
                  {
                     int ij = i * dim + j;
                     int kl = k * dim + l;
                     real_t n_i = eigenvectors(i, m);
                     real_t n_j = eigenvectors(j, m);
                     real_t n_k = eigenvectors(k, m);
                     real_t n_l = eigenvectors(l, m);
                     tangent(ij, kl) += 2.0 * mu * n_i * n_j * n_k * n_l;
                  }
               }
            }
         }
      }
   }

   // Add cross-derivative terms
   const real_t tol = 1e-10;
   for (int a = 0; a < dim; a++)
   {
      for (int b = a + 1; b < dim; b++)
      {
         real_t diff = eigenvalues[a] - eigenvalues[b];
         if (std::abs(diff) > tol)
         {
            real_t lambda_a_neg = MacaulayNeg(eigenvalues[a]);
            real_t lambda_b_neg = MacaulayNeg(eigenvalues[b]);
            real_t coeff = 2.0 * mu * (lambda_a_neg - lambda_b_neg) / diff;

            for (int i = 0; i < dim; i++)
            {
               for (int j = 0; j < dim; j++)
               {
                  for (int k = 0; k < dim; k++)
                  {
                     for (int l = 0; l < dim; l++)
                     {
                        int ij = i * dim + j;
                        int kl = k * dim + l;
                        real_t n_ai = eigenvectors(i, a);
                        real_t n_aj = eigenvectors(j, a);
                        real_t n_bi = eigenvectors(i, b);
                        real_t n_bj = eigenvectors(j, b);
                        real_t n_ak = eigenvectors(k, a);
                        real_t n_al = eigenvectors(l, a);
                        real_t n_bk = eigenvectors(k, b);
                        real_t n_bl = eigenvectors(l, b);

                        real_t term = 0.5 * (n_ai * n_bj + n_bi * n_aj) *
                                      0.5 * (n_ak * n_bl + n_bk * n_al);
                        tangent(ij, kl) += coeff * term;
                     }
                  }
               }
            }
         }
      }
   }
}

} // namespace pff

} // namespace mfem

#endif // MFEM_PFF_SPECTRAL_DECOMPOSITION_HPP
