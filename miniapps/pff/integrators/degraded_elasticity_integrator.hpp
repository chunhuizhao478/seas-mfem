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

#ifndef MFEM_PFF_DEGRADED_ELASTICITY_INTEGRATOR_HPP
#define MFEM_PFF_DEGRADED_ELASTICITY_INTEGRATOR_HPP

#include "mfem.hpp"
#include "../materials/pff_material.hpp"
#include "../materials/degradation_function.hpp"
#include "../materials/spectral_decomposition.hpp"

namespace mfem
{

namespace pff
{

/** @brief Nonlinear integrator for degraded elasticity.
 *
 * Implements the weak form:
 *   ∫_Ω σ(u,d) : ε(v) dΩ
 *
 * where the stress is computed with tension-compression split:
 *   σ = g(d) σ⁺ + σ⁻
 *   σ⁺ = λ⟨tr(ε)⟩₊ I + 2μ ε⁺
 *   σ⁻ = λ⟨tr(ε)⟩₋ I + 2μ ε⁻
 *
 * The positive (tensile) part is degraded by g(d), while the negative
 * (compressive) part remains undegraded.
 */
class DegradedElasticityIntegrator : public NonlinearFormIntegrator
{
public:
   /** @brief Construct degraded elasticity integrator.
    *
    * @param[in] mat Material parameters
    * @param[in] damage Damage field (scalar GridFunction)
    */
   DegradedElasticityIntegrator(const PFFMaterialParameters &mat,
                                const GridFunction &damage)
      : mat_(mat), damage_(damage), g_(mat.eta, mat.p) {}

   /** @brief Assemble element residual vector.
    *
    * Computes: r_i = ∫_Ω σ : ε(φ_i) dΩ
    *
    * @param[in] el Finite element
    * @param[in] Tr Element transformation
    * @param[in] elfun Element DOF values for u (vector field)
    * @param[out] elvect Element residual vector
    */
   void AssembleElementVector(const FiniteElement &el,
                              ElementTransformation &Tr,
                              const Vector &elfun,
                              Vector &elvect) override
   {
      int nd = el.GetDof();
      int dim = el.GetDim();

      elvect.SetSize(nd * dim);
      elvect = 0.0;

      const IntegrationRule *ir = IntRule;
      if (ir == nullptr)
      {
         int order = 2 * el.GetOrder() + Tr.OrderGrad(&el);
         ir = &IntRules.Get(el.GetGeomType(), order);
      }

      DenseMatrix dshape(nd, dim);
      DenseMatrix dshapedxt(nd, dim);
      DenseMatrix grad_u(dim, dim);
      DenseMatrix strain(dim, dim);
      DenseMatrix stress(dim, dim);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);
         Tr.SetIntPoint(&ip);

         el.CalcDShape(ip, dshape);
         Mult(dshape, Tr.InverseJacobian(), dshapedxt);

         // Compute gradient of u
         ComputeGradient(elfun, dshapedxt, nd, dim, grad_u);

         // Compute strain: ε = 0.5 * (∇u + ∇uᵀ)
         ComputeStrain(grad_u, strain);

         // Evaluate damage at quadrature point
         real_t d_val = damage_.GetValue(Tr, ip);
         d_val = std::max(0.0, std::min(1.0, d_val));

         // Compute degraded stress
         ComputeStress(strain, d_val, dim, stress);

         real_t w = ip.weight * Tr.Weight();

         // Add contribution: ∫ σ : ε(φ_i) dΩ
         // For each test function i and component k
         for (int j = 0; j < nd; j++)
         {
            for (int k = 0; k < dim; k++)
            {
               real_t contrib = 0.0;
               for (int l = 0; l < dim; l++)
               {
                  // σ_kl * ∂φ_j/∂x_l (symmetric stress)
                  contrib += stress(k, l) * dshapedxt(j, l);
               }
               elvect(j + k * nd) += w * contrib;
            }
         }
      }
   }

   /** @brief Assemble element Jacobian matrix.
    *
    * Computes: K_ij = ∫_Ω C : ε(φ_i) : ε(φ_j) dΩ
    * where C is the tangent stiffness tensor.
    *
    * @param[in] el Finite element
    * @param[in] Tr Element transformation
    * @param[in] elfun Element DOF values for u
    * @param[out] elmat Element Jacobian matrix
    */
   void AssembleElementGrad(const FiniteElement &el,
                            ElementTransformation &Tr,
                            const Vector &elfun,
                            DenseMatrix &elmat) override
   {
      int nd = el.GetDof();
      int dim = el.GetDim();
      int ndim = nd * dim;

      elmat.SetSize(ndim, ndim);
      elmat = 0.0;

      const IntegrationRule *ir = IntRule;
      if (ir == nullptr)
      {
         int order = 2 * el.GetOrder() + Tr.OrderGrad(&el);
         ir = &IntRules.Get(el.GetGeomType(), order);
      }

      DenseMatrix dshape(nd, dim);
      DenseMatrix dshapedxt(nd, dim);
      DenseMatrix grad_u(dim, dim);
      DenseMatrix strain(dim, dim);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);
         Tr.SetIntPoint(&ip);

         el.CalcDShape(ip, dshape);
         Mult(dshape, Tr.InverseJacobian(), dshapedxt);

         // Compute gradient of u
         ComputeGradient(elfun, dshapedxt, nd, dim, grad_u);

         // Compute strain
         ComputeStrain(grad_u, strain);

         // Evaluate damage at quadrature point
         real_t d_val = damage_.GetValue(Tr, ip);
         d_val = std::max(0.0, std::min(1.0, d_val));

         // Compute tangent stiffness coefficients
         real_t g = g_.Eval(d_val);

         real_t w = ip.weight * Tr.Weight();

         // Simplified: use isotropic elasticity with degradation
         // C_ijkl = g * (λ δ_ij δ_kl + μ (δ_ik δ_jl + δ_il δ_jk))
         // For simplicity, we use full degradation (no spectral split in Jacobian)
         // This gives a consistent but simplified tangent

         real_t lambda_eff = g * mat_.lambda;
         real_t mu_eff = g * mat_.mu;

         // Add contribution to element matrix
         // K_{(j,k),(m,n)} = ∫ C_klmn * ∂φ_j/∂x_l * ∂φ_m/∂x_n dΩ
         for (int j = 0; j < nd; j++)
         {
            for (int k = 0; k < dim; k++)
            {
               int row = j + k * nd;
               for (int m = 0; m < nd; m++)
               {
                  for (int n = 0; n < dim; n++)
                  {
                     int col = m + n * nd;

                     real_t val = 0.0;

                     // Lambda term: λ * ∂φ_j/∂x_k * ∂φ_m/∂x_n
                     val += lambda_eff * dshapedxt(j, k) * dshapedxt(m, n);

                     // Mu terms: μ * (∂φ_j/∂x_l * ∂φ_m/∂x_l * δ_kn + ∂φ_j/∂x_n * ∂φ_m/∂x_k)
                     for (int l = 0; l < dim; l++)
                     {
                        if (k == n)
                        {
                           val += mu_eff * dshapedxt(j, l) * dshapedxt(m, l);
                        }
                     }
                     val += mu_eff * dshapedxt(j, n) * dshapedxt(m, k);

                     elmat(row, col) += w * val;
                  }
               }
            }
         }
      }
   }

   /** @brief Compute element energy.
    *
    * Energy: ∫_Ω [g(d) ψ⁺(ε) + ψ⁻(ε)] dΩ
    *
    * @param[in] el Finite element
    * @param[in] Tr Element transformation
    * @param[in] elfun Element DOF values for u
    * @return Element energy
    */
   real_t GetElementEnergy(const FiniteElement &el,
                           ElementTransformation &Tr,
                           const Vector &elfun) override
   {
      int nd = el.GetDof();
      int dim = el.GetDim();
      real_t energy = 0.0;

      const IntegrationRule *ir = IntRule;
      if (ir == nullptr)
      {
         int order = 2 * el.GetOrder() + Tr.OrderGrad(&el);
         ir = &IntRules.Get(el.GetGeomType(), order);
      }

      DenseMatrix dshape(nd, dim);
      DenseMatrix dshapedxt(nd, dim);
      DenseMatrix grad_u(dim, dim);
      DenseMatrix strain(dim, dim);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);
         Tr.SetIntPoint(&ip);

         el.CalcDShape(ip, dshape);
         Mult(dshape, Tr.InverseJacobian(), dshapedxt);

         // Compute gradient of u
         ComputeGradient(elfun, dshapedxt, nd, dim, grad_u);

         // Compute strain
         ComputeStrain(grad_u, strain);

         // Evaluate damage
         real_t d_val = damage_.GetValue(Tr, ip);
         d_val = std::max(0.0, std::min(1.0, d_val));
         real_t g = g_.Eval(d_val);

         // Compute strain energy density
         real_t psi_pos, psi_neg;
         ComputeStrainEnergy(strain, dim, psi_pos, psi_neg);

         real_t w = ip.weight * Tr.Weight();
         energy += w * (g * psi_pos + psi_neg);
      }

      return energy;
   }

   /** @brief Compute positive strain energy at quadrature points.
    *
    * This is used as the driving force for damage evolution.
    *
    * @param[in] el Finite element
    * @param[in] Tr Element transformation
    * @param[in] elfun Element DOF values for u
    * @param[out] psi_active Positive strain energy at each quadrature point
    */
   void ComputePositiveStrainEnergy(const FiniteElement &el,
                                    ElementTransformation &Tr,
                                    const Vector &elfun,
                                    Vector &psi_active)
   {
      int nd = el.GetDof();
      int dim = el.GetDim();

      const IntegrationRule *ir = IntRule;
      if (ir == nullptr)
      {
         int order = 2 * el.GetOrder() + Tr.OrderGrad(&el);
         ir = &IntRules.Get(el.GetGeomType(), order);
      }

      int nqp = ir->GetNPoints();
      psi_active.SetSize(nqp);

      DenseMatrix dshape(nd, dim);
      DenseMatrix dshapedxt(nd, dim);
      DenseMatrix grad_u(dim, dim);
      DenseMatrix strain(dim, dim);

      for (int i = 0; i < nqp; i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);
         Tr.SetIntPoint(&ip);

         el.CalcDShape(ip, dshape);
         Mult(dshape, Tr.InverseJacobian(), dshapedxt);

         // Compute gradient of u
         ComputeGradient(elfun, dshapedxt, nd, dim, grad_u);

         // Compute strain
         ComputeStrain(grad_u, strain);

         // Compute positive strain energy
         real_t psi_pos, psi_neg;
         ComputeStrainEnergy(strain, dim, psi_pos, psi_neg);

         psi_active(i) = psi_pos;
      }
   }

private:
   const PFFMaterialParameters &mat_;
   const GridFunction &damage_;
   DegradationFunction g_;

   /// Compute gradient of u from element DOFs
   void ComputeGradient(const Vector &elfun, const DenseMatrix &dshapedxt,
                        int nd, int dim, DenseMatrix &grad_u)
   {
      grad_u.SetSize(dim, dim);
      grad_u = 0.0;

      for (int i = 0; i < dim; i++)  // component of u
      {
         for (int j = 0; j < dim; j++)  // derivative direction
         {
            for (int k = 0; k < nd; k++)
            {
               grad_u(i, j) += elfun(k + i * nd) * dshapedxt(k, j);
            }
         }
      }
   }

   /// Compute symmetric strain from gradient
   void ComputeStrain(const DenseMatrix &grad_u, DenseMatrix &strain)
   {
      int dim = grad_u.Height();
      strain.SetSize(dim, dim);

      for (int i = 0; i < dim; i++)
      {
         for (int j = 0; j < dim; j++)
         {
            strain(i, j) = 0.5 * (grad_u(i, j) + grad_u(j, i));
         }
      }
   }

   /// Compute eigenvalues of 2D symmetric matrix
   void ComputeEigenvalues2D(const DenseMatrix &A, real_t &lam1, real_t &lam2)
   {
      real_t a11 = A(0, 0), a22 = A(1, 1), a12 = A(0, 1);
      real_t trace = a11 + a22;
      real_t det = a11 * a22 - a12 * a12;
      real_t disc = std::sqrt(std::max(0.0, 0.25 * trace * trace - det));
      lam1 = 0.5 * trace - disc;  // smaller eigenvalue
      lam2 = 0.5 * trace + disc;  // larger eigenvalue
   }

   /// Compute stress from strain and damage using SPECTRAL decomposition
   /// σ = g(d) * σ⁺ + σ⁻
   /// Only tensile part is degraded; compressive part retains full stiffness
   void ComputeStress(const DenseMatrix &strain, real_t d, int dim,
                      DenseMatrix &stress)
   {
      stress.SetSize(dim, dim);
      stress = 0.0;

      // Degradation factor
      real_t g = g_.Eval(d);

      real_t tr_strain = 0.0;
      for (int i = 0; i < dim; i++)
      {
         tr_strain += strain(i, i);
      }

      if (dim == 2)
      {
         // 2D spectral decomposition
         real_t lam1, lam2;
         ComputeEigenvalues2D(strain, lam1, lam2);

         // Positive and negative eigenvalues
         real_t lam1_pos = std::max(lam1, 0.0);
         real_t lam2_pos = std::max(lam2, 0.0);
         real_t lam1_neg = std::min(lam1, 0.0);
         real_t lam2_neg = std::min(lam2, 0.0);

         // Positive and negative trace
         real_t tr_pos = std::max(tr_strain, 0.0);
         real_t tr_neg = std::min(tr_strain, 0.0);

         // Compute eigenvectors for reconstruction
         real_t a11 = strain(0, 0), a22 = strain(1, 1), a12 = strain(0, 1);

         // Eigenvector for lam2 (larger eigenvalue)
         real_t n2x, n2y;
         if (std::abs(a12) > 1e-14)
         {
            n2x = lam2 - a22;
            n2y = a12;
            real_t norm = std::sqrt(n2x * n2x + n2y * n2y);
            n2x /= norm;
            n2y /= norm;
         }
         else
         {
            n2x = (a11 >= a22) ? 1.0 : 0.0;
            n2y = (a11 >= a22) ? 0.0 : 1.0;
         }
         // Eigenvector for lam1 (perpendicular)
         real_t n1x = -n2y, n1y = n2x;

         // Reconstruct positive strain tensor: ε⁺ = Σ <λᵢ>₊ nᵢ⊗nᵢ
         DenseMatrix strain_pos(2, 2), strain_neg(2, 2);
         strain_pos(0, 0) = lam1_pos * n1x * n1x + lam2_pos * n2x * n2x;
         strain_pos(0, 1) = lam1_pos * n1x * n1y + lam2_pos * n2x * n2y;
         strain_pos(1, 0) = strain_pos(0, 1);
         strain_pos(1, 1) = lam1_pos * n1y * n1y + lam2_pos * n2y * n2y;

         strain_neg(0, 0) = lam1_neg * n1x * n1x + lam2_neg * n2x * n2x;
         strain_neg(0, 1) = lam1_neg * n1x * n1y + lam2_neg * n2x * n2y;
         strain_neg(1, 0) = strain_neg(0, 1);
         strain_neg(1, 1) = lam1_neg * n1y * n1y + lam2_neg * n2y * n2y;

         // σ⁺ = λ <tr(ε)>₊ I + 2μ ε⁺  (DEGRADED by g)
         // σ⁻ = λ <tr(ε)>₋ I + 2μ ε⁻  (NOT degraded)
         for (int i = 0; i < 2; i++)
         {
            for (int j = 0; j < 2; j++)
            {
               real_t sigma_pos = 2.0 * mat_.mu * strain_pos(i, j);
               real_t sigma_neg = 2.0 * mat_.mu * strain_neg(i, j);
               if (i == j)
               {
                  sigma_pos += mat_.lambda * tr_pos;
                  sigma_neg += mat_.lambda * tr_neg;
               }
               stress(i, j) = g * sigma_pos + sigma_neg;
            }
         }
      }
      else // dim == 3
      {
         // For 3D, use eigenvalue decomposition
         DenseMatrix strain_copy(strain);
         Vector eig_vals(3);
         DenseMatrix eig_vecs(3, 3);
         strain_copy.Eigenvalues(eig_vals, eig_vecs);

         // Positive and negative eigenvalues
         Vector lam_pos(3), lam_neg(3);
         for (int i = 0; i < 3; i++)
         {
            lam_pos(i) = std::max(eig_vals(i), 0.0);
            lam_neg(i) = std::min(eig_vals(i), 0.0);
         }

         real_t tr_pos = std::max(tr_strain, 0.0);
         real_t tr_neg = std::min(tr_strain, 0.0);

         // Reconstruct strain tensors
         DenseMatrix strain_pos(3, 3), strain_neg(3, 3);
         strain_pos = 0.0;
         strain_neg = 0.0;
         for (int k = 0; k < 3; k++)
         {
            for (int i = 0; i < 3; i++)
            {
               for (int j = 0; j < 3; j++)
               {
                  strain_pos(i, j) += lam_pos(k) * eig_vecs(i, k) * eig_vecs(j, k);
                  strain_neg(i, j) += lam_neg(k) * eig_vecs(i, k) * eig_vecs(j, k);
               }
            }
         }

         // Compute stress with spectral split
         for (int i = 0; i < 3; i++)
         {
            for (int j = 0; j < 3; j++)
            {
               real_t sigma_pos = 2.0 * mat_.mu * strain_pos(i, j);
               real_t sigma_neg = 2.0 * mat_.mu * strain_neg(i, j);
               if (i == j)
               {
                  sigma_pos += mat_.lambda * tr_pos;
                  sigma_neg += mat_.lambda * tr_neg;
               }
               stress(i, j) = g * sigma_pos + sigma_neg;
            }
         }
      }
   }

   /// Compute positive and negative strain energy densities using SPECTRAL decomposition
   /// ψ⁺ = λ/2 <tr(ε)>₊² + μ ε⁺:ε⁺
   /// ψ⁻ = λ/2 <tr(ε)>₋² + μ ε⁻:ε⁻
   void ComputeStrainEnergy(const DenseMatrix &strain, int dim,
                            real_t &psi_pos, real_t &psi_neg)
   {
      real_t tr_strain = 0.0;
      for (int i = 0; i < dim; i++)
      {
         tr_strain += strain(i, i);
      }

      real_t tr_pos = std::max(tr_strain, 0.0);
      real_t tr_neg = std::min(tr_strain, 0.0);

      if (dim == 2)
      {
         real_t lam1, lam2;
         ComputeEigenvalues2D(strain, lam1, lam2);

         real_t lam1_pos = std::max(lam1, 0.0);
         real_t lam2_pos = std::max(lam2, 0.0);
         real_t lam1_neg = std::min(lam1, 0.0);
         real_t lam2_neg = std::min(lam2, 0.0);

         // ε⁺:ε⁺ = Σ <λᵢ>₊²
         real_t eps_pos_sq = lam1_pos * lam1_pos + lam2_pos * lam2_pos;
         real_t eps_neg_sq = lam1_neg * lam1_neg + lam2_neg * lam2_neg;

         psi_pos = 0.5 * mat_.lambda * tr_pos * tr_pos + mat_.mu * eps_pos_sq;
         psi_neg = 0.5 * mat_.lambda * tr_neg * tr_neg + mat_.mu * eps_neg_sq;
      }
      else // dim == 3
      {
         DenseMatrix strain_copy(strain);
         Vector eig_vals(3);
         DenseMatrix eig_vecs(3, 3);
         strain_copy.Eigenvalues(eig_vals, eig_vecs);

         real_t eps_pos_sq = 0.0, eps_neg_sq = 0.0;
         for (int i = 0; i < 3; i++)
         {
            real_t lam_pos = std::max(eig_vals(i), 0.0);
            real_t lam_neg = std::min(eig_vals(i), 0.0);
            eps_pos_sq += lam_pos * lam_pos;
            eps_neg_sq += lam_neg * lam_neg;
         }

         psi_pos = 0.5 * mat_.lambda * tr_pos * tr_pos + mat_.mu * eps_pos_sq;
         psi_neg = 0.5 * mat_.lambda * tr_neg * tr_neg + mat_.mu * eps_neg_sq;
      }
   }
};

} // namespace pff

} // namespace mfem

#endif // MFEM_PFF_DEGRADED_ELASTICITY_INTEGRATOR_HPP
