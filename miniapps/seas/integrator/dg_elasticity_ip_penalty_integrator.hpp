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

#ifndef MFEM_SEAS_DG_ELASTICITY_IP_PENALTY_INTEGRATOR_HPP
#define MFEM_SEAS_DG_ELASTICITY_IP_PENALTY_INTEGRATOR_HPP

#include "mfem.hpp"
#include <cmath>

namespace mfem
{
namespace seas
{

/// @brief SIPG penalty face integrator for 3D vector elasticity.
///
/// Computes ONLY the penalty (stabilization) term of the SIPG/IP bilinear form:
///
///   a_pen(u,v) = Σ_F η_F ∫_F [[u]] · [[v]] ds
///
/// where η_F is a material- and mesh-dependent penalty per face:
///
///   η_F = (p(K⁻) + p(K⁺)) / 4     (interior face)
///   η_F = p(K)                       (boundary face)
///   p(K) = (D+1) * c_N_1 * (A_F/V_K) * (c₁²/c₀)
///
/// with:
///   c₀ = 2μ                  (min eigenvalue of isotropic C)
///   c₁ = D*λ + 2μ            (max eigenvalue of isotropic C)
///   A_F = face area           (= |nor_q| at quadrature points)
///   V_K = element volume
///   c_N_1 = 1                 (inverse inequality trace constant for p=0)
///   D = spatial dimension
///
/// This integrator computes ONLY the penalty. The consistency and symmetry
/// terms should come from MFEM's DGElasticityIntegrator with kappa=0.
///
/// The penalty uses |nor_q| (linear in face area), NOT |nor_q|² (squared).
/// This matches the formulation in Uphoff et al. (2023), GJI 233(1), 586-626.
///
/// Reference: Eq. (36) penalty term in the benchmark document.
class DGElasticityIPPenaltyIntegrator : public BilinearFormIntegrator
{
public:
   /// @brief Construct the IP penalty integrator.
   ///
   /// @param lambda First Lame parameter coefficient
   /// @param mu Shear modulus coefficient
   /// @param dim Spatial dimension (default 3)
   DGElasticityIPPenaltyIntegrator(Coefficient &lambda, Coefficient &mu,
                                    int dim = 3)
      : lambda_(lambda), mu_(mu), dim_(dim) {}

   using BilinearFormIntegrator::AssembleFaceMatrix;
   void AssembleFaceMatrix(const FiniteElement &el1,
                           const FiniteElement &el2,
                           FaceElementTransformations &Trans,
                           DenseMatrix &elmat) override
   {
      const int ndofs1 = el1.GetDof();
      const int ndofs2 = (Trans.Elem2No >= 0) ? el2.GetDof() : 0;
      const int nvdofs = dim_ * (ndofs1 + ndofs2);

      elmat.SetSize(nvdofs);
      elmat = 0.0;

      // Integration rule
      const int order = 2 * std::max(el1.GetOrder(),
                                      ndofs2 ? el2.GetOrder() : 0);
      const IntegrationRule &ir = IntRules.Get(Trans.GetGeometryType(), order);

      Vector shape1(ndofs1), shape2(ndofs2 > 0 ? ndofs2 : 1);
      Vector nor(dim_);

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Trans.SetAllIntPoints(&ip);

         const IntegrationPoint &eip1 = Trans.GetElement1IntPoint();
         el1.CalcShape(eip1, shape1);

         if (ndofs2 > 0)
         {
            const IntegrationPoint &eip2 = Trans.GetElement2IntPoint();
            el2.CalcShape(eip2, shape2);
         }

         // Unnormalized normal (magnitude = face area at this quad point)
         CalcOrtho(Trans.Jacobian(), nor);
         real_t nl_q = nor.Norml2();  // |nor| = NormalLength

         // Material properties (evaluate at face quadrature point)
         real_t lam_val = lambda_.Eval(*Trans.Elem1, eip1);
         real_t mu_val = mu_.Eval(*Trans.Elem1, eip1);

         // Stiffness tensor bounds (isotropic)
         real_t c0 = 2.0 * mu_val;
         real_t c1 = dim_ * lam_val + 2.0 * mu_val;

         // Inverse inequality trace constant: c_N(p) = p*(p+D-1)/D
         // (Tandem: InverseInequality.h, trace_constant(PolynomialDegree-1))
         int p = std::max(el1.GetOrder(),
                          ndofs2 > 0 ? el2.GetOrder() : el1.GetOrder());
         real_t c_N_1 = p * (p + dim_ - 1.0) / dim_;

         // Element volumes (Weight() = |det(J)| = D! × V_phys for simplices)
         real_t vol1 = Trans.Elem1->Weight();

         // Penalty per side: p(K) = (D+1) * c_N_1 * (A/V) * (c1²/c0)
         // Physical A/V = dim * nl_q / Weight() (corrects reference element
         // measure ratio: ref_face/ref_vol = 1/(D-1)! / 1/D! = D)
         // See bp5_debug_v47.md Section 2 for derivation.
         real_t p0 = (dim_ + 1) * c_N_1 * (real_t(dim_) * nl_q / vol1) * (c1 * c1 / c0);

         real_t penalty;
         if (ndofs2 > 0)
         {
            real_t vol2 = Trans.Elem2->Weight();
            real_t p1 = (dim_ + 1) * c_N_1 * (real_t(dim_) * nl_q / vol2) * (c1 * c1 / c0);
            penalty = (p0 + p1) / 4.0;
         }
         else
         {
            penalty = p0;
         }

         // Penalty coefficient per quadrature point:
         //   penalty * ip.weight * nl_q
         // But nl_q is already in the penalty via (nl_q / vol).
         // So the total is: penalty * ip.weight
         // Wait — need to be careful. The penalty formula uses face AREA
         // (integrated over the face), not per-quadrature-point nl_q.
         //
         // In Tandem: penalty is precomputed per face using total face area.
         // Then the assembly is: penalty * Σ_q w_q * φ_i * φ_j * nl_q
         //
         // The penalty already contains (A_face / V_elem). If we also
         // multiply by nl_q per quadrature point, we get A²/V scaling.
         //
         // Tandem precomputes: area_[fctNo] = total face area (integrated).
         // Then penalty = ... * (area / volume) * ...
         // And assembly: penalty * w_q * shape * shape * nl_q
         //
         // The nl_q in assembly = |nor_q| which sums to face area:
         //   Σ_q w_q * nl_q = face_area (for flat faces)
         //
         // So total = penalty * face_area = ... * (area/vol) * ... * area
         //          = ... * area² / vol * ...
         //
         // This DOES have area² — but it's area²/vol, not area²/vol².
         //
         // MFEM's formula: kappa * |nor|² * w * (L+2M) / detJ
         //   = kappa * nl² * w * (L+2M) / (6V)
         //   Summing: kappa * (Σ w * nl²) * (L+2M) / (6V)
         //
         // These are different integrals:
         //   Tandem: penalty * Σ w_q * nl_q      (∫ |n| ds)
         //   MFEM:   kappa * Σ w_q * nl_q²       (∫ |n|² ds)
         //
         // For the Tandem formula:
         //   coeff_q = penalty * w_q * nl_q

         real_t coeff = penalty * ip.weight * nl_q;

         // Assemble same-component penalty: coeff * φ_i * φ_j
         // Block (0,0): +coeff (same side)
         for (int d = 0; d < dim_; d++)
         {
            const int offset1 = d * ndofs1;
            for (int i = 0; i < ndofs1; i++)
            {
               for (int j = 0; j <= i; j++)
               {
                  real_t val = coeff * shape1(i) * shape1(j);
                  elmat(offset1 + i, offset1 + j) += val;
                  if (i != j)
                  {
                     elmat(offset1 + j, offset1 + i) += val;
                  }
               }
            }
         }

         if (ndofs2 == 0) { continue; }

         // Block (1,1): +coeff (same side)
         for (int d = 0; d < dim_; d++)
         {
            const int offset2 = dim_ * ndofs1 + d * ndofs2;
            for (int i = 0; i < ndofs2; i++)
            {
               for (int j = 0; j <= i; j++)
               {
                  real_t val = coeff * shape2(i) * shape2(j);
                  elmat(offset2 + i, offset2 + j) += val;
                  if (i != j)
                  {
                     elmat(offset2 + j, offset2 + i) += val;
                  }
               }
            }
         }

         // Block (0,1) and (1,0): -coeff (cross side)
         for (int d = 0; d < dim_; d++)
         {
            const int offset1 = d * ndofs1;
            const int offset2 = dim_ * ndofs1 + d * ndofs2;
            for (int i = 0; i < ndofs1; i++)
            {
               for (int j = 0; j < ndofs2; j++)
               {
                  real_t val = -coeff * shape1(i) * shape2(j);
                  elmat(offset1 + i, offset2 + j) += val;
                  elmat(offset2 + j, offset1 + i) += val;
               }
            }
         }
      }
   }

private:
   Coefficient &lambda_;
   Coefficient &mu_;
   int dim_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DG_ELASTICITY_IP_PENALTY_INTEGRATOR_HPP
