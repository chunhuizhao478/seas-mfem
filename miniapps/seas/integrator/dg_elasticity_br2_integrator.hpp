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

#ifndef MFEM_SEAS_DG_ELASTICITY_BR2_INTEGRATOR_HPP
#define MFEM_SEAS_DG_ELASTICITY_BR2_INTEGRATOR_HPP

#include "mfem.hpp"
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief BR2 Interior Face Integrator for 3D Linear Elasticity
///
/// Generalizes BR2InteriorFaceIntegrator from scalar Laplace to vector
/// elasticity, following Tandem's Elasticity.cpp lift_skeleton kernel.
///
/// The DG bilinear form on each interior face has three terms:
///   a(u,v) = -∫_F {{σ(u)·n}} · [[v]] ds              (consistency)
///          - ∫_F {{σ(v)·n}} · [[u]] ds              (symmetry, SIPG)
///          + σ_BR2 ∫_F C:R_h(u) : R_h(v) ds          (BR2 lifting penalty)
///   where [[·]] = (·)⁻ − (·)⁺ (normal points from K⁻ to K⁺)
///
/// The elasticity tensor coupling in the BR2 lifting is encoded via
/// the test_normal operator (Tandem's elasticity.py lines 118-120):
///   test_normal(x)_{iu,sq} = λ·δ_us·n_i + μ·(δ_iu·n_s + δ_is·n_u)
///
/// For scalar Laplace, test_normal reduces to K·n_i, recovering the
/// scalar BR2InteriorFaceIntegrator.
///
/// The BR2 penalty parameter is σ = num_faces_per_element:
///   6 for hexahedra (2·dim), 4 for tetrahedra (dim+1).
///
class DGElasticityBR2Integrator : public BilinearFormIntegrator
{
public:
   /// @brief Construct BR2 interior face integrator for elasticity
   ///
   /// @param lambda First Lame parameter coefficient
   /// @param mu Shear modulus coefficient
   /// @param epsilon SIPG sign (-1 for symmetric)
   /// @param elem_mass_inv Precomputed scalar element mass matrix inverses
   /// @param dim Spatial dimension (3 for 3D elasticity)
   DGElasticityBR2Integrator(Coefficient &lambda, Coefficient &mu,
                              real_t epsilon,
                              const std::vector<DenseMatrix> &elem_mass_inv,
                              int dim = 3)
      : lambda_(lambda), mu_(mu), epsilon_(epsilon),
        elem_mass_inv_(elem_mass_inv), dim_(dim)
   {}

   void AssembleFaceMatrix(const FiniteElement &el1,
                           const FiniteElement &el2,
                           FaceElementTransformations &Trans,
                           DenseMatrix &elmat) override;

protected:
   Coefficient &lambda_, &mu_;
   real_t epsilon_;
   const std::vector<DenseMatrix> &elem_mass_inv_;
   int dim_;

   /// Compute test_normal tensor at a quadrature point.
   /// test_normal_{iu,sq} = λ·δ_us·n_i + μ·(δ_iu·n_s + δ_is·n_u)
   /// Stored as flat array [i*dim*dim*dim + u*dim*dim + s*dim + q_unused]
   /// but we compute it on the fly for each (i,u,s) given normal n.
   inline real_t TestNormal(real_t lam, real_t mu_val,
                            const Vector &n,
                            int i, int u, int s) const
   {
      real_t val = lam * (u == s ? 1.0 : 0.0) * n(i);
      val += mu_val * ((i == u ? 1.0 : 0.0) * n(s) +
                       (i == s ? 1.0 : 0.0) * n(u));
      return val;
   }
};

/// @brief BR2 Boundary Face Integrator for 3D Linear Elasticity (Dirichlet)
///
/// Inherits from DGElasticityBR2Integrator to share TestNormal and members.
/// Full factor (no 0.5 averaging) on lifting and traction terms.
class DGElasticityBR2BoundaryIntegrator : public DGElasticityBR2Integrator
{
public:
   DGElasticityBR2BoundaryIntegrator(Coefficient &lambda, Coefficient &mu,
                                      real_t epsilon,
                                      const std::vector<DenseMatrix> &elem_mass_inv,
                                      int dim = 3)
      : DGElasticityBR2Integrator(lambda, mu, epsilon, elem_mass_inv, dim)
   {}

   void AssembleFaceMatrix(const FiniteElement &el1,
                           const FiniteElement &el2,
                           FaceElementTransformations &Trans,
                           DenseMatrix &elmat) override;
};

// ============================================================================
// Implementation: Interior Face
// ============================================================================

inline void DGElasticityBR2Integrator::AssembleFaceMatrix(
   const FiniteElement &el1,
   const FiniteElement &el2,
   FaceElementTransformations &Trans,
   DenseMatrix &elmat)
{
   int ndof1 = el1.GetDof();
   int ndof2 = el2.GetDof();
   int vdof1 = ndof1 * dim_;
   int vdof2 = ndof2 * dim_;
   int vdofs = vdof1 + vdof2;

   elmat.SetSize(vdofs);
   elmat = 0.0;

   // Detect element type for penalty (local variable for thread safety)
   Geometry::Type geom = el1.GetGeomType();
   const real_t sigma = (geom == Geometry::TETRAHEDRON)
                            ? real_t(dim_ + 1)   // 4 faces for tet
                            : real_t(2 * dim_);  // 6 faces for hex

   int elem1 = Trans.Elem1No;
   int elem2 = Trans.Elem2No;

   const DenseMatrix &Minv1 = elem_mass_inv_[elem1];
   const DenseMatrix &Minv2 = elem_mass_inv_[elem2];

   // Integration order
   int order = 2 * std::max(el1.GetOrder(), el2.GetOrder()) + 1;
   const IntegrationRule &ir = IntRules.Get(Trans.FaceGeom, order);
   int nqp = ir.GetNPoints();

   // Precompute shapes, gradients, normals at all quadrature points
   DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
   // dshapes stored as (ndof * dim, nqp) — physical gradients
   DenseMatrix dshapes1_phys(ndof1 * dim_, nqp), dshapes2_phys(ndof2 * dim_, nqp);
   Vector nor_all(dim_ * nqp);
   Vector w_all(nqp);
   Vector lam_all(nqp), mu_all(nqp);
   // Store 1/detJ for each element at each q
   Vector invdetJ1_all(nqp), invdetJ2_all(nqp);

   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      Trans.SetAllIntPoints(&ip);

      const IntegrationPoint &eip1 = Trans.GetElement1IntPoint();
      const IntegrationPoint &eip2 = Trans.GetElement2IntPoint();

      // Shape functions
      Vector shape1_q(shapes1.GetColumn(q), ndof1);
      Vector shape2_q(shapes2.GetColumn(q), ndof2);
      el1.CalcShape(eip1, shape1_q);
      el2.CalcShape(eip2, shape2_q);

      // Gradients in reference, then transform via adjugate
      DenseMatrix dshape1_ref(ndof1, dim_), dshape2_ref(ndof2, dim_);
      el1.CalcDShape(eip1, dshape1_ref);
      el2.CalcDShape(eip2, dshape2_ref);

      DenseMatrix adjJ1(dim_), adjJ2(dim_);
      CalcAdjugate(Trans.Elem1->Jacobian(), adjJ1);
      CalcAdjugate(Trans.Elem2->Jacobian(), adjJ2);

      DenseMatrix dshape1_phys(ndof1, dim_), dshape2_phys(ndof2, dim_);
      Mult(dshape1_ref, adjJ1, dshape1_phys);
      Mult(dshape2_ref, adjJ2, dshape2_phys);

      // Store flattened physical gradients
      for (int k = 0; k < ndof1; k++)
         for (int d = 0; d < dim_; d++)
            dshapes1_phys(k * dim_ + d, q) = dshape1_phys(k, d);
      for (int k = 0; k < ndof2; k++)
         for (int d = 0; d < dim_; d++)
            dshapes2_phys(k * dim_ + d, q) = dshape2_phys(k, d);

      // Normal vector (unnormalized, |nor| = det(J_face))
      Vector nor_q(&nor_all[q * dim_], dim_);
      CalcOrtho(Trans.Jacobian(), nor_q);

      // Quadrature weight
      w_all[q] = ip.weight;

      // Material coefficients
      lam_all[q] = lambda_.Eval(*Trans.Face, ip);
      mu_all[q] = mu_.Eval(*Trans.Face, ip);

      // 1/detJ for each element
      invdetJ1_all[q] = 1.0 / Trans.Elem1->Weight();
      invdetJ2_all[q] = 1.0 / Trans.Elem2->Weight();
   }

   // =========================================================================
   // Compute BR2 Lifting with elasticity tensor coupling
   //
   // For each source element y, source DOF l, component u:
   //   Lift[elem][l,s,m] = 0.5 * Minv[elem][m,o] * sum_q E[elem][o,q] * E[y][l,q] * n[s,q] * w[q]
   //   L_q[y][l,i,u,q]  = 0.5 * sum_{elem,s,m} test_normal(elem)_{iu,sq_norm} *
   //                        E[elem][m,q] * Lift[elem][l,s,m]
   //
   // For assembly, we compute the contribution to elmat block (x_elem, y_elem)
   // component (i, u) as:
   //   a_br2[(x,k,i),(y,l,u)] = sigma * sum_q w[q] * E[x][k,q] * L_q[y][l,i,u,q]
   //
   // To keep memory manageable, we compute L_q for each source DOF l on the fly.
   // =========================================================================

   // For each source element y={0,1} and source DOF l, compute face integrals
   // IntFace[elem_dest][l*dim+s, m_dest] = sum_q E_dest[m,q] * E_src[l,q] * n[s,q] * w[q]

   // Then Lift[elem][l,s,m] = 0.5 * Minv[elem] * IntFace[elem]
   // Then L_q evaluated at each q

   // We use the normalized normal for test_normal but unnormalized for lifting.
   // Actually, following the scalar BR2 exactly: the lifting uses the unnormalized
   // normal (from CalcOrtho), and L_q contracts with unnormalized normal again.
   // This gives the correct face-area scaling. The test_normal also uses the
   // unnormalized normal.

   // Helper: compute L_q for all (l, i, u, q) for a given source element.
   // Returns L_q as a (ndof_src * dim_ * dim_, nqp) matrix,
   // indexed as L_q[(l*dim+i)*dim+u, q].

   // Source = element 1 (ndof1 source DOFs)
   DenseMatrix L_q_y0(ndof1 * dim_ * dim_, nqp);
   L_q_y0 = 0.0;

   // Source = element 2 (ndof2 source DOFs)
   DenseMatrix L_q_y1(ndof2 * dim_ * dim_, nqp);
   L_q_y1 = 0.0;

   // Compute L_q for source y=0 (element 1)
   {
      // IntFace[dest_elem][l*dim+s, m] = sum_q E_dest[m,q] * E_src[l,q] * n[s,q] * w[q]
      // For dest=elem1:
      DenseMatrix IntFace_11(ndof1 * dim_, ndof1);
      IntFace_11 = 0.0;
      // For dest=elem2:
      DenseMatrix IntFace_21(ndof1 * dim_, ndof2);
      IntFace_21 = 0.0;

      for (int q = 0; q < nqp; q++)
      {
         real_t wq = w_all[q];
         for (int l = 0; l < ndof1; l++)
         {
            real_t E_src = shapes1(l, q);
            for (int s = 0; s < dim_; s++)
            {
               real_t n_s = nor_all[q * dim_ + s];
               real_t factor = E_src * n_s * wq;
               for (int m = 0; m < ndof1; m++)
               {
                  IntFace_11(l * dim_ + s, m) += shapes1(m, q) * factor;
               }
               for (int m = 0; m < ndof2; m++)
               {
                  IntFace_21(l * dim_ + s, m) += shapes2(m, q) * factor;
               }
            }
         }
      }

      // Lift[elem][l,s,m] = 0.5 * Minv[elem][m,o] * IntFace[elem][l*dim+s, o]
      DenseMatrix Lift_1(ndof1 * dim_, ndof1); // [l*dim+s, m]
      MultABt(IntFace_11, Minv1, Lift_1);
      Lift_1 *= 0.5;

      DenseMatrix Lift_2(ndof1 * dim_, ndof2); // [l*dim+s, m]
      MultABt(IntFace_21, Minv2, Lift_2);
      Lift_2 *= 0.5;

      // L_q[y=0][l,i,u,q] = 0.5 * sum_{elem,s,m} test_normal(elem)_{iu,sq} *
      //                       E[elem][m,q] * Lift[elem][l,s,m]
      for (int q = 0; q < nqp; q++)
      {
         real_t lam = lam_all[q];
         real_t mu_val = mu_all[q];
         Vector n_q(dim_);
         for (int d = 0; d < dim_; d++) { n_q(d) = nor_all[q * dim_ + d]; }

         for (int l = 0; l < ndof1; l++)
         {
            for (int i = 0; i < dim_; i++)
            {
               for (int u = 0; u < dim_; u++)
               {
                  real_t sum = 0.0;
                  // Contribution from elem1
                  for (int s = 0; s < dim_; s++)
                  {
                     real_t tn = TestNormal(lam, mu_val, n_q, i, u, s);
                     real_t eval_lift = 0.0;
                     for (int m = 0; m < ndof1; m++)
                     {
                        eval_lift += shapes1(m, q) * Lift_1(l * dim_ + s, m);
                     }
                     sum += tn * eval_lift;
                  }
                  // Contribution from elem2
                  for (int s = 0; s < dim_; s++)
                  {
                     real_t tn = TestNormal(lam, mu_val, n_q, i, u, s);
                     real_t eval_lift = 0.0;
                     for (int m = 0; m < ndof2; m++)
                     {
                        eval_lift += shapes2(m, q) * Lift_2(l * dim_ + s, m);
                     }
                     sum += tn * eval_lift;
                  }
                  L_q_y0((l * dim_ + i) * dim_ + u, q) = 0.5 * sum;
               }
            }
         }
      }
   }

   // Compute L_q for source y=1 (element 2)
   {
      DenseMatrix IntFace_12(ndof2 * dim_, ndof1);
      IntFace_12 = 0.0;
      DenseMatrix IntFace_22(ndof2 * dim_, ndof2);
      IntFace_22 = 0.0;

      for (int q = 0; q < nqp; q++)
      {
         real_t wq = w_all[q];
         for (int l = 0; l < ndof2; l++)
         {
            real_t E_src = shapes2(l, q);
            for (int s = 0; s < dim_; s++)
            {
               real_t n_s = nor_all[q * dim_ + s];
               real_t factor = E_src * n_s * wq;
               for (int m = 0; m < ndof1; m++)
               {
                  IntFace_12(l * dim_ + s, m) += shapes1(m, q) * factor;
               }
               for (int m = 0; m < ndof2; m++)
               {
                  IntFace_22(l * dim_ + s, m) += shapes2(m, q) * factor;
               }
            }
         }
      }

      DenseMatrix Lift_1(ndof2 * dim_, ndof1);
      MultABt(IntFace_12, Minv1, Lift_1);
      Lift_1 *= 0.5;

      DenseMatrix Lift_2(ndof2 * dim_, ndof2);
      MultABt(IntFace_22, Minv2, Lift_2);
      Lift_2 *= 0.5;

      for (int q = 0; q < nqp; q++)
      {
         real_t lam = lam_all[q];
         real_t mu_val = mu_all[q];
         Vector n_q(dim_);
         for (int d = 0; d < dim_; d++) { n_q(d) = nor_all[q * dim_ + d]; }

         for (int l = 0; l < ndof2; l++)
         {
            for (int i = 0; i < dim_; i++)
            {
               for (int u = 0; u < dim_; u++)
               {
                  real_t sum = 0.0;
                  for (int s = 0; s < dim_; s++)
                  {
                     real_t tn = TestNormal(lam, mu_val, n_q, i, u, s);
                     real_t eval_lift = 0.0;
                     for (int m = 0; m < ndof1; m++)
                     {
                        eval_lift += shapes1(m, q) * Lift_1(l * dim_ + s, m);
                     }
                     sum += tn * eval_lift;
                  }
                  for (int s = 0; s < dim_; s++)
                  {
                     real_t tn = TestNormal(lam, mu_val, n_q, i, u, s);
                     real_t eval_lift = 0.0;
                     for (int m = 0; m < ndof2; m++)
                     {
                        eval_lift += shapes2(m, q) * Lift_2(l * dim_ + s, m);
                     }
                     sum += tn * eval_lift;
                  }
                  L_q_y1((l * dim_ + i) * dim_ + u, q) = 0.5 * sum;
               }
            }
         }
      }
   }

   // =========================================================================
   // Assemble the face matrix
   //
   // The face matrix has block structure:
   //   (dim*ndof1 + dim*ndof2) x (dim*ndof1 + dim*ndof2)
   //
   // For DOF (x_elem, k, i) and (y_elem, l, u):
   //   a[(x,k,i),(y,l,u)] = c0 * traction_consistency
   //                       + c1 * traction_symmetry
   //                       + c2 * BR2_lifting
   //
   // Consistency: σ_i(φ_k^{test})·n contracted with φ_l (component u)
   //   σ_j(φ_k^a)·n_s = λ·(∂φ_k/∂x_a)·n_j + μ·(δ_ja·∂φ_k/∂x_s·n_s + ∂φ_k/∂x_j·n_a)
   //   → For component i of test function with basis k:
   //     traction_k_i = λ·(∂φ_k/∂x_i)·(n·n) ... no, actually
   //     σ_{js}(φ_k e_i) = λ·(∂φ_k/∂x_i)·δ_{js} + μ·(δ_{ji}·∂φ_k/∂x_s + δ_{si}·∂φ_k/∂x_j)
   //     σ(φ_k e_i)·n_s = λ·(∂φ_k/∂x_i)·n_j + μ·(δ_{ji}·(∇φ_k·n) + ∂φ_k/∂x_j·n_i)
   //   Contracted with e_u (trial component):
   //     [σ(φ_k e_i)·n]_u = λ·(∂φ_k/∂x_i)·n_u + μ·(δ_{ui}·(∇φ_k·n) + ∂φ_k/∂x_u·n_i)
   //
   // This is exactly test_normal_{iu,sq} contracted with ∂φ_k/∂x_s:
   //     [σ(φ_k e_i)·n]_u = sum_s (∂φ_k/∂x_s) * test_normal_{iu,s}
   //                       = sum_s (∂φ_k/∂x_s) * [λ·δ_us·n_i + μ·(δ_iu·n_s + δ_is·n_u)]
   //
   // Actually let's compute it directly:
   //   [σ(φ_k e_i)·n]_u = λ·(∂φ_k/∂x_i)·n_u + μ·(δ_{iu}·(∂φ_k/∂x_s·n_s) + ∂φ_k/∂x_u·n_i)
   //
   // Sign conventions (with [[·]] = (·)⁻ − (·)⁺):
   //   c0 = coeff for [σ(φ_test)·n] · φ_trial   (symmetry term)
   //   c1 = coeff for [σ(φ_trial)·n] · φ_test   (consistency term)
   //   Block (0,0): c0=-0.5, c1=-0.5,  c2=+σ
   //   Block (0,1): c0=+0.5, c1=-0.5,  c2=-σ
   //   Block (1,0): c0=-0.5, c1=+0.5,  c2=-σ
   //   Block (1,1): c0=+0.5, c1=+0.5,  c2=+σ
   //
   // Note: c0 and c1 are named after code variables, not the standard DG terms.
   // "c0" (test traction) implements the symmetry term −∫{{σ(v)·n}}·[[u]].
   // "c1" (trial traction) implements the consistency term −∫{{σ(u)·n}}·[[v]].
   // =========================================================================

   for (int q = 0; q < nqp; q++)
   {
      real_t wq = w_all[q];
      real_t lam = lam_all[q];
      real_t mu_val = mu_all[q];
      real_t invdetJ1 = invdetJ1_all[q];
      real_t invdetJ2 = invdetJ2_all[q];

      // Normal
      Vector n_q(dim_);
      for (int d = 0; d < dim_; d++) { n_q(d) = nor_all[q * dim_ + d]; }

      // Precompute dn1[k] = ∇φ_k · n for each DOF k (element 1)
      // and dn2[k] for element 2
      Vector dn1(ndof1), dn2(ndof2);
      for (int k = 0; k < ndof1; k++)
      {
         real_t val = 0.0;
         for (int d = 0; d < dim_; d++)
         {
            val += dshapes1_phys(k * dim_ + d, q) * n_q(d);
         }
         dn1(k) = val;
      }
      for (int k = 0; k < ndof2; k++)
      {
         real_t val = 0.0;
         for (int d = 0; d < dim_; d++)
         {
            val += dshapes2_phys(k * dim_ + d, q) * n_q(d);
         }
         dn2(k) = val;
      }

      // Traction operator: [σ(φ_k e_i)·n]_u for element e
      // = λ * dφ_k/dx_i * n_u + μ * (δ_{iu} * dn_k + dφ_k/dx_u * n_i)
      // where dn_k = ∇φ_k · n

      // Block (0,0): x=elem1, y=elem1
      for (int k = 0; k < ndof1; k++)
      {
         for (int l = 0; l < ndof1; l++)
         {
            for (int i = 0; i < dim_; i++)
            {
               for (int u = 0; u < dim_; u++)
               {
                  // byNODES ordering: component i, DOF k
                  int row = i * ndof1 + k;
                  int col = u * ndof1 + l;

                  // Traction of test: [σ(φ_k e_i)·n]_u for elem1
                  real_t trac_test = lam * dshapes1_phys(k * dim_ + i, q) * n_q(u)
                     + mu_val * ((i == u ? 1.0 : 0.0) * dn1(k)
                                 + dshapes1_phys(k * dim_ + u, q) * n_q(i));

                  // Traction of trial: [σ(φ_l e_u)·n]_i for elem1
                  real_t trac_trial = lam * dshapes1_phys(l * dim_ + u, q) * n_q(i)
                     + mu_val * ((u == i ? 1.0 : 0.0) * dn1(l)
                                 + dshapes1_phys(l * dim_ + i, q) * n_q(u));

                  real_t val = 0.0;
                  // Symmetry (c0): -0.5 * trac_test * φ_l / detJ1
                  val += (-0.5) * trac_test * shapes1(l, q) * invdetJ1;
                  // Consistency (c1): ε * 0.5 * trac_trial * φ_k / detJ1
                  val += epsilon_ * 0.5 * trac_trial * shapes1(k, q) * invdetJ1;
                  // BR2 lifting: +σ * φ_k * L_q[y=0][(l,i,u),q]
                  val += sigma * shapes1(k, q)
                     * L_q_y0((l * dim_ + i) * dim_ + u, q);

                  elmat(row, col) += wq * val;
               }
            }
         }
      }

      // Block (0,1): x=elem1, y=elem2
      for (int k = 0; k < ndof1; k++)
      {
         for (int l = 0; l < ndof2; l++)
         {
            for (int i = 0; i < dim_; i++)
            {
               for (int u = 0; u < dim_; u++)
               {
                  // byNODES ordering: component i, DOF k
                  int row = i * ndof1 + k;
                  int col = vdof1 + u * ndof2 + l;

                  // Traction of test (from elem1)
                  real_t trac_test = lam * dshapes1_phys(k * dim_ + i, q) * n_q(u)
                     + mu_val * ((i == u ? 1.0 : 0.0) * dn1(k)
                                 + dshapes1_phys(k * dim_ + u, q) * n_q(i));

                  // Traction of trial (from elem2)
                  real_t trac_trial = lam * dshapes2_phys(l * dim_ + u, q) * n_q(i)
                     + mu_val * ((u == i ? 1.0 : 0.0) * dn2(l)
                                 + dshapes2_phys(l * dim_ + i, q) * n_q(u));

                  real_t val = 0.0;
                  // Symmetry (c0): +0.5 * trac_test * φ_l / detJ1
                  val += 0.5 * trac_test * shapes2(l, q) * invdetJ1;
                  // Consistency (c1): ε * 0.5 * trac_trial * φ_k / detJ2
                  val += epsilon_ * 0.5 * trac_trial * shapes1(k, q) * invdetJ2;
                  // BR2 lifting: -σ * φ_k * L_q[y=1][(l,i,u),q]
                  val += (-sigma) * shapes1(k, q)
                     * L_q_y1((l * dim_ + i) * dim_ + u, q);

                  elmat(row, col) += wq * val;
               }
            }
         }
      }

      // Block (1,0): x=elem2, y=elem1
      for (int k = 0; k < ndof2; k++)
      {
         for (int l = 0; l < ndof1; l++)
         {
            for (int i = 0; i < dim_; i++)
            {
               for (int u = 0; u < dim_; u++)
               {
                  // byNODES ordering: component i, DOF k
                  int row = vdof1 + i * ndof2 + k;
                  int col = u * ndof1 + l;

                  // Traction of test (from elem2)
                  real_t trac_test = lam * dshapes2_phys(k * dim_ + i, q) * n_q(u)
                     + mu_val * ((i == u ? 1.0 : 0.0) * dn2(k)
                                 + dshapes2_phys(k * dim_ + u, q) * n_q(i));

                  // Traction of trial (from elem1)
                  real_t trac_trial = lam * dshapes1_phys(l * dim_ + u, q) * n_q(i)
                     + mu_val * ((u == i ? 1.0 : 0.0) * dn1(l)
                                 + dshapes1_phys(l * dim_ + i, q) * n_q(u));

                  real_t val = 0.0;
                  // Symmetry (c0): -0.5 * trac_test * φ_l / detJ2
                  val += (-0.5) * trac_test * shapes1(l, q) * invdetJ2;
                  // Consistency (c1): -ε * 0.5 * trac_trial * φ_k / detJ1
                  val += (-epsilon_) * 0.5 * trac_trial * shapes2(k, q) * invdetJ1;
                  // BR2 lifting: -σ * φ_k * L_q[y=0][(l,i,u),q]
                  val += (-sigma) * shapes2(k, q)
                     * L_q_y0((l * dim_ + i) * dim_ + u, q);

                  elmat(row, col) += wq * val;
               }
            }
         }
      }

      // Block (1,1): x=elem2, y=elem2
      for (int k = 0; k < ndof2; k++)
      {
         for (int l = 0; l < ndof2; l++)
         {
            for (int i = 0; i < dim_; i++)
            {
               for (int u = 0; u < dim_; u++)
               {
                  // byNODES ordering: component i, DOF k
                  int row = vdof1 + i * ndof2 + k;
                  int col = vdof1 + u * ndof2 + l;

                  real_t trac_test = lam * dshapes2_phys(k * dim_ + i, q) * n_q(u)
                     + mu_val * ((i == u ? 1.0 : 0.0) * dn2(k)
                                 + dshapes2_phys(k * dim_ + u, q) * n_q(i));

                  real_t trac_trial = lam * dshapes2_phys(l * dim_ + u, q) * n_q(i)
                     + mu_val * ((u == i ? 1.0 : 0.0) * dn2(l)
                                 + dshapes2_phys(l * dim_ + i, q) * n_q(u));

                  real_t val = 0.0;
                  // Symmetry (c0): +0.5 * trac_test * φ_l / detJ2
                  val += 0.5 * trac_test * shapes2(l, q) * invdetJ2;
                  // Consistency (c1): -ε * 0.5 * trac_trial * φ_k / detJ2
                  val += (-epsilon_) * 0.5 * trac_trial * shapes2(k, q) * invdetJ2;
                  // BR2 lifting: +σ * φ_k * L_q[y=1][(l,i,u),q]
                  val += sigma * shapes2(k, q)
                     * L_q_y1((l * dim_ + i) * dim_ + u, q);

                  elmat(row, col) += wq * val;
               }
            }
         }
      }
   }
}

// ============================================================================
// Implementation: Boundary Face
// ============================================================================

inline void DGElasticityBR2BoundaryIntegrator::AssembleFaceMatrix(
   const FiniteElement &el1,
   const FiniteElement &el2,
   FaceElementTransformations &Trans,
   DenseMatrix &elmat)
{
   int ndof = el1.GetDof();
   int vdof = ndof * dim_;

   elmat.SetSize(vdof);
   elmat = 0.0;

   // Detect element type for penalty (local variable for thread safety)
   Geometry::Type geom = el1.GetGeomType();
   const real_t sigma = (geom == Geometry::TETRAHEDRON)
                            ? real_t(dim_ + 1)   // 4 faces for tet
                            : real_t(2 * dim_);  // 6 faces for hex

   int elem = Trans.Elem1No;
   const DenseMatrix &Minv = elem_mass_inv_[elem];

   int order = 2 * el1.GetOrder() + 1;
   const IntegrationRule &ir = IntRules.Get(Trans.FaceGeom, order);
   int nqp = ir.GetNPoints();

   // Precompute shapes, gradients, normals
   DenseMatrix shapes(ndof, nqp);
   DenseMatrix dshapes_phys(ndof * dim_, nqp);
   Vector nor_all(dim_ * nqp);
   Vector w_all(nqp);
   Vector lam_all(nqp), mu_all(nqp);
   Vector invdetJ_all(nqp);

   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      Trans.SetAllIntPoints(&ip);
      const IntegrationPoint &eip = Trans.GetElement1IntPoint();

      Vector shape_q(shapes.GetColumn(q), ndof);
      el1.CalcShape(eip, shape_q);

      DenseMatrix dshape_ref(ndof, dim_);
      el1.CalcDShape(eip, dshape_ref);

      DenseMatrix adjJ(dim_);
      CalcAdjugate(Trans.Elem1->Jacobian(), adjJ);
      DenseMatrix dshape_phys(ndof, dim_);
      Mult(dshape_ref, adjJ, dshape_phys);

      for (int k = 0; k < ndof; k++)
         for (int d = 0; d < dim_; d++)
            dshapes_phys(k * dim_ + d, q) = dshape_phys(k, d);

      Vector nor_q(&nor_all[q * dim_], dim_);
      CalcOrtho(Trans.Jacobian(), nor_q);

      w_all[q] = ip.weight;
      lam_all[q] = lambda_.Eval(*Trans.Face, ip);
      mu_all[q] = mu_.Eval(*Trans.Face, ip);
      invdetJ_all[q] = 1.0 / Trans.Elem1->Weight();
   }

   // Compute boundary BR2 lifting
   // IntFace[l*dim+s, m] = sum_q E[m,q] * E[l,q] * n[s,q] * w[q]
   DenseMatrix IntFace(ndof * dim_, ndof);
   IntFace = 0.0;

   for (int q = 0; q < nqp; q++)
   {
      real_t wq = w_all[q];
      for (int l = 0; l < ndof; l++)
      {
         real_t shape_l = shapes(l, q);
         for (int s = 0; s < dim_; s++)
         {
            real_t n_s = nor_all[q * dim_ + s];
            real_t factor = shape_l * n_s * wq;
            for (int m = 0; m < ndof; m++)
            {
               IntFace(l * dim_ + s, m) += shapes(m, q) * factor;
            }
         }
      }
   }

   // Lift[l,s,m] = Minv[m,o] * IntFace[l*dim+s, o]
   DenseMatrix Lift(ndof * dim_, ndof);
   MultABt(IntFace, Minv, Lift);

   // L_q[l,i,u,q] = sum_{s,m} test_normal_{iu,s} * E[m,q] * Lift[l,s,m]
   DenseMatrix L_q(ndof * dim_ * dim_, nqp);
   L_q = 0.0;

   for (int q = 0; q < nqp; q++)
   {
      real_t lam = lam_all[q];
      real_t mu_val = mu_all[q];
      Vector n_q(dim_);
      for (int d = 0; d < dim_; d++) { n_q(d) = nor_all[q * dim_ + d]; }

      for (int l = 0; l < ndof; l++)
      {
         for (int i = 0; i < dim_; i++)
         {
            for (int u = 0; u < dim_; u++)
            {
               real_t sum = 0.0;
               for (int s = 0; s < dim_; s++)
               {
                  real_t tn = TestNormal(lam, mu_val, n_q, i, u, s);
                  real_t eval_lift = 0.0;
                  for (int m = 0; m < ndof; m++)
                  {
                     eval_lift += shapes(m, q) * Lift(l * dim_ + s, m);
                  }
                  sum += tn * eval_lift;
               }
               L_q((l * dim_ + i) * dim_ + u, q) = sum;
            }
         }
      }
   }

   // Assemble boundary face matrix
   // For boundary: c0=-1, c1=ε, c2=σ (no 0.5 averaging)
   for (int q = 0; q < nqp; q++)
   {
      real_t wq = w_all[q];
      real_t lam = lam_all[q];
      real_t mu_val = mu_all[q];
      real_t invdetJ = invdetJ_all[q];

      Vector n_q(dim_);
      for (int d = 0; d < dim_; d++) { n_q(d) = nor_all[q * dim_ + d]; }

      Vector dn(ndof);
      for (int k = 0; k < ndof; k++)
      {
         real_t val = 0.0;
         for (int d = 0; d < dim_; d++)
         {
            val += dshapes_phys(k * dim_ + d, q) * n_q(d);
         }
         dn(k) = val;
      }

      for (int k = 0; k < ndof; k++)
      {
         for (int l = 0; l < ndof; l++)
         {
            for (int i = 0; i < dim_; i++)
            {
               for (int u = 0; u < dim_; u++)
               {
                  // byNODES ordering: component i, DOF k
                  int row = i * ndof + k;
                  int col = u * ndof + l;

                  real_t trac_test = lam * dshapes_phys(k * dim_ + i, q) * n_q(u)
                     + mu_val * ((i == u ? 1.0 : 0.0) * dn(k)
                                 + dshapes_phys(k * dim_ + u, q) * n_q(i));

                  real_t trac_trial = lam * dshapes_phys(l * dim_ + u, q) * n_q(i)
                     + mu_val * ((u == i ? 1.0 : 0.0) * dn(l)
                                 + dshapes_phys(l * dim_ + i, q) * n_q(u));

                  real_t val = 0.0;
                  val += (-1.0) * trac_test * shapes(l, q) * invdetJ;
                  val += epsilon_ * trac_trial * shapes(k, q) * invdetJ;
                  val += sigma * shapes(k, q)
                     * L_q((l * dim_ + i) * dim_ + u, q);

                  elmat(row, col) += wq * val;
               }
            }
         }
      }
   }
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DG_ELASTICITY_BR2_INTEGRATOR_HPP
