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

#ifndef MFEM_SEAS_DG_BR2_INTEGRATOR_HPP
#define MFEM_SEAS_DG_BR2_INTEGRATOR_HPP

#include "mfem.hpp"
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief BR2 (Bassi-Rebay 2) Interior Face Integrator for DG methods
///
/// This integrator implements the BR2 stabilization method for the Laplace
/// equation using lifting operators, following Tandem's implementation.
///
/// The BR2 bilinear form for interior faces consists of three terms:
///   a(u,v) = -∫_e {{K∇u·n}} [[v]] ds     (consistency)
///          - ε∫_e {{K∇v·n}} [[u]] ds     (symmetry, ε=-1 for SIPG)
///          + σ ∫_e E_q[x] · L_q[y] ds    (BR2 lifting)
///
/// where L_q[y] is the lifted flux computed following Tandem's formula:
///   Lift[i] = 0.5 * M_i^{-1} * ∫_e φ_i ⊗ n * φ_y ds
///   L_q[y] = 0.5 * n · (K * E_0 * Lift[0] + K * E_1 * Lift[1])
///
/// The BR2 penalty parameter is σ = D + 1 (= 3 for 2D).
///
/// Reference: Tandem's Poisson.cpp and poisson.py kernels
///
class BR2InteriorFaceIntegrator : public BilinearFormIntegrator
{
public:
   /// @brief Construct BR2 interior face integrator
   ///
   /// @param K Diffusion coefficient (usually constant = 1)
   /// @param epsilon SIPG sign (-1 for symmetric, +1 for non-symmetric)
   /// @param elem_mass_inv Precomputed element mass matrix inverses
   /// @param dim Spatial dimension (2 for antiplane)
   BR2InteriorFaceIntegrator(Coefficient &K, real_t epsilon,
                             const std::vector<DenseMatrix> &elem_mass_inv,
                             int dim = 2)
      : K_(K), epsilon_(epsilon), elem_mass_inv_(elem_mass_inv), dim_(dim)
   {
      // BR2 penalty: σ = D + 1
      sigma_ = dim_ + 1;
   }

   /// @brief Assemble face matrix for interior face
   ///
   /// Follows Tandem's assembleSurface:
   /// a[x][y] = c0[y]*∫K∇φ_x·n φ_y + c1[x]*∫K∇φ_y·n φ_x + c2[|x-y|]*∫φ_x L_q[y]
   void AssembleFaceMatrix(const FiniteElement &el1,
                           const FiniteElement &el2,
                           FaceElementTransformations &Trans,
                           DenseMatrix &elmat) override;

protected:
   Coefficient &K_;
   real_t epsilon_;   // SIPG sign (-1 for symmetric)
   real_t sigma_;     // BR2 penalty: D + 1
   const std::vector<DenseMatrix> &elem_mass_inv_;
   int dim_;
};

/// @brief BR2 Boundary Face Integrator for Dirichlet BC
///
/// Similar to BR2InteriorFaceIntegrator but for boundary faces.
/// Uses the boundary lifting formula from Tandem's lift_boundary.
///
class BR2BoundaryFaceIntegrator : public BilinearFormIntegrator
{
public:
   BR2BoundaryFaceIntegrator(Coefficient &K, real_t epsilon,
                             const std::vector<DenseMatrix> &elem_mass_inv,
                             int dim = 2)
      : K_(K), epsilon_(epsilon), elem_mass_inv_(elem_mass_inv), dim_(dim)
   {
      sigma_ = dim_ + 1;
   }

   void AssembleFaceMatrix(const FiniteElement &el1,
                           const FiniteElement &el2,
                           FaceElementTransformations &Trans,
                           DenseMatrix &elmat) override;

protected:
   Coefficient &K_;
   real_t epsilon_;
   real_t sigma_;
   const std::vector<DenseMatrix> &elem_mass_inv_;
   int dim_;
};

// ============================================================================
// Implementation
// ============================================================================

inline void BR2InteriorFaceIntegrator::AssembleFaceMatrix(
   const FiniteElement &el1,
   const FiniteElement &el2,
   FaceElementTransformations &Trans,
   DenseMatrix &elmat)
{
   int ndof1 = el1.GetDof();
   int ndof2 = el2.GetDof();
   int ndofs = ndof1 + ndof2;

   elmat.SetSize(ndofs);
   elmat = 0.0;

   // Get element indices
   int elem1 = Trans.Elem1No;
   int elem2 = Trans.Elem2No;

   // Get mass matrix inverses
   const DenseMatrix &Minv1 = elem_mass_inv_[elem1];
   const DenseMatrix &Minv2 = elem_mass_inv_[elem2];

   // Integration order
   int order = 2 * std::max(el1.GetOrder(), el2.GetOrder()) + 1;
   const IntegrationRule &ir = IntRules.Get(Trans.FaceGeom, order);
   int nqp = ir.GetNPoints();

   // Precompute shape functions at all quadrature points
   DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
   DenseMatrix dshapes1(ndof1 * dim_, nqp), dshapes2(ndof2 * dim_, nqp);
   Vector nor_all(dim_ * nqp);
   Vector w_all(nqp);

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

      // Gradients in reference coordinates
      DenseMatrix dshape1_ref(ndof1, dim_), dshape2_ref(ndof2, dim_);
      el1.CalcDShape(eip1, dshape1_ref);
      el2.CalcDShape(eip2, dshape2_ref);

      // Store flattened
      for (int i = 0; i < ndof1; i++)
         for (int d = 0; d < dim_; d++)
            dshapes1(i * dim_ + d, q) = dshape1_ref(i, d);
      for (int i = 0; i < ndof2; i++)
         for (int d = 0; d < dim_; d++)
            dshapes2(i * dim_ + d, q) = dshape2_ref(i, d);

      // Normal vector
      Vector nor_q(&nor_all[q * dim_], dim_);
      CalcOrtho(Trans.Jacobian(), nor_q);

      // Quadrature weight
      w_all[q] = ip.weight;
   }

   // =========================================================================
   // Compute BR2 Lift matrices following Tandem's lift_skeleton
   // Lift[elem][l,i,u] = 0.5 * Minv[elem][u,s] * Σ_q E[elem][s,q] * E[y][l,q] * n[i,q] * w[q]
   //
   // For each source element y ∈ {0,1}, we compute Lift[0] and Lift[1]
   // then L_q[y][l,q] = 0.5 * n[i,q] * (K * E[0][u,q] * Lift[0][l,i,u] + K * E[1][v,q] * Lift[1][l,i,v])
   // =========================================================================

   // L_q[y=0] and L_q[y=1] matrices: (ndof_y, nqp)
   DenseMatrix L_q0(ndof1, nqp), L_q1(ndof2, nqp);
   L_q0 = 0.0;
   L_q1 = 0.0;

   // For K=1 (constant diffusion), the formula simplifies
   // We compute the lifting for each source element separately

   // Lift contribution from element 1 basis (y=0) to L_q[0]
   // and Lift contribution from element 2 basis (y=1) to L_q[1]

   // Step 1: Compute Lift matrices
   // Lift[elem][l,i,u] = 0.5 * Minv[elem][u,s] * Σ_q E[elem][s,q] * E[y][l,q] * n[i,q] * w[q]

   // For y=0 (element 1 source):
   // Lift[0] contributions: 0.5 * Minv1[u,s] * Σ_q shapes1[s,q] * shapes1[l,q] * n[i,q] * w[q]
   // Lift[1] contributions: 0.5 * Minv2[v,t] * Σ_q shapes2[t,q] * shapes1[l,q] * n[i,q] * w[q]

   // Then L_q[0][l,q] = 0.5 * n[i,q] * (shapes1[u,q] * Lift[0][l,i,u] + shapes2[v,q] * Lift[1][l,i,v])

   // This is O(ndof^3 * nqp * dim) which can be expensive.
   // For efficiency, we reformulate:

   // Let Face1[l,i] = Σ_q shapes1[l,q] * n[i,q] * w[q]  (face integral for element 1 source)
   // Let Face2[l,i] = Σ_q shapes2[l,q] * n[i,q] * w[q]  (face integral for element 2 source)

   // For y=0 (source from elem1):
   // IntFace1_y0[l,i,s] = Σ_q shapes1[s,q] * shapes1[l,q] * n[i,q] * w[q]
   // Lift1_y0[l,i,u] = 0.5 * Minv1[u,s] * IntFace1_y0[l,i,s]
   //                 = 0.5 * Σ_s Minv1[u,s] * IntFace1_y0[l,i,s]

   // Then at each q:
   // L_q[0][l,q] = 0.5 * Σ_i n[i,q] * (Σ_u shapes1[u,q] * Lift1_y0[l,i,u] + Σ_v shapes2[v,q] * Lift2_y0[l,i,v])

   // For computational efficiency, precompute the face integrals

   // Face integral tensor: IntFace[elem][y][l,i,s] = Σ_q E[elem][s,q] * E[y][l,q] * n[i,q] * w[q]

   // For y=0 (ndof1 source DOFs):
   DenseMatrix IntFace11_y0(ndof1 * dim_, ndof1);  // [l*dim+i, s] for elem1
   DenseMatrix IntFace21_y0(ndof1 * dim_, ndof2);  // [l*dim+i, t] for elem2
   IntFace11_y0 = 0.0;
   IntFace21_y0 = 0.0;

   for (int q = 0; q < nqp; q++)
   {
      real_t wq = w_all[q];
      for (int l = 0; l < ndof1; l++)
      {
         real_t shape1_l = shapes1(l, q);
         for (int i = 0; i < dim_; i++)
         {
            real_t n_i = nor_all[q * dim_ + i];
            real_t factor = shape1_l * n_i * wq;
            for (int s = 0; s < ndof1; s++)
            {
               IntFace11_y0(l * dim_ + i, s) += shapes1(s, q) * factor;
            }
            for (int t = 0; t < ndof2; t++)
            {
               IntFace21_y0(l * dim_ + i, t) += shapes2(t, q) * factor;
            }
         }
      }
   }

   // For y=1 (ndof2 source DOFs):
   DenseMatrix IntFace12_y1(ndof2 * dim_, ndof1);  // [l*dim+i, s] for elem1
   DenseMatrix IntFace22_y1(ndof2 * dim_, ndof2);  // [l*dim+i, t] for elem2
   IntFace12_y1 = 0.0;
   IntFace22_y1 = 0.0;

   for (int q = 0; q < nqp; q++)
   {
      real_t wq = w_all[q];
      for (int l = 0; l < ndof2; l++)
      {
         real_t shape2_l = shapes2(l, q);
         for (int i = 0; i < dim_; i++)
         {
            real_t n_i = nor_all[q * dim_ + i];
            real_t factor = shape2_l * n_i * wq;
            for (int s = 0; s < ndof1; s++)
            {
               IntFace12_y1(l * dim_ + i, s) += shapes1(s, q) * factor;
            }
            for (int t = 0; t < ndof2; t++)
            {
               IntFace22_y1(l * dim_ + i, t) += shapes2(t, q) * factor;
            }
         }
      }
   }

   // Now compute Lift matrices: Lift[elem][l,i,u] = 0.5 * Minv[elem][u,s] * IntFace[elem][l*dim+i,s]
   // Lift1_y0[l,i,u] = 0.5 * Σ_s Minv1[u,s] * IntFace11_y0[l*dim+i,s]
   // Lift2_y0[l,i,v] = 0.5 * Σ_t Minv2[v,t] * IntFace21_y0[l*dim+i,t]

   // For y=0:
   DenseMatrix Lift1_y0(ndof1 * dim_, ndof1);  // [l*dim+i, u]
   DenseMatrix Lift2_y0(ndof1 * dim_, ndof2);  // [l*dim+i, v]
   // Lift1_y0 = 0.5 * IntFace11_y0 * Minv1^T
   MultABt(IntFace11_y0, Minv1, Lift1_y0);
   Lift1_y0 *= 0.5;
   MultABt(IntFace21_y0, Minv2, Lift2_y0);
   Lift2_y0 *= 0.5;

   // For y=1:
   DenseMatrix Lift1_y1(ndof2 * dim_, ndof1);  // [l*dim+i, u]
   DenseMatrix Lift2_y1(ndof2 * dim_, ndof2);  // [l*dim+i, v]
   MultABt(IntFace12_y1, Minv1, Lift1_y1);
   Lift1_y1 *= 0.5;
   MultABt(IntFace22_y1, Minv2, Lift2_y1);
   Lift2_y1 *= 0.5;

   // Compute L_q at each quadrature point
   // L_q[y][l,q] = 0.5 * Σ_i n[i,q] * (Σ_u E[0][u,q] * Lift[0]_y[l,i,u] + Σ_v E[1][v,q] * Lift[1]_y[l,i,v])

   for (int q = 0; q < nqp; q++)
   {
      // For y=0 (element 1 source)
      for (int l = 0; l < ndof1; l++)
      {
         real_t sum = 0.0;
         for (int i = 0; i < dim_; i++)
         {
            real_t n_i = nor_all[q * dim_ + i];
            real_t contrib1 = 0.0, contrib2 = 0.0;
            for (int u = 0; u < ndof1; u++)
            {
               contrib1 += shapes1(u, q) * Lift1_y0(l * dim_ + i, u);
            }
            for (int v = 0; v < ndof2; v++)
            {
               contrib2 += shapes2(v, q) * Lift2_y0(l * dim_ + i, v);
            }
            sum += n_i * (contrib1 + contrib2);
         }
         L_q0(l, q) = 0.5 * sum;
      }

      // For y=1 (element 2 source)
      for (int l = 0; l < ndof2; l++)
      {
         real_t sum = 0.0;
         for (int i = 0; i < dim_; i++)
         {
            real_t n_i = nor_all[q * dim_ + i];
            real_t contrib1 = 0.0, contrib2 = 0.0;
            for (int u = 0; u < ndof1; u++)
            {
               contrib1 += shapes1(u, q) * Lift1_y1(l * dim_ + i, u);
            }
            for (int v = 0; v < ndof2; v++)
            {
               contrib2 += shapes2(v, q) * Lift2_y1(l * dim_ + i, v);
            }
            sum += n_i * (contrib1 + contrib2);
         }
         L_q1(l, q) = 0.5 * sum;
      }
   }

   // =========================================================================
   // Assemble the face matrix following Tandem's assembleSurface
   // a[x][y] = c0[y]*∫K∇φ_x·n φ_y ds + c1[x]*∫K∇φ_y·n φ_x ds + c2[|x-y|]*∫φ_x L_q[y] ds
   //
   // c0 = -0.5, c1 = epsilon * 0.5
   // c2[0] = sigma (same element), c2[1] = -sigma (cross element)
   // =========================================================================

   real_t c0 = -0.5;
   real_t c1 = epsilon_ * 0.5;
   real_t c20 = sigma_;    // same element
   real_t c21 = -sigma_;   // cross element

   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      Trans.SetAllIntPoints(&ip);

      const IntegrationPoint &eip1 = Trans.GetElement1IntPoint();
      const IntegrationPoint &eip2 = Trans.GetElement2IntPoint();

      // Shape functions at this point
      Vector shape1(ndof1), shape2(ndof2);
      el1.CalcShape(eip1, shape1);
      el2.CalcShape(eip2, shape2);

      // Gradient of shape functions (physical)
      DenseMatrix dshape1_ref(ndof1, dim_), dshape2_ref(ndof2, dim_);
      el1.CalcDShape(eip1, dshape1_ref);
      el2.CalcDShape(eip2, dshape2_ref);

      // Transform gradients to physical space
      DenseMatrix adjJ1(dim_), adjJ2(dim_);
      CalcAdjugate(Trans.Elem1->Jacobian(), adjJ1);
      CalcAdjugate(Trans.Elem2->Jacobian(), adjJ2);

      DenseMatrix dshape1_phys(ndof1, dim_), dshape2_phys(ndof2, dim_);
      Mult(dshape1_ref, adjJ1, dshape1_phys);
      Mult(dshape2_ref, adjJ2, dshape2_phys);

      // Normal
      Vector nor(dim_);
      CalcOrtho(Trans.Jacobian(), nor);

      // ∇φ · n
      Vector dn1(ndof1), dn2(ndof2);
      dshape1_phys.Mult(nor, dn1);
      dshape2_phys.Mult(nor, dn2);

      // Jacobian determinants for scaling
      real_t detJ1 = Trans.Elem1->Weight();
      real_t detJ2 = Trans.Elem2->Weight();

      real_t wq = ip.weight;

      // Assembly of 4 blocks: a[0][0], a[0][1], a[1][0], a[1][1]

      // Block (0,0): x=0, y=0
      // c0[0] * ∫ K∇φ_0·n φ_0 + c1[0] * ∫ K∇φ_0·n φ_0 + c2[0] * ∫ φ_0 L_q[0]
      // c0[0] = c0, c1[0] = c1, c2[0] = c20
      for (int k = 0; k < ndof1; k++)
      {
         for (int l = 0; l < ndof1; l++)
         {
            real_t val = 0.0;
            // Consistency: c0 * K∇φ_x·n * φ_y / detJ
            val += c0 * dn1(k) * shape1(l) / detJ1;
            // Symmetry: c1 * K∇φ_y·n * φ_x / detJ
            val += c1 * dn1(l) * shape1(k) / detJ1;
            // BR2 lifting: c20 * φ_x * L_q[y]
            val += c20 * shape1(k) * L_q0(l, q);
            elmat(k, l) += wq * val;
         }
      }

      // Block (0,1): x=0, y=1
      // c0[1] = -c0 (note: in Tandem c0[y] changes sign for y=1)
      // c1[0] = c1, c2[1] = c21
      for (int k = 0; k < ndof1; k++)
      {
         for (int l = 0; l < ndof2; l++)
         {
            real_t val = 0.0;
            // Consistency: -c0 * K∇φ_x·n * φ_y / detJ
            // Note: for y from element 2, we use -c0 (opposite sign)
            val += (-c0) * dn1(k) * shape2(l) / detJ1;
            // Symmetry: c1 * K∇φ_y·n * φ_x / detJ
            // ∇φ_y is from element 2
            val += c1 * dn2(l) * shape1(k) / detJ2;
            // BR2 lifting: c21 * φ_x * L_q[y]
            val += c21 * shape1(k) * L_q1(l, q);
            elmat(k, ndof1 + l) += wq * val;
         }
      }

      // Block (1,0): x=1, y=0
      // c0[0] = c0, c1[1] = -c1, c2[1] = c21
      for (int k = 0; k < ndof2; k++)
      {
         for (int l = 0; l < ndof1; l++)
         {
            real_t val = 0.0;
            // Consistency: c0 * K∇φ_x·n * φ_y / detJ
            val += c0 * dn2(k) * shape1(l) / detJ2;
            // Symmetry: -c1 * K∇φ_y·n * φ_x / detJ
            val += (-c1) * dn1(l) * shape2(k) / detJ1;
            // BR2 lifting: c21 * φ_x * L_q[y]
            val += c21 * shape2(k) * L_q0(l, q);
            elmat(ndof1 + k, l) += wq * val;
         }
      }

      // Block (1,1): x=1, y=1
      // c0[1] = -c0, c1[1] = -c1, c2[0] = c20
      for (int k = 0; k < ndof2; k++)
      {
         for (int l = 0; l < ndof2; l++)
         {
            real_t val = 0.0;
            // Consistency: -c0 * K∇φ_x·n * φ_y / detJ
            val += (-c0) * dn2(k) * shape2(l) / detJ2;
            // Symmetry: -c1 * K∇φ_y·n * φ_x / detJ
            val += (-c1) * dn2(l) * shape2(k) / detJ2;
            // BR2 lifting: c20 * φ_x * L_q[y]
            val += c20 * shape2(k) * L_q1(l, q);
            elmat(ndof1 + k, ndof1 + l) += wq * val;
         }
      }
   }
}

inline void BR2BoundaryFaceIntegrator::AssembleFaceMatrix(
   const FiniteElement &el1,
   const FiniteElement &el2,
   FaceElementTransformations &Trans,
   DenseMatrix &elmat)
{
   int ndof = el1.GetDof();
   elmat.SetSize(ndof);
   elmat = 0.0;

   int elem = Trans.Elem1No;
   const DenseMatrix &Minv = elem_mass_inv_[elem];

   int order = 2 * el1.GetOrder() + 1;
   const IntegrationRule &ir = IntRules.Get(Trans.FaceGeom, order);
   int nqp = ir.GetNPoints();

   // Compute boundary lift: Lift[l,i,u] = Minv[u,s] * Σ_q E[s,q] * E[l,q] * n[i,q] * w[q]
   // Then L_q[l,q] = n[i,q] * E[u,q] * Lift[l,i,u]

   // Precompute shape functions
   DenseMatrix shapes(ndof, nqp);
   Vector nor_all(dim_ * nqp);
   Vector w_all(nqp);

   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      Trans.SetAllIntPoints(&ip);

      const IntegrationPoint &eip = Trans.GetElement1IntPoint();
      Vector shape_q(shapes.GetColumn(q), ndof);
      el1.CalcShape(eip, shape_q);

      Vector nor_q(&nor_all[q * dim_], dim_);
      CalcOrtho(Trans.Jacobian(), nor_q);

      w_all[q] = ip.weight;
   }

   // Face integral: IntFace[l,i,s] = Σ_q E[s,q] * E[l,q] * n[i,q] * w[q]
   DenseMatrix IntFace(ndof * dim_, ndof);
   IntFace = 0.0;

   for (int q = 0; q < nqp; q++)
   {
      real_t wq = w_all[q];
      for (int l = 0; l < ndof; l++)
      {
         real_t shape_l = shapes(l, q);
         for (int i = 0; i < dim_; i++)
         {
            real_t n_i = nor_all[q * dim_ + i];
            real_t factor = shape_l * n_i * wq;
            for (int s = 0; s < ndof; s++)
            {
               IntFace(l * dim_ + i, s) += shapes(s, q) * factor;
            }
         }
      }
   }

   // Lift[l,i,u] = Minv[u,s] * IntFace[l*dim+i,s]
   DenseMatrix Lift(ndof * dim_, ndof);
   MultABt(IntFace, Minv, Lift);

   // L_q[l,q] = Σ_i n[i,q] * Σ_u E[u,q] * Lift[l,i,u]
   DenseMatrix L_q(ndof, nqp);
   for (int q = 0; q < nqp; q++)
   {
      for (int l = 0; l < ndof; l++)
      {
         real_t sum = 0.0;
         for (int i = 0; i < dim_; i++)
         {
            real_t n_i = nor_all[q * dim_ + i];
            real_t contrib = 0.0;
            for (int u = 0; u < ndof; u++)
            {
               contrib += shapes(u, q) * Lift(l * dim_ + i, u);
            }
            sum += n_i * contrib;
         }
         L_q(l, q) = sum;
      }
   }

   // Assembly: a[k,l] = c0 * ∫ K∇φ·n φ + c1 * ∫ K∇φ·n φ + c2 * ∫ φ L_q
   real_t c0 = -1.0;
   real_t c1 = epsilon_;
   real_t c2 = sigma_;

   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      Trans.SetAllIntPoints(&ip);

      const IntegrationPoint &eip = Trans.GetElement1IntPoint();

      Vector shape(ndof);
      el1.CalcShape(eip, shape);

      DenseMatrix dshape_ref(ndof, dim_);
      el1.CalcDShape(eip, dshape_ref);

      DenseMatrix adjJ(dim_);
      CalcAdjugate(Trans.Elem1->Jacobian(), adjJ);

      DenseMatrix dshape_phys(ndof, dim_);
      Mult(dshape_ref, adjJ, dshape_phys);

      Vector nor(dim_);
      CalcOrtho(Trans.Jacobian(), nor);

      Vector dn(ndof);
      dshape_phys.Mult(nor, dn);

      real_t detJ = Trans.Elem1->Weight();
      real_t wq = ip.weight;

      for (int k = 0; k < ndof; k++)
      {
         for (int l = 0; l < ndof; l++)
         {
            real_t val = 0.0;
            // Consistency: c0 * K∇φ·n * φ / detJ
            val += c0 * dn(k) * shape(l) / detJ;
            // Symmetry: c1 * K∇φ·n * φ / detJ
            val += c1 * dn(l) * shape(k) / detJ;
            // BR2 lifting: c2 * φ * L_q
            val += c2 * shape(k) * L_q(l, q);
            elmat(k, l) += wq * val;
         }
      }
   }
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DG_BR2_INTEGRATOR_HPP
