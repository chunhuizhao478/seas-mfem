// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-443211.
// All rights reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_FACE_QUADRATURE_HPP
#define MFEM_SEAS_FACE_QUADRATURE_HPP

#include "mfem.hpp"
#include <cmath>

namespace mfem
{
namespace seas
{

/// @brief Face quadrature and L2 projection for multi-DOF fault discretization.
///
/// Provides basis functions, reference mass matrix inverse, and quadrature
/// for projecting traction from quadrature points to fault DOFs and
/// interpolating slip from DOFs to quadrature points.
///
/// Follows Tandem's AdapterBase architecture (AdapterBase.cpp lines 28-96):
///   - e_q:       basis functions at quadrature points [nbf × nq]
///   - M_ref_inv: inverse of reference mass matrix [nbf × nbf]
///   - ir:        quadrature rule matching DG face integrals
///
/// Design:
///   - face_order = 0  → nbf = 1 (constant, equivalent to face averaging)
///   - face_order ≥ 1  → nbf = (p+1)(p+2)/2 via H1_TriangleElement
///
/// Usage in SEAS:
///   - p=1 (vol_order=1): face_order=0, nbf=1 → proven stable, backward compatible
///   - p≥2 (vol_order≥2): face_order=vol_order, nbf=(p+1)(p+2)/2 → matches Tandem
///
/// Math: For flat faces, the physical mass matrix M_phys = |J_face| × M_ref.
/// Since |J_face| also appears in the projection integral, it cancels:
///   τ_k = M_ref_inv × Σ_q w_q × φ_k(q) × T_q
/// At nbf=1: this reduces to face_average(T_q). See v44 debug doc.
class FaceQuadrature
{
public:
   /// Construct face quadrature for given polynomial orders.
   /// @param face_order Basis polynomial order (0 = constant/face-averaged, ≥ 1 = H1)
   /// @param vol_order Volume polynomial order (controls quadrature rule to match
   ///                  the DG face integral rule used in AssembleSlipContribution
   ///                  and ComputeTraction)
   /// @param face_geom Face geometry type (default: TRIANGLE for tet meshes)
   /// @param basis_type BasisType for face DOF nodes (default: GaussLobatto).
   ///                   v50g: ClosedUniform has 50× better mass matrix conditioning
   ///                   at p=4 on triangles (cond=58 vs 2901 for GaussLobatto).
   FaceQuadrature(int face_order, int vol_order,
                  Geometry::Type face_geom = Geometry::TRIANGLE,
                  int basis_type = BasisType::GaussLobatto)
       : order_(face_order),
         face_fe_(std::max(face_order, 1), basis_type)
   {
      MFEM_VERIFY(face_order >= 0, "FaceQuadrature requires face_order >= 0");
      MFEM_VERIFY(vol_order >= 1, "FaceQuadrature requires vol_order >= 1");
      MFEM_VERIFY(face_geom == Geometry::TRIANGLE,
                  "FaceQuadrature currently supports TRIANGLE faces only");

      // Quadrature rule: 2*vol_order+1 matches the DG face integral rule
      // used in both ComputeTraction and AssembleSlipContribution.
      const IntegrationRule &ir_ref = IntRules.Get(face_geom,
                                                    2 * vol_order + 1);
      nq_ = ir_ref.GetNPoints();

      // Copy the quadrature rule (we need to own it)
      ir_.SetSize(nq_);
      for (int q = 0; q < nq_; q++)
      {
         ir_.IntPoint(q) = ir_ref.IntPoint(q);
      }

      if (face_order == 0)
      {
         // ================================================================
         // Order 0: single constant basis function φ₀ = 1
         // This gives nbf=1, equivalent to face-averaged (current v42 behavior).
         // ================================================================
         nbf_ = 1;

         // Basis = constant 1 at all quad points
         e_q_.SetSize(1, nq_);
         for (int q = 0; q < nq_; q++)
         {
            e_q_(0, q) = 1.0;
         }

         // Reference mass matrix: M_ref = Σ_q w_q × 1 × 1 = Σ w_q = area_ref
         // For the reference triangle with vertices (0,0),(1,0),(0,1): area = 1/2
         // M_ref_inv = 1/area_ref
         M_ref_inv_.SetSize(1);
         real_t M_ref_00 = 0.0;
         for (int q = 0; q < nq_; q++)
         {
            M_ref_00 += ir_.IntPoint(q).weight;
         }
         M_ref_inv_(0, 0) = 1.0 / M_ref_00;

         // Nodal rule: single point at face centroid
         nodal_ir_order0_.SetSize(1);
         nodal_ir_order0_.IntPoint(0).x = 1.0 / 3.0;
         nodal_ir_order0_.IntPoint(0).y = 1.0 / 3.0;
         nodal_ir_order0_.IntPoint(0).weight = 0.5;
      }
      else
      {
         // ================================================================
         // Order ≥ 1: use H1 triangle element (GaussLobatto nodes)
         // nbf = (p+1)(p+2)/2 for triangles
         // Matches Tandem's NodalRefElement<2>(PolynomialDegree)
         // ================================================================
         nbf_ = face_fe_.GetDof();

         // Evaluate basis functions at quadrature points: e_q_[nbf × nq]
         e_q_.SetSize(nbf_, nq_);
         Vector shape(nbf_);
         for (int q = 0; q < nq_; q++)
         {
            face_fe_.CalcShape(ir_.IntPoint(q), shape);
            for (int k = 0; k < nbf_; k++)
            {
               e_q_(k, q) = shape(k);
            }
         }

         // Compute reference mass matrix: M_ref[i,j] = Σ_q w_q × φ_i(q) × φ_j(q)
         DenseMatrix M_ref(nbf_);
         M_ref = 0.0;
         for (int i = 0; i < nbf_; i++)
         {
            for (int j = 0; j < nbf_; j++)
            {
               real_t val = 0.0;
               for (int q = 0; q < nq_; q++)
               {
                  val += ir_.IntPoint(q).weight * e_q_(i, q) * e_q_(j, q);
               }
               M_ref(i, j) = val;
            }
         }

         // Invert: M_ref_inv = M_ref^{-1}
         DenseMatrixInverse M_ref_inv_solver(M_ref);
         M_ref_inv_.SetSize(nbf_);
         M_ref_inv_solver.GetInverseMatrix(M_ref_inv_);
      }
   }

   /// Number of basis functions per face.
   /// Order 0: 1 (constant).
   /// Order p ≥ 1: (p+1)(p+2)/2 for triangles.
   int NumBasisFunctions() const { return nbf_; }

   /// Number of quadrature points (determined by 2*vol_order+1 rule).
   int NumQuadPoints() const { return nq_; }

   /// Basis functions evaluated at quadrature points [nbf × nq].
   /// e_q(k, q) = φ_k(x_q) where x_q is the q-th quadrature point.
   const DenseMatrix &BasisAtQuadPoints() const { return e_q_; }

   /// Reference mass matrix inverse [nbf × nbf].
   /// M_ref_inv = (Σ_q w_q φ_i(q) φ_j(q))^{-1}
   const DenseMatrix &RefMassInverse() const { return M_ref_inv_; }

   /// Face quadrature rule (2*vol_order+1 points on reference triangle).
   const IntegrationRule &GetQuadRule() const { return ir_; }

   /// Nodal integration rule (node positions on reference face).
   /// For order 0: single centroid point (1/3, 1/3).
   /// For order ≥ 1: GaussLobatto nodes from H1_TriangleElement.
   const IntegrationRule &GetNodalRule() const
   {
      return (order_ == 0) ? nodal_ir_order0_ : face_fe_.GetNodes();
   }

   /// @brief Interpolate nodal values to quadrature points.
   ///
   /// For each component c:
   ///   quad_vals[c*nq + q] = Σ_k e_q(k, q) × nodal_vals[c*nbf + k]
   ///
   /// @param ncomp Number of components (e.g. 2 for dip+strike)
   /// @param nodal_vals Input nodal values [ncomp × nbf]
   /// @param quad_vals  Output quad-point values [ncomp × nq]
   void InterpolateToQuadPoints(int ncomp, const Vector &nodal_vals,
                                Vector &quad_vals) const
   {
      MFEM_ASSERT(nodal_vals.Size() == ncomp * nbf_,
                  "nodal_vals size mismatch: got " << nodal_vals.Size()
                  << ", expected " << ncomp * nbf_);
      quad_vals.SetSize(ncomp * nq_);
      quad_vals = 0.0;

      for (int c = 0; c < ncomp; c++)
      {
         for (int q = 0; q < nq_; q++)
         {
            real_t val = 0.0;
            for (int k = 0; k < nbf_; k++)
            {
               val += e_q_(k, q) * nodal_vals(c * nbf_ + k);
            }
            quad_vals(c * nq_ + q) = val;
         }
      }
   }

   /// @brief L2 (Galerkin) project quad-point values to nodal DOFs.
   ///
   /// For each component c:
   ///   nodal_vals[c*nbf + k] = Σ_l M_ref_inv(k,l) × Σ_q w_q × φ_l(q) × quad_vals[c*nq + q]
   ///
   /// This is the standard L2 projection onto the face polynomial space.
   /// For flat faces, the physical mass matrix factors cancel (see class doc).
   ///
   /// At order 0 (nbf=1): reduces to weighted average = face average.
   /// At order ≥ 1: full L2 projection preserving spatial variation.
   ///
   /// @param ncomp Number of components
   /// @param quad_vals  Input quad-point values [ncomp × nq]
   /// @param nodal_vals Output nodal values [ncomp × nbf]
   void GalerkinProject(int ncomp, const Vector &quad_vals,
                        Vector &nodal_vals) const
   {
      MFEM_ASSERT(quad_vals.Size() == ncomp * nq_,
                  "quad_vals size mismatch: got " << quad_vals.Size()
                  << ", expected " << ncomp * nq_);
      nodal_vals.SetSize(ncomp * nbf_);
      nodal_vals = 0.0;

      // Step 1: Compute RHS = Σ_q w_q × φ_l(q) × f_q
      Vector rhs(ncomp * nbf_);
      rhs = 0.0;
      for (int c = 0; c < ncomp; c++)
      {
         for (int l = 0; l < nbf_; l++)
         {
            real_t val = 0.0;
            for (int q = 0; q < nq_; q++)
            {
               val += ir_.IntPoint(q).weight * e_q_(l, q)
                      * quad_vals(c * nq_ + q);
            }
            rhs(c * nbf_ + l) = val;
         }
      }

      // Step 2: Apply M_ref_inv: result = M_ref_inv × rhs
      for (int c = 0; c < ncomp; c++)
      {
         for (int k = 0; k < nbf_; k++)
         {
            real_t val = 0.0;
            for (int l = 0; l < nbf_; l++)
            {
               val += M_ref_inv_(k, l) * rhs(c * nbf_ + l);
            }
            nodal_vals(c * nbf_ + k) = val;
         }
      }
   }

private:
   int order_;  ///< Face polynomial order (0 = constant, ≥1 = H1)
   int nbf_;    ///< Number of basis functions: 1 or (p+1)(p+2)/2
   int nq_;     ///< Number of quadrature points

   H1_TriangleElement face_fe_;       ///< Face finite element (used for order >= 1)
   DenseMatrix e_q_;                  ///< Basis at quad points [nbf × nq]
   DenseMatrix M_ref_inv_;            ///< Reference mass matrix inverse [nbf × nbf]
   IntegrationRule ir_;               ///< Face quadrature rule
   IntegrationRule nodal_ir_order0_;  ///< Nodal rule for order 0 (centroid)
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FACE_QUADRATURE_HPP
