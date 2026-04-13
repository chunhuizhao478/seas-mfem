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

#include "wave_operator.hpp"
#include <cmath>
#include <limits>

namespace mfem
{
namespace seas
{

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
WaveOperator::WaveOperator(Mesh &mesh, int order,
                           real_t lambda, real_t mu, real_t rho,
                           const BoundaryConfig &bc)
   : TimeDependentOperator(0),  // size set below
     mesh_(mesh),
     order_(order),
     flux_(lambda, mu, rho),
     bc_(bc)
{
   int dim = mesh_.Dimension();
   MFEM_VERIFY(dim == 3, "WaveOperator requires 3D mesh, got dim=" << dim);

   // Create L2 DG finite element space (scalar)
   fec_ = std::make_unique<L2_FECollection>(order, dim, BasisType::GaussLobatto);
   fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());

   ne_ = mesh_.GetNE();
   ndof_per_el_ = fes_->GetFE(0)->GetDof();
   ndof_total_ = ne_ * ndof_per_el_;

   // Set TimeDependentOperator size: 9 components * ndof_total
   height = width = NUM_STATE * ndof_total_;

   // Precompute Jacobian matrices (constant for homogeneous material)
   flux_.BuildJacobian(0, Ax_);
   flux_.BuildJacobian(1, Ay_);
   flux_.BuildJacobian(2, Az_);

   // Assemble per-element inverse mass matrices
   AssembleElementMassInverse();

   // Compute minimum characteristic element size for CFL.
   // h = V^{1/dim} is the correct CFL length for both hex and tet elements.
   h_min_ = std::numeric_limits<real_t>::max();
   for (int e = 0; e < ne_; e++)
   {
      real_t h = std::pow(mesh_.GetElementVolume(e), 1.0 / mesh_.Dimension());
      h_min_ = std::min(h_min_, h);
   }

   // Build face → boundary attribute map (R-001 fix).
   // face_bdr_attr_[f] = 0 for interior faces, bdr attribute for boundary faces.
   face_bdr_attr_.assign(mesh_.GetNumFaces(), 0);
   for (int b = 0; b < mesh_.GetNBE(); b++)
   {
      int face_idx = mesh_.GetBdrFace(b);
      face_bdr_attr_[face_idx] = mesh_.GetBdrAttribute(b);
   }

   // Fault detection deferred to Phase 3 (FaultFaceFlux integration).
   num_fault_dofs_ = 0;
}

WaveOperator::~WaveOperator() = default;

// ---------------------------------------------------------------------------
// Mult: dQ/dt = M^{-1} * (-Face + Vol)   [Eq. (3)]
// ---------------------------------------------------------------------------
void WaveOperator::Mult(const Vector &Q, Vector &dQdt) const
{
   MFEM_VERIFY(Q.Size() == height,
               "Q size mismatch: " << Q.Size() << " vs " << height);

   dQdt.SetSize(height);
   dQdt = 0.0;

   // Accumulate volume integral (positive contribution)
   ComputeVolumeRHS(Q, dQdt);

   // Accumulate face flux (negative contribution, sign handled inside)
   ComputeFaceFluxRHS(Q, dQdt);

   // Apply per-element inverse mass matrix
   ApplyMassInverse(dQdt);
}

// ---------------------------------------------------------------------------
// Volume integral: Eq. (2c)
// ---------------------------------------------------------------------------
void WaveOperator::ComputeVolumeRHS(const Vector &Q, Vector &rhs) const
{
   const real_t *Q_data = Q.GetData();

   for (int e = 0; e < ne_; e++)
   {
      const FiniteElement *fe = fes_->GetFE(e);
      ElementTransformation *Tr = fes_->GetElementTransformation(e);
      int ndof = fe->GetDof();

      // Integration rule (order 2*p for exact integration of polynomial flux)
      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2*order_);
      int nqp = ir.GetNPoints();

      // Element DOF offset
      int dof_offset = e * ndof_per_el_;

      // Temporary arrays for quadrature point evaluation
      Vector shape(ndof);
      DenseMatrix dshape(ndof, 3);  // physical gradients

      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetIntPoint(&ip);
         real_t w = ip.weight * Tr->Weight();

         // Evaluate basis functions and their physical gradients
         fe->CalcShape(ip, shape);
         fe->CalcPhysDShape(*Tr, dshape);

         // Interpolate Q at this quadrature point: Q_qp[c] = sum_i shape(i) * Q_e[c*ndof+i]
         real_t Q_qp[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_qp[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               Q_qp[c] += shape(i) * Q_data[c * ndof_total_ + dof_offset + i];
            }
         }

         // Compute fluxes: F[j][c] = (A_j * Q_qp)[c] for j = 0,1,2
         real_t F[3][NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            F[0][c] = 0.0; F[1][c] = 0.0; F[2][c] = 0.0;
            for (int k = 0; k < NUM_STATE; k++)
            {
               F[0][c] += Ax_(c, k) * Q_qp[k];
               F[1][c] += Ay_(c, k) * Q_qp[k];
               F[2][c] += Az_(c, k) * Q_qp[k];
            }
         }

         // Accumulate: rhs[c*ndof_total+dof_offset+i] += w * dshape(i,j) * F[j][c]
         for (int c = 0; c < NUM_STATE; c++)
         {
            for (int i = 0; i < ndof; i++)
            {
               real_t val = 0.0;
               for (int j = 0; j < 3; j++)
               {
                  val += dshape(i, j) * F[j][c];
               }
               rhs[c * ndof_total_ + dof_offset + i] += w * val;
            }
         }
      }
   }
}

// ---------------------------------------------------------------------------
// Face flux: Eq. (2d)
// ---------------------------------------------------------------------------
void WaveOperator::ComputeFaceFluxRHS(const Vector &Q, Vector &rhs) const
{
   const real_t *Q_data = Q.GetData();

   // Process all faces
   for (int f = 0; f < mesh_.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh_.GetFaceElementTransformations(f);
      if (!ftr) { continue; }

      int e1 = ftr->Elem1No;
      int e2 = ftr->Elem2No;

      // Classify face type using precomputed boundary attribute map (R-001/R-002 fix).
      // face_bdr_attr_[f] > 0 for boundary faces, 0 for interior/shared faces.
      int bdr_attr = face_bdr_attr_[f];
      bool is_boundary = (e2 < 0) && (bdr_attr > 0);

      // In parallel (Phase 4b), shared faces have e2 < 0 but bdr_attr == 0.
      // Skip them here; they need ParMesh::GetSharedFaceTransformations.
      if (e2 < 0 && bdr_attr == 0) { continue; }

      const FiniteElement *fe1 = fes_->GetFE(e1);
      int ndof = fe1->GetDof();
      int dof_offset1 = e1 * ndof_per_el_;

      // Face integration rule: order 2p (R-004 fix)
      const IntegrationRule &ir = IntRules.Get(
         ftr->GetGeometryType(), 2*order_);
      int nqp = ir.GetNPoints();

      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         // Outward normal from Elem1 (unnormalized, |n| = face Jacobian determinant)
         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0) { nor_vec /= nor_len; }

         // Face weight: ip.weight * |CalcOrtho| = ip.weight * face_area_at_qp
         real_t w = ip.weight * nor_len;
         real_t nor[3] = {nor_vec(0), nor_vec(1), nor_vec(2)};

         // Evaluate Q on Elem1 side using reference coordinates
         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);

         real_t Q_self[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_self[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               Q_self[c] += shape1(i) * Q_data[c * ndof_total_ + dof_offset1 + i];
            }
         }

         real_t F_h[NUM_STATE];

         if (is_boundary)
         {
            // Boundary face: apply appropriate BC using precomputed attribute
            FaceBC bc_type = ClassifyBoundaryFace(bdr_attr);

            switch (bc_type)
            {
               case FaceBC::Absorbing:
                  flux_.Absorbing(nor, Q_self, F_h);
                  break;
               case FaceBC::FreeSurface:
                  flux_.FreeSurface(nor, Q_self, F_h);
                  break;
               case FaceBC::Fault:
                  // Fault faces handled separately by FaultFaceFlux (Phase 3).
                  // For now, treat as interior (welded) or absorbing.
                  flux_.Absorbing(nor, Q_self, F_h);
                  break;
               default:
                  // Default to absorbing for unclassified boundaries
                  flux_.Absorbing(nor, Q_self, F_h);
                  break;
            }

            // Accumulate into Elem1 only (subtract: face term has negative sign)
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  rhs[c * ndof_total_ + dof_offset1 + i] -= w * shape1(i) * F_h[c];
               }
            }
         }
         else
         {
            // Interior face: compute neighbor state
            const FiniteElement *fe2 = fes_->GetFE(e2);
            int dof_offset2 = e2 * ndof_per_el_;

            IntegrationPoint ip2;
            ftr->Loc2.Transform(ip, ip2);
            Vector shape2(ndof);
            fe2->CalcShape(ip2, shape2);

            real_t Q_nbr[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               Q_nbr[c] = 0.0;
               for (int i = 0; i < ndof; i++)
               {
                  Q_nbr[c] += shape2(i) * Q_data[c * ndof_total_ + dof_offset2 + i];
               }
            }

            // Compute Godunov flux from Elem1's perspective
            flux_.Interior(nor, Q_self, Q_nbr, F_h);

            // Subtract from Elem1 RHS
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  rhs[c * ndof_total_ + dof_offset1 + i] -= w * shape1(i) * F_h[c];
               }
            }

            // Add to Elem2 RHS (flux with reversed normal = -F_h for Elem2)
            // From Elem2's perspective, the normal is -nor, so the flux is:
            // F_h_2 = A_{-n}^+ Q_nbr + A_{-n}^- Q_self = -(A_n^- Q_nbr + A_n^+ Q_self)
            // But the face integral for Elem2 also has the opposite sign due to
            // the outward normal convention. So the contribution to Elem2 is:
            // rhs_2 += w * shape2 * F_h  (positive, not negative)
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  rhs[c * ndof_total_ + dof_offset2 + i] += w * shape2(i) * F_h[c];
               }
            }
         }
      }
   }
}

// ---------------------------------------------------------------------------
// Apply per-element inverse mass matrix
// ---------------------------------------------------------------------------
void WaveOperator::ApplyMassInverse(Vector &dQdt) const
{
   Vector elem_rhs(ndof_per_el_);
   Vector elem_result(ndof_per_el_);

   for (int e = 0; e < ne_; e++)
   {
      int dof_offset = e * ndof_per_el_;
      const DenseMatrix &Minv = elem_mass_inv_[e];

      for (int c = 0; c < NUM_STATE; c++)
      {
         // Extract element RHS for component c
         for (int i = 0; i < ndof_per_el_; i++)
         {
            elem_rhs(i) = dQdt[c * ndof_total_ + dof_offset + i];
         }

         // Apply inverse mass: result = M^{-1} * rhs
         Minv.Mult(elem_rhs, elem_result);

         // Write back
         for (int i = 0; i < ndof_per_el_; i++)
         {
            dQdt[c * ndof_total_ + dof_offset + i] = elem_result(i);
         }
      }
   }
}

// ---------------------------------------------------------------------------
// Assemble per-element inverse mass matrices
// ---------------------------------------------------------------------------
void WaveOperator::AssembleElementMassInverse()
{
   elem_mass_inv_.resize(ne_);

   for (int e = 0; e < ne_; e++)
   {
      const FiniteElement *fe = fes_->GetFE(e);
      ElementTransformation *Tr = fes_->GetElementTransformation(e);
      int ndof = fe->GetDof();

      // Integration rule for mass matrix (order 2*p)
      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2*order_);

      DenseMatrix M(ndof);
      M = 0.0;

      Vector shape(ndof);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetIntPoint(&ip);
         real_t w = ip.weight * Tr->Weight();

         fe->CalcShape(ip, shape);

         // M_kl += w * shape(k) * shape(l)
         AddMult_a_VVt(w, shape, M);
      }

      // Invert: M^{-1}
      elem_mass_inv_[e].SetSize(ndof);
      DenseMatrixInverse M_solver(M);
      M_solver.GetInverseMatrix(elem_mass_inv_[e]);
   }
}

// ---------------------------------------------------------------------------
// CFL time step
// ---------------------------------------------------------------------------
real_t WaveOperator::ComputeMaxDt(real_t cfl) const
{
   // dt <= cfl * h_min / c_p, where cfl ~ 1/(2N+1) for DG of order N
   return cfl * h_min_ / flux_.GetCp();
}

// ---------------------------------------------------------------------------
// Classify boundary face by attribute
// ---------------------------------------------------------------------------
WaveOperator::FaceBC WaveOperator::ClassifyBoundaryFace(int bdr_attr) const
{
   if (bc_.absorbing_attrs.count(bdr_attr)) { return FaceBC::Absorbing; }
   if (bc_.natural_attrs.count(bdr_attr))   { return FaceBC::FreeSurface; }
   if (bdr_attr == bc_.fault_attr)          { return FaceBC::Fault; }
   // Default: absorbing (safe fallback)
   return FaceBC::Absorbing;
}

} // namespace seas
} // namespace mfem
