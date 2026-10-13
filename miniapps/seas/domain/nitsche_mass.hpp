// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_mass.hpp — Nitsche dynamic rupture, plan item C-02.
//
// Density-weighted inverse of the DG vector mass matrix,
//
//    M = diag_e( rho_e M_e (x) I_3 ),   M_e(k,l) = int_e phi_k phi_l dx,
//
// applied component-wise:  a_c|_e = (1/rho_e) * (M_e^{-1} r_c|_e),  c = 0,1,2.
//
// Why not ElasticityDomainOperator::PrecomputeMassInverse
// (elasticity_operator_setup.inl): it carries no density and is only built
// on the BR2 path.  This class mirrors it (same scalar DG space, same rule
// IntRules.Get(geom, 2*order), weight ip.weight * T->Weight(),
// DenseMatrixInverse) and adds rho_e.
//
// rho is passed per LOCAL element (a Vector of size GetNE()), so "rho is
// constant on each element" is part of the interface.  Build it from the
// same per-element sample the stiffness uses with RhoPerElement() below.
//
// Layout: the operator's vector DG space is Ordering::byNODES, so component
// c of scalar dof i sits at c * ndofs + i (ndofs = scalar NDofs).
// Element-local: no face-neighbour data, no communication.
//
// The damping-block variant (M + dt/2 C)^{-1} (plan rev 7) is not part of
// P2; it arrives with the absorbing boundary in P4.
//
// Scope (Chunhui Zhao, 2026-10-13): tetrahedral meshes only.

#ifndef MFEM_SEAS_NITSCHE_MASS_HPP
#define MFEM_SEAS_NITSCHE_MASS_HPP

#include "mfem.hpp"
#include "nitsche_wave_operator.hpp"

#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

class NitscheMass
{
public:
   NitscheMass(ParMesh &mesh, int order, const Vector &rho_per_elem)
      : fec_(order, mesh.Dimension(), BasisType::GaussLobatto),
        fes_(&mesh, &fec_)
   {
      const int ne = mesh.GetNE();
      MFEM_VERIFY(rho_per_elem.Size() == ne,
                  "NitscheMass: rho_per_elem has size " << rho_per_elem.Size()
                  << " but the mesh has " << ne << " local elements.");
      for (int e = 0; e < ne; ++e)
      {
         MFEM_VERIFY(mesh.GetElementBaseGeometry(e) == Geometry::TETRAHEDRON,
                     "NitscheMass: element " << e << " is a "
                     << Geometry::Name[mesh.GetElementBaseGeometry(e)]
                     << "; the Nitsche path supports tetrahedra only.");
         MFEM_VERIFY(rho_per_elem(e) > 0.0,
                     "NitscheMass: rho (" << rho_per_elem(e)
                     << ") must be > 0 at element " << e);
      }

      ndofs_ = fes_.GetNDofs();
      inv_rho_.SetSize(ne);
      minv_.resize(ne);
      elem_dofs_.resize(ne);

      for (int e = 0; e < ne; ++e)
      {
         const FiniteElement *fe = fes_.GetFE(e);
         ElementTransformation *T = mesh.GetElementTransformation(e);
         const int nd = fe->GetDof();
         DenseMatrix M(nd);
         M = 0.0;
         Vector shape(nd);

         const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(),
                                                  2 * fe->GetOrder());
         for (int j = 0; j < ir.GetNPoints(); ++j)
         {
            const IntegrationPoint &ip = ir.IntPoint(j);
            T->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);
            const real_t w = ip.weight * T->Weight();
            for (int k = 0; k < nd; ++k)
            {
               for (int l = 0; l < nd; ++l)
               {
                  M(k, l) += w * shape(k) * shape(l);
               }
            }
         }

         minv_[e].SetSize(nd);
         DenseMatrixInverse M_inv(M);
         M_inv.GetInverseMatrix(minv_[e]);

         fes_.GetElementDofs(e, elem_dofs_[e]);
         inv_rho_(e) = 1.0 / rho_per_elem(e);
      }
   }

   /// a = M^{-1} r, component-wise.  Fixed operation order per element and
   /// component: y = M_e^{-1} * r_c|_e first, then a_c|_e = y * (1/rho_e),
   /// with 1/rho_e computed once in the constructor.
   void Apply(const Vector &r, Vector &a) const
   {
      MFEM_VERIFY(r.Size() == 3 * ndofs_,
                  "NitscheMass::Apply: input size " << r.Size()
                  << " != 3 * ndofs = " << 3 * ndofs_);
      a.SetSize(r.Size());
      Vector x, y;
      for (size_t e = 0; e < minv_.size(); ++e)
      {
         const Array<int> &dofs = elem_dofs_[e];
         const int nd = dofs.Size();
         x.SetSize(nd);
         y.SetSize(nd);
         for (int c = 0; c < 3; ++c)
         {
            const int off = c * ndofs_;
            for (int k = 0; k < nd; ++k) { x(k) = r(off + dofs[k]); }
            minv_[e].Mult(x, y);
            for (int k = 0; k < nd; ++k) { a(off + dofs[k]) = y(k) * inv_rho_(e); }
         }
      }
   }

   int NumScalarDofs() const { return ndofs_; }
   int VSize() const { return 3 * ndofs_; }
   const Array<int> &ElementScalarDofs(int e) const { return elem_dofs_[e]; }
   const DenseMatrix &ElementMassInverse(int e) const { return minv_[e]; }

private:
   DG_FECollection fec_;
   ParFiniteElementSpace fes_;
   int ndofs_ = 0;
   Vector inv_rho_;
   std::vector<DenseMatrix> minv_;
   std::vector<Array<int>> elem_dofs_;
};

/// Per-local-element density from the stiffness's material sample, so M and
/// K share one sampling.
inline void RhoPerElement(const NitscheMaterialSample &sample, Vector &rho)
{
   const ParGridFunction &gf = sample.RhoGF();
   const ParFiniteElementSpace *fes = gf.ParFESpace();
   const int ne = fes->GetNE();
   rho.SetSize(ne);
   Array<int> dofs;
   for (int e = 0; e < ne; ++e)
   {
      fes->GetElementDofs(e, dofs);
      MFEM_ASSERT(dofs.Size() == 1, "order-0 L2 has one DOF per element");
      rho(e) = gf(dofs[0]);
   }
}

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_MASS_HPP
