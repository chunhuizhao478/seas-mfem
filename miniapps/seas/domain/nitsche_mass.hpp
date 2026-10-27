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
// Damping blocks (C-02, P4): for elements that own an absorbing-boundary
// dashpot block C_e (3 nd x 3 nd, element-local byNODES, see
// integrator/dg_absorbing_bc_integrator.hpp), the constructor optionally
// inverts the dense SPD block  rho_e M_e (x) I_3 + half_dt * C_e  once, and
// Apply uses it there.  Every other element keeps the scalar path; with no
// blocks the result is bit-identical to P2/P3.
//
// Device (GPU dynamic-rupture plan D2): with an enabled mfem::Device and r or
// a using it, Apply runs one mfem::forall over the elements with the same
// stored inverses, copied once to flat device arrays on the first such call,
// and the same operation order as the host loop (kernels::Mult is the
// product of DenseMatrix::Mult); the debug device therefore gives the host
// bits.  Without a device the host loop runs, unchanged.
//
// Local time stepping (profiling plan item 7, lts-kernels-design-rev1 §3.4
// step 6): half_dt_per_elem, when given, replaces half_dt per element, so
// each damped element uses half of its own cluster's step.  Without it every
// bit is as before.
//
// Scope (Chunhui Zhao, 2026-10-13): tetrahedral meshes only.

#ifndef MFEM_SEAS_NITSCHE_MASS_HPP
#define MFEM_SEAS_NITSCHE_MASS_HPP

#include "mfem.hpp"
#include "linalg/kernels.hpp"   // kernels::Mult in the device kernels (plan D2)
#include "nitsche_material_sample.hpp"

#include <map>
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

class NitscheMass
{
public:
   NitscheMass(ParMesh &mesh, int order, const Vector &rho_per_elem,
               const std::map<int, DenseMatrix> *damping_blocks = nullptr,
               real_t half_dt = 0.0, const Vector *half_dt_per_elem = nullptr)
      : fec_(order, mesh.Dimension(), BasisType::GaussLobatto),
        fes_(&mesh, &fec_), half_dt_(half_dt)
   {
      const int ne = mesh.GetNE();
      MFEM_VERIFY(rho_per_elem.Size() == ne,
                  "NitscheMass: rho_per_elem has size " << rho_per_elem.Size()
                  << " but the mesh has " << ne << " local elements.");
      MFEM_VERIFY(!half_dt_per_elem || half_dt_per_elem->Size() == ne,
                  "NitscheMass: half_dt_per_elem has size " << half_dt_per_elem->Size()
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
         if (e == 0) { nd_ = nd; }
         MFEM_VERIFY(nd == nd_, "NitscheMass: element " << e << " has " << nd
                     << " DOFs, element 0 " << nd_ << " (one order on the mesh)");
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

         if (damping_blocks)
         {
            auto it = damping_blocks->find(e);
            if (it != damping_blocks->end())
            {
               const DenseMatrix &Ce = it->second;
               MFEM_VERIFY(Ce.Height() == 3 * nd && Ce.Width() == 3 * nd,
                           "NitscheMass: damping block of element " << e
                           << " is " << Ce.Height() << " x " << Ce.Width()
                           << ", expected " << 3 * nd << " x " << 3 * nd);
               DenseMatrix B(3 * nd);
               B = 0.0;
               for (int c = 0; c < 3; ++c)
                  for (int k = 0; k < nd; ++k)
                     for (int l = 0; l < nd; ++l)
                     {
                        B(c * nd + k, c * nd + l) = rho_per_elem(e) * M(k, l);
                     }
               B.Add(half_dt_per_elem ? (*half_dt_per_elem)(e) : half_dt, Ce);
               DenseMatrix &Binv = block_inv_[e];
               Binv.SetSize(3 * nd);
               DenseMatrixInverse B_inv(B);
               B_inv.GetInverseMatrix(Binv);
            }
         }

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
      if (Device::IsEnabled() && (r.UseDevice() || a.UseDevice()))
      {
         if (minv_.empty()) { return; }   // a rank without elements
         BuildDevice_();
         switch (nd_)
         {
            case 4: ApplyDevice_<4>(r, a); return;     // p = 1
            case 10: ApplyDevice_<10>(r, a); return;   // p = 2
            case 20: ApplyDevice_<20>(r, a); return;   // p = 3
            default:
               MFEM_ABORT("NitscheMass::Apply: no device kernel for " << nd_
                          << " DOFs per element (p = 1, 2, 3 only)");
         }
      }
      Vector x, y;
      for (size_t e = 0; e < minv_.size(); ++e)
      {
         const Array<int> &dofs = elem_dofs_[e];
         const int nd = dofs.Size();
         if (!block_inv_.empty())
         {
            auto it = block_inv_.find(static_cast<int>(e));
            if (it != block_inv_.end())
            {
               // Damped element: (rho_e M_e (x) I_3 + half_dt C_e)^{-1}.
               x.SetSize(3 * nd);
               y.SetSize(3 * nd);
               for (int c = 0; c < 3; ++c)
                  for (int k = 0; k < nd; ++k)
                  {
                     x(c * nd + k) = r(c * ndofs_ + dofs[k]);
                  }
               it->second.Mult(x, y);
               for (int c = 0; c < 3; ++c)
                  for (int k = 0; k < nd; ++k)
                  {
                     a(c * ndofs_ + dofs[k]) = y(c * nd + k);
                  }
               continue;
            }
         }
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
   /// half_dt the damping blocks were built with (0 without blocks).
   real_t HalfDt() const { return half_dt_; }   // the scalar, not half_dt_per_elem
   int NumDampedElements() const { return static_cast<int>(block_inv_.size()); }
   /// (rho_e M_e (x) I_3 + half_dt C_e)^{-1} of damped element e (3 nd x
   /// 3 nd, element-local byNODES), or null for an undamped element.
   const DenseMatrix *ElementBlockInverse(int e) const
   {
      auto it = block_inv_.find(e);
      return it == block_inv_.end() ? nullptr : &it->second;
   }

   /// Scalar DOFs per element (one order and geometry on the mesh).
   int ElementDofs() const { return nd_; }
   /// ElementScalarDofs of every element, element-major (GetNE() x nd), as
   /// a device array (plan D2); built on the first call.
   const Array<int> &DeviceElementDofs() const { BuildDevice_(); return d_dofs_; }

private:
   /// Flat copies of elem_dofs_, minv_ (column-major, as DenseMatrix
   /// stores it), 1/rho and the damped blocks for the device kernel; d_blk_
   /// holds each element's damped block (-1: the scalar path).
   void BuildDevice_() const
   {
      if (d_built_) { return; }
      const int ne = static_cast<int>(minv_.size()), nd = nd_;
      d_dofs_.SetSize(ne * nd);
      d_minv_.SetSize(ne * nd * nd);
      d_inv_rho_.SetSize(ne);
      d_blk_.SetSize(ne);
      d_binv_.SetSize(static_cast<int>(block_inv_.size()) * 9 * nd * nd);
      int *hd = d_dofs_.HostWrite(), *hb = d_blk_.HostWrite();
      real_t *hm = d_minv_.HostWrite(), *hr = d_inv_rho_.HostWrite(),
              *hbi = d_binv_.HostWrite();
      for (int e = 0; e < ne; ++e)
      {
         for (int k = 0; k < nd; ++k) { hd[e * nd + k] = elem_dofs_[e][k]; }
         const real_t *m = minv_[e].Data();
         for (int i = 0; i < nd * nd; ++i) { hm[static_cast<size_t>(e) * nd * nd + i] = m[i]; }
         hr[e] = inv_rho_(e);
         hb[e] = -1;
      }
      int b = 0;
      for (const auto &kv : block_inv_)
      {
         hb[kv.first] = b;
         const real_t *m = kv.second.Data();
         for (int i = 0; i < 9 * nd * nd; ++i) { hbi[static_cast<size_t>(b) * 9 * nd * nd + i] = m[i]; }
         ++b;
      }
      d_built_ = true;
   }

public:
   // The function holding the kernel is public: nvcc refuses an extended
   // __host__ __device__ lambda inside a private or protected member.
   template <int ND>
   void ApplyDevice_(const Vector &r, Vector &a) const
   {
      const int ne = static_cast<int>(minv_.size()), n = ndofs_;
      const real_t *dr = r.Read();
      real_t *da = a.Write();
      const int *dofs = d_dofs_.Read(), *blk = d_blk_.Read();
      const real_t *mi = d_minv_.Read(), *ir = d_inv_rho_.Read(), *bi = d_binv_.Read();
      mfem::forall(ne, [=] MFEM_HOST_DEVICE (int e)
      {
         const int *de = dofs + e * ND;
         const int b = blk[e];
         if (b >= 0)
         {
            // damped element: (rho_e M_e (x) I_3 + half_dt C_e)^{-1}
            real_t x[3 * ND], y[3 * ND];
            for (int c = 0; c < 3; ++c)
               for (int k = 0; k < ND; ++k) { x[c * ND + k] = dr[c * n + de[k]]; }
            kernels::Mult(3 * ND, 3 * ND, bi + static_cast<size_t>(b) * 9 * ND * ND, x, y);
            for (int c = 0; c < 3; ++c)
               for (int k = 0; k < ND; ++k) { da[c * n + de[k]] = y[c * ND + k]; }
            return;
         }
         const real_t *m = mi + static_cast<size_t>(e) * ND * ND;
         const real_t s = ir[e];
         for (int c = 0; c < 3; ++c)
         {
            real_t x[ND], y[ND];
            for (int k = 0; k < ND; ++k) { x[k] = dr[c * n + de[k]]; }
            kernels::Mult(ND, ND, m, x, y);
            for (int k = 0; k < ND; ++k) { da[c * n + de[k]] = y[k] * s; }
         }
      });
   }

private:
   DG_FECollection fec_;
   ParFiniteElementSpace fes_;
   int ndofs_ = 0;
   int nd_ = 0;   ///< scalar DOFs per element
   real_t half_dt_ = 0.0;
   Vector inv_rho_;
   std::map<int, DenseMatrix> block_inv_;
   std::vector<DenseMatrix> minv_;
   std::vector<Array<int>> elem_dofs_;
   mutable bool d_built_ = false;
   mutable Array<int> d_dofs_, d_blk_;
   mutable Vector d_minv_, d_inv_rho_, d_binv_;
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
