// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// rheology_field.hpp — where the rheology state lives (rheology plan rev 1,
// sections 4.3 and 6.4).
//
// The state is stored at the nodes of the discontinuous space of degree
// p - 1 on each tetrahedron (Gauss-Legendre nodes, all inside the element:
// one at p = 1, 4 at p = 2, 10 at p = 3).  On an affine tetrahedron the
// strain of a degree-p displacement is a polynomial of degree p - 1, so
// sampling it at these nodes is exact and the stress correction
//
//    Delta sigma_h(x) = sum_i psi_i(x) Delta sigma_i
//
// (psi_i the nodal basis) is its interpolant.  Per node the field keeps the
// context a pointwise model may read (position, sigma_0', C_0 of the
// element), the model state q and Delta sigma; per element the active flag
// (Delta sigma nonzero at any node).  Delta sigma is a ParGridFunction with
// 6 components so that ExchangeFaceNbrData gives the neighbour values the
// shared faces of P and Q need.
//
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_RHEOLOGY_FIELD_HPP
#define MFEM_SEAS_RHEOLOGY_FIELD_HPP

#include "mfem.hpp"
#include "rheology_stack.hpp"
#include "../domain/nitsche_material_sample.hpp"

#include <functional>
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

class RheologyField
{
public:
   /// sigma_0'(x) as a row-major 3 x 3 tensor (Pa, tension positive), the
   /// function that also gives the fault prestress (plan section 6.3).
   using StressFn = std::function<void(const real_t x[3], real_t s[9])>;

   /// Collective (face-neighbour exchange of the space).
   RheologyField(ParMesh &pm, int order, const NitscheMaterialSample &sample,
                 const RheologyStack &stack, const StressFn &sigma0_eff)
      : pm_(pm), stack_(stack),
        fec_(order - 1, pm.Dimension(), BasisType::GaussLegendre),
        fes_(&pm, &fec_, 6, Ordering::byVDIM),
        dsig_(&fes_)
   {
      MFEM_VERIFY(pm.Dimension() == 3, "RheologyField: 3-D meshes only");
      MFEM_VERIFY(order >= 1, "RheologyField: order " << order << " < 1");
      MFEM_VERIFY(!stack.Empty(), "RheologyField: empty rheology stack");
      MFEM_VERIFY(static_cast<bool>(sigma0_eff), "RheologyField: no sigma_0'");
      ne_ = pm.GetNE();
      nr_ = fec_.FiniteElementForGeometry(Geometry::TETRAHEDRON)->GetDof();
      ns_ = stack.NumStateVars();
      const IntegrationRule &nodes =
         fec_.FiniteElementForGeometry(Geometry::TETRAHEDRON)->GetNodes();
      ctx_.resize(static_cast<std::size_t>(ne_) * nr_);
      q_.assign(static_cast<std::size_t>(ne_) * nr_ * ns_, 0.0);
      active_.assign(ne_, 0);
      Array<int> dofs;
      Vector x(3);
      for (int e = 0; e < ne_; e++)
      {
         MFEM_VERIFY(pm.GetElementBaseGeometry(e) == Geometry::TETRAHEDRON,
                     "RheologyField: element " << e << " is not a tetrahedron");
         // node (e, i) is DOF e nr + i, its 6 components at 6 (e nr + i) + c
         fes_.GetElementDofs(e, dofs);
         for (int i = 0; i < nr_; i++)
         {
            MFEM_VERIFY(dofs[i] == e * nr_ + i, "RheologyField: unexpected DOF numbering");
         }
         ElementTransformation *T = pm.GetElementTransformation(e);
         const IntegrationPoint &c = Geometries.GetCenter(Geometry::TETRAHEDRON);
         const real_t lam = sample.LambdaGF().GetValue(e, c);
         const real_t mu = sample.MuGF().GetValue(e, c);
         for (int i = 0; i < nr_; i++)
         {
            RheologyPointContext &ctx = ctx_[Node(e, i)];
            T->Transform(nodes.IntPoint(i), x);
            real_t s9[9];
            for (int a = 0; a < 3; a++) { ctx.x[a] = x(a); }
            sigma0_eff(ctx.x, s9);
            sym6::FromRowMajor(s9, ctx.sigma0_eff);
            ctx.lambda0 = lam;
            ctx.mu0 = mu;
            stack.InitState(ctx, State(e, i));
         }
      }
      dsig_ = 0.0;
      pm.ExchangeFaceNbrData();
      fes_.ExchangeFaceNbrData();
      const int nnbr = pm.GetNFaceNeighborElements();
      nbr_active_.assign(nnbr, 0);
      nbr_vdofs_.resize(static_cast<std::size_t>(nnbr) * 6 * nr_);
      Array<int> vd;
      for (int k = 0; k < nnbr; k++)
      {
         fes_.GetFaceNbrElementVDofs(k, vd);
         MFEM_VERIFY(vd.Size() == 6 * nr_, "RheologyField: face neighbour " << k
                     << " has " << vd.Size() << " DOFs, expected " << 6 * nr_);
         for (int a = 0; a < 6 * nr_; a++) { nbr_vdofs_[static_cast<std::size_t>(k) * 6 * nr_ + a] = vd[a]; }
      }
   }

   int NumElements() const { return ne_; }
   int NodesPerElement() const { return nr_; }
   int StateVarsPerNode() const { return ns_; }
   const RheologyStack &Stack() const { return stack_; }

   /// The scalar element of the rheology space (nodes and basis psi_i).
   const FiniteElement &NodeElement() const
   {
      return *fec_.FiniteElementForGeometry(Geometry::TETRAHEDRON);
   }

   std::size_t Node(int e, int i) const { return static_cast<std::size_t>(e) * nr_ + i; }
   RheologyPointContext &Context(int e, int i) { return ctx_[Node(e, i)]; }
   const RheologyPointContext &Context(int e, int i) const { return ctx_[Node(e, i)]; }
   real_t *State(int e, int i) { return q_.data() + Node(e, i) * ns_; }
   const real_t *State(int e, int i) const { return q_.data() + Node(e, i) * ns_; }

   /// Delta sigma of node (e, i), [xx, yy, zz, xy, yz, xz].
   real_t *Dsig(int e, int i) { return dsig_.GetData() + 6 * Node(e, i); }
   const real_t *Dsig(int e, int i) const { return dsig_.GetData() + 6 * Node(e, i); }

   bool Active(int e) const { return active_[e] != 0; }
   void SetActive(int e, bool on) { active_[e] = on ? 1 : 0; }

   /// Collective: neighbour Delta sigma for the shared faces.
   void ExchangeNeighbours()
   {
      dsig_.ExchangeFaceNbrData();
      const Vector &nd = dsig_.FaceNbrData();
      for (std::size_t k = 0; k < nbr_active_.size(); k++)
      {
         char on = 0;
         for (int a = 0; a < 6 * nr_ && !on; a++)
         {
            on = nd(nbr_vdofs_[k * 6 * nr_ + a]) != 0.0;
         }
         nbr_active_[k] = on;
      }
   }

   bool NeighbourActive(int k) const { return nbr_active_[k] != 0; }

   /// Delta sigma component c of node i of face-neighbour element k (after
   /// ExchangeNeighbours).
   real_t NeighbourDsig(int k, int i, int c) const
   {
      return dsig_.FaceNbrData()(nbr_vdofs_[static_cast<std::size_t>(k) * 6 * nr_ + c * nr_ + i]);
   }

private:
   ParMesh &pm_;
   const RheologyStack &stack_;
   L2_FECollection fec_;
   ParFiniteElementSpace fes_;
   ParGridFunction dsig_;
   int ne_ = 0, nr_ = 0, ns_ = 0;
   std::vector<RheologyPointContext> ctx_;
   std::vector<real_t> q_;
   std::vector<char> active_, nbr_active_;
   std::vector<int> nbr_vdofs_;   // face-neighbour element vdofs, [c nr + i]
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_RHEOLOGY_FIELD_HPP
