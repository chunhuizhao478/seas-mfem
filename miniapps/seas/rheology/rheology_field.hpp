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

#include <algorithm>
#include <cmath>
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
      MFEM_VERIFY(pm.Dimension() == 2 || pm.Dimension() == 3,
                  "RheologyField: 2-D (plane strain) or 3-D meshes only");
      geom_ = (pm.Dimension() == 3) ? Geometry::TETRAHEDRON : Geometry::TRIANGLE;
      MFEM_VERIFY(order >= 1, "RheologyField: order " << order << " < 1");
      MFEM_VERIFY(!stack.Empty(), "RheologyField: empty rheology stack");
      MFEM_VERIFY(static_cast<bool>(sigma0_eff), "RheologyField: no sigma_0'");
      ne_ = pm.GetNE();
      nr_ = fec_.FiniteElementForGeometry(geom_)->GetDof();
      ns_ = stack.NumStateVars();
      const IntegrationRule &nodes = fec_.FiniteElementForGeometry(geom_)->GetNodes();
      ctx_.resize(static_cast<std::size_t>(ne_) * nr_);
      q_.assign(static_cast<std::size_t>(ne_) * nr_ * ns_, 0.0);
      active_.assign(ne_, 0);
      Array<int> dofs;
      Vector x(pm.SpaceDimension());
      for (int e = 0; e < ne_; e++)
      {
         MFEM_VERIFY(pm.GetElementBaseGeometry(e) == geom_,
                     "RheologyField: element " << e << " is not a "
                     << (geom_ == Geometry::TETRAHEDRON ? "tetrahedron" : "triangle"));
         // node (e, i) is DOF e nr + i, its 6 components at 6 (e nr + i) + c
         fes_.GetElementDofs(e, dofs);
         for (int i = 0; i < nr_; i++)
         {
            MFEM_VERIFY(dofs[i] == e * nr_ + i, "RheologyField: unexpected DOF numbering");
         }
         ElementTransformation *T = pm.GetElementTransformation(e);
         const IntegrationPoint &c = Geometries.GetCenter(geom_);
         const real_t lam = sample.LambdaGF().GetValue(e, c);
         const real_t mu = sample.MuGF().GetValue(e, c);
         const real_t rho = sample.RhoGF().GetValue(e, c);
         for (int i = 0; i < nr_; i++)
         {
            RheologyPointContext &ctx = ctx_[Node(e, i)];
            T->Transform(nodes.IntPoint(i), x);
            real_t s9[9];
            for (int a = 0; a < 3; a++) { ctx.x[a] = (a < x.Size()) ? x(a) : 0.0; }
            sigma0_eff(ctx.x, s9);
            sym6::FromRowMajor(s9, ctx.sigma0_eff);
            ctx.lambda0 = lam;
            ctx.mu0 = mu;
            ctx.rho0 = rho;
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
      return *fec_.FiniteElementForGeometry(geom_);
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
   void SetNeighbourActive(int k, bool on) { nbr_active_[k] = on ? 1 : 0; }
   int NumNeighbours() const { return static_cast<int>(nbr_active_.size()); }

   // Ghost memory (attenuation plan rev6, 5.8; local time stepping only):
   // this rank's copy of the state of each face-neighbour element, stepped
   // from the displacement halo.  BuildGhosts exchanges the neighbours'
   // node contexts and exclusion flags once (collective); the ghost dsig is
   // written into the face-neighbour data that NeighbourDsig reads.
   void BuildGhosts(const std::vector<char> &excluded)
   {
      MFEM_VERIFY(static_cast<int>(excluded.size()) == ne_, "RheologyField::BuildGhosts: "
                  << excluded.size() << " flags for " << ne_ << " elements");
      constexpr int nv = 13;   // x (3), sigma0_eff (6), lambda0, mu0, rho0, excluded
      ParFiniteElementSpace cfes(&pm_, &fec_, nv, Ordering::byVDIM);
      ParGridFunction cgf(&cfes);
      Array<int> vd;
      for (int e = 0; e < ne_; e++)
      {
         cfes.GetElementVDofs(e, vd);
         for (int i = 0; i < nr_; i++)
         {
            const RheologyPointContext &ctx = Context(e, i);
            real_t v[nv];
            for (int a = 0; a < 3; a++) { v[a] = ctx.x[a]; }
            for (int a = 0; a < 6; a++) { v[3 + a] = ctx.sigma0_eff[a]; }
            v[9] = ctx.lambda0;
            v[10] = ctx.mu0;
            v[11] = ctx.rho0;
            v[12] = excluded[e] ? 1.0 : 0.0;
            for (int a = 0; a < nv; a++) { cgf(vd[a * nr_ + i]) = v[a]; }
         }
      }
      cfes.ExchangeFaceNbrData();
      cgf.ExchangeFaceNbrData();
      const Vector &nd = cgf.FaceNbrData();
      const int nnbr = NumNeighbours();
      nbr_ctx_.assign(static_cast<std::size_t>(nnbr) * nr_, RheologyPointContext());
      nbr_excluded_.assign(nnbr, 0);
      nbr_q_.assign(static_cast<std::size_t>(nnbr) * nr_ * ns_, 0.0);
      for (int k = 0; k < nnbr; k++)
      {
         cfes.GetFaceNbrElementVDofs(k, vd);
         MFEM_VERIFY(vd.Size() == nv * nr_, "RheologyField::BuildGhosts: neighbour " << k);
         for (int i = 0; i < nr_; i++)
         {
            RheologyPointContext &ctx = nbr_ctx_[static_cast<std::size_t>(k) * nr_ + i];
            for (int a = 0; a < 3; a++) { ctx.x[a] = nd(vd[a * nr_ + i]); }
            for (int a = 0; a < 6; a++) { ctx.sigma0_eff[a] = nd(vd[(3 + a) * nr_ + i]); }
            ctx.lambda0 = nd(vd[9 * nr_ + i]);
            ctx.mu0 = nd(vd[10 * nr_ + i]);
            ctx.rho0 = nd(vd[11 * nr_ + i]);
            if (nd(vd[12 * nr_ + i]) != 0.0) { nbr_excluded_[k] = 1; }
         }
      }
      ghosts_ = true;
   }
   bool HasGhosts() const { return ghosts_; }
   bool NeighbourExcluded(int k) const { return nbr_excluded_[k] != 0; }
   RheologyPointContext &NeighbourContext(int k, int i)
   {
      return nbr_ctx_[static_cast<std::size_t>(k) * nr_ + i];
   }
   real_t *NeighbourState(int k, int i)
   {
      return nbr_q_.data() + (static_cast<std::size_t>(k) * nr_ + i) * ns_;
   }
   /// Writable ghost dsig of node i of face-neighbour element k, component c.
   real_t &NeighbourDsigRef(int k, int i, int c)
   {
      return dsig_.FaceNbrData()(nbr_vdofs_[static_cast<std::size_t>(k) * 6 * nr_ + c * nr_ + i]);
   }
   /// Test of the ghosts (collective): the largest |ghost - owner| over
   /// every ghost dsig component after a face-neighbour exchange of the
   /// owners' values, which then replace the ghosts (equal when they agree).
   real_t GhostMismatch()
   {
      Vector g(dsig_.FaceNbrData());
      dsig_.ExchangeFaceNbrData();
      const Vector &o = dsig_.FaceNbrData();
      real_t m = 0.0;
      for (int a = 0; a < g.Size(); a++) { m = std::max(m, std::abs(g(a) - o(a))); }
      MPI_Allreduce(MPI_IN_PLACE, &m, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, pm_.GetComm());
      return m;
   }

   /// Delta sigma component c of node i of face-neighbour element k (after
   /// ExchangeNeighbours).
   real_t NeighbourDsig(int k, int i, int c) const
   {
      return dsig_.FaceNbrData()(nbr_vdofs_[static_cast<std::size_t>(k) * 6 * nr_ + c * nr_ + i]);
   }

private:
   ParMesh &pm_;
   const RheologyStack &stack_;
   Geometry::Type geom_ = Geometry::TETRAHEDRON;
   L2_FECollection fec_;
   ParFiniteElementSpace fes_;
   ParGridFunction dsig_;
   int ne_ = 0, nr_ = 0, ns_ = 0;
   std::vector<RheologyPointContext> ctx_;
   std::vector<real_t> q_;
   std::vector<char> active_, nbr_active_;
   std::vector<int> nbr_vdofs_;   // face-neighbour element vdofs, [c nr + i]
   bool ghosts_ = false;
   std::vector<RheologyPointContext> nbr_ctx_;   // ghost node contexts
   std::vector<char> nbr_excluded_;
   std::vector<real_t> nbr_q_;                   // ghost state
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_RHEOLOGY_FIELD_HPP
