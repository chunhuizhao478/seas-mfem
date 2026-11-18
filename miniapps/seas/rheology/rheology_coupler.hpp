// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// rheology_coupler.hpp — box R of the time step and the two loads the
// Nitsche operator calls (rheology plan rev 1, sections 5 and 6.6).
//
// ImplicitSolve(x = U^{n+1}, ...), t = t_{n+1}, dt = h:
//
//    R   Advance(U^{n+1}, t, dt):
//          eps = S U^{n+1} at every rheology node
//          per node: stack.Update(ctx, eps, dt, q) -> q^{n+1}, dsig^{n+1}
//          active elements: dsig nonzero at any node; one reduction says
//          whether any element is active on any rank, and only then the
//          neighbour dsig is exchanged (collective)
//   3-11 the fault, with T += Q dsig      (NitscheWaveOperator::AdvanceFault_)
//   12-13 M a = ... + P dsig              (AddLoad, NitscheWaveOperator::Acceleration_)
//
// While nothing has yielded AnyActive() is false on every rank, the
// operator skips both loads, and every bit is that of the run without the
// coupler (test O1).  The DG operators are in dg_rheology_operators.hpp, the
// state in rheology_field.hpp.
//
// Attenuation plan rev6, section 5.5 (2026-11-17):
//  - Reset(U_entry): a new run or a cycle's dynamic phase.  Every node's
//    state is InitState'd again, the entry strain stored (SetEntryStrain),
//    dsig = 0 and the time check restarted (DynPhase restarts the operator's
//    clock at 0 in every phase).
//  - SetExcluded(mask): elements that carry no rheology (the PML band,
//    decision D4); their dsig stays 0.
//  - A stack that is linear in the strain (attenuation alone) has every
//    element outside the mask active from its first step: the flags are set
//    once (one reduction at the first Refresh) and the per-step reduction is
//    skipped; the neighbour exchange stays.  YieldingNodes() is then this
//    rank's count.
//
// Local time stepping (attenuation plan rev6, section 5.8, phase V6; a
// stack linear in the strain only): SetLts(cl, cl_nbr) gives each local and
// face-neighbour element its cluster.  At a tick where clusters 0..top end,
// AdvanceLts steps the memory of every element of those clusters with its
// cluster's dt_c, and the ghost copy of each face-neighbour element of those
// clusters from the displacement halo (no message: both ranks do the same
// arithmetic on the same values).  AddLoadCluster(c, F) adds the owner
// split (5.6) of the elements of cluster c (DGRheologyOperators::AddOwned).
// Global stepping does not call these.
//
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_RHEOLOGY_COUPLER_HPP
#define MFEM_SEAS_RHEOLOGY_COUPLER_HPP

#include "mfem.hpp"
#include "dg_rheology_operators.hpp"
#include "rheology_field.hpp"
#include "rheology_stack.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

class RheologyCoupler
{
public:
   /// ufes: the displacement space of the operator; dirichlet: the
   /// boundary attributes with the SIPG boundary integrator; sigma0_eff:
   /// the function that also gives the fault prestress.  Collective.
   RheologyCoupler(ParFiniteElementSpace &ufes, const NitscheMaterialSample &sample,
                   const std::set<int> &dirichlet, std::unique_ptr<RheologyStack> stack,
                   const RheologyField::StressFn &sigma0_eff)
      : stack_(std::move(stack)), comm_(ufes.GetComm()),
        field_(*ufes.GetParMesh(), ufes.FEColl()->GetOrder(), sample, *stack_, sigma0_eff),
        ops_(ufes, field_, dirichlet)
   {
      eps_.resize(6 * static_cast<std::size_t>(field_.NodesPerElement()));
      linear_ = stack_->IsLinearInStrain();
      excluded_.assign(field_.NumElements(), 0);
   }

   /// Elements without rheology (mask per local element, 1 = excluded).
   /// Before the first Advance or right after Reset.
   void SetExcluded(const std::vector<char> &mask)
   {
      MFEM_VERIFY(static_cast<int>(mask.size()) == field_.NumElements(),
                  "RheologyCoupler::SetExcluded: " << mask.size() << " flags for "
                  << field_.NumElements() << " elements");
      excluded_ = mask;
      linear_flags_set_ = false;
   }
   bool Excluded(int e) const { return excluded_[e] != 0; }

   /// New run or new dynamic phase from the state U (plan (2.42d)).
   /// Collective.
   void Reset(const Vector &U)
   {
      last_t_ = -std::numeric_limits<real_t>::infinity();
      const int nr = field_.NodesPerElement();
      for (int e = 0; e < field_.NumElements(); e++)
      {
         ops_.Strain(e, U, eps_.data());
         for (int i = 0; i < nr; i++)
         {
            const RheologyPointContext &ctx = field_.Context(e, i);
            stack_->InitState(ctx, field_.State(e, i));
            stack_->SetEntryStrain(ctx, &eps_[6 * i], field_.State(e, i));
            real_t *ds = field_.Dsig(e, i);
            for (int c = 0; c < 6; c++) { ds[c] = 0.0; }
         }
         field_.SetActive(e, false);
      }
      yielding_ = 0;
      yielding_global_ = 0;
      active_global_ = 0;
      linear_flags_set_ = false;
      field_.ExchangeNeighbours();
   }

   // ---------------- local time stepping (file comment)

   /// Clusters of the local and face-neighbour elements; after SetExcluded.
   /// Collective (one exchange of the neighbours' contexts and flags).
   void SetLts(const std::vector<int> &cl, const std::vector<int> &cl_nbr)
   {
      MFEM_VERIFY(linear_, "RheologyCoupler::SetLts: local time stepping needs a stack "
                  "linear in the strain (attenuation alone); \"" << stack_->Describe()
                  << "\" is not (plasticity under LTS is not implemented)");
      MFEM_VERIFY(static_cast<int>(cl.size()) == field_.NumElements()
                  && static_cast<int>(cl_nbr.size()) == field_.NumNeighbours(),
                  "RheologyCoupler::SetLts: " << cl.size() << " / " << cl_nbr.size()
                  << " clusters for " << field_.NumElements() << " / "
                  << field_.NumNeighbours() << " elements");
      cl_ = cl;
      cl_nbr_ = cl_nbr;
      field_.BuildGhosts(excluded_);
      int nc = 0;
      for (int c : cl_) { nc = std::max(nc, c + 1); }
      for (int c : cl_nbr_) { nc = std::max(nc, c + 1); }
      MPI_Allreduce(MPI_IN_PLACE, &nc, 1, MPI_INT, MPI_MAX, comm_);
      own_.assign(nc, std::vector<char>(field_.NumElements(), 0));
      own_nbr_.assign(nc, std::vector<char>(field_.NumNeighbours(), 0));
      for (int e = 0; e < field_.NumElements(); e++) { own_[cl_[e]][e] = 1; }
      for (int k = 0; k < field_.NumNeighbours(); k++) { own_nbr_[cl_nbr_[k]][k] = 1; }
      lts_ = true;
   }
   bool Lts() const { return lts_; }

   /// Entry state of the ghosts (after Reset on every rank): halo(k) points
   /// at a copy of face-neighbour element k's displacement DOFs, stride s
   /// (DGRheologyOperators::NeighbourStrain).
   void ResetGhosts(const std::function<const real_t *(int)> &halo, int s)
   {
      MFEM_VERIFY(lts_, "RheologyCoupler::ResetGhosts: SetLts first");
      const int nr = field_.NodesPerElement();
      for (int k = 0; k < field_.NumNeighbours(); k++)
      {
         ops_.NeighbourStrain(k, halo(k), s, eps_.data());
         for (int i = 0; i < nr; i++)
         {
            const RheologyPointContext &ctx = field_.NeighbourContext(k, i);
            stack_->InitState(ctx, field_.NeighbourState(k, i));
            stack_->SetEntryStrain(ctx, &eps_[6 * i], field_.NeighbourState(k, i));
            for (int c = 0; c < 6; c++) { field_.NeighbourDsigRef(k, i, c) = 0.0; }
         }
         field_.SetNeighbourActive(k, !field_.NeighbourExcluded(k));
      }
      SetLinearFlags_();   // collective; no dsig exchange (the ghosts hold it)
   }

   /// One tick (clusters 0..top end at t): the memory of every local element
   /// of those clusters and the ghost of every face-neighbour element of
   /// those clusters steps with dt(c).  U: the local displacement at t;
   /// halo(k): face-neighbour element k's displacement at t (received at
   /// this tick, stride s).  Not collective.
   void AdvanceLts(int top, const Vector &U, const std::function<const real_t *(int)> &halo,
                   int s, real_t t, const std::function<real_t(int)> &dt)
   {
      MFEM_VERIFY(lts_, "RheologyCoupler::AdvanceLts: SetLts first");
      const int nr = field_.NodesPerElement();
      long long nyield = 0;
      for (int e = 0; e < field_.NumElements(); e++)
      {
         if (excluded_[e] || cl_[e] > top) { continue; }
         const real_t h = dt(cl_[e]);
         ops_.Strain(e, U, eps_.data());
         for (int i = 0; i < nr; i++)
         {
            RheologyPointContext ctx = field_.Context(e, i);
            ctx.t = t;
            if (stack_->Update(ctx, &eps_[6 * i], h, field_.State(e, i), field_.Dsig(e, i)))
            {
               nyield++;
            }
         }
      }
      real_t ds[6];
      for (int k = 0; k < field_.NumNeighbours(); k++)
      {
         if (field_.NeighbourExcluded(k) || cl_nbr_[k] > top) { continue; }
         const real_t h = dt(cl_nbr_[k]);
         ops_.NeighbourStrain(k, halo(k), s, eps_.data());
         for (int i = 0; i < nr; i++)
         {
            RheologyPointContext ctx = field_.NeighbourContext(k, i);
            ctx.t = t;
            stack_->Update(ctx, &eps_[6 * i], h, field_.NeighbourState(k, i), ds);
            for (int c = 0; c < 6; c++) { field_.NeighbourDsigRef(k, i, c) = ds[c]; }
         }
      }
      yielding_ = nyield;
      yielding_global_ = nyield;
   }

   /// F += the owner split of P dsig of the elements of cluster c (5.6).
   void AddLoadCluster(int c, Vector &F) const
   {
      MFEM_VERIFY(lts_, "RheologyCoupler::AddLoadCluster: SetLts first");
      if (!AnyActive() || c >= static_cast<int>(own_.size())) { return; }
      ops_.AddOwned(field_, own_[c], own_nbr_[c], F);
   }

   /// Test of the ghosts (collective): largest |ghost dsig - owner dsig|.
   real_t GhostMismatch() { return field_.GhostMismatch(); }

   const RheologyStack &Stack() const { return *stack_; }
   RheologyField &Field() { return field_; }
   const RheologyField &Field() const { return field_; }
   const DGRheologyOperators &Operators() const { return ops_; }

   /// Box R.  U = U^{n+1} (local vector of the displacement space).
   void Advance(const Vector &U, real_t t, real_t dt)
   {
      MFEM_VERIFY(!(t <= last_t_), "RheologyCoupler::Advance: t = " << t
                  << " does not follow the last step's t = " << last_t_);
      last_t_ = t;
      const int nr = field_.NodesPerElement();
      long long nyield = 0;
      for (int e = 0; e < field_.NumElements(); e++)
      {
         if (excluded_[e]) { continue; }
         ops_.Strain(e, U, eps_.data());
         for (int i = 0; i < nr; i++)
         {
            RheologyPointContext ctx = field_.Context(e, i);
            ctx.t = t;
            if (stack_->Update(ctx, &eps_[6 * i], dt, field_.State(e, i), field_.Dsig(e, i)))
            {
               nyield++;
            }
         }
      }
      yielding_ = nyield;
      Refresh();
   }

   /// Active flags from dsig, the global flag, the yield count and the
   /// neighbour exchange.  Collective.  Advance calls it; a test that sets
   /// dsig directly calls it too.
   void Refresh()
   {
      if (linear_)
      {
         if (!linear_flags_set_) { SetLinearFlags_(); }
         yielding_global_ = yielding_;
         if (active_global_ > 0) { field_.ExchangeNeighbours(); }
         return;
      }
      const int nr = field_.NodesPerElement();
      long long red[2] = {yielding_, 0};
      for (int e = 0; e < field_.NumElements(); e++)
      {
         bool on = false;
         for (int i = 0; i < nr && !on; i++)
         {
            const real_t *ds = field_.Dsig(e, i);
            for (int c = 0; c < 6; c++) { on = on || ds[c] != 0.0; }
         }
         field_.SetActive(e, on);
         red[1] += on ? 1 : 0;
      }
      MPI_Allreduce(MPI_IN_PLACE, red, 2, MPI_LONG_LONG, MPI_SUM, comm_);
      yielding_global_ = red[0];
      active_global_ = red[1];
      if (active_global_ > 0) { field_.ExchangeNeighbours(); }
   }

   /// Whether any element on any rank has a nonzero dsig (same on every rank).
   bool AnyActive() const { return active_global_ > 0; }
   /// Nodes inelastic in the last step, all ranks.
   long long YieldingNodes() const { return yielding_global_; }
   /// Elements with a nonzero dsig, all ranks.
   long long ActiveElements() const { return active_global_; }

   /// F += P dsig on the active elements and their faces (after Advance).
   void AddLoad(Vector &F) const
   {
      if (!AnyActive()) { return; }
      for (int e = 0; e < field_.NumElements(); e++)
      {
         if (field_.Active(e)) { ops_.AddVolume(e, field_, F); }
      }
      ops_.AddFaces(field_, F);
   }

   /// {dsig_h} at a fault face point (integration points of Tr set).
   void FaultStress(FaceElementTransformations &Tr, real_t s[9]) const
   {
      ops_.FaceAverage(field_, Tr, s);
   }

   /// Global largest value of the diagnostic named `name` over all nodes,
   /// or -inf when no model of the stack has it.  Collective.
   real_t MaxDiagnostic(const std::string &name) const
   {
      real_t m = -std::numeric_limits<real_t>::infinity();
      std::vector<real_t> out;
      for (int k = 0; k < stack_->Size(); k++)
      {
         const IRheologyModel &M = stack_->Model(k);
         int idx = -1;
         for (int a = 0; a < M.NumDiagnostics(); a++)
         {
            if (name == M.DiagnosticName(a)) { idx = a; }
         }
         if (idx < 0) { continue; }
         out.resize(M.NumDiagnostics());
         for (int e = 0; e < field_.NumElements(); e++)
         {
            for (int i = 0; i < field_.NodesPerElement(); i++)
            {
               M.Diagnostics(field_.Context(e, i), field_.State(e, i) + stack_->StateOffset(k),
                             out.data());
               m = std::max(m, out[idx]);
            }
         }
      }
      MPI_Allreduce(MPI_IN_PLACE, &m, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm_);
      return m;
   }

   /// Names of every diagnostic of the stack, in stack order.
   std::vector<std::string> DiagnosticNames() const
   {
      std::vector<std::string> n;
      for (int k = 0; k < stack_->Size(); k++)
      {
         const IRheologyModel &M = stack_->Model(k);
         for (int a = 0; a < M.NumDiagnostics(); a++) { n.push_back(M.DiagnosticName(a)); }
      }
      return n;
   }

   /// All diagnostics at every node: out[(node) * N + a], N = the number of
   /// names, node = e nr + i.
   void Diagnostics(std::vector<real_t> &out) const
   {
      const std::size_t N = DiagnosticNames().size();
      const int nr = field_.NodesPerElement();
      out.assign(static_cast<std::size_t>(field_.NumElements()) * nr * N, 0.0);
      std::vector<real_t> tmp;
      for (int e = 0; e < field_.NumElements(); e++)
      {
         for (int i = 0; i < nr; i++)
         {
            std::size_t a0 = field_.Node(e, i) * N;
            for (int k = 0; k < stack_->Size(); k++)
            {
               const IRheologyModel &M = stack_->Model(k);
               tmp.resize(M.NumDiagnostics());
               M.Diagnostics(field_.Context(e, i), field_.State(e, i) + stack_->StateOffset(k),
                             tmp.data());
               for (int a = 0; a < M.NumDiagnostics(); a++) { out[a0 + a] = tmp[a]; }
               a0 += M.NumDiagnostics();
            }
         }
      }
   }

private:
   /// Linear stack: every element outside the mask active.  Collective.
   void SetLinearFlags_()
   {
      long long on = 0;
      for (int e = 0; e < field_.NumElements(); e++)
      {
         field_.SetActive(e, !excluded_[e]);
         on += excluded_[e] ? 0 : 1;
      }
      MPI_Allreduce(MPI_IN_PLACE, &on, 1, MPI_LONG_LONG, MPI_SUM, comm_);
      active_global_ = on;
      linear_flags_set_ = true;
   }

   std::unique_ptr<RheologyStack> stack_;
   MPI_Comm comm_;
   RheologyField field_;
   DGRheologyOperators ops_;
   std::vector<real_t> eps_;
   real_t last_t_ = -std::numeric_limits<real_t>::infinity();
   long long yielding_ = 0, yielding_global_ = 0, active_global_ = 0;
   bool linear_ = false, linear_flags_set_ = false;
   std::vector<char> excluded_;
   bool lts_ = false;
   std::vector<int> cl_, cl_nbr_;
   std::vector<std::vector<char>> own_, own_nbr_;   ///< per cluster
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_RHEOLOGY_COUPLER_HPP
