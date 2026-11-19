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
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_RHEOLOGY_COUPLER_HPP
#define MFEM_SEAS_RHEOLOGY_COUPLER_HPP

#include "mfem.hpp"
#include "dg_rheology_operators.hpp"
#include "rheology_field.hpp"
#include "rheology_stack.hpp"

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
   }

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
   std::unique_ptr<RheologyStack> stack_;
   MPI_Comm comm_;
   RheologyField field_;
   DGRheologyOperators ops_;
   std::vector<real_t> eps_;
   real_t last_t_ = -std::numeric_limits<real_t>::infinity();
   long long yielding_ = 0, yielding_global_ = 0, active_global_ = 0;
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_RHEOLOGY_COUPLER_HPP
