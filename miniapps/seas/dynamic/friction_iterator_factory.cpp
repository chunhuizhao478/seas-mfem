// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Implementation of MakeFrictionIterator (Phase 2).

#include "friction_iterator_factory.hpp"
#include "friction_solver.hpp"   // FrictionSolver::V0 (R-009 defensive guard)

#include <cmath>
#include <memory>

namespace mfem
{
namespace seas
{

std::unique_ptr<IFrictionIterator> MakeFrictionIterator(
   const spatial::SpatialFrictionConfig &cfg,
   FaultFaceFlux &flux,
   const spatial::RateStatePerDOFParams *rs /*nullptr for LSW*/)
{
   // `rs` (resolved per-DOF params) is consumed by the slip-law-SRW branch
   // below — it carries the per-DOF V_w injected into the SRW iterator.  The
   // aging adapter reads only the scalar RateStateBlock from cfg.rate_state, so
   // on the aging / LSW paths `rs` is legitimately unread (no (void) needed: it
   // is referenced in the SRW branch, so there is no unused-parameter warning).

   switch (cfg.law)
   {
      case spatial::FrictionLawKind::SlipWeakening:
         return std::make_unique<LswFrictionIterator>(flux);

      case spatial::FrictionLawKind::RateState:
      {
         MFEM_VERIFY(cfg.rate_state.has_value(),
                     "MakeFrictionIterator: cfg.law == rate_state requires a "
                     "[friction.rate_state] block.");

         // R-009 defensive guard: the config V_0 and the force solve's
         // hard-coded FrictionSolver::V0 must agree, else the fault is not
         // at V_init at t=0 (same invariant as SeedEquilibriumPsi_RS).
         MFEM_VERIFY(std::abs(cfg.rate_state->V_0_default - FrictionSolver::V0)
                     <= 1e-30 + 1e-12 * FrictionSolver::V0,
                     "MakeFrictionIterator: [friction.rate_state].V_0 ("
                     << cfg.rate_state->V_0_default
                     << ") must equal FrictionSolver::V0 ("
                     << FrictionSolver::V0
                     << "); the force solve hardcodes V0.");

         // State-evolution dispatch (slip-law-SRW plan §4.5).  The V_0 guard
         // above applies to BOTH laws (the force solve hardcodes V0 either way).
         switch (cfg.rate_state->state_evolution)
         {
            case spatial::StateEvolutionKind::Aging:
               return std::make_unique<RateStateAgingFrictionIterator>(
                  flux, *cfg.rate_state);

            case spatial::StateEvolutionKind::SlipSRW:
            {
               // The SRW iterator needs the resolved per-DOF V_w (ResolveRateState
               // fills rs->V_w only when state_evolution==SlipSRW).
               MFEM_VERIFY(rs != nullptr,
                           "MakeFrictionIterator: state_evolution='slip_srw' "
                           "requires resolved per-DOF params (rs); got nullptr.");
               // R-001: do NOT require Size() > 0 — ranks that own no fault DOFs
               // legitimately have N == 0 (rs.a.Size() == rs.V_w.Size() == 0), and
               // this factory is called on EVERY rank (the iterator is needed for
               // SetSubSteps even where Advance is later skipped).  Requiring > 0
               // would MFEM_VERIFY-abort the whole MPI job on far-field ranks.  The
               // "==" check still catches a SlipSRW config whose V_w was left
               // unfilled (Size()==0 while a.Size()==N>0).
               MFEM_VERIFY(rs->V_w.Size() == rs->a.Size(),
                           "MakeFrictionIterator: rs->V_w (size " << rs->V_w.Size()
                           << ") must equal the per-DOF count (rs->a size "
                           << rs->a.Size() << ") for slip_srw; ResolveRateState "
                           "fills V_w to N (possibly 0) when state_evolution==SlipSRW.");
               auto iter = std::make_unique<SlipLawSRWFrictionIterator>(
                  flux, *cfg.rate_state);
               iter->SetVw(std::vector<real_t>(
                  rs->V_w.GetData(), rs->V_w.GetData() + rs->V_w.Size()));
               return iter;
            }
         }
         MFEM_ABORT("MakeFrictionIterator: unhandled state_evolution = "
                    << static_cast<int>(cfg.rate_state->state_evolution));
         return nullptr;
      }
   }

   MFEM_ABORT("MakeFrictionIterator: unhandled cfg.law = "
              << static_cast<int>(cfg.law));
   return nullptr;
}

} // namespace seas
} // namespace mfem
