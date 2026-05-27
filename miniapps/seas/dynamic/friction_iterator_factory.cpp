// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Implementation of MakeFrictionIterator (Phase 2).

#include "friction_iterator_factory.hpp"
#include "friction_substep_iterator.hpp"   // Phase 5: unified iterators
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
   // `rs` is reserved for the Phase-3 driver wiring and the future
   // slip-law-SRW path; the aging adapter reads only the scalar
   // RateStateBlock from cfg.rate_state, so rs is unused here today.
   (void) rs;

   switch (cfg.law)
   {
      case spatial::FrictionLawKind::SlipWeakening:
         // Phase 5: the unified LinearSlipWeakeningIterator (reproduces the
         // standalone Tpv205SubStepIterator bit-for-bit — see
         // test_friction_substep_iterator_parity).
         return std::make_unique<LinearSlipWeakeningIterator>(flux);

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

         // Phase 5: the unified RateStateAgingIterator owns its AgingLawPsi
         // by value (b, V0, f0 from the scalar RateStateBlock) and uses the
         // production RS friction method Brent (CLAUDE.md).  Reproduces the
         // standalone Tpv102SubStepIterator bit-for-bit (parity test).
         //
         // Aging remains the ONLY RS law dispatched here: the slip-law
         // strong-rate-weakening (SRW) iterator exists (Phase 5,
         // RateStateSlipLawSrwIterator) but its per-DOF V_w side-channel and
         // the [friction.rate_state] state_evolution selector land in
         // Phase 6.  When that field exists, branch on it here to return a
         // RateStateSlipLawSrwIterator (constructed with &rs->V_w); until
         // then there is no SRW config to dispatch on.
         return std::make_unique<RateStateAgingIterator>(
            flux,
            AgingLawPsi(cfg.rate_state->b_default,
                        cfg.rate_state->V_0_default,
                        cfg.rate_state->f_0_default),
            FrictionSolver::Method::Brent);
      }
   }

   MFEM_ABORT("MakeFrictionIterator: unhandled cfg.law = "
              << static_cast<int>(cfg.law));
   return nullptr;
}

} // namespace seas
} // namespace mfem
