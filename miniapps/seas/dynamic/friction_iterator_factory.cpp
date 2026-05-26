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
   // `rs` is reserved for the Phase-3 driver wiring and the future
   // slip-law-SRW path; the aging adapter reads only the scalar
   // RateStateBlock from cfg.rate_state, so rs is unused here today.
   (void) rs;

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

         // Aging is the only RS law wired in Phases 1-3.  The slip-law
         // strong-rate-weakening (SRW) variant lands in Phase 5, gated on
         // a future [friction.rate_state] state_evolution field; until that
         // field exists there is no SRW config to abort on here.  When it
         // is added, branch on it and MFEM_ABORT("slip-law SRW lands in
         // Phase 5") before this return.
         return std::make_unique<RateStateAgingFrictionIterator>(
            flux, *cfg.rate_state);
      }
   }

   MFEM_ABORT("MakeFrictionIterator: unhandled cfg.law = "
              << static_cast<int>(cfg.law));
   return nullptr;
}

} // namespace seas
} // namespace mfem
