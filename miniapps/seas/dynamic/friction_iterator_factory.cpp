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
   // `rs` (resolved per-DOF params) is consumed only by the SRW branch below
   // (for the per-DOF V_w side-channel, rs->V_w); the LSW + aging paths build
   // from the scalar config and ignore it.

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

         // Phase 6 req 4: dispatch on the state-evolution selector.
         if (cfg.rate_state->state_evolution
             == spatial::StateEvolutionKind::SlipLawStrongRateWeakening)
         {
            // SRW (TPV104).  The law carries the scalar globals (b, V0, f0,
            // muW=f_w_default); the per-DOF weakening velocity comes from the
            // resolver's rs->V_w (an mfem::Vector — matches the policy Extra).
            // The scalar `a` in the law is unused (the policy reads per-DOF
            // d.a); pass a_default as a valid placeholder.
            MFEM_VERIFY(rs != nullptr,
                        "MakeFrictionIterator: state_evolution=slip_law_strong_"
                        "rate_weakening requires resolved per-DOF params (rs) "
                        "for the V_w side-channel.");
            // V_w must be resolved to the SAME per-DOF fault count as the other
            // rate-state arrays (ResolveRateState SetSize(N)s them together).
            // Do NOT require Size() > 0: at np>1 the fault is partitioned across
            // ranks, so a rank that owns no local fault face has N==0 and a
            // legitimately empty V_w (== empty rs->a).  The SRW iterator already
            // tolerates this — RateStateSlipLawSrwPolicy::ValidateExtra checks
            // Size()==dof_data.size() (0==0) and the empty QP loop never indexes
            // (*V_w)(i).  A `> 0` guard here aborted TPV104 (SRW) on every
            // fault-less rank at np=800 while aging (TPV102, no guard) ran fine.
            MFEM_VERIFY(rs->V_w.Size() == rs->a.Size(),
                        "MakeFrictionIterator: SRW per-DOF V_w (size "
                        << rs->V_w.Size() << ") must be resolved to the same "
                        "fault-DOF count as the other rate-state arrays (rs->a "
                        "size " << rs->a.Size() << "); thread V_w through "
                        "ResolveRateState.  (An empty V_w is valid only when "
                        "rs->a is also empty — a rank that owns no fault DOFs.)");
            // REVIEW R-007: put the law in production mode before wiring it
            // into the iterator.  The R-001 guard in SlipLawSRWPsi aborts if a
            // base virtual is invoked with the scalar V_w_default_/a_ instead
            // of per-QP values; without SetProductionMode() that guard is dead,
            // so a future call to law_.Rate(...)/SteadyState(...) would
            // silently use the scalar (VW-core) V_w for every QP.
            SlipLawSRWPsi srw(cfg.rate_state->a_default,
                              cfg.rate_state->b_default,
                              cfg.rate_state->V_0_default,
                              cfg.rate_state->f_0_default,
                              cfg.rate_state->f_w_default,
                              cfg.rate_state->V_w_default);
            srw.SetProductionMode();
            return std::make_unique<RateStateSlipLawSrwIterator>(
               flux, std::move(srw), FrictionSolver::Method::Brent, &rs->V_w);
         }

         // Aging (default; TPV102 / SAFS).  The unified RateStateAgingIterator
         // owns its AgingLawPsi by value (b, V0, f0 from the scalar block) and
         // uses the production RS friction method Brent (CLAUDE.md); reproduces
         // the standalone Tpv102SubStepIterator bit-for-bit (parity test).
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
