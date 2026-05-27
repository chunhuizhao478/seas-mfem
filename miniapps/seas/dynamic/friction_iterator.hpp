// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/friction_iterator.hpp — Phase 2 of
// PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.
//
// A minimal runtime-dispatch strategy interface for the SAFS spatial
// driver's per-sub-step friction iterator, with thin adapters over the
// existing, proven standalone iterators:
//
//   IFrictionIterator               <- strategy interface
//     LswFrictionIterator           -> owns a Tpv205SubStepIterator (LSW)
//     RateStateAgingFrictionIterator-> owns a Tpv102SubStepIterator (aging RS)
//
// The driver (`AdvanceADERWithSubStep_Spatial`) takes an
// `IFrictionIterator&` and calls only the interface methods; the LSW path
// is a transparent forwarder, so an existing SAFS LSW run is byte-
// identical to the pre-Phase-2 direct-Tpv205 path.  The full method-
// oriented unification (SubStepIteratorBase + SRW) is Phase 5; the
// standalone tpv{205,102}_substep_iterator stay the byte-exact oracle.
//
// `FaultFrictionLaw` (the WaveOpLaw() return) lives in wave_operator.hpp;
// the driver and the factory already include it, so pulling it here is no
// new dependency on the wired path.  wave_operator.inl is purely template
// member functions, so a translation unit that includes this header but
// never instantiates WaveOperator<...> does not link wave_operator.o.

#ifndef MFEM_SEAS_FRICTION_ITERATOR_HPP
#define MFEM_SEAS_FRICTION_ITERATOR_HPP

#include "mfem.hpp"

#include "fault_face_flux.hpp"
#include "friction_solver.hpp"
#include "wave_operator.hpp"                      // FaultFrictionLaw
#include "tpv205_substep_iterator.hpp"
#include "tpv102_substep_iterator.hpp"
#include "../friction/state_evolution.hpp"        // AgingLawPsi
#include "../spatial/code/spatial_friction.hpp"   // spatial::RateStateBlock

#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Per-sub-step friction-iterator strategy used by the SAFS
/// spatial dynamic-rupture driver.  Implementations adapt the existing
/// benchmark sub-step iterators; the driver depends only on this
/// interface so adding a friction law is a new adapter + one factory
/// line, never a driver edit.
class IFrictionIterator
{
public:
   /// Configure the ADER sub-step quadrature (forwarded to the owned
   /// iterator; same constraints — see Tpv*SubStepIterator::SetSubSteps).
   virtual void SetSubSteps(std::vector<real_t> deltaT,
                            std::vector<real_t> time_weights) = 0;

   /// R-010: `AdvanceADERWithSubStep_Spatial`'s body also reads these two
   /// (driver scales deltaT and re-reads it), so they MUST be on the
   /// interface or the parameter retype will not compile.
   virtual const std::vector<real_t> &GetDeltaT()      const = 0;
   virtual const std::vector<real_t> &GetTimeWeights() const = 0;

   /// Advance every fault QP through one macro-step.  `nuc_callback`
   /// fires once per ADER sub-step before the per-QP friction pipeline
   /// (the Phase-1 callback overload on the owned iterator); pass a no-op
   /// `[](real_t,real_t){}` to opt out.
   virtual void Advance(
      std::vector<DOFData> &dof_data,
      const std::vector<Vector> &fault_coords,
      const std::vector<std::vector<real_t>> &Q_pointwise_plus,
      const std::vector<std::vector<real_t>> &Q_pointwise_minus,
      real_t dt_macro,
      real_t t_macro_start,
      real_t *I_imp_plus_flat,
      real_t *I_imp_minus_flat,
      const std::function<void(real_t, real_t)> &nuc_callback) = 0;

   /// Diagnostic [SLIP] is_shared hook (forwarded; no-op for TPV102).
   virtual void SetDiagNumLocalFaultQPs(int n) = 0;

   /// The wave-operator fault friction law this iterator drives.  The
   /// driver passes it to `wave.SetFaultFrictionLaw(...)`.
   virtual FaultFrictionLaw WaveOpLaw() const = 0;

   virtual ~IFrictionIterator() = default;
};

/// DEPRECATED (Phase 5): `MakeFrictionIterator` now returns the unified
/// `LinearSlipWeakeningIterator` (dynamic/friction_substep_iterator.hpp),
/// which reproduces `Tpv205SubStepIterator` bit-for-bit
/// (test_friction_substep_iterator_parity).  This adapter is retained only
/// for the Phase-2 compile/factory tests; do NOT use it in new code — it
/// still forwards to the standalone oracle, not the production unified path.
///
/// LSW adapter: a transparent forwarder over the existing
/// `Tpv205SubStepIterator`.  Behaviour is byte-identical to the pre-
/// Phase-2 driver, which constructed a `Tpv205SubStepIterator` directly.
class LswFrictionIterator : public IFrictionIterator
{
public:
   explicit LswFrictionIterator(FaultFaceFlux &flux) : it_(flux) {}

   void SetSubSteps(std::vector<real_t> deltaT,
                    std::vector<real_t> time_weights) override
   { it_.SetSubSteps(std::move(deltaT), std::move(time_weights)); }

   const std::vector<real_t> &GetDeltaT() const override
   { return it_.GetDeltaT(); }
   const std::vector<real_t> &GetTimeWeights() const override
   { return it_.GetTimeWeights(); }

   void Advance(std::vector<DOFData> &dof_data,
                const std::vector<Vector> &fault_coords,
                const std::vector<std::vector<real_t>> &Q_pointwise_plus,
                const std::vector<std::vector<real_t>> &Q_pointwise_minus,
                real_t dt_macro,
                real_t t_macro_start,
                real_t *I_imp_plus_flat,
                real_t *I_imp_minus_flat,
                const std::function<void(real_t, real_t)> &nuc_callback)
      override
   {
      // R-016: the owned Tpv205SubStepIterator callback overload calls
      // nuc_callback unconditionally (no empty-callback guard, unlike the
      // Tpv102 overload's R-013 guard), so reject an empty callback here —
      // symmetric with the RS adapter — instead of letting an opaque
      // std::bad_function_call escape from the sub-step loop.  Pass a no-op
      // [](real_t,real_t){} to opt out.  (Does not touch the byte-exact
      // Tpv205 oracle.)
      if (!nuc_callback)
      {
         throw std::runtime_error(
            "LswFrictionIterator::Advance: nuc_callback is empty; pass a "
            "no-op [](real_t,real_t){} to opt out of nucleation.");
      }
      it_.AdvanceWithSubStepStates(dof_data, fault_coords,
                                   Q_pointwise_plus, Q_pointwise_minus,
                                   dt_macro, t_macro_start,
                                   I_imp_plus_flat, I_imp_minus_flat,
                                   nuc_callback);
   }

   void SetDiagNumLocalFaultQPs(int n) override
   { it_.SetDiagNumLocalFaultQPs(n); }

   FaultFrictionLaw WaveOpLaw() const override
   { return FaultFrictionLaw::LSW; }

private:
   Tpv205SubStepIterator it_;
};

/// DEPRECATED (Phase 5): `MakeFrictionIterator` now returns the unified
/// `RateStateAgingIterator` (= RateStateSubStepIterator<RateStateAgingPolicy>,
/// dynamic/friction_substep_iterator.hpp), which reproduces
/// `Tpv102SubStepIterator` bit-for-bit.  This adapter is retained only for the
/// Phase-2 compile/factory tests; do NOT use it in new code.
///
/// Aging rate-and-state adapter over the existing `Tpv102SubStepIterator`.
/// Uses the Phase-1 nucleation-callback overload and the production RS
/// friction method `FrictionSolver::Method::Brent` (CLAUDE.md).
class RateStateAgingFrictionIterator : public IFrictionIterator
{
public:
   /// `blk` supplies the scalar aging-law globals (b, V0, f0).  R-009:
   /// `law_(b, V_0, f_0)` makes `AgingLawPsi` read the config V_0; the
   /// force solve hard-codes `FrictionSolver::V0`, and `MakeFrictionIterator`
   /// (and the Phase-1 seed) assert the two agree.
   RateStateAgingFrictionIterator(FaultFaceFlux &flux,
                                  const spatial::RateStateBlock &blk)
      : law_(blk.b_default, blk.V_0_default, blk.f_0_default),
        it_(flux, law_) {}

   // R-015: it_ holds a `const AgingLawPsi&` bound to law_; a member-wise
   // copy/move would leave the copy's reference dangling at the source's
   // law_.  This adapter is only ever held via unique_ptr (make_unique
   // constructs it in place), so forbid copy/move to make the footgun a
   // compile error rather than a silent use-after-free.
   RateStateAgingFrictionIterator(const RateStateAgingFrictionIterator&) = delete;
   RateStateAgingFrictionIterator&
      operator=(const RateStateAgingFrictionIterator&) = delete;
   RateStateAgingFrictionIterator(RateStateAgingFrictionIterator&&) = delete;
   RateStateAgingFrictionIterator&
      operator=(RateStateAgingFrictionIterator&&) = delete;

   void SetSubSteps(std::vector<real_t> deltaT,
                    std::vector<real_t> time_weights) override
   { it_.SetSubSteps(std::move(deltaT), std::move(time_weights)); }

   const std::vector<real_t> &GetDeltaT() const override
   { return it_.GetDeltaT(); }
   const std::vector<real_t> &GetTimeWeights() const override
   { return it_.GetTimeWeights(); }

   void Advance(std::vector<DOFData> &dof_data,
                const std::vector<Vector> &fault_coords,
                const std::vector<std::vector<real_t>> &Q_pointwise_plus,
                const std::vector<std::vector<real_t>> &Q_pointwise_minus,
                real_t dt_macro,
                real_t t_macro_start,
                real_t *I_imp_plus_flat,
                real_t *I_imp_minus_flat,
                const std::function<void(real_t, real_t)> &nuc_callback)
      override
   {
      it_.AdvanceWithSubStepStates(dof_data, fault_coords,
                                   Q_pointwise_plus, Q_pointwise_minus,
                                   dt_macro, t_macro_start,
                                   I_imp_plus_flat, I_imp_minus_flat,
                                   FrictionSolver::Method::Brent,
                                   nuc_callback);
   }

   void SetDiagNumLocalFaultQPs(int n) override
   { it_.SetDiagNumLocalFaultQPs(n); }

   FaultFrictionLaw WaveOpLaw() const override
   { return FaultFrictionLaw::RateAndState; }

private:
   // Member-init order matters: it_ holds a `const AgingLawPsi&`, so law_
   // MUST be declared (and thus constructed) BEFORE it_.
   AgingLawPsi           law_;
   Tpv102SubStepIterator it_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FRICTION_ITERATOR_HPP
