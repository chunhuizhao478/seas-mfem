// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/friction_iterator_factory.hpp — Phase 2 of
// PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.
//
// One-line factory that selects the per-sub-step friction iterator from a
// SpatialFrictionConfig.  Adding a friction law is a new IFrictionIterator
// adapter + one case here — never a driver edit.

#ifndef MFEM_SEAS_FRICTION_ITERATOR_FACTORY_HPP
#define MFEM_SEAS_FRICTION_ITERATOR_FACTORY_HPP

#include "friction_iterator.hpp"
#include "fault_face_flux.hpp"
#include "../spatial/code/spatial_friction.hpp"   // SpatialFrictionConfig

#include <memory>

namespace mfem
{
namespace seas
{

/// @brief Build the friction sub-step iterator selected by `cfg.law`.
///
/// Returns the unified Phase-5 iterators (dynamic/friction_substep_iterator.hpp),
/// which honour `SetFaultResample` (the deprecated `*FrictionIterator` adapters
/// in friction_iterator.hpp do NOT — do not wire those on the production path):
///   cfg.law == SlipWeakening -> LinearSlipWeakeningIterator (WaveOpLaw()==LSW)
///   cfg.law == RateState     -> RateStateSlipLawSrwIterator  (slip-law SRW), or
///                               RateStateAgingIterator       (aging; default)
///                               (WaveOpLaw()==RateAndState)
///
/// Aborts (MFEM_VERIFY) on the RateState branch if `cfg.rate_state` is
/// unset, or if its `V_0_default` disagrees with the force solve's
/// hard-coded `FrictionSolver::V0` (R-009 — same invariant as the
/// Phase-1 SeedEquilibriumPsi_RS guard, asserted up front so a mis-set
/// V_0 is caught before the iterator is built).
///
/// @param cfg  Parsed spatial-friction config (selects the law).
/// @param flux Fault-face Riemann solver bound to the owned iterator.
/// @param rs   Per-DOF RS params (nullptr for LSW).  Reserved for the
///             Phase-3 driver wiring / the future slip-law-SRW path; the
///             aging adapter reads only the scalar RateStateBlock today.
std::unique_ptr<IFrictionIterator> MakeFrictionIterator(
   const spatial::SpatialFrictionConfig &cfg,
   FaultFaceFlux &flux,
   const spatial::RateStatePerDOFParams *rs /*nullptr for LSW*/);

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FRICTION_ITERATOR_FACTORY_HPP
