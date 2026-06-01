// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// friction/state_policies.hpp — Phase 5 of
// PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.
//
// State-evolution policies for the unified `RateStateSubStepIterator<Policy>`
// (dynamic/friction_substep_iterator.hpp).  Each policy bundles:
//   - `Law`   : the state-evolution law type held (by const ref) by the iterator,
//   - `Extra` : a per-iterator side-channel carried alongside (V_w for SRW;
//               none for aging), and
//   - `UpdatePsi(...)`: the single per-QP ψ-update that is the ONLY point of
//               variation between the aging (TPV102) and slip-law-SRW (TPV104)
//               rate-and-state iterators.
//
// The policies live next to their laws (friction/) so the friction-method
// coupling sits with the friction code, per plan §5.1.
//
// === Deviations from the plan's §5.1 snippet (recorded; verified at impl time) ===
//  (a) BOTH `UpdatePsi` policies read the PER-DOF `d.b`, NOT `Law::GetB()`.
//      Phase 11a made `b` per-DOF on the aging path (live oracle
//      `tpv102_substep_iterator.cpp:189` passes `d.b`); the 2026-06-01
//      depth-profile-SRW fix extends the SAME per-DOF `d.b` to the SRW policy
//      so a depth-varying `b` evolves ψ consistently with the equilibrium seed
//      (`SeedEquilibriumPsi_RS` already uses per-DOF `rs.b`).  For scalar TPV104
//      `d.b == b_default == L.GetB()`, so both policies stay byte-identical to
//      the standalone `tpv102_substep_iterator.cpp:189` / `tpv104_..._.cpp:416`
//      (the native-TPV104-only SRW standalone keeps the scalar `GetB()`; it is
//      always scalar `b`).
//  (b) `Extra` for SRW is `const mfem::Vector*` (non-owning), NOT
//      `std::span<const real_t>`: the toolchain is pre-C++20 and `std::span` is
//      unavailable, and an `mfem::Vector` matches the resolver's
//      `RateStatePerDOFParams.V_w` (Phase 6 req 4) so the factory passes
//      `&rs->V_w` directly with no copy.  Aging carries no side-channel
//      (`std::nullptr_t`).

#ifndef MFEM_SEAS_STATE_POLICIES_HPP
#define MFEM_SEAS_STATE_POLICIES_HPP

#include "mfem.hpp"                          // real_t

#include "../dynamic/fault_face_flux.hpp"    // DOFData
#include "state_evolution.hpp"               // AgingLawPsi, UpdateStateAnalytic
#include "slip_law_srw_psi.hpp"              // SlipLawSRWPsi, UpdateStateAnalyticSlipLawSRW

#include <cstddef>                           // std::nullptr_t
#include <vector>

namespace mfem
{
namespace seas
{

/// Aging law — TPV102 / SAFS rate-and-state.  The per-QP ψ-update is the
/// analytic aging-law step.  Phase 11a: `b` is per-DOF (`d.b`), so this
/// reproduces `tpv102_substep_iterator.cpp:186-190` bit-for-bit (the iterator
/// passes `d.b`, not `state_evo_.GetB()`).  `f0`/`V0` remain scalar globals.
struct RateStateAgingPolicy
{
   using Law   = AgingLawPsi;
   using Extra = std::nullptr_t;    ///< no side-channel on the aging path

   static real_t UpdatePsi(const Law &L, const DOFData &d, real_t V,
                           real_t dt, const Extra & /*unused*/, int /*i*/)
   {
      return UpdateStateAnalytic(d.psi, V, d.Dc, dt,
                                 L.GetF0(), d.b, L.GetV0());
   }

   /// No per-QP side-channel on the aging path — nothing to validate.
   static void ValidateExtra(const Extra & /*unused*/, int /*n*/) {}
};

/// Slip-law with strong rate weakening — TPV104 / SAFS.  The per-QP ψ-update
/// uses the per-QP weakening velocity `V_w[i]` (the `Extra` side-channel), the
/// per-QP direct-effect `d.a`, AND the per-QP state-evolution `d.b` (depth
/// profile), plus the scalar globals `V0/f0/muW` carried by the law.  The
/// per-DOF `d.b` matches the equilibrium seed (`SeedEquilibriumPsi_RS` uses
/// `rs.b`), so a depth-varying `b` seeds the fault IN equilibrium at t=0.  For
/// scalar TPV104 `d.b == b_default == L.GetB()`, so this stays byte-identical to
/// `tpv104_substep_iterator.cpp:413-419` (which uses the scalar `GetB()`).
struct RateStateSlipLawSrwPolicy
{
   using Law   = SlipLawSRWPsi;
   // V_w side-channel, non-owning.  An mfem::Vector (NOT std::vector) to match
   // the resolver's RateStatePerDOFParams.V_w (Phase 6 req 4), so the factory
   // can pass &rs->V_w directly with no copy/dangle.
   using Extra = const mfem::Vector *;

   static real_t UpdatePsi(const Law &L, const DOFData &d, real_t V,
                           real_t dt, const Extra &Vw, int i)
   {
      return UpdateStateAnalyticSlipLawSRW(d.psi, V, d.Dc, dt,
                                           (*Vw)(i), d.a,
                                           d.b, L.GetV0(),
                                           L.GetF0(), L.GetMuW());
   }

   /// The SRW path indexes `(*Vw)(i)` for every fault QP, so the V_w
   /// side-channel must be non-null and sized to the fault-QP count.
   static void ValidateExtra(const Extra &Vw, int n)
   {
      MFEM_VERIFY(Vw != nullptr,
                  "RateStateSlipLawSrwPolicy: V_w side-channel pointer is null; "
                  "the SRW iterator requires a per-QP V_w vector.");
      MFEM_VERIFY(Vw->Size() == n,
                  "RateStateSlipLawSrwPolicy: V_w size (" << Vw->Size()
                  << ") must equal the fault-QP count (" << n << ").");
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_STATE_POLICIES_HPP
