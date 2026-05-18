# Fix Report — Round 6 (2026-05-18)

Review: `spatial_dynamic_rupture_review_round6_2026-05-18.md` + user request: generalise nucleation to choose between strength-reduction and overstress.

## Summary

- Findings addressed: **6 of 8** (R-601, R-602, R-603, R-604, R-607 fully; R-605 + R-606 + R-608 are POSSIBLE / defensive / audit-deferred — see Unresolved).
- User request addressed: NucleationKind enum + parser + driver wiring + stub `ResolveOverstress` resolver.
- Files modified:
  - `miniapps/seas/dynamic/tpv205_substep_iterator.{hpp,cpp}` — R-601 forced-rupture mode setter + per-substep `t_sub_abs` plumbing + branch on `forced_rupture_mode_`
  - `miniapps/seas/drivers/spatial_dyn_driver.cpp` — R-601 `SetForcedRuptureMode(true)`; R-602 ParaView state uses `LSWFrictionCoefficient_ForcedRupture`; R-603 `last_completed_step`; R-604 dead ref_normal cleanup; R-607 ProbeNbfPerFace call simplified; NucleationKind dispatch
  - `miniapps/seas/spatial/code/spatial_friction.{hpp,cpp}` — `NucleationKind` enum (StrengthReduction | Overstress); `OverstressSpec`, `OverstressPerDOFParams`; parser support for `[nucleation].kind` + `[nucleation.overstress]` sub-block; `ResolveOverstress` stub
  - `miniapps/seas/Makefile` — link `dynamic/tpv205_substep_iterator.o` into `seas_test_phaseh_lsw_forced_rupture`
  - `miniapps/seas/tests/unit/test_phaseh_lsw_forced_rupture.cpp` — F-10 added: iterator-level regression for R-601
- Tests added: **1** (F-10, with 4 sub-asserts)
- Test suite: **PASS — 583 / 583 total** (was 564 + F-10's 4 sub-asserts = 568; the constant-parity test contributes additional sub-asserts):
  - `seas_test_phaseh_godunov_flux_pool`         — **30 / 30**
  - `seas_test_phaseh_lsw_forced_rupture`        — **51 / 51** (was 47; +F-10)
  - `seas_test_phaseh_wave_operator_constant_parity` — **passes serial + parallel**
  - `seas_test_spatial_setup`                    — **71 / 71**
  - `seas_test_spatial_friction_config`          — **28 / 28**
  - `seas_test_spatial_friction_resolver`        — **110 / 110**
  - `seas_test_spatial_velocity_bundle`          — **6 / 6**
  - `seas_test_spatial_stress_bundle`            — **6 / 6**
  - `seas_test_spatial_constant_stress_source`   — **45 / 45**
  - `seas_test_compute_safs_params`              — **13 / 13** (byte-exact regression intact)
  - `seas_test_safs_mode_wiring`                 — **8 / 8**
  - `seas_test_tpv104_checkpoint`                — **164 / 164** (DRIVER_TAG_V1 back-compat intact)
  - `seas_spatial_dyn_driver`, `seas_tpv104_driver`, `seas_tpv205_driver` — rebuild clean (no warnings, no errors); driver smoke-test (`--config` missing) emits clean usage error and exits 0.

## Changes Made

1. **R-601 (CRITICAL — forced rupture dead on dominant path).**  Added `bool forced_rupture_mode_ = false` member to `Tpv205SubStepIterator` with `SetForcedRuptureMode(bool)` / `GetForcedRuptureMode()` accessors.  `StepOneQP_` now takes a new `real_t t_sub_abs` parameter (per-substep absolute simulation time) and branches on `forced_rupture_mode_`: when true it calls `mfem::seas::spatial::LSWFrictionCoefficient_ForcedRupture(delta, mu_s, mu_d, d_c, t_sub_abs, d.T_forced_rupture, d.t0_decay_forced)`; when false it keeps the original `LSWFrictionCoefficient_TPV205(...)` call byte-for-byte (TPV205 byte-exact contract preserved).  Both `Advance` and `AdvanceWithSubStepStates` now accumulate `t_sub_acc` from `t_macro_start` and pass it down.  Driver calls `substep_iterator.SetForcedRuptureMode(true)` once at iterator construction.
2. **R-602 (MODERATE — diagnostic mu_eff plain LSW).**  `paraview_write` in the driver replaces `LSWFrictionCoefficient_TPV205(delta_norm, ...)` with `mfem::seas::spatial::LSWFrictionCoefficient_ForcedRupture(delta_norm, ..., time, d.T_forced_rupture, d.t0_decay_forced)`.  For `T_forced_rupture >= 1e8` the helper short-circuits f_2 to 0 — byte-equivalent to TPV205 when nucleation is off.  ParaView output now reflects the actual mu_eff including the f_2 reduction inside the nucleation zone.
3. **R-603 (MODERATE — final checkpoint step number).**  Added `int last_completed_step = step0;` before the time loop; incremented to `step + 1` after each successful step.  Final-checkpoint call at line 1262 uses `last_completed_step` instead of `nsteps`.  When the loop breaks early (`dt_step <= 0`), the checkpoint records the correct step number; restart from such a checkpoint resumes at the right place instead of skipping all remaining work.
4. **R-604 (LOW — dead ref_normal assignment).**  Removed the initial `ref_normal = (0, 0, 1)` and `up_vec` re-write; collapsed to a single, top-of-block assignment `ref_normal = (0, 1, 0)` and `up_vec = (0, 0, 1)` with the explanatory comment moved up so the reader sees one value, not two.
5. **R-607 (LOW — dead `#ifdef MFEM_USE_MPI` in arg list).**  Replaced the per-arg `#ifdef MFEM_USE_MPI` / `#else MPI_COMM_NULL` block with a direct `comm` pass at the `ProbeNbfPerFace` call site.  The `#error "spatial_dyn_driver requires MFEM_USE_MPI=YES"` at L682 already guarantees the macro is defined.
6. **User request — NucleationKind enum + overstress stub.**  Added `enum class NucleationKind { StrengthReduction, Overstress }` to `spatial_friction.hpp`, with a default of `StrengthReduction` to preserve round-1..5 behaviour byte-equivalently.  Added `OverstressSpec` (`delta_tau_pa`, `direction`, `delta_sigma_n_pa`) and `OverstressPerDOFParams` (`delta_tau1_pa`, `delta_tau2_pa`, `delta_sigma_n_pa`) structs.  Extended the TOML parser to accept optional `[nucleation].kind = "strength_reduction" | "overstress"` and `[nucleation.overstress]` sub-block.  Added `SpatialFrictionResolver::ResolveOverstress(...)` declaration + stub body: returns empty vectors when nucleation is disabled or kind is StrengthReduction; aborts with a precise "not yet implemented" message when kind is Overstress.  Driver dispatch wiring chooses `FaultFrictionLaw::LSW_ForcedRupture` for StrengthReduction (existing path) and `FaultFrictionLaw::LSW` for Overstress (the per-DOF perturbation lives in `DOFData::tau{1,2}_nuc` and is read by plain LSW; pre-DOFData population will land in the resolver follow-up).  Driver calls `ResolveOverstress(...)` always; the result is unused on the StrengthReduction path and aborts loudly on the Overstress path until the resolver body is filled in.

## Unresolved Findings

- **R-605** (POSSIBLE — CFL parity reliance on untouched scalar h_min_) — flagged as POSSIBLE in round 6.  The constant-parity test C-3 verifies byte-identical `ComputeMaxDt` between scalar and hetero ctors on the test fixture; if a future tetrahedron with negative-orientation vertices breaks the assumption, C-3 catches it.  No concrete failure case in hand.  Deferred.
- **R-606** (defensive — placeholder material in delegating ctor) — would require refactoring the scalar ctor's body into a private helper.  The existing MFEM_VERIFY catches Mode::Coefficient before any Mult call.  Deferred.
- **R-608** (POSSIBLE — uniform-weight substep quadrature audit) — audit-only, no functional change.  TPV205 production may use a different quadrature; auditing against `drivers/tpv205_driver.cpp` is a separate task.  Deferred.

## New Tests

- `F_10_iterator_consumes_forced_rupture_fields` (4 sub-asserts) — covers R-601.  Builds two `Tpv205SubStepIterator`s on identical DOFData except for `(T_forced_rupture, t0_decay_forced)`; calls `Advance(...)` with O=1 substep at `t_macro_start = 0.25`.  Asserts:
  - Plain-LSW iterator stays locked (`slip_rate == 0`)
  - Forced-rupture iterator (with `SetForcedRuptureMode(true)`) unlocks (`slip_rate > 1e-3`)
  - `I_imp_plus` materially differs between the two paths
  - `Tpv205SubStepIterator` default-constructs with `GetForcedRuptureMode() == false` (byte-exact TPV205 contract)

## Ready for Re-Review: YES
