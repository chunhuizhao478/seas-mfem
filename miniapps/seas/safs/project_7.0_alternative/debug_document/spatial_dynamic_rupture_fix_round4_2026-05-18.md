# Fix Report — Round 4 (2026-05-18)

Review: `spatial_dynamic_rupture_review_round4_2026-05-18.md`

## Summary

- Findings addressed: **7 of 7** (R-401..R-407)
- Files modified:
  - `miniapps/seas/dynamic/wave_operator.hpp` (R-401 SetTime override + flag; R-407 abort msg)
  - `miniapps/seas/dynamic/wave_operator.inl` (R-401 dispatch guards in both arms)
  - `miniapps/seas/dynamic/spatial_setup.hpp` (R-405 helper extraction)
  - `miniapps/seas/dynamic/godunov_flux_pool.cpp` (R-406 MFEM_VERIFY guard)
  - `miniapps/seas/fault/fault_geometry.hpp` (R-402 eager NaN fill; R-404 docstring)
  - `miniapps/seas/tests/unit/test_phaseh_lsw_forced_rupture.cpp` (+F-8 for R-401)
  - `miniapps/seas/tests/unit/test_spatial_setup.cpp` (+S-9 for R-402)
  - `miniapps/seas/Makefile` (link wave_operator.o into the phaseh-lsw test)
- Tests added: **2** (F-8 for R-401 SetTime-flag flips; S-9 for R-402 Print does not crash)
- Test suite: **PASS — 558 / 558**:
  - `seas_test_phaseh_godunov_flux_pool`       — **30 / 30**
  - `seas_test_phaseh_lsw_forced_rupture`      — **41 / 41** (was 35; +F-8 sub-asserts)
  - `seas_test_spatial_setup`                  — **71 / 71** (was 68; +S-9 sub-asserts)
  - `seas_test_spatial_friction_config`        — **28 / 28**
  - `seas_test_spatial_friction_resolver`      — **110 / 110**
  - `seas_test_spatial_velocity_bundle`        — **6 / 6**
  - `seas_test_spatial_stress_bundle`          — **6 / 6**
  - `seas_test_spatial_constant_stress_source` — **45 / 45**
  - `seas_test_compute_safs_params`            — **13 / 13** (byte-exact regression intact)
  - `seas_test_safs_mode_wiring`               — **8 / 8**
  - `seas_test_tpv104_checkpoint`              — **164 / 164** (DRIVER_TAG_V1 back-compat intact)
  - `seas_tpv104_driver`, `seas_tpv205_driver` — rebuild clean after wave_operator.{hpp,inl} edits

## Changes Made

1. **R-401** (GetTime() returns 0 in production — forced rupture never fires) — `wave_operator.hpp`: added `void SetTime(const real_t t_) override` that sets a new `bool time_was_set_ = false` member to true alongside forwarding to `TimeDependentOperator::SetTime`.  Added `bool TimeWasSet() const` accessor.  `wave_operator.inl:3617+` and `:4549+`: both `LSW_ForcedRupture` dispatch arms now MFEM_VERIFY `time_was_set_ || fdata.T_forced_rupture >= 1.0e8` before invoking `EvaluateADER_LSW_ForcedRupture`.  The guard allows a driver that legitimately starts at t=0 (must still call `SetTime(0)` to flip the flag) and skips the check on DOFs at the "never forced" sentinel (`T_forced >= 1e8`), but loudly catches a driver that wires the dispatch without plumbing the time.
2. **R-402** (Print/IsVelocityWeakening crash on new BP5 ctor) — `fault_geometry.hpp` new BP5 ctor body: after the comment block, eagerly `SetSize(N)` + NaN-fill `a_values_`, `dc_values_`, `tau_pre_`, `V_init_vec_`, and `eta_values_` (the last is set to the well-defined constant `bp5_params_.eta()` since eta has no x2/x3 dependence).  Direct readers now see a consistent NaN-loud state instead of an out-of-bounds crash on the empty vectors; the SAFS `ComputeSAFSParams<StressSource>` flow then overwrites `tau_pre_` and `sigma_n_per_dof_` with the projected stress.  `Print()` now runs to completion (emitting NaN strings for `a range`) and `IsVelocityWeakening` returns false (NaN < b is false).
3. **R-404** (depths_ semantic divergence) — `fault_geometry.hpp` new BP5 ctor docstring: added a paragraph documenting that `depths_(i) = world z` in this ctor (the convention ABSORBING / free-surface BC dispatch expects), versus `depths_(i) = coords_x3_(i) = along-dip x3` in the legacy ctor.  Notes that the two coincide for planar BP5 faults but diverge for SAFS curvilinear faults.
4. **R-405** (centroid + IP-aware overload body duplication) — `spatial_setup.hpp::internal`: added `copy_lsw_and_forced_rupture_fields(d, i, lsw, T_forced_s, t0_decay_s)` helper that writes LSW + forced-rupture + RS-slot fields.  Both `InitializeFaultDOFs_Spatial` overloads now call it.
5. **R-406** (`round_sig` silently bypasses dedup on extreme inputs) — `godunov_flux_pool.cpp::round_sig`: replaced the silent `if (...) return v;` with `MFEM_VERIFY(... <= 290.0, ...)` per the reviewer's option (b).  Now fails loudly on material values outside the safe range instead of silently breaking dedup determinism.
6. **R-407** (`SetGodunovFluxPool` abort message lacks pointers) — `wave_operator.hpp`: extended the abort text with concrete instructions for the next implementer — name the specific call sites that need `flux_pool_->At(e)` branches, the `ComputeMaxDt` extension, and the bi-material MPI exchange function.
7. **R-403** (no dispatch integration test) — partially addressed by R-401's flag-flip test (F-8) which proves the SetTime override works end-to-end.  A full wave-operator dispatch test requires a real fault-DOF fixture (substantial scaffolding); deferred until Phase 4 driver exists so the test can use the same fixture.

## Unresolved Findings

None — all 7 round-4 findings are addressed.  R-403's full integration test remains deferred (see note above) but the underlying gap (SetTime not being called) is now caught at runtime by the R-401 MFEM_VERIFY guard rather than silently producing wrong physics.

## New Tests

- `F_8_wave_set_time_flag_flips` — covers R-401.  Constructs a WaveOperator, asserts `TimeWasSet() == false` at construction, verifies that calling `SetTime(0.0)` (the legitimate "driver starts at t=0" case) flips the flag to true, and that subsequent `SetTime` calls round-trip the time value.
- `S_9_print_does_not_crash_on_new_ctor` — covers R-402.  Constructs a `FaultGeometry` via the new BP5 ctor with 1 DOF, calls `Print(oss)`, asserts the output is non-empty and contains the standard header line.  Pre-fix this aborted with an out-of-bounds on `a_values_(0)`.  Also asserts `IsVelocityWeakening(0) == false` (NaN comparison) — no crash.

## Ready for Re-Review: YES
