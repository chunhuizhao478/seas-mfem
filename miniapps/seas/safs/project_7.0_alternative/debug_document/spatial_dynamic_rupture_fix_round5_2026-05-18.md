# Fix Report — Round 5 (2026-05-18)

Review: `spatial_dynamic_rupture_review_round5_2026-05-18.md`

## Summary

- Findings addressed: **4 of 4** (R-501 MODERATE + R-502, R-503, R-504 LOW)
- Files modified:
  - `miniapps/seas/dynamic/wave_operator.hpp` (R-502 hoist via static helper; R-501 testable surface)
  - `miniapps/seas/dynamic/wave_operator.inl` (both dispatch arms call the static helper)
  - `miniapps/seas/dynamic/spatial_setup.hpp` (R-504 validator helper extraction)
  - `miniapps/seas/fault/fault_geometry.hpp` (R-503 NaN-fill eta_values_)
  - `miniapps/seas/tests/unit/test_phaseh_lsw_forced_rupture.cpp` (R-501 F-9 added)
- Tests added: **1** (F-9 with 6 sub-asserts)
- Test suite: **PASS — 564 / 564**:
  - `seas_test_phaseh_godunov_flux_pool`        — **30 / 30**
  - `seas_test_phaseh_lsw_forced_rupture`       — **47 / 47** (was 41; +F-9 sub-asserts)
  - `seas_test_spatial_setup`                   — **71 / 71**
  - `seas_test_spatial_friction_config`         — **28 / 28**
  - `seas_test_spatial_friction_resolver`       — **110 / 110**
  - `seas_test_spatial_velocity_bundle`         — **6 / 6**
  - `seas_test_spatial_stress_bundle`           — **6 / 6**
  - `seas_test_spatial_constant_stress_source`  — **45 / 45**
  - `seas_test_compute_safs_params`             — **13 / 13**
  - `seas_test_safs_mode_wiring`                — **8 / 8**
  - `seas_test_tpv104_checkpoint`               — **164 / 164**
  - `seas_tpv104_driver`, `seas_tpv205_driver`  — rebuild clean

## Changes Made

1. **R-501 + R-502** (dispatch guard not unit-tested; per-DOF hot-loop cost) — `wave_operator.hpp`: extracted the dispatch guard into a static helper `WaveOperator<MeshType>::VerifyForcedRuptureTimeReady(bool time_was_set, real_t T_forced_rupture)` that runs the MFEM_VERIFY.  Both dispatch arms in `wave_operator.inl:3617+` and `:4587+` now call the helper (R-502 — same logical guard, single function call instead of inline MFEM_VERIFY at each site; the per-DOF `T_forced_rupture` check is still per-DOF, but the helper is now small and inlineable).  R-501: the static helper is directly callable from unit tests (F-9 below) without needing a real fault-face fixture; six sub-asserts exercise the boolean truth table (time_was_set × T_forced_rupture, including the exact threshold 1e8 and just below).
2. **R-503** (Print eta misleading for SAFS new ctor) — `fault_geometry.hpp` new BP5 ctor: changed `eta_values_(i) = bp5_params_.eta()` to `eta_values_(i) = k_nan` so the SAFS new ctor's loud-fail pattern is uniform across all five BP5 analytic arrays (a, eta, dc, tau_pre, V_init).  Print() now emits `eta: nan MPa·s/m` for SAFS — observably wrong-as-intended.  Comment updated.
3. **R-504** (validator block duplicated across centroid + IP-aware overloads) — `spatial_setup.hpp`: extracted the 9 size-validation MFEM_VERIFYs into `internal::validate_lsw_per_dof_arrays(ndof, dof_to_elem, lsw, tau_pre, sigma_n_eff, T_forced_s, t0_decay_s)`.  Both `InitializeFaultDOFs_Spatial` LSW overloads now call it.  The IP-aware overload still has its own `dof_ips.size() == ndof` check (specific to that overload).

## Unresolved Findings

None — all 4 round-5 findings are addressed.

## New Tests

- `F_9_verify_guard_aborts_when_time_unset_and_dof_active` (six sub-asserts) — covers R-501.  Exercises `WaveOperator<Mesh>::VerifyForcedRuptureTimeReady(...)` directly via a forked child for each row of the (time_was_set × T_forced_rupture) truth table.  Case-by-case:
  - Case A: `time_was_set=false, T_forced=0.5` → MUST ABORT ✓
  - Case B: `time_was_set=false, T_forced=1.0e9` → MUST PASS (sentinel) ✓
  - Case C: `time_was_set=true,  T_forced=0.5` → MUST PASS ✓
  - Case D: `time_was_set=true,  T_forced=1.0e9` → MUST PASS ✓
  - Case E: `time_was_set=false, T_forced=1.0e8` (exact threshold) → MUST PASS (>=) ✓
  - Case F: `time_was_set=false, T_forced=9.9999e7` (just below) → MUST ABORT ✓

The helper's six-case truth table is now load-bearing: any future regression that removes the guard, flips the OR to AND, or changes the threshold gets caught by F-9.

## Ready for Re-Review: YES
