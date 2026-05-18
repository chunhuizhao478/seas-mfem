# Fix Report — Round 2 (2026-05-18)

Review: `spatial_dynamic_rupture_review_round2_2026-05-18.md`

## Summary

- Findings addressed: **8 of 8** (R-201..R-208)
- Files modified:
  - `miniapps/seas/spatial/code/spatial_friction.cpp`
  - `miniapps/seas/io/field_coefficient.hpp`
  - `miniapps/seas/tests/unit/test_spatial_friction_config.cpp`
  - `miniapps/seas/tests/unit/test_spatial_friction_resolver.cpp`
  - `miniapps/seas/tests/unit/test_spatial_stress_bundle.cpp`
- Tests added: **4** (T-17/18/19/20 for R-201/202/203/205)
- Test suite: **PASS** — all 5 spatial test binaries plus the byte-exact regression:
  - `seas_test_spatial_friction_config`        — **28 / 28** (was 24/24 — added T-17/18/19/20)
  - `seas_test_spatial_friction_resolver`      — **110 / 110**
  - `seas_test_spatial_velocity_bundle`        — **6 / 6**
  - `seas_test_spatial_stress_bundle`          — **6 / 6** (now ~0.3 s wall-time vs ~30 s pre-cache)
  - `seas_test_spatial_constant_stress_source` — **45 / 45**
  - `seas_test_compute_safs_params`            — **13 / 13** (byte-exact regression intact)

## Changes Made

1. **R-201** (`[pore_pressure]` missing-block silence) — `spatial_friction.cpp:parse_root`: added `MFEM_VERIFY(root.contains("pore_pressure"), ...)` to the new top-level mandatory-block guard list.  Test added (T-17).
2. **R-202** (dt_initial sentinel collision) — `spatial_friction.cpp`: introduced a local `parse_dt_initial(const toml::value&)` helper that returns the -1 sentinel **only** when the value is the literal string `"auto"` (or the key is absent).  Any numeric value (including `-1.0`) must be `> 0` or the helper aborts.  Replaced the call site in `parse_root` and dropped the redundant downstream validator.  Test added (T-18).
3. **R-203** (paraview_*_dt missing > 0 validators) — `spatial_friction.cpp:parse_root`: added three `MFEM_VERIFY(... > 0.0, ...)` checks for `paraview_volume_dt`, `paraview_bulk_dt`, `paraview_fault_dt`.  Test added (T-19).
4. **R-204** (R-1b misleadingly named) — `test_spatial_friction_resolver.cpp`: renamed `R_1b_rs_eta_auto_per_dof` → `R_1b_rs_eta_auto_smoke_const_material`; updated docstring to point at `R_9_eta_auto_uses_centroid` as the genuine R-003 regression.  Updated `main()` to call the new name.
5. **R-205** (silent `stress.kind` default) — `spatial_friction.cpp:parse_root`: added `MFEM_VERIFY(s.contains("kind"), ...)` before the `parse_stress_kind(...)` call.  Test added (T-20).
6. **R-206** (BuildFixture triplicated) — `test_spatial_stress_bundle.cpp`: replaced `std::unique_ptr<Fixture> BuildFixture()` with `Fixture* GetSharedFixture()` that uses a function-static `unique_ptr` to memoise the constructed fixture.  S-2/S-3/S-4 now call `GetSharedFixture()`.  Test wall-time dropped from ~30 s to ~0.33 s.
7. **R-207** (unused locals in F-3) — `test_spatial_friction_resolver.cpp:F_3_T_uses_centroid_not_corner`: removed the unused `const int N = 1;` and `Array<int> attr(1)`; updated comment to note `ResolveForcedRupture` does not take `dof_to_attr`.  Same `const int N = 1` removed from `R_9_eta_auto_uses_centroid` (attr is still used there by `ResolveRateState`).
8. **R-208** (`AbortRangeFailure` unnecessarily public) — `io/field_coefficient.hpp`: moved `AbortRangeFailure` back to `private:` with a comment explaining the R-208 reversion.  `AbortContainmentFailure` stays `public:` because `spatial_velocity.cpp` consumes it.

## Unresolved Findings

None — all 8 round-2 findings are fixed.

## New Tests

- `T-17_missing_pore_pressure_block_aborts` — covers R-201
- `T-18_dt_initial_literal_negative_aborts` — covers R-202
- `T-19_paraview_fault_dt_zero_aborts`      — covers R-203
- `T-20_stress_kind_missing_aborts`         — covers R-205

(R-204 / R-207 / R-208 are quality-only and were verified by the existing suite continuing to pass with the renames / cleanups in place.  R-206 is a performance fix; visible via wall-time only.)

## Ready for Re-Review: YES
