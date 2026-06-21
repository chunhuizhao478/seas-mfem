# Fix Report: Phase 0+1 review findings (REVIEW_phase01_code_2026-06-03.md) — 2026-06-03

## Summary
- Findings addressed: **4 of 5 fixed**; 1 deliberately deferred (R-104, per the review's own "not a fix-agent edit" instruction).
- Files modified: `spatial/code/spatial_friction.cpp`, `drivers/spatial_seas_driver.cpp`, `tests/unit/test_spatial_seas_config.cpp`, `Makefile`.
- Tests added: 2 assertions (T-6) + 1 Makefile smoke target (`test-spatial-seas-warning`).
- Test suite: `seas_test_spatial_seas_config` **40/40 pass**; `test-spatial-seas-warning` **PASS**; no-regression intact (dyn driver builds, 0 warning leak; friction-config test links clean, 29/0, same pre-existing NUM-1 abort).

## Changes Made
1. **R-101 (MODERATE) — QD `double` knobs underflow on single-precision builds.** Added a `double`-preserving helper `toml_double` in `spatial_friction.cpp` (next to `toml_real`; accepts float **and** int literals) and switched all five QD `[time]` knob reads (`rk45_atol`, `rk45_rtol`, `dt_init`, `dt_max_years`, `plate_rate_vp`) from `toml_real` to `toml_double`. The `double` defaults (notably `rk45_rtol=1e-50`) no longer narrow through a `real_t` parameter, so on a single-precision MFEM build the default no longer underflows to 0 and the `rk45_rtol>0` guard no longer spuriously aborts the shared parser. **Deviation from the reviewer's suggested fix (noted):** the review proposed an inline lambda using `as_floating()` only, which would `MFEM_ABORT` on a TOML *integer* literal; `toml_double` handles the int branch like `toml_real`, so integer-valued knobs (`dt_init = 3`) still parse.
2. **R-102 (LOW) — QD `solver.type` validation in the shared parser.** Kept the parse-time validation (the plan's Edge Case mandates the hard error *at parse*) and added the review-requested clarifying comment documenting that it is intentional, dyn-safe shared validation (no spatial-schema config carries `[solver]`; BP5 `[solver]` is the `seas_driver` schema with key `solver_type`, parsed elsewhere). No behavioral change.
3. **R-103 (LOW) — warning untested.** Added the `test-spatial-seas-warning` Makefile target (runs the QD driver `--dry-run` on the rate-state SAFS config, asserts the `ignored by spatial_seas` warning reaches stderr). Chose the driver-level smoke (review option A) over extracting a shared helper, to avoid re-introducing the R-102 coupling.
4. **R-105 (LOW) — warning printed before the banner.** Moved the dynamic-key warning block from section 2b (before the banner) to section 3b (immediately after the banner, still before the dry-run exit). Verified ordering: banner → WARNING → dry-run message.

## Verification
- [x] R-101: `toml_double` added + used; `seas_test_spatial_seas_config` T-2 still asserts `rk45_rtol==1e-50` on absent key (the single-precision guard); new T-6 asserts integer `dt_init=3`→3.0 and `dt_max_years=1`→1.0 (guards the int branch). 40/40 pass.
- [x] R-102: clarifying comment added; validation unchanged (still aborts unknown `solver.type` at parse — T-3 passes).
- [x] R-103: `test-spatial-seas-warning` PASS.
- [ ] R-104: **not fixed — by design.** The review explicitly states "not a code change for Phase 1 (plan-faithful)… Track as a plan-doc decision, not a fix-agent edit." The `dt_init`/`dt_max_years` vs `dt_initial`/`dt_max` overlap is left as-is; recommend the plan reconcile R-004 before Phase 5.
- [x] R-105: warning relocated after the banner; ordering verified (banner line 2, WARNING line 25, dry-run line 26).

## No-regression (shared `spatial_friction.o` changed again with `toml_double`)
- `spatial_dyn_driver` rebuilds clean; QD warning leak count = **0**.
- `seas_test_spatial_friction_config` links clean (no undefined symbols), **29 PASSED / 0 FAILED**, exit 134 at the *same* pre-existing NUM-1 `material_profile` abort (unrelated, untouched).
- No debug prints / TODOs / commented-out code left behind.

## Unresolved Findings
- R-104 (dt overlap) — deferred to a plan-level R-004 reconciliation per the review's instruction (not a fix-agent edit).

## Ready for Re-Review: YES
