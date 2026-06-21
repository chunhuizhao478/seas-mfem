# Code Review: Phase 0 + Phase 1 code changes — spatial_seas QD driver (2026-06-03)

## Review Scope
- Plan: `document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md` (Phases 0 + 1, incl. manifest §C / R-001 / R-008 supersessions).
- Files reviewed (this changeset):
  - `drivers/spatial_seas_driver.cpp` — new (Phase 0 skeleton + Phase 1 dynamic-key warning).
  - `spatial/code/spatial_friction.hpp` — modified (fwd-decl `SolverType`; `SolverSpec`; `TimeSpec` QD knobs; `solver` member; `ParseQDSolverType` decl).
  - `spatial/code/spatial_friction.cpp` — modified (`elasticity_operator.hpp` include; `[solver]`/QD-`[time]` parse + validation; `ParseQDSolverType` def).
  - `tests/unit/test_spatial_seas_config.cpp` — new (38 assertions).
  - `Makefile` — modified (test target/obj/vars).
- Domain context: `CLAUDE.md` (root + worktree), project memory (no-touch shared-parser/dyn-driver constraint, dt-overlap R-004 trail), plan manifest §C, the implementer completion report, prior reviews `REVIEW_phase0_code_2026-06-03*.md`.
- Method: three adversarial passes; build + run exercised; cross-checked every existing `*.toml` against the shared parser for regression.

## Build / regression evidence
- `seas_test_spatial_seas_config`: **38/38 pass**, exit 0.
- `spatial_friction.o` compiles clean with the new heavy include (no `-DSEAS_USE_MPI`); `seas_test_spatial_friction_config` **links clean** (ABI/minimal-link-set preserved) — 29 pass, then the **same pre-existing** NUM-1 `material_profile` abort (exit 134), unchanged by this work.
- `spatial_dyn_driver` rebuilds clean against the modified shared header+object; the QD dynamic-key warning does **not** leak into it (0 occurrences); banner prints. (Its dry-run exit 1 is an absent 500 m mesh file — environmental.)
- Tree scan: **no spatial-schema config carries a `[solver]` table** — all 7 `[solver]` tables are BP5/`seas_driver` schema (`solver_type`, not parsed by `ParseSpatialFrictionConfigString`).

## Findings

### [R-101] MODERATE [POSSIBLE] [spatial_friction.cpp:parse_root] — QD `double` time knobs read through `toml_real` (real_t) underflow on a single-precision MFEM build → default config aborts in the SHARED parser

**Category:** BUG (latent / precision)

**Description:**
The new QD time fields are `double` (deliberately — `rk45_rtol=1e-50` cannot be represented in `float`). But they are populated via `toml_real(...)`, whose default-value parameter is `real_t`:
```cpp
cfg.time.rk45_rtol = toml_real(t, "rk45_rtol", cfg.time.rk45_rtol);  // default_val is real_t
```
If MFEM is built single-precision (`real_t == float`), the default `1e-50` (passed as the `real_t default_val`) underflows to `0.0f` **before** it ever reaches the field. When the key is absent (the common case — no QD config sets `rk45_rtol`), the field becomes `0.0`, and the new guard
```cpp
MFEM_VERIFY(cfg.time.rk45_rtol > 0.0, "[time].rk45_rtol must be > 0; …");
```
**aborts the parse**. Because this VERIFY lives in the *shared* `parse_root`, on a single-precision build it would abort **every** spatial config — including all dynamic-rupture (TPV/SAFS) configs — i.e. a CRITICAL dynamic-driver regression. On the current double build it is harmless (hence MODERATE/POSSIBLE), but it silently defeats the plan's deliberate `double` choice.

**Trigger:** MFEM built with `MFEM_USE_SINGLE=YES` (or `-DMFEM_USE_SINGLE`); parse any spatial config that omits `[time].rk45_rtol`.

**Actual behavior:** `rk45_rtol` defaults to `0.0` (underflow) → `MFEM_VERIFY(... > 0.0)` aborts.

**Expected behavior:** the `double` default `1e-50` survives regardless of `real_t`, so the default config parses.

**Suggested fix:** read the QD `double` fields as true doubles (bypass the `real_t` default-arg narrowing). Minimal, targeted for the underflow-prone fields:
```diff
-      cfg.time.rk45_atol     = toml_real(t, "rk45_atol",     cfg.time.rk45_atol);
-      cfg.time.rk45_rtol     = toml_real(t, "rk45_rtol",     cfg.time.rk45_rtol);
-      cfg.time.dt_init       = toml_real(t, "dt_init",       cfg.time.dt_init);
-      cfg.time.dt_max_years  = toml_real(t, "dt_max_years",  cfg.time.dt_max_years);
-      cfg.time.plate_rate_vp = toml_real(t, "plate_rate_vp", cfg.time.plate_rate_vp);
+      // Read as true double (the fields are double on purpose: rk45_rtol=1e-50
+      // underflows float).  contains-guard keeps the struct default otherwise.
+      auto rd = [](const toml::value& tb, const char* k, double d) -> double
+      { return tb.contains(k) ? tb.at(k).as_floating() : d; };
+      cfg.time.rk45_atol     = rd(t, "rk45_atol",     cfg.time.rk45_atol);
+      cfg.time.rk45_rtol     = rd(t, "rk45_rtol",     cfg.time.rk45_rtol);
+      cfg.time.dt_init       = rd(t, "dt_init",       cfg.time.dt_init);
+      cfg.time.dt_max_years  = rd(t, "dt_max_years",  cfg.time.dt_max_years);
+      cfg.time.plate_rate_vp = rd(t, "plate_rate_vp", cfg.time.plate_rate_vp);
```
(`as_floating()` returns `double` in toml11; an integer-typed TOML value would need `is_integer()` handling like `toml_real` if integer literals must be accepted — acceptable to require float syntax for these tolerances, or replicate `toml_real`'s int branch in the lambda.)

**Test case:**
```cpp
// test_R101_rk45_rtol_default_survives_precision:
//   Parse a minimal config WITHOUT [time].rk45_rtol.
//   Assert cfg.time.rk45_rtol == 1e-50 (not 0).
//   On a single-precision build this currently aborts; with the fix it passes.
```

---

### [R-102] LOW [POSSIBLE] [spatial_friction.cpp:parse_root] — QD-specific `[solver]` validation runs in the SHARED parser for the dynamic driver too (latent coupling)

**Category:** ASSUMPTION / DEVIATION (dyn-safety philosophy)

**Description:**
`(void)ParseQDSolverType(cfg.solver.type);` plus the four `[solver].*` `MFEM_VERIFY`s execute unconditionally in `parse_root`, i.e. for **every** config the shared parser handles, including the dynamic driver's. This is the same class of coupling that R-001 deliberately moved *out* of the parser (the dynamic-key warning). It is **currently harmless**: no spatial-schema config has a `[solver]` table, so `cfg.solver.type` is always the default `"cg_amg"` (valid) and the knob VERIFYs see valid defaults. But it is a latent dyn-safety regression surface: if a future dynamic config ever sets `[solver].type` to a value valid for the dynamic path but outside the QD set, the **shared** parser would abort it. The plan's Edge Case ("unknown `solver.type` ⇒ hard error at parse") motivates parse-time validation, so this is plan-sanctioned — flagged so the coupling is explicit, not silently assumed safe.

**Trigger:** a future `[solver].type="<non-QD-value>"` in a config parsed by `ParseSpatialFrictionConfigString` for the dynamic driver.

**Actual behavior:** shared parser aborts on any `[solver].type` not in the QD set.

**Expected behavior (per R-001 philosophy):** QD-only validation should not gate the dynamic driver's parse.

**Suggested fix (optional; only if the coupling is judged unacceptable):** validate `solver.type` in the QD driver after load (alongside the dynamic-key warning) instead of inside `parse_root`, OR gate the `[solver]` validation on a QD marker. If kept in the parser, add a one-line comment that this is intentional shared validation. No code change required for current correctness.

**Test case:**
```cpp
// test_R102_dyn_config_with_foreign_solver_type_still_parses (only meaningful
// if the validation is relocated): a config with [solver].type="petsc_gamg"
// parses without abort under the dynamic driver's path.
```

---

### [R-103] LOW [spatial_seas_driver.cpp:main] — dynamic-key "ignored" warning has no automated test

**Category:** QUALITY (coverage)

**Description:**
The Phase 1 deliverable "ignores dynamic-only keys with a clear warning" is implemented in the driver (correct per §C) and was verified only by manual smoke (it printed `WARNING: [numerics] key(s) 'ader_order' 'mixed_flux' 'use_pml' …`). The unit test (`test_spatial_seas_config.cpp` T-5) only asserts the parser does **not** throw — it does not exercise the driver-side warning emission, which could regress silently (e.g., a typo'd key name in `dyn_keys[]`, or the block being skipped).

**Trigger:** any future edit to the warning block.

**Suggested fix:** add a small driver-level smoke target (run `seas_spatial_seas_driver --config <cfg-with-ader_order> --dry-run`, grep stderr for `ignored by spatial_seas`), or refactor the detection into a tiny pure helper (e.g. `std::vector<std::string> DetectIgnoredNumericsKeys(const toml::value&)`) that a unit test can call directly. Low priority; the logic is simple and currently correct.

**Test case:**
```cpp
// test_R103_warning_lists_present_dynamic_keys:
//   feed a toml::value with [numerics].ader_order + use_pml to the helper;
//   assert it returns {"ader_order","use_pml"} and NOT mixed_flux/interior_flux.
```

---

### [R-104] LOW [spatial_friction.hpp:TimeSpec] — `dt_init`/`dt_max_years` duplicate `dt_initial`/`dt_max` (R-004 overlap), wrong-field hazard

**Category:** QUALITY / DEVIATION (known, plan-faithful)

**Description:**
`TimeSpec` now carries both the dynamic `dt_initial`(s)/`dt_max`(s) and the QD `dt_init`(s)/`dt_max_years`(yr). The implementer followed the plan/Phase-5 appendix literally (which reads `cfg.time.dt_init` / `cfg.time.dt_max_years`) and disambiguated with comments, but four near-homonym fields in one struct invite a future Phase-2+ author to read the wrong one (e.g. `dt_max` seconds vs `dt_max_years` years — a 10⁶–10⁹× error in the QD dt ceiling). This is the unresolved R-004 from `REVIEW_dyn_impact_2026-06-02.md`.

**Trigger:** Phase 2+ wiring reading `cfg.time.dt_max` (seconds) where it meant `dt_max_years`.

**Suggested fix:** not a code change for Phase 1 (plan-faithful). Recommend the plan reconcile R-004 before Phase 5 — e.g. drop `dt_init` and reuse `dt_initial` (same `<0 ⇒ derive` semantics), keeping only `dt_max_years` as the genuinely new (years) field. Track as a plan-doc decision, not a fix-agent edit.

**Test case:** N/A (design/plan decision).

---

### [R-105] LOW [spatial_seas_driver.cpp:main] — dynamic-key warning prints before the banner

**Category:** QUALITY

**Description:** the warning block (section 2b) runs before the rank-0 banner (section 3), so stderr shows the `WARNING:` line above the config echo. Cosmetic; some users expect the banner first for context. Both go to different streams (warning→stderr, banner→stdout), so interleaving depends on the terminal.

**Suggested fix (optional):** move the warning block to just after the banner, or leave as-is (stderr/stdout separation makes ordering moot in piped output). No correctness impact.

**Test case:** N/A.

---

## Summary
- Critical issues: **0** (no bug fires on the current double-precision build with any existing config).
- Moderate issues: 1 (R-101 — latent single-precision underflow that would abort the shared parser; POSSIBLE).
- Low issues: 4 (R-102 shared-parser coupling; R-103 warning untested; R-104 dt overlap; R-105 ordering).
- Plan compliance: **FULL** for Phase 0 + Phase 1 normative requirements —
  - `SolverSpec` (8 fields, exact defaults incl. `residual_check=true`), `TimeSpec` QD knobs, `SolverSpec solver;` member ✔
  - `[solver]` + QD `[time]` parsed; unknown `solver.type` is a hard parse error ✔ (Edge Case)
  - `ParseQDSolverType` returns `mfem::seas::SolverType` (R-008), all 7 mappings correct, aborts on unknown ✔
  - dynamic-key warning in the **driver**, not the shared parser (§C / R-001) ✔ — verified it does not leak into `spatial_dyn_driver`
  - new test `seas_test_spatial_seas_config` (38/38) ✔; existing config test unaffected (same pre-existing failure) ✔
  - `SolverType` forward-decl in header + heavy include only in `.cpp`: ABI verified (dyn driver + config test link clean) ✔
- Verdict: **PASS WITH FIXES** — no CRITICAL bug on the current build; fix R-101 (precision/underflow) to harden the shared parser, and consider R-102 (coupling). R-104 is a plan-level reconciliation, not a fix-agent edit.

## Unreviewed Areas
- Full `make all && make test` not run end-to-end (heavy; one pre-existing unrelated failure). Non-interference verified structurally + by rebuilding the dynamic driver and the existing config test.
- Single-precision (`MFEM_USE_SINGLE`) build not actually exercised — R-101 is reasoned from the type flow, flagged POSSIBLE.
- Phase 2+ consumption of the new config fields (operator wiring, RK45) — out of scope.
