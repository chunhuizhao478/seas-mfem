# Fix Report: Phase 2 review findings (REVIEW_phase2_code_2026-06-03.md) — 2026-06-03

## Summary
- Findings addressed: **5 of 6 fixed**; 1 deliberately deferred (R-205 — the review states "none required for Phase 2; track for Phase 7").
- Files modified: `drivers/spatial_seas_driver.cpp` (only).
- Tests added: 0 (all findings LOW). R-201's new abort and R-203's fail-fast ordering are verified by driver smokes (below); no unit-test harness exists at the driver level for Phase 2.
- Build: clean. Smokes + no-regression: all PASS.

## Changes Made
1. **R-201** (BP5 `dirichlet_attrs` fallback foot-gun) → the bc now distinguishes "[boundary] present" (`fault_attr > 0`) from "absent": when present, `dirichlet_attrs` MUST be non-empty (clear QD-specific abort); the BP5 `{5}` fallback applies only when `[boundary]` is entirely absent. Verified: a config with `[boundary]` but no `dirichlet_attrs` now aborts with "…requires far-field Dirichlet plate-loading walls…" instead of a downstream "Dirichlet attr 5 not in mesh".
2. **R-203** (hetero abort after mesh load) → moved the `hetero_requested` `MFEM_VERIFY` to a new section 4.0, **before** the mesh load (the check is config-only). Verified: a heterogeneous config aborts with **0** mesh-load diagnostics printed before the abort (fail-fast).
3. **R-202** (dead `MaterialField` + vacuous mode assert) → removed; section 4.4 now builds `LinearElastic` directly from `[material_constant_fallback]`. The real heterogeneous gate is R-203's 4.0 check.
4. **R-204** (rank-0 local counts mislabelled) → the getter print header is now `"…constructed (rank 0 local counts):"`.
5. **R-206** (absorbing in QD) → added a clarifying comment that `absorbing_attrs` is a dynamic-wave concept inert in the quasi-static operator (unmarked → natural/traction-free). No behavioral change (the review prescribed "document" for Phase 2).

## Verification
- [x] R-201: `[boundary]` without `dirichlet_attrs` → clear QD abort (exit 1); smoke config (which sets `dirichlet_attrs=[5]`) still constructs.
- [x] R-202: `MaterialField` removed; build clean; construction unchanged (np=1 getters `27936/27936/9312/3`).
- [x] R-203: hetero gate fires before mesh load (0 mesh diagnostics before abort).
- [x] R-204: "rank 0 local counts" label present.
- [x] R-206: documenting comment added.
- [ ] R-205: **not fixed — by design.** The review says "none required for Phase 2 … Track for Phase 7: confirm BP5 uses the DomainConfig defaults for penalty_factor/face_basis_type, or add the knobs to SolverSpec." Adding `penalty_factor`/`face_basis_type`/`match_quad_order` to `SolverSpec` is a Phase-1 config change deferred to the Phase-7 parity work; the Phase 2 defaults are correct for construction.

## Smoke + no-regression evidence
- np=1 dry-run on `bp5_tandem_exact.msh`: constructs, getters printed, **exit 0**.
- np=4 dry-run: **exit 0** (construction path unchanged for valid configs).
- R-201 abort: exit 1 with the QD-specific message.
- R-203 fail-fast: heterogeneous config aborts before any mesh-load output.
- `spatial_dyn_driver` rebuilds clean (no dyn-safety regression).
- `seas_test_spatial_seas_config` 40/40 (exit 0); `seas_test_spatial_friction_config` 29 PASSED + same pre-existing NUM-1 abort (exit 134, unrelated).
- No debug prints / TODOs / commented-out code left behind.

## Unresolved Findings
- R-205 — deferred to Phase 7 per the review's instruction (config-schema decision, not a fix-agent edit).

## Ready for Re-Review: YES
