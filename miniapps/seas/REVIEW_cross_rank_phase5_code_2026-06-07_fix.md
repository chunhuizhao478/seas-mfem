# Fix Report — Cross-rank Phase 5 (against REVIEW_cross_rank_phase5_code_2026-06-07.md)

**Date:** 2026-06-07. Review verdict: **0 critical**, 1 moderate, 2 low. The reviewer confirmed the
deprecation's core claim — `seam_continuous_` has ZERO live readers (only the setter assignment, the
member decl, and a comment; the central-build abort consumer was removed in Phase 2). All findings
addressed. No critical ⇒ phase complete.

## Changes made
- **P5-001 (MODERATE) + P5-003 (LOW):** `tests/parallel/test_bimaterial_mixed_flux_shared.cpp` —
  rewrote the file-header 3.4a/3.4b/3.4c descriptions and the inline comments that still narrated
  `SetSeamContinuous(true)` as a LIVE central-flux gate / abort trigger. Now: 3.4b is OBSOLETE
  (the abort was removed in cross-rank Phase 2); 3.4a/3.4c build for ANY material via the true peer;
  the `SetSeamContinuous(true)` call is annotated as a deprecated no-op; test labels no longer say
  "seam_continuous=true". The test still passes 16/16 (the setter is a harmless store).
- **P5-002 (LOW, noted):** `--partition-fault-locality` is not exposed by `spatial_dyn_driver.cpp`
  at all (it is an opt-in `else if` only in the native `tpv*_driver.cpp`) and the spatial driver
  never required it. So "make it optional" is moot for the spatial path; the cross-rank exchange's
  correctness without fault-locality is demonstrated by the np=2 gates, which deliberately CUT the
  fault across the seam (PartitionByYSign / explicit seam partition) and pass.

## Phase-5 deprecation sites (the actual diff — all accurate + consistent, verified by the reviewer)
- `dynamic/bimaterial_wave_operator.hpp` — `SetSeamContinuous` + `seam_continuous_` docs: "DEPRECATED,
  no effect, stored for config back-compat, never read".
- `drivers/spatial_dyn_driver.cpp` — the propagation-call comment + the `[mixed-flux] (...)` banner
  ("DEPRECATED, no effect — cross-rank seam material uses the TRUE peer").
- `spatial/code/spatial_friction.hpp` + `.cpp` — the `MaterialSpec::seam_continuous` field + TOML
  parse comments ("DEPRECATED, parsed for back-compat, no effect; key may be omitted").

## Verification
- Config-parse gate `seas_test_tpv_config_parse` (R-208): **184/184** — `seam_continuous` absent ⇒
  defaults false; `=true` ⇒ round-trips; both parse with NO abort (back-compat preserved). (Built
  with the main-repo toml11 override, the worktree submodule being absent.)
- `seas_test_bimaterial_mixed_flux_shared` (np=2): **16/16** after the comment/label edits.
- Full suite green (see the overall summary): seam_material 135/135, seam_fault 24/24, perside 11/11,
  parity 9/9, central 6/6, dispatch 10/10, constant-parity 19/19.

## np=8 TPV6 dry-run smoke (NON-gating): deferred to Frontera (the TPV6 mesh is gitignored; this is
a production smoke, not a local gate, per the plan).

## Phase complete: YES (0 critical; deprecation accurate + back-compat preserved; tests green).
