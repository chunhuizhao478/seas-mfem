# Plan Fix Report: PLAN_cgal_corefine_multifault.md

Source review: `PLAN_cgal_corefine_multifault_check.md` (2026-05-08)
Target document: `PLAN_cgal_corefine_multifault.md`

## Summary
- Findings addressed: **10 of 10** (R-001 .. R-010)
- Files modified: `PLAN_cgal_corefine_multifault.md` (only this file; the review document is read-only per project conventions)
- Tests added: spec-level only — Phase 1 acceptance criteria added for R-007 (deterministic sort) and R-004 (OFF round-trip); the existing Phase 1 STL round-trip criterion was tightened for R-005/R-009. No code exists to test (this is a planning document).
- Test suite: N/A (no executable code in this artifact)

## Changes Made (in priority order)

### CRITICAL

1. **R-001 — `tangential_relaxation` removed from `quality_repair::run`**
   - Phase 1 algorithm step 9 (line 151): rewrote to explicitly state the function does NOT call `tangential_relaxation`, with citation to `tangential_relaxation.h:96` showing the default `vertex_is_constrained_map` is all-`false`.
   - Phase 1 §Dependencies (line 259): removed `tangential_relaxation` from the CGAL feature list, with a back-reference to step 9.

2. **R-002 — `corefine_set` re-entrancy fixed via separate `workdir/in/` and `workdir/out/`**
   - Phase 2 CLI table (line 285): `--workdir` description rewritten to specify the two-subdirectory layout (`in/` and `out/`).
   - Phase 2 algorithm (lines 292–311): all path references updated. Step 3 writes inputs to `workdir / "in"`; step 5's subprocess invocation uses `workdir/in` as IN_DIR and `workdir/out` as OUT_DIR; steps 6–8 read from `workdir/out`. An explicit comment on step 5 cites R-002 as the rationale.

3. **R-003 — Phase 2 `--mesh-edge-size` default aligned to 1500.0**
   - Phase 2 CLI table (line 281): default changed from `--res (e.g. 2000)` to `1500.0`, with a justification noting the cross-phase consistency (matches Phase 1 default and Phase 3 `LC_NEAR`).
   - Phase 2 CLI table (line 280): `--res` description clarified to "input fixture resolution suffix to scan", explicitly noting independence from `--mesh-edge-size`.

4. **R-004 — High-precision OFF writer specified (`%.15g`)**
   - Phase 1 §Files to Create (line 107): renamed `write_off` → `write_off_high_precision`, specified `%.15g`, added explicit "Both writers" requirement and the `out.precision(15)` implementation note.
   - Phase 1 algorithm step 12 (line 154): updated to call `write_off_high_precision`.
   - Phase 1 acceptance criteria (line 254): added new criterion "OFF ASCII round-trip preserves conformality ... < 1e-9 m absolute".

### MODERATE

5. **R-005 — STL writer precision tightened from `%.10g` to `%.15g`**
   - Phase 1 §Files to Create (line 107): writer signature changed to `%.15g`; rationale section discusses the LSB-at-UTM-Y reasoning and the 6-order-of-magnitude margin vs. `Geometry.Tolerance`.
   - §Numerical constraints (line 72): rewrote the precision requirement from "≥ 10 significant decimal digits" to "≥ 15 significant decimal digits", with a worked example showing why 10 digits is non-deterministic.
   - Phase 1 acceptance criterion 4 (line 253): tightened test threshold from `< 1e-6 · bbox_diag` to `< 1e-6 m absolute`, with cross-reference to R-009.
   - Phase 3 §Conformality preservation through STL Merge (line 478): rewrote to remove the now-stale claim about `%.10g` producing sub-mm precision.
   - §Hard constraints conformality bullet (line 52): tightened to absolute units (`< 1e-9 m` for OFF, `< 1e-6 m` for STL).

6. **R-006 — `protect_constraints=true` rationale rewritten**
   - §Defence section (corresponding bullet): replaced the (incorrect) `≤ 4/3·target` precondition justification with the correct rationale (cross-mesh polyline vertex sequence consistency) and an explicit note that the cited precondition belongs to `isotropic_remeshing`, not `surface_Delaunay_remeshing`.
   - Phase 1 algorithm step 7: same fix applied to the inline rationale, with a back-reference to §Defence.

7. **R-007 — Tie-breaker for `sort_by_overlap_volume` specified**
   - Phase 1 §Files to Create entry for `intersection_graph.h` (line 106): added the deterministic comparator spec — equality test `|vol_a - vol_b| > 1e-12 · max(...)`, lexicographic `(i, j)` tie-breaker.
   - Phase 1 acceptance criteria (line 256): added a new criterion verifying determinism across rebuilds and the lexicographic tie order.

8. **R-008 — Generator propagates `LC_NEAR` and `LC_MIN` from manifest**
   - Phase 3 `.geo` template §2 size-controls block (lines 388–390): `LC_MIN`, `LC_NEAR`, `LC_FAR` are now placeholder substitutions (`<LC_MIN>`, `<LC_NEAR>`, `<LC_FAR>`) with comments naming their source fields.
   - Phase 3 `generate_multifault_geo.py` requirements (line 468): added a new bullet specifying that the generator reads `mesh_edge_size` and `min_edge` from `manifest.json`. Added `--lc-near`, `--lc-min`, `--lc-far` to the CLI flag list with their manifest-derived defaults.

### LOW

9. **R-009 — Phase 5 conformality test threshold tightened to absolute units**
   - Phase 5 test 5 (`test_pairwise_conformality`, line 590): threshold changed from `1e-6 · bbox_diag` (≈ 0.5 m for SAFS) to `1e-6 m` absolute, with a worked-out justification that the absolute threshold is 1000× tighter than `Geometry.Tolerance = 1e-3` so a passing test guarantees gmsh dedup.

10. **R-010 — Generator emits Physical Surface lines from manifest**
    - Phase 3 `.geo` template §7 Physical groups block (line 445): replaced the hardcoded six-line `Physical Surface` block with a single `<PHYSICAL_FAULT_SURFACES>` placeholder; added a comment block specifying the name-derivation rule and showing what the generator would emit for the canonical 6-fault SAFS 2000 m manifest.
    - Phase 3 `generate_multifault_geo.py` requirements (line 472): added a "Generator-emits" bullet describing the per-fault Physical Surface lines and tag assignment `100 + k + 1`.
    - Phase 3 §Acceptance Criteria (manifest-driven count): updated the "1 volume + 6 fault surface physical groups (IDs 101..106)" criterion to "1 volume + N fault surface physical groups (IDs 101..(100+N))" with the canonical 6-fault example called out.
    - The post-template prose (formerly line 446: "must be updated by hand") was rewritten to reflect that ordering is data-driven, removing the old recommendation to hand-edit the IDs.

## Unresolved Findings
- None. All 10 findings from the review (R-001 through R-010) have a corresponding fix applied to the plan.

## New Tests / Acceptance Criteria
The plan is a planning document, not executable code, so "tests" here means new acceptance criteria and proposed test functions. The check document already proposed pytest-style tests for each finding; these are summarized below for the implementer to write when Phase 1 / Phase 2 / Phase 3 are built:

| Finding | Test name (proposed in check document) | Now traced to plan section |
|---|---|---|
| R-001 | `test_R001_polyline_vertices_unchanged_by_quality_repair` | Phase 1 acceptance criterion: synthetic two-rectangles test (line 251) — already verifies polyline vertex preservation; the R-001 fix removes the source of drift |
| R-002 | `test_R002_idempotent_under_repeated_runs` | Phase 2 acceptance criterion (line 322) — "Re-running with the same args is idempotent" |
| R-003 | `test_R003_default_mesh_edge_size_consistent` | Phase 1 + Phase 2 + Phase 3 default consistency now enforced by table values |
| R-004 | `test_R004_off_writer_precision_at_utm_scale` | Phase 1 acceptance criterion: OFF round-trip preserves conformality (line 254, NEW) |
| R-005 | `test_R005_writer_lsb_below_geometry_tolerance` | Phase 1 acceptance criterion: STL round-trip < 1e-6 m absolute (line 253, tightened) |
| R-006 | N/A (rationale-only) | §Defence, Phase 1 step 7 |
| R-007 | `test_R007_sort_by_overlap_volume_deterministic` | Phase 1 acceptance criterion: deterministic sort (line 256, NEW) |
| R-008 | `test_R008_generator_propagates_mesh_edge_size` | Phase 3 generator requirements (line 468), Phase 3 acceptance: "MD5 stable across re-runs" implies `LC_NEAR` matches manifest |
| R-009 | `test_R009_conformality_tolerance_below_geometry_tolerance` | Phase 5 test 5 (tightened, line 590) |
| R-010 | `test_R010_generator_emits_physical_surfaces_from_manifest` | Phase 3 acceptance: "N fault surface physical groups (IDs 101..(100+N))" |

## Verification

Each finding traced to the corresponding edit:

- [x] R-001: tangential_relaxation removed from spec — verified at lines 151 & 259
- [x] R-002: separate workdir/in and workdir/out — verified at lines 285, 292, 303, 306, 309, 310, 311
- [x] R-003: Phase 2 default → 1500.0 — verified at line 281
- [x] R-004: write_off_high_precision (%.15g) added — verified at lines 107, 154, 254
- [x] R-005: STL writer → %.15g, conformality bullets tightened — verified at lines 52, 72, 107, 154, 253, 478
- [x] R-006: §Defence rationale rewritten + Phase 1 step 7 in-line note — verified at the §Defence section and line 149
- [x] R-007: sort_by_overlap_volume tie-breaker specified — verified at lines 106, 256
- [x] R-008: generator reads mesh_edge_size/min_edge from manifest — verified at lines 388–390, 468
- [x] R-009: 1e-6 · bbox_diag → 1e-6 m absolute — verified at line 590
- [x] R-010: <PHYSICAL_FAULT_SURFACES> placeholder + generator emits per-fault lines — verified at lines 445, 472

## Ready for Re-Review: YES

All review findings addressed. The plan is internally consistent: defaults aligned across Phases 1/2/3/4; precision contract upgraded to `%.15g` end-to-end; conformality risk from `tangential_relaxation` eliminated; re-entrancy bug eliminated; generator is data-driven from the manifest. The remaining open questions (multi-mesh corefine API in CGAL 6.1.1; `data_corefined/` version-control policy) are pre-existing flags for the implementer, not new issues introduced by these fixes.
