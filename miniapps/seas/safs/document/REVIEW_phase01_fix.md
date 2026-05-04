# Fix Report: REVIEW_phase01_code.md (2026-04-29)

## Summary

- Findings addressed: **11 of 11**
- Files modified: 6
- Tests added: 3 sub-tests inside `test_corefine_smoke` (R-002, R-004, expanded coincident)
- Test suite:
  - **C++ ctest** (`tools/build/`): 1/1 pass — 6 sub-tests inside the smoke binary all pass.
  - **Python pytest** (`miniapps/seas/safs/mesh/tests/`): 42/42 pass.
  - **End-to-end**: two_squares fixture + Mill Creek × SBMT-SAF 2000 m fixture both run cleanly; all four gates report PASS in `intersection_report.json`.
  - **Validator can fire**: corrupting the on-disk JSON's `n_total_triangles` or a `range` entry now triggers exit code 3 (verified manually with two corruptions).

## Changes Made

### Critical

1. **R-001 — tautological schema check replaced with on-disk re-read.** `tools/corefine_faults/io.{hpp,cpp}` now exports `validate_triangle_to_fault_json(path, per_fault)` that reads the just-written JSON, parses `n_total_triangles` and per-fault `range` entries, and verifies the partition is contiguous and sums match the in-memory vector. `main.cpp` calls this validator after writing the JSON; old tautological loop deleted.
2. **R-002 — plan acceptance criterion reconciled.** `PLAN_cgal_corefine.md` Phase 1 acceptance line (originally `>200` / `≥500 floor`) updated to pin `143 ± 10` constrained edges per side on the 2000 m fixture, plus a `> 50` floor for "corefine found nothing", with rationale documenting why the original 500-floor was an over-extrapolation from Python's snap-induced sub-edge over-count. The smoke test gains a `test_mill_creek_pin` that locks the count and skips gracefully when the fixture is absent (so a fresh checkout doesn't fail CI before the Python pipeline has run).

### Moderate

3. **R-003 — manifold-gate exit propagation.** Both `manifold_gate(meshes[i])` and `manifold_gate(meshes[j])` now return exit 2 with a fault-named diagnostic if they report anything other than `"PASS"`. Same wiring added for `gate_interior_only`. The JSON still records the exact FAIL string for diagnostic continuity, but the binary cannot exit 0 with a FAIL gate anymore.
4. **R-004 — real polyline-edge-coincidence gate.** Added `safs::corefine::polyline_edge_coincidence_gate` in `corefine.{hpp,cpp}`: collects canonical-ordered constrained-edge endpoint pairs from each mesh, returns true iff the sets are equal. Replaces the unconditional `"PASS"` with the actual check + a fatal exit on disagreement. Smoke test `test_polyline_edge_coincidence_gate` covers both the positive case (perpendicular squares) and the negative control (mis-paired ECMs).
5. **R-005 — independent cross-check between `meshes[*].number_of_faces()` and `per_fault[*].n_triangles`.** Added before the on-disk validator, so a stale-snapshot bug between mesh mutation and per_fault population is caught with exit 3.
6. **R-006 — `find_coincident_triangle` uses `f.idx()`.** Switched both A-loop and B-loop from enumeration counters to the stable `Face_index::idx()` value, so the reported triangle indices remain correct under any future `remove_face` / `garbage_collect` cycle.
7. **R-007 — `stl_path` is filename-only.** `write_intersection_report_json` now emits `std::filesystem::path(f.stl_path).filename().string()`, so the JSON remains correct after the output directory is moved. Verified: the new Mill Creek run shows `"stl_path": "safs_sbmt_millcreek.stl"` (was `"/tmp/.../safs_sbmt_millcreek.stl"`).

### Low

8. **R-008 — `<cstring>` added explicitly** to `corefine.cpp` so `std::memcpy` no longer relies on the transitive include from `<unordered_map>`.
9. **R-009 — `+0.0` / `-0.0` canonicalised** in both `make_tri_key` and `find_coincident_triangle`'s lex/equality comparators via a `canon_zero(v)` helper. Hash and explicit comparison now agree on IEEE-754's `+0.0 == -0.0` semantics.
10. **R-010 — `repair_polygon_soup(require_same_orientation(true))`.** Switched from `false` so opposite-winding duplicate triangles (a buggy `ts_to_stl.py` emitting both sides of a closed surface) are no longer silently merged into a non-orientable mesh.
11. **R-011 — `corefine_pair` takes ECMs by const reference.** Header signature and implementation now match the plan; underlying handle semantics are unchanged but the "shared across pairs (P-002)" intent is no longer ambiguous.

## Files Modified

- `miniapps/seas/safs/tools/corefine_faults/io.hpp` — added `validate_triangle_to_fault_json` declaration.
- `miniapps/seas/safs/tools/corefine_faults/io.cpp` — added validator implementation; changed `repair_polygon_soup` flag (R-010); filename-only `stl_path` (R-007).
- `miniapps/seas/safs/tools/corefine_faults/corefine.hpp` — `corefine_pair` ECMs now `const&` (R-011); declared `polyline_edge_coincidence_gate` (R-004).
- `miniapps/seas/safs/tools/corefine_faults/corefine.cpp` — `<cstring>` include (R-008); `canon_zero` for `+0/-0` (R-009); `f.idx()` (R-006); `polyline_edge_coincidence_gate` impl (R-004); ECMs by `const&` in definition (R-011).
- `miniapps/seas/safs/tools/corefine_faults/main.cpp` — gate exit propagation (R-003); R-005 cross-check; R-001 on-disk validator call.
- `miniapps/seas/safs/tools/tests/test_corefine_smoke.cpp` — `test_polyline_edge_coincidence_gate` (R-004); `test_mill_creek_pin` (R-002); now also includes `io.hpp` for the pin test.
- `miniapps/seas/safs/tools/CMakeLists.txt` — `test_corefine_smoke` now links `corefine_faults/io.cpp` (needed by the Mill Creek pin test).
- `miniapps/seas/safs/PLAN_cgal_corefine.md` — Phase 1 acceptance line for Mill Creek constrained-edge count (R-002).

## New Tests

- `test_polyline_edge_coincidence_gate` — covers R-004 (positive: perpendicular squares' constrained edges match; negative: empty B-side ECM disagrees with populated A-side).
- `test_mill_creek_pin` — covers R-002 (locks `ce ∈ [133, 153]` on the real fixture; skips if fixture absent).
- (Indirect) Validator test was performed manually against two intentional JSON corruptions; both correctly returned exit 3. The smoke test does not yet inject corruption — adding such a test would require a separate validator-level test fixture; left out of scope of this fix pass to avoid leaking debug code into the smoke binary.

## Verification Checklist

- [x] R-001 — Validator re-reads the on-disk JSON; manually verified it catches a corrupted `n_total_triangles` (rc=3) and a corrupted `range` entry (rc=3).
- [x] R-002 — Plan updated to pin 143 ± 10; smoke test enforces the bound.
- [x] R-003 — Manifold gate FAIL → exit 2; interior gate FAIL → exit 2.
- [x] R-004 — Real coincidence gate replaces tautological PASS; smoke test covers positive + negative cases.
- [x] R-005 — Cross-check between live mesh face counts and `per_fault` sums.
- [x] R-006 — `find_coincident_triangle` returns `Face_index::idx()`.
- [x] R-007 — `stl_path` is filename-only (verified: Mill Creek JSON shows `"safs_sbmt_millcreek.stl"`).
- [x] R-008 — `<cstring>` explicitly included.
- [x] R-009 — `canon_zero` applied in `make_tri_key`, hash, lex sort, equality compare.
- [x] R-010 — `require_same_orientation(true)`.
- [x] R-011 — `corefine_pair` ECMs now `const&` in both header and implementation.

## Unresolved Findings

None. All 11 findings addressed.

## Notes on Reviewer's Suggested Fixes

- **R-001 suggested fix**: applied with one small refactor — extracted the validator into `safs::io::validate_triangle_to_fault_json` so it can be unit-tested in isolation rather than living as an inline block in `main.cpp`. Functionally equivalent.
- **R-002 fix mode**: applied **option A** (pin the empirical count). The reviewer said this is preferred unless evidence suggests CGAL is under-counting; the geometric reasoning at 2 km mesh resolution × ~286 km crossing length matches 143 segments, so option A is correct.
- **R-003 suggested fix**: applied verbatim plus the symmetric `gate_interior_only` propagation (review explicitly asked for both).
- **R-004 suggested fix**: applied. Helper lives in `corefine` namespace alongside the ECM type rather than in `main.cpp` (the review's pseudocode placed it in `main.cpp` as a free function; namespace-scoping is more idiomatic and lets the smoke test reuse it).
- **R-007 suggested fix**: applied to the writer only; `main.cpp:per_fault[i].stl_path = out` left as the absolute path because the writer canonicalises to filename anyway, and keeping the absolute path in the in-memory struct is useful for the verbose log line.

## Ready for Re-Review: YES

Build green, all C++ smoke tests pass, Python regression 100%, end-to-end on real fixture clean, validator demonstrably catches schema corruption.
