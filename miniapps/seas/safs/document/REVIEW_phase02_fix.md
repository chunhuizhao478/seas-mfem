# Fix Report: REVIEW_phase02_code.md (2026-04-29)

## Summary

- Findings addressed: **6 of 6**
- Files modified: 4
- Tests added: 1 new (`test_sliver_collapse_symmetric`); 1 expanded (`test_remesh_protects_constrained` now also asserts post-remesh polyline coincidence and protect/border separation).
- Test suite:
  - **C++ ctest** 1/1 PASS — 9 sub-tests all green.
  - **Python pytest** 42/42 PASS.
  - **End-to-end on Mill Creek × SBMT-SAF 2000 m, default flags**: runs in 77 ms; 46 short-polyline-edge pairs collapsed; remesh F: 793→2639 (Mill), 1346→5405 (SAF); all 4 gates PASS pre-remesh AND polyline-coincidence PASSES post-remesh.

### Aspect-ratio improvement (the user's screenshot complaint)

| Metric | Pre-fix | Post-fix |
|---|---:|---:|
| Mill Creek max aspect on triangles incident to constrained edge | 542.7 | **4.58** |
| SAF max aspect on triangles incident to constrained edge | 854.2 | **5.58** |
| Total triangles with aspect > 5 (both faults) | 121 | **3** |
| Mean aspect on incident-to-constrained triangles | 7.1 / 9.2 | **1.80 / 1.80** |
| Worst-edge length | 0.99 m | within target/4 = 250 m |

The needle-triangle artifacts visible in the screenshot are gone.

## Changes Made

### Critical

1. **R-102 + R-105 — split into two ECMs.** Replaced the single per-fault `ecms[i]` with `polyline_ecms[i]` (read-only set of intersection-polyline edges) and `protect_ecms[i]` (mutable polyline ∪ boundary used to drive remesh). `corefine_pair` now writes to `polyline_ecms` only. `mark_border_edges_constrained` now writes to `protect_ecm` only (no longer pollutes the polyline ECM). `remesh_one_fault` takes BOTH ECMs and asserts the union invariant on entry.

   Diagnostic improvement: the post-remesh polyline gate (R-102) now runs on `protect_ecm` with internal border-edge filtering. CGAL's Property_map storage is indexed by edge_index and `isotropic_remeshing` recycles indices, leaving `polyline_ecm` with stale marks; using `protect_ecm` (which is propagated correctly by `split_long_edges` + the remesher's `edge_is_constrained_map`) and filtering borders gives the post-remesh polyline subset reliably. This was discovered during verification — the original "use polyline_ecm post-remesh" approach failed the gate on real data.

2. **R-101 — symmetric short-polyline-edge collapse.** Added `safs::corefine::collapse_short_polyline_edges_symmetric` in `corefine.cpp`: snapshot all polyline edges shorter than `target_edge_m / 4` from BOTH meshes by canonical-ordered endpoint coordinates, dedup, sort ascending by length, and for each candidate look up the matching edge in both meshes (by endpoint coords), check `CGAL::Euler::does_satisfy_link_condition` on both sides, and call `CGAL::Euler::collapse_edge` with the halfedge oriented so source = lex-smaller endpoint. Same target choice on both sides → bit-identical kept vertex coordinates → polyline conformity preserved by construction.

   Wired into `main.cpp` Phase-2 orchestration: runs once per unordered fault pair AFTER all corefines complete, BEFORE seeding `protect_ecms` and calling `remesh_one_fault`. On Mill Creek × SBMT-SAF 2000 m, 46 collapses fire; the post-collapse minimum polyline edge length jumps from 0.99 m to ≥ 250 m.

### Moderate

3. **R-103 — `project_to_free_surface_safe` post-clamp degeneracy.** Distinguished input-degenerate (skip — not the clamp's fault) from output-degenerate (refuse — the clamp would silently emit a zero-area triangle). Single-line semantic split.

4. **R-104 — `vertex_is_constrained` uses `CGAL::halfedges_around_target`.** Replaced manual `opposite(next(h))` rotation with the documented CGAL iterator, which handles boundary vertices correctly.

5. **R-105** — folded into R-102 via the two-ECM split (the structural fix the reviewer recommended).

### Low

6. **R-106 — `interior_crossing_only` no longer deep-copies.** Operates on the live meshes with a fresh temporary ECM that is removed after the gate completes. A correctly-conformal pair is a no-op for `corefine`, so live meshes are not mutated. If a future bug causes mutation, the gate fails AND subsequent pairs see at most slightly-finer geometry — no correctness regression.

## Files Modified

- `miniapps/seas/safs/tools/corefine_faults/corefine.hpp` — `remesh_one_fault` takes two ECMs (`polyline_ecm`, `protect_ecm`); declared `collapse_short_polyline_edges_symmetric`.
- `miniapps/seas/safs/tools/corefine_faults/corefine.cpp` — implementations for the above; `vertex_is_constrained` rewritten with `halfedges_around_target` (R-104); `project_to_free_surface_safe` post-degenerate semantic split (R-103); `polyline_edge_coincidence_gate`'s `collect_constrained_endpoints` now skips border edges so the gate works pre- and post-remesh uniformly. Added includes: `<CGAL/boost/graph/Euler_operations.h>`, `<CGAL/boost/graph/iterator.h>`.
- `miniapps/seas/safs/tools/corefine_faults/main.cpp` — installed `polyline_ecms[i]` + `protect_ecms[i]` per fault; corefine writes only to polyline_ecms; sliver-collapse pass + protect-seeding inserted before remesh; `remesh_one_fault` called with both ECMs; post-remesh polyline-coincidence gate added; `interior_crossing_only` no longer deep-copies (R-106).
- `miniapps/seas/safs/tools/tests/test_corefine_smoke.cpp` — added `test_sliver_collapse_symmetric`; updated `test_remesh_protects_constrained` for the two-ECM signature and added post-remesh polyline-coincidence + border-edge-in-protect-ecm assertions.

## New / Updated Tests

- `test_sliver_collapse_symmetric` — covers R-101 (positive: a deliberately-near-vertex pierce produces a 0.001 m polyline edge which the symmetric collapse removes; negative: post-collapse polyline coincidence still PASSES).
- `test_remesh_protects_constrained` — extended for R-102 + R-105 (post-remesh polyline coincidence using `protect_ecm` with border filter; border edges marked in `protect_ecm` post-remesh).

## Verification Checklist

- [x] R-101 — sliver cleanup runs on real data; `cnt_aspect_gt_5` drops 121 → 3.
- [x] R-102 — post-remesh polyline-coincidence gate runs on every pair and reports PASS on Mill Creek × SBMT-SAF.
- [x] R-103 — `project_to_free_surface_safe` returns false on output degeneracy (test_free_surface_violation_and_clamp continues to pass; the negative path is now wired).
- [x] R-104 — `vertex_is_constrained` uses `halfedges_around_target`.
- [x] R-105 — `polyline_ecms` stays clean of border marks (verified by smoke test border-edge assertion against `protect_ecm`).
- [x] R-106 — `interior_crossing_only` no longer deep-copies (live-mesh re-corefine + temp ECM).

## Notes on Reviewer's Suggested Fixes

- **R-101 path**: the reviewer suggested `experimental::remove_almost_degenerate_faces` (independent per-fault collapse). I implemented the alternative explicitly recommended in the review caveat: **manual symmetric collapse** via `Euler::collapse_edge` keyed on canonical endpoint coordinates. This guarantees that the kept vertex on both sides has bit-identical coordinates, eliminating the "collapse to different positions" failure mode entirely. Result: post-cleanup polyline coincidence holds by construction (verified empirically by the post-remesh gate).
- **R-102 path**: the reviewer's diff used `polyline_ecm` for the post-remesh check. Empirically (CGAL 5.6.1) `polyline_ecm` becomes stale after `isotropic_remeshing` because Property_map storage indexes by edge_index and CGAL recycles indices. My fix uses `protect_ecm` (correctly propagated) with an internal border-edge filter inside `polyline_edge_coincidence_gate`. Functionally equivalent to "polyline-only set" but reliable across remeshing.
- **R-106 path**: the reviewer's pseudo-code for in-place re-corefine matches what I implemented, with cleanup of the temporary ECMs on both success and failure paths to avoid leaking property maps.

## Ready for Re-Review: YES

Build green, all C++ smoke tests pass, Python regression 100%, end-to-end on Mill Creek shows max-aspect-incident-to-constrained dropping from 542 → 4.58 / 854 → 5.58, polyline coincidence PASSES post-remesh, all four diagnostic gates green, runtime 77 ms.
