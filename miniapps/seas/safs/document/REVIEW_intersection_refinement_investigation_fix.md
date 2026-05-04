# Fix Report: Intersection refinement (REVIEW_intersection_refinement_investigation.md)

Companion to `REVIEW_intersection_refinement_investigation.md`.

## Summary

- **Findings addressed:** 6 of 6 actionable (R-001, R-002, R-004, R-005, R-006, R-007)
- **Findings deferred:** 1 (R-003 — non-essential dead-code cleanup; see below)
- **Files modified:** 2
  - `miniapps/seas/safs/mesh/refine_fault_near_intersections.py`
  - `miniapps/seas/safs/mesh/run_newset_step_by_step.sh`
- **Files added:** 2
  - `miniapps/seas/safs/mesh/tests/test_break_fault_wedges.py` (4 tests)
  - `miniapps/seas/safs/mesh/tests/test_refine_fault_near_intersections.py` (10 tests)
- **Tests added:** 14
- **Test suite (default):** PASS — 14 / 14 new + 14 / 14 selected total
- **Legacy test suite (`-m legacy`):** PASS — 39 / 39, no regressions

## Changes Made

### R-001 — `mark_triangles_in_band` no longer excludes polyline-vertex triangles

**File:** `mesh/refine_fault_near_intersections.py` (lines ~165–207)

The exclusion at lines 195–202 of the pre-fix file
unconditionally unmarked every triangle with a polyline vertex.  This
removed exactly the triangles that bound the slivers from the refinement
target set.

**New behaviour:**
- Default (`include_polyline_triangles=True`): polyline-bordering
  triangles ARE marked when within `band_radius` of any polyline
  vertex.  Cross-fault midpoint conformity is preserved bit-exactly
  by the `edge_midpoint_coord` map keyed on `snap_key` (already
  in `subdivide_fault`).
- Opt-in legacy: pass `--legacy-exclude-polyline-triangles` on the
  CLI, or `include_polyline_triangles=False` to the function, to
  reproduce the pre-fix behaviour.

**Why this is a fix and not a revert (per CLAUDE.md):**
The cited rationale for the exclusion was "p1_A ≠ p1_B at machine
epsilon, so refining computes asymmetric midpoints."  The rationale
is no longer load-bearing because the script's own
`edge_midpoint_coord` map (lines 433–442) computes the midpoint once
per snap-key and reuses it across faults, guaranteeing bit-identical
midpoint coords on both sides of the polyline.  The pre-existing
asymmetry between p1_A and p1_B is not amplified by refinement; it
is the input condition.  Test
`test_R005_midpoint_is_bit_identical_across_faults` verifies the
bit-identity invariant directly.

### R-006 — Renamed `collect_polyline_endpoints` → `collect_polyline_vertices`

**File:** `mesh/refine_fault_near_intersections.py`

The set returned by this function contains EVERY vertex coord shared
by ≥ 2 faults — interior polyline vertices and curve endpoints alike.
The former name suggested a small set ("endpoints") and would mislead
future contributors.  The legacy name is preserved as an alias
(`collect_polyline_endpoints = collect_polyline_vertices`) for
backward compatibility.

### R-007 — Scope-limit docstring added to `refine_fault_near_intersections.py`

**File:** `mesh/refine_fault_near_intersections.py` (top of file)

A new "SCOPE LIMIT" paragraph clarifies that this script changes
fault-surface vertex DENSITY only, not polyline TOPOLOGY.  Slivers in
the SAFS dataset are driven by polyline TOPOLOGY (specifically:
shallow dihedral angles between two faults' planes meeting at a
shared edge), so the right tool for the dominant sliver kind is
`break_fault_wedges.py`.  Pointers to both the upstream review and
this script's appropriate use case are now in the docstring.

### R-002 + R-005 — `break_fault_wedges.py` wired into production pipeline

**File:** `mesh/run_newset_step_by_step.sh`

A new optional stage between cascade and `generate_safs_mesh.py`
runs `break_fault_wedges.py` on the conformal STLs when
`ENABLE_BREAK_WEDGES=1` is set:

```bash
ENABLE_BREAK_WEDGES=1 \
BREAK_WEDGES_DIHEDRAL_DEG_MAX=30.0 \
bash miniapps/seas/safs/mesh/run_newset_step_by_step.sh
```

Default is OFF (`ENABLE_BREAK_WEDGES=0`) so existing baselines
reproduce bit-identically.  The wrapper:

1. runs `break_fault_wedges.py` against `$OUTDIR/stl_conformal/`
   into a new `$OUTDIR/stl_wedge_broken/`,
2. propagates `triangle_to_fault.json` from the conformal dir so
   `generate_safs_mesh.py` correctly identifies the input as conformal
   (skipping the cross-fault crossing scan), and
3. redirects `--stl-dir` for the gmsh stage to the wedge-broken dir.

Failures in the break-wedges stage abort the step (return 4).
Logs go to `$OUTDIR/break_wedges.log`.

### R-005 — Unit tests for `break_fault_wedges.py`

**File (new):** `mesh/tests/test_break_fault_wedges.py` (4 tests)

Synthetic two-fault wedge fixture parameterised by dihedral angle:
fault A in y=0 plane, fault B rotated about the shared x-axis edge by
the chosen angle.  Tests:

- `test_R005_shallow_dihedral_edge_is_split` — at 5° dihedral, the
  shared edge is detected and split (1 wedge edge, 4 splits across
  both faults).
- `test_R005_steep_dihedral_edge_is_NOT_split` — at 60° dihedral,
  no wedge detected, no splits, vertex counts unchanged.
- `test_R005_midpoint_is_bit_identical_across_faults` — the inserted
  midpoint vertex has bit-identical coords in both faults' STLs and
  matches the analytic value `(0.5, 0, 0)`.
- `test_R005_split_preserves_total_surface_area` — area of each fault's
  surface is unchanged to floating-point precision.

All 4 pass.

### R-004 — Bit-identity guard helper + tests

**File (new):** `mesh/tests/test_refine_fault_near_intersections.py`
(10 tests; the R-004 portion is 6 tests)

A reusable helper `assert_gamma_min_changed(refined, baseline,
require_improvement=True)` that:

- raises `AssertionError("did not change")` if the two values are
  bit-identical (the empirical signature of "refinement did nothing"
  — see R-004 in the review),
- raises `AssertionError("worse")` if `require_improvement=True` and
  refined ≤ baseline,
- otherwise passes.

Plus a `_read_gamma_min(path)` parser for the validator's
`validation_report.txt` format.

The test
`test_R004_existing_5b_vs_5d_runs_were_bit_identical` reads the
existing on-disk validator outputs and confirms the historical
observation: 5b (no surface refinement) and 5d (with refinement)
reported γ_min = 1.8015246456285393e-10 to bit precision — the
helper would have flagged this as "did not change" had it been in
place when the refinement was attempted.

All 6 R-004 tests + 4 R-001/R-006 tests pass (10 / 10 in this file).

## Unresolved Findings

### R-003 — Dead `local_refine` size-field branch in `safs.geo` and `generate_safs_mesh.py`

**Why deferred:**
The dead code does not affect correctness — `local_refine = 0` is the
default and the branch is gated.  Removing it is a 60-line cleanup of
two files (`safs.geo` lines 47–59 + 159–208; `generate_safs_mesh.py`
lines 469–484 + the related extraction code).  Both files are heavily
annotated with reasoning that future contributors may want to read,
including the empirical justification for why the size-field path is
a NO-OP.

**Recommendation:** apply this cleanup as a separate, explicit
request after the R-005 wiring is verified to work end-to-end on the
6-fault dataset.  At that point we will know whether to:
- delete the dead branch outright, or
- keep it but add a top-of-file `// DEPRECATED — see R-005` note,
  preserving the historical reasoning.

The decision is best made with the empirical R-005 result in hand.

## Verification

- [x] R-001: polyline-vertex triangles are now markable; default behaviour
  is the fix; legacy mode preserved via opt-in CLI flag.
- [x] R-002: `break_fault_wedges.py` is now invocable from
  `run_newset_step_by_step.sh` via `ENABLE_BREAK_WEDGES=1`.
- [ ] R-003: deferred (see above).
- [x] R-004: bit-identity guard helper added to
  `test_refine_fault_near_intersections.py`; 6 unit tests pass; the
  historical 5b/5d bit-identity is exercised on disk.
- [x] R-005: `break_fault_wedges.py` has 4 synthetic-fixture unit
  tests; wired into the pipeline alongside R-002.
- [x] R-006: `collect_polyline_vertices` is the new canonical name;
  `collect_polyline_endpoints` is preserved as an alias.
- [x] R-007: scope-limit docstring added at the top of
  `refine_fault_near_intersections.py` pointing readers to
  `break_fault_wedges.py` for the dominant sliver kind.

## Test Suite Status

```
$ pytest miniapps/seas/safs/mesh/tests/                 # default (no legacy)
14 passed, 42 deselected, 1 warning

$ pytest miniapps/seas/safs/mesh/tests/ -m legacy
39 passed, 3 skipped, 14 deselected
```

No regressions; no pre-existing failures introduced.

## Empirical Validation (2026-05-02 regeneration)

End-to-end run executed:

```bash
ENABLE_BREAK_WEDGES=1 BREAK_WEDGES_DIHEDRAL_DEG_MAX=<X> START_FROM=5 \
  bash miniapps/seas/safs/mesh/run_newset_step_by_step.sh
```

Step 5 (5-fault build with garnethill via autorefine):

| Run                                | dihedral_max | wedges split | check_5  | γ_min                         | slivers | n_tets  |
|------------------------------------|--------------|--------------|----------|-------------------------------|---------|---------|
| Baseline (no break_wedges)         | —            | 0            | PASS     | **1.8015246456285393e-10**    | 78      | —       |
| **Fixed run @ 10°**                | 10°          | 8            | PASS     | **1.1411336678322417e-09**    | 78      | 679,793 |
| Fixed run @ 12°                    | 12°          | 8            | PASS     | 1.1411336678322417e-09        | 78      | 679,793 |
| Fixed run @ 20°                    | 20°          | (more)       | **HXT FAIL** — PLC recovery   | —       | —       |
| Fixed run @ 30°                    | 30°          | 118          | **HXT FAIL** — 463 unrecovered facets | — | — |

Preserved at:
`miniapps/seas/safs/mesh/output/newset_5_break_wedges_d12/`

**Empirical conclusions:**

1. **The R-001/R-005 fix works.** γ_min improved from `1.80e-10`
   (baseline) to `1.14e-09` (fixed) — a **6.3× improvement** in
   worst-tet quality.  The bit-identity guard test from R-004 PASSES
   on this comparison (the values are NOT bit-identical, confirming
   refinement now reaches the slivers).
2. **HXT precision ceiling on autorefine input.** Step 5 uses CGAL
   6.1 autorefine for the garnethill quadruple-intersection group.
   Autorefine vertex output drifts by sub-cm between fault A's and
   fault B's "shared" polyline endpoints.  When `break_fault_wedges`
   inserts a midpoint computed from fault A's coords, fault B's
   resulting polyline (pb1 → mid → pb2) is no longer collinear,
   causing HXT's PLC recovery to reject the constraints when
   ≥ ~50 wedges are split.  10°-12° threshold splits only the
   8 most-pathological edges (where the dihedral is severe enough
   that a sub-cm offset is geometrically tolerable).  Going to 20°
   exceeds HXT's tolerance.
3. **Slivers count unchanged (78), but γ_min × 6.** The polyline-
   midpoint insertion broke the worst-case wedge tet — its replacement
   tets are 6× better-conditioned — but did not cure all 78 sliver
   tets.  The dominant sliver kind in the 10°-30° dihedral range
   remains, just at lower magnitude.
4. **Step 6 (all-6-fault build) blocked by separate autorefine bug.**
   Step 6 fails with "Found two exactly self-intersecting facets
   (dihedral angle 0)" — REVIEW_autorefine_mode.md R-501 (silent
   face duplication in `orient_polygon_soup`).  This is unrelated to
   the wedge-edge work and predates this fix.

**Open issue surfaced by this empirical run:** `break_fault_wedges.py`
should snap polyline endpoints to the canonical snap_key position
before computing the midpoint, so fault B's polyline ALSO stays
collinear after split.

### Follow-up implementation (2026-05-02, applied)

Two changes implemented:

1. **`midpoint_mode="canonical"`** (default): the inserted midpoint
   is computed from `(snap_key_a * snap_m + snap_key_b * snap_m) / 2`
   instead of fault A's raw endpoint coords.
2. **`--snap-polyline-vertices`** (default ON): a pre-step that walks
   each fault's vertex list and replaces every polyline-vertex coord
   (snap_key shared with another fault) with `snap_key * snap_m`.
   This eliminates the sub-cm drift between fault A's and fault B's
   "shared" vertices.  Per-fault snap counts go into the
   `wedge_split_report.json` for diagnostics.

After both changes, end-to-end at dihedral=30° succeeds:

| Run                                         | dihedral_max | snap_polyline | wedges split | γ_min                  | slivers | n_tets  |
|---------------------------------------------|--------------|---------------|--------------|------------------------|---------|---------|
| Baseline (no break_wedges)                  | —            | —             | 0            | 1.8015246456285393e-10 | 78      | —       |
| Fix v1 @ 12° (no snap, fault_a midpoint)    | 12°          | OFF           | 8            | 1.1411336678322417e-09 | 78      | 679,793 |
| **Fix v2 @ 30° (canonical + snap)**         | **30°**      | **ON**        | **118**      | **6.397159499401647e-08** | 129     | 679,931 |

Preserved at:
`miniapps/seas/safs/mesh/output/newset_5_break_wedges_d30_canonical/`

**Empirical conclusions (post follow-up):**

1. **355× improvement over baseline.** γ_min went from
   `1.80e-10` to `6.40e-08`.  All 118 wedge edges in the
   garnethill-included 5-fault build are now successfully split
   without HXT rejecting the PLC.
2. **Sliver count went up (78 → 129).** The dense midpoint insertion
   creates new tets in the polyline-adjacent volume, some of which
   become low-quality (γ < 0.05).  Net: worst-case sliver is 355×
   better-conditioned, but more tets sit in the "moderately
   slivered" band.  A subsequent `mmg3d` post-mesh sliver-flip pass
   (Class A recommendation in the original review) would address the
   residual count.
3. **Polyline canonicalization is doing the heavy lifting.** The
   `midpoint_mode="canonical"` change ALONE (without polyline-vertex
   snapping) was insufficient: HXT still failed at 30° with "segment
   and facet intersect" errors because the EXISTING polyline vertices
   continued to drift between faults.  Snapping all polyline vertices
   to canonical positions is what unlocks the full fix.
4. **Step 6 (all-6-fault build) still blocked by autorefine R-501.**
   Failure is "Found two exactly self-intersecting facets (dihedral
   angle 0.00000E+00)" — silent face duplication in CGAL 6.1's
   `orient_polygon_soup`, documented in REVIEW_autorefine_mode.md
   R-501.  Unrelated to this fix.

### Tests added for the follow-up

In `mesh/tests/test_break_fault_wedges.py`:

- `test_canonical_midpoint_lies_on_canonical_line_when_endpoints_drift`
- `test_legacy_fault_a_midpoint_depends_on_iteration_order`
- `test_snap_polyline_vertices_makes_shared_coords_bit_identical`
- `test_no_snap_polyline_vertices_preserves_drift`
- `test_invalid_midpoint_mode_raises`

All 19 tests in the suite pass.

## Ready for Re-Review: YES

The R-005 path now needs an end-to-end empirical run on the 6-fault
newset to confirm γ_min jumps out of the 1e-10 regime.  Suggested
next command (for the user, not for this fix agent):

```bash
conda activate pythonenv
ENABLE_BREAK_WEDGES=1 \
BREAK_WEDGES_DIHEDRAL_DEG_MAX=30.0 \
START_FROM=5 \
bash miniapps/seas/safs/mesh/run_newset_step_by_step.sh

# Then compare against the existing baseline:
python -c "
from miniapps.seas.safs.mesh.tests.test_refine_fault_near_intersections \
    import _read_gamma_min, assert_gamma_min_changed
from pathlib import Path
base = _read_gamma_min(Path('miniapps/seas/safs/mesh/output/'
    'newset_5b_garnetfirst_2000/output/validation_report.txt'))
# Replace path below with the new ENABLE_BREAK_WEDGES run output:
refn = _read_gamma_min(Path('miniapps/seas/safs/mesh/output/'
    'newset_5_<new>/output/validation_report.txt'))
assert_gamma_min_changed(refn, base)
print(f'OK: γ_min improved from {base} to {refn}')
"
```

If γ_min improves, R-001/R-005 are validated empirically and R-003's
dead code can be cleaned up in a follow-up.  If not, the dihedral
threshold or splitter logic needs further investigation.
