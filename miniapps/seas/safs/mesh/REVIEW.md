# Code Review: SAFS Multi-Fault Conforming Mesh Pipeline (2026-04-29)

## Review Scope

- **Plan documents**:
  - `miniapps/seas/safs/PLAN.md` (Phases 3-5: STL import → fault intersection
    → 3D meshing)
  - `miniapps/seas/safs/PLAN_domain.md` (size-field rules)
- **Files reviewed**:
  - `miniapps/seas/safs/mesh/generate_safs_mesh.py` (driver, vertex weld)
  - `miniapps/seas/safs/mesh/safs.geo` (Gmsh model: Box + Embed + size field)
  - `miniapps/seas/safs/mesh/audit_ts_quality.py` (`find_overlap_pairs` —
    the only geometric-intersection primitive that survives in the repo)
  - `miniapps/seas/safs/mesh/ts_to_stl.py` (per-fault cleanup, no awareness
    of other faults)
  - `miniapps/seas/safs/mesh/validate_msh.py` (check 5 = SEAS interface
    invariant)
  - `miniapps/seas/safs/mesh/output/all8_2000m/output/validation_report.txt`
    (the all-8 failure record)
- **Domain context consulted**: `seas-mfem/CLAUDE.md`,
  `miniapps/seas/CLAUDE.md`, `miniapps/seas/safs/mesh/README.md` (the
  "Caveats" and "Roadmap" sections explicitly enumerate the known gap).

The user's prompt frames the review as a *design audit* of the existing
multi-fault path with the goal of producing a doable plan for two
intersecting traces. Findings below are the gating items the plan must
fix; the plan itself is appended after the findings.

---

## Smoking-Gun Evidence

`output/all8_2000m/output/validation_report.txt` (line 17–21):

```
[FAIL] 5_internal_interface: 4917 fault triangle(s) have wrong tet adjacency
        n_internal_with_2_tets: 0
        n_trace_with_1_tet:    0
        n_offending:           4917
```

Every fault triangle in the all-8 mesh is *detached* from the tet mesh —
not 99% wrong, 100% wrong. This is not "slivers near intersection knots";
the entire fault surface is a phantom mesh entity that no tet shares a
face with. The tet mesh was generated as if the fault triangulation
weren't there.

---

## Findings

### [R-001] [CRITICAL] [generate_safs_mesh.py:_combine_stls] — Vertex weld is a no-op for true cross-fault intersections

**Category:** ASSUMPTION / DEVIATION

**Description:**
`_combine_stls` (lines 169–271) snaps every per-fault STL vertex to a
1 cm grid (`snap_m=0.01`) and dedups via `np.unique`. The docstring
(line 174–179) claims this welds "shared CFM vertices at fault-fault
junctions" so that gmsh's STL Merge does not duplicate them. The premise
is wrong: CFM fault surfaces are independently sampled and triangulated
*point clouds*, with no a-priori shared vertices between faults. Where
fault A and fault B cross in 3-D, neither carries a vertex on the
intersection curve unless the CFM authors happened to place one there
(they didn't, per the cross-fault crossing pair counts in
`README.md`'s caveat 1).

The 1 cm dedup therefore catches *only* accidentally-coincident vertices,
not actual geometric intersections. The combined STL shipped to gmsh is
still topologically two non-conforming surfaces that pass through each
other.

**Trigger:** Any `--include-fault` set with two or more faults that
cross in 3-D (10 of 28 CFM pairs at 2000 m per the README).

**Actual behavior:** Combined STL is welded only at coincidental matches;
intersection curves remain unrepresented; gmsh's bulk mesher silently
ignores fault triangles → all 4917 fault triangles in the all-8 mesh
have 0 adjacent tets.

**Expected behavior:** Either (a) the function actually computes
fault-fault intersection segments and inserts them into both
triangulations, or (b) the function refuses to operate on inputs that
have unprocessed crossings and points the caller at the new
intersection-resolving step.

**Suggested fix:** Replace `_combine_stls` with a `_combine_stls_conformal`
that delegates to a new `fault_intersect.py` (see plan below). For now,
the docstring at lines 174–179 must stop claiming this function fixes
multi-fault conformity:

```diff
-    """Combine per-fault ASCII STLs into a single ASCII STL with vertex
-    deduplication at `snap_m` precision so shared CFM vertices at
-    fault-fault junctions become a single point in the combined mesh.
-
-    Without this, Gmsh's STL Merge sees two faults each carrying a copy
-    of the shared vertex; HXT then "filters" the duplicates and the
-    constraint topology breaks.
-    """
+    """Concatenate per-fault STLs and snap-dedup at `snap_m` to remove
+    coincidental vertex duplicates only.  This DOES NOT make the
+    combined surface conformal across fault-fault crossings — true
+    intersections need a geometric splitter (see fault_intersect.py
+    Phase 2 of PLAN_multifault_intersections.md).  Use this only after
+    each pair of input STLs has been verified non-crossing (e.g., the
+    disjoint subsets in run_subsets_2000m.sh)."""
```

**Test case:**
```python
def test_R001_two_crossing_planes_remain_non_conforming():
    """Two perpendicular square planes that cross at x=0; combine_stls
    must NOT silently produce a conforming surface, because no shared
    edge vertices exist along the intersection line."""
    import tempfile, json, numpy as np
    from pathlib import Path
    # Plane A: y=0, x in [-1,1], z in [-1,0], 2 triangles
    # Plane B: x=0, y in [-1,1], z in [-1,0], 2 triangles
    # The two planes cross along the segment {x=0, y=0, z in [-1,0]}.
    # Neither input has any vertex on that segment.
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        _write_test_stl(tmp / "A.stl", "A", _square_triangles(axis="y"))
        _write_test_stl(tmp / "B.stl", "B", _square_triangles(axis="x"))
        meta = _combine_stls(tmp / "out.stl", ["A", "B"], tmp, snap_m=0.01)
        # The combined surface must NOT contain a shared edge along the
        # intersection line — the function should not pretend it does.
        v, t = _read_ascii_stl(tmp / "out.stl")
        edge_count = _count_shared_edges_between(v, t,
                                                  fault_a_range=meta["fault_ranges"]["A"],
                                                  fault_b_range=meta["fault_ranges"]["B"])
        assert edge_count == 0, (
            "vertex-weld must not invent shared edges; got "
            f"{edge_count} (would mean false conformity)")
```

---

### [R-002] [CRITICAL] [safs.geo:Algorithm3D=4 + Surface{...} In Volume{1}] — Frontal-Delaunay does not enforce internal-face conformity for non-conforming embedded constraint surfaces

**Category:** DEVIATION

**Description:**
Line 118 sets `Mesh.Algorithm3D = 4` (Frontal-Delaunay) with an
explanatory comment that HXT (=10) rejects multi-fault inputs and FD-3D
"is more permissive about non-conformal Steiner constraints and produces
a valid multi-fault mesh." The "valid" claim is wrong. FD-3D produces a
mesh in which the bulk tets do not even attempt to share faces with the
embedded fault triangulation when those triangles cross another embedded
surface. The validator's check 5 catches this every time
(`n_internal_with_2_tets: 0` for all 4917 fault tris in all-8).

PLAN.md Phase 5 acceptance criterion explicitly requires "every fault
triangle has 2 adjacent tets" (the SEAS DG friction solver depends on
this). The current `safs.geo` violates this acceptance criterion silently
in multi-fault mode and only the validator stops the bad mesh from
shipping.

**Trigger:** Multi-fault `--include-fault` list where at least two
included faults intersect in 3-D.

**Actual behavior:** Mesh generation completes; the `.msh` file contains
4917 fault triangles that the bulk tet mesh ignores. Validator check 5
exits non-zero. Caller has no signal earlier in the pipeline.

**Expected behavior:** Either (a) `generate_safs_mesh.py` runs an
intersection scan first and rejects any input set with unresolved
crossings (with a clear "run conformalize_faults.py first" message), or
(b) the pipeline includes a conformalization step that makes the input
crossing-free before reaching gmsh.

**Suggested fix (gating change, until R-001's intersection pipeline
exists):** Add an early-fail in `generate_safs_mesh.py` when the included
fault set has any cross-fault crossings:

```diff
     # 5. Compute domain box
     box = _compute_box(bbox, args.buf_x, args.buf_y,
                        args.depth, args.round_to)
     _validate_clearance(box, bbox, args.buf_x, args.buf_y, args.depth)

+    # 5b. Reject multi-fault input that has unresolved 3-D crossings.
+    if len(included) > 1:
+        from fault_intersect import scan_cross_fault_crossings
+        crossings = scan_cross_fault_crossings(args.stl_dir, included)
+        if crossings:
+            n_pairs = sum(len(v) for v in crossings.values())
+            details = "; ".join(
+                f"{a}×{b}: {len(crossings[(a,b)])} pairs"
+                for (a, b) in crossings
+            )
+            raise SystemExit(
+                f"refusing to mesh: {len(crossings)} cross-fault "
+                f"intersection cluster(s), {n_pairs} crossing triangle "
+                f"pair(s):\n  {details}\n"
+                f"run conformalize_faults.py first to insert the "
+                f"intersection curves into both triangulations, or "
+                f"split into disjoint subsets (see "
+                f"run_subsets_2000m.sh).")
+
     # 6. Resolve geo + includes paths
```

**Test case:**
```python
def test_R002_multi_fault_with_crossings_fails_loudly():
    """generate_safs_mesh.py must refuse to run on the all-8 set
    until conformalization exists."""
    import subprocess, sys
    # Use the existing all-8 STL dir, expect non-zero exit and a
    # "refusing to mesh" message.
    out = subprocess.run(
        [sys.executable, "generate_safs_mesh.py",
         "--stl-dir", "output/all8_2000m/stl",
         "--bbox-json", "output/all8_2000m/bbox.json",
         "--transform-json", "output/all8_2000m/transform.json",
         "-o", "/tmp/should_not_be_written.msh"],
        capture_output=True, text=True)
    assert out.returncode != 0
    assert "refusing to mesh" in (out.stdout + out.stderr)
    assert "intersection" in (out.stdout + out.stderr)
```

---

### [R-003] [CRITICAL] [missing module] — No fault-fault intersection geometry exists in the repo

**Category:** DEVIATION (PLAN.md Phase 4 not implemented)

**Description:**
PLAN.md Phase 4 specifies "Fault–Fault Intersection Resolution" with
either OCC `BooleanFragments` (primary) or HXT discrete embedding
(fallback). Neither is operative for the multi-fault use case:

- The OCC path requires `ClassifySurfaces` + `CreateGeometry` to promote
  STL meshes to NURBS. README.md line 5–7 documents that this fails on
  CFM .ts triangulations in gmsh 4.13: "Gmsh's CreateGeometry /
  BooleanFragments do not promote discrete STL surfaces to OCC NURBS".
- The HXT discrete-embedding fallback rejects multi-fault input with PLC
  errors (README.md "Caveats" §1).

Therefore the only multi-fault path that runs is the
Algorithm3D=4 / non-conforming-embed combination flagged in R-002.

The README's "Path A — Triangle-splitting preprocessor" mentions four
geometric primitives (`tri_tri_intersect_3d`, `edge_pierces_tri`,
`_point_to_tri_distance`, `seg_seg_cross_2d`) that "would have to be
rewritten" — they are not in the repo today. `find_overlap_pairs` exists
in `audit_ts_quality.py` but only handles within-fault coplanar overlaps
(triangles sharing an edge in the *same* surface), not cross-fault
3-D crossings.

**Trigger:** Any multi-fault build with at least one crossing pair.

**Actual behavior:** The pipeline silently produces a non-conforming
mesh.

**Expected behavior:** A `fault_intersect.py` module that, given two (or
more) per-fault triangulations, returns a list of polylines along which
the surfaces intersect, with the index of every triangle each polyline
crosses. A `conformalize_faults.py` driver that consumes those polylines
and re-triangulates each crossed triangle via constrained Delaunay so
the intersection polylines become exact edge sequences in both surfaces.

**Suggested fix:** This is the substance of the new plan
(`PLAN_multifault_intersections.md`, Phase 1 + 2). The fix-agent should
implement those phases in order; do not attempt to use OCC
`BooleanFragments` or HXT — both paths are dead-ends for CFM input
per README §Caveats.

**Test case:** see Phase 1 acceptance criteria in the plan below.

---

### [R-004] [MODERATE] [safs.geo:Field[2] Threshold lines 100-105] — Size field has no uniform near-fault region; ramps immediately

**Category:** DEVIATION (deviates from user requirement, not from
written plan — PLAN_domain.md does not address this either, so this is
also a plan gap)

**Description:**
The threshold-distance field is configured with `DistMin = 0`,
`SizeMin = res_f`, `DistMax = ramp_dist`, `SizeMax = res_ff`. This means
*the moment* a tet leaves the fault triangle, the target size starts
ramping toward `res_ff`. There is no "uniform near-fault tube" zone of
constant `res_f`, even though the user's expected behavior (and most
SEAS solver requirements for resolving the cohesive zone) is uniform
elements within some tube radius around the fault, then graded coarsening
beyond it.

The user explicitly asked: "I would like to keep near the fault region
mesh size to be uniform … this means we might need to define subdomain
partitions that surrounds the fault traces, within which the mesh is
uniform, and coarsen as move away from the fault traces."

**Trigger:** Any mesh build (single-fault smoke test included).

**Actual behavior:** Tets one element away from the fault are already
~`res_f * 1.4` (linear ramp). The expected uniform near-fault layer
does not exist.

**Expected behavior:** Add a `tube_radius` parameter; size = `res_f`
for distance ≤ `tube_radius`, then ramp from `tube_radius` to
`tube_radius + ramp_dist`. With Gmsh's Threshold field the change is a
one-liner: set `DistMin = tube_radius` (not 0).

**Suggested fix:**

```diff
 // --- Mesh sizing (driver overrides via -setnumber) ---
 DefineConstant[ res_f     = { 1000, Name "On-fault target edge length (m)"} ];
 DefineConstant[ res_ff    = {20000, Name "Far-field edge length (m)"} ];
+DefineConstant[ tube_radius = { 5000, Name "Uniform near-fault tube radius (m)"} ];
 DefineConstant[ ramp_dist = {30000, Name "Linear ramp width (m)"} ];
 ...
 Field[2] = Threshold;
 Field[2].InField  = 1;
 Field[2].SizeMin  = res_f;
 Field[2].SizeMax  = res_ff;
-Field[2].DistMin  = 0;
-Field[2].DistMax  = ramp_dist;
+Field[2].DistMin  = tube_radius;
+Field[2].DistMax  = tube_radius + ramp_dist;
```

And expose it in `generate_safs_mesh.py`:

```diff
+    parser.add_argument("--tube-radius", type=float, default=5_000.0,
+                        metavar="METRES",
+                        help="Distance from fault inside which the mesh "
+                             "is uniform at --res-f.  Default: 5 km "
+                             "(= 5 * default res_f).")
     ...
     cmd = [
         "gmsh", "-3", str(run_geo),
         ...
         "-setnumber", "res_f",       repr(args.res_f),
         "-setnumber", "res_ff",      repr(args.res_ff),
+        "-setnumber", "tube_radius", repr(args.tube_radius),
         "-setnumber", "ramp_dist",   repr(args.ramp_dist),
         "-o", str(args.output),
     ]
```

The validator must learn about the tube as well: extend
`check_6_fault_edge_length` or add `check_11_tube_uniformity` that
samples tets within `tube_radius` of the fault and checks
`std(edge_len) / mean(edge_len) < 0.2`.

**Test case:**
```python
def test_R004_tube_uniformity_within_radius():
    """A 2-fault mesh with --tube-radius 4000 must have all tets whose
    centroid is within 4000 m of any fault triangle exhibit a tight
    edge-length distribution (mean ≈ res_f, p95 < 1.3 * res_f)."""
    import meshio, numpy as np
    from scipy.spatial import cKDTree
    msh = meshio.read("output/two_faults_test/safs_two_faults.msh")
    pts = msh.points
    tets = msh.cells_dict["tetra"]
    fault_tris_idx = ...  # tag-100 from cell_data_dict
    fault_centroids = pts[fault_tris_idx].mean(axis=1)
    tree = cKDTree(fault_centroids)
    tet_centroids = pts[tets].mean(axis=1)
    d, _ = tree.query(tet_centroids, k=1)
    near = tets[d < 4000.0]
    a, b, c, d_ = (pts[near[:, k]] for k in range(4))
    edges = np.concatenate([
        np.linalg.norm(b-a, axis=1), np.linalg.norm(c-a, axis=1),
        np.linalg.norm(d_-a, axis=1), np.linalg.norm(c-b, axis=1),
        np.linalg.norm(d_-b, axis=1), np.linalg.norm(d_-c, axis=1)])
    assert 800 < edges.mean() < 1300, edges.mean()       # res_f=1000
    assert np.percentile(edges, 95) < 1300, np.percentile(edges, 95)
```

---

### [R-005] [MODERATE] [missing test fixture] — No 2-fault crossing test exists

**Category:** EDGE_CASE / TESTING

**Description:**
The test inventory:
- `run_smoke_millcreek.sh`: 1 fault, no crossings.
- `run_subsets_2000m.sh`: 4 *disjoint* subsets, no crossings.
- `run_all8_2000m.sh`: 8 faults, ~10 crossing pairs, fails check 5.

There is no minimal `run_two_crossing_2000m.sh` (or equivalent) exercising
the simplest non-trivial multi-fault configuration. The user's prompt
specifically calls for "just two fault traces (that are intersecting at
some places)" as the proving ground for the new approach.

**Trigger:** Any work on the conformalization pipeline; without a
dedicated 2-fault target we cannot decouple intersection-handling bugs
from cascading multi-pair bugs.

**Actual behavior:** The conformalization fix has to be debugged against
the all-8 input or against ad-hoc local tests.

**Expected behavior:** A 2-fault wrapper using the smallest CFM crossing
pair (Pinto × Mill Creek = 42 pairs, per README) — but Pinto is
CFM-defective standalone, so use **Mill Creek × SBMT-SAF (163 pairs)**
or **Mission Creek SBMT × SAF (135 pairs)**. Both members are individually
SEAS-grade. Pair 1 has the highest pair count (good signal); pair 2 has
geometric simplicity (Mission Creek SBMT is 449 tris, smaller).

**Suggested fix:** Create `run_two_crossing_2000m.sh` that drives the
new pipeline through `audit → ts_to_stl → conformalize_faults →
generate_safs_mesh → write_fault_provenance → validate_msh`. Acceptance
is `10/10 checks passed` *with* multi-fault enabled, including check 5.

```bash
#!/usr/bin/env bash
# Minimal multi-fault crossing test: Mill Creek strand × SBMT San Andreas.
# 163 known crossing triangle pairs at 2000m → exercises the
# conformalization pipeline end-to-end while staying small enough to
# debug interactively.
set -euo pipefail
cd "$(dirname "$0")"

OUTDIR="output/two_crossing_2000m"
RES=2000; RES_F=2000; RES_FF=25000; RAMP=20000; TUBE=4000
CFM_DIR="${SAFS_CFM_DIR:-$HOME/Documents/.../CFM_data}"

rm -rf "$OUTDIR"; mkdir -p "$OUTDIR/stl" "$OUTDIR/output"

INCL=(--include-fault safs_sbmt_millcreek --include-fault safs_sbmt_saf)

python audit_ts_quality.py --cfm-dir "$CFM_DIR" --res $RES "${INCL[@]}" \
    --out "$OUTDIR/cfm_audit_${RES}m.csv"

python ts_to_stl.py --cfm-dir "$CFM_DIR" --res $RES "${INCL[@]}" \
    --out-dir "$OUTDIR/stl" \
    --transform-json "$OUTDIR/transform.json" \
    --bbox-json     "$OUTDIR/bbox.json" \
    --log-csv       "$OUTDIR/cleanup_log_${RES}m.csv" \
    --free-surface-clearance 100 --repair-overlaps

# NEW: insert intersection curves before meshing.
python conformalize_faults.py "${INCL[@]}" \
    --in-stl-dir  "$OUTDIR/stl" \
    --out-stl-dir "$OUTDIR/stl_conformal" \
    --report      "$OUTDIR/intersection_report.json"

python generate_safs_mesh.py "${INCL[@]}" \
    --stl-dir "$OUTDIR/stl_conformal" \
    --bbox-json "$OUTDIR/bbox.json" \
    --transform-json "$OUTDIR/transform.json" \
    --res-f $RES_F --res-ff $RES_FF \
    --tube-radius $TUBE --ramp-dist $RAMP \
    -o "$OUTDIR/output/safs_two_crossing_${RES}m.msh"

python write_fault_provenance.py "${INCL[@]}" \
    --msh "$OUTDIR/output/safs_two_crossing_${RES}m.msh" \
    --stl-dir "$OUTDIR/stl_conformal" \
    --transform-json "$OUTDIR/transform.json" \
    --out "$OUTDIR/output/fault_provenance.json"

python validate_msh.py \
    --msh "$OUTDIR/output/safs_two_crossing_${RES}m.msh" \
    --transform-json "$OUTDIR/transform.json" \
    --bbox-json     "$OUTDIR/bbox.json" \
    --provenance-json "$OUTDIR/output/fault_provenance.json" \
    --report "$OUTDIR/output/validation_report.txt"
```

**Test case:** `bash run_two_crossing_2000m.sh` exits 0 and the report
shows `10/10 checks passed`.

---

### [R-006] [LOW] [generate_safs_mesh.py:_combine_stls] — Per-fault provenance is destroyed by the merge

**Category:** QUALITY / ROBUSTNESS

**Description:**
`_combine_stls` writes a single `solid SAFS:combined` STL with no
per-fault tagging in the file. Downstream, `write_fault_provenance.py`
recovers per-fault attribution via KDTree against the original
per-fault STL centroids. When the new conformalization pipeline splits
triangles, the post-split centroids of the new sub-triangles are still
on the original fault surface, so KDTree attribution still works — but
this is a fragile coincidence, not a guarantee. A triangle split exactly
at an intersection curve has a centroid arbitrarily close to triangles
of the *other* fault; KDTree can mis-attribute under round-off.

**Trigger:** Conformalization produces a sub-triangle whose centroid is
within `min_inter_fault_distance` of an unrelated fault.

**Actual behavior:** Risk of mis-attribution after splits at intersection
curves; nothing currently guards against it.

**Expected behavior:** Carry per-fault IDs through the conformalization
step (each post-split sub-triangle inherits its parent's fault ID) and
write a per-triangle `fault_id` array next to the combined STL. The
`write_fault_provenance.py` step becomes a trivial passthrough instead
of a KDTree reconstruction.

**Suggested fix:** Defer to PLAN_multifault_intersections.md Phase 2,
which writes `<conformal>/triangle_to_fault.json` (per-triangle fault
short-name) alongside the conformal STL. This obsoletes the KDTree path
in `write_fault_provenance.py`.

**Test case:** When the new pipeline is in place, attribute every fault
triangle in the 2-crossing test mesh and confirm 100% of triangles
agree with the conformalizer's `triangle_to_fault.json` (no KDTree mis-
attribution at the intersection curve).

---

## Summary

- Critical issues: 3 (R-001, R-002, R-003)
- Moderate issues: 2 (R-004, R-005)
- Low issues: 1 (R-006)
- Plan compliance: **INCOMPLETE** — PLAN.md Phase 4 (fault-fault
  intersection resolution) is not implemented; the workaround in the
  current `safs.geo` produces meshes that fail Phase 5 acceptance
  criteria for multi-fault input.
- Verdict: **FAIL** — must implement Phase 1 (intersection geometry) +
  Phase 2 (constrained re-triangulation) of the plan below before any
  multi-fault SEAS run. R-004 (uniform tube) and R-005 (2-fault test)
  are required for the user's stated scientific workflow.

## Unreviewed Areas

- `write_fault_provenance.py` — only inspected at the interface level
  via R-006. Its KDTree implementation is correct for non-conformalized
  inputs and out of scope for this review.
- `safs_origin.py` — coordinate-frame module, unchanged since the smoke
  test passed. Not relevant to multi-fault conformity.
- `convert_msh.py` — output-format converter, post-mesh.
- `audit_ts_quality.py` parser path (only `find_overlap_pairs` was
  reviewed for relevance to the intersection problem).

---

# Recommended Plan: Two-Fault Conformalization + Tube-Uniform Sizing

This is the actionable plan the user asked for. It is split out into
`PLAN_multifault_intersections.md` (sibling of `PLAN.md`); the summary
below is what the fix-agent should follow. Each phase has a concrete
acceptance gate so we can stop and inspect after each phase.

## Phase 1 — Cross-fault intersection geometry (`fault_intersect.py`)

**Goal:** given two cleaned per-fault STLs, return the polylines along
which they intersect plus the per-triangle pierce list.

**Files to create:**
- `mesh/fault_intersect.py` — pure-Python module (numpy only; no gmsh,
  no OCC). Public API:
  ```python
  def tri_tri_intersect_3d(V_A, T_A, V_B, T_B, eps=1e-6) -> list[Segment]:
      """For each (i, j) where T_A[i] crosses T_B[j], return Segment(
          p0, p1, tri_a=i, tri_b=j, fault_a=name_a, fault_b=name_b)."""

  def chain_segments(segments, snap_m=0.01) -> list[Polyline]:
      """Connect end-to-end-coincident segments into open/closed
      polylines.  Polyline.points is dedup'd at snap_m."""

  def scan_cross_fault_crossings(stl_dir, included) -> dict:
      """One-shot: read every pair, return {(short_a, short_b):
      [Segment, ...]} for all crossing pairs.  Used by R-002 gate."""
  ```
- `mesh/tests/test_fault_intersect.py` — unit tests on synthetic inputs.

**Robustness:** orientation predicate via plain float64; if a candidate
crossing has |orient3d| < eps_relative, fall back to `gmpy2.mpfr` exact
arithmetic (optional; declare the dependency, do not vendor it).

**Acceptance:**
- Two perpendicular squares (90° crossing) → exactly one polyline,
  segment count = number of crossed triangles in the finer surface.
- Two near-coplanar surfaces (dihedral angle 1°) → no false positives
  (the predicate must declare them non-crossing because the
  intersection is a fat strip, not a curve, by construction the input
  faults are non-coplanar in CFM).
- `scan_cross_fault_crossings("output/all8_2000m/stl", all_8)` returns
  the same crossing-pair counts the README enumerates (Mill × SBMT-SAF
  ≈ 163, etc., to ±10%).

## Phase 2 — Constrained re-triangulation (`conformalize_faults.py`)

**Goal:** insert each intersection polyline as a fixed edge sequence in
both the surfaces it crosses, preserving fault provenance.

**Files to create:**
- `mesh/conformalize_faults.py` — driver. For each fault and each
  polyline that crosses it: project the polyline onto the fault (it
  already lies on it within `eps`); for every triangle the polyline
  enters: insert the polyline's pierce points as Steiner vertices and
  re-triangulate the triangle using **2-D constrained Delaunay** in the
  triangle's local plane (use the `triangle` Python wrapper around
  Shewchuk's library, already a `pythonenv` dependency, or `scipy`-only
  fallback). Carry the parent triangle's fault short-name to all
  children.

- Output: `<conformal>/safs_<short>.stl` — re-triangulated per-fault STLs
  with intersection edges baked in. Plus `triangle_to_fault.json` —
  list of `{short_name, n_triangles}` matching the post-split STLs.
  Plus `intersection_report.json` listing every polyline, the two
  faults it belongs to, the # of new vertices inserted, and the # of
  new sub-triangles per side.

**Edge cases the implementation must handle:**
- A polyline endpoint sits exactly on an existing edge of one fault
  (T-junction): both adjacent triangles must split, not just one.
- Three faults crossing at a single point (the San Gorgonio knot — even
  in the 2-fault test it does not occur, but the cascading logic must
  handle it: insert each polyline against the *current* triangulation,
  not the original).
- Numerical near-tangency: an intersection segment shorter than
  `0.05 * mean_edge_len` is dropped (it would create degenerate
  sub-triangles); the report records the drop.
- The polyline crosses a triangle near a vertex: snap the pierce point
  to the existing vertex when within `0.01 * edge_len` of it; do not
  insert a Steiner vertex.

**Acceptance:**
- Synthetic 90°-cross test: each square goes from 2 tris to 4 tris
  (one new vertex on each side at the centerline crossing).
- Mill Creek × SBMT-SAF at 2000m: post-split triangle count grows by
  the # of pierce points (= # of polyline vertices) in each fault, to
  within ±5%.
- `find_overlap_pairs` (the existing within-fault detector) still
  reports 0 same-side overlaps in each post-split fault (i.e., the
  re-triangulation produced a clean 2-manifold).
- Cross-fault intersection scan run on the post-split STLs returns
  zero crossing pairs (= true conformality at the polyline edges; HXT
  will accept the input).

## Phase 3 — Tube-uniform size field (`safs.geo` + driver)

**Goal:** within `tube_radius` of any fault, mesh size is uniform at
`res_f`; coarsens to `res_ff` over `ramp_dist` beyond the tube.

**Files to modify:**
- `mesh/safs.geo` — line 100-105 (R-004 fix above).
- `mesh/generate_safs_mesh.py` — add `--tube-radius` CLI flag, pass via
  `-setnumber tube_radius`.
- `mesh/validate_msh.py` — new `check_11_tube_uniformity` (R-004 test
  above), turned on when `sizing.json` carries a non-zero `tube_radius_m`.

**Implementation note:** The Distance field already samples each fault
surface (`Field[1].Sampling = 100`). For the tube to be a true 3-D
shell rather than a 2-D ring, `Sampling` may need to scale with fault
area; defer that micro-tuning unless `check_11` fails on the 2-fault
target.

**Acceptance:**
- 2-fault mesh built with `--tube-radius 4000 --res-f 1000 --res-ff
  20000 --ramp-dist 16000` (= 4 km uniform + 16 km ramp = 20 km total
  influence): mean fault-edge length 1000 ± 100 m; mean within-tube tet
  edge length 1000 ± 100 m; far-field (>20 km) mean edge length within
  ±10% of `res_ff`.
- Visual check in ParaView: distinct uniform "shell" around each fault,
  smooth gradation outside.

## Phase 4 — Wire it together; the 2-fault target

**Goal:** prove the 2-fault test (Mill Creek × SBMT San Andreas) passes
all 10 validator checks, including `check_5` (4917-style failure → 0).

**Files to create:**
- `mesh/run_two_crossing_2000m.sh` (R-005 fix above; full pipeline:
  audit → cleanup → **conformalize** → mesh → provenance → validate).

**Files to modify:**
- `mesh/generate_safs_mesh.py` — early-fail gate (R-002 fix); when
  `--stl-dir` already contains a `triangle_to_fault.json`, treat it as
  conformal-input mode and skip the cross-fault scan; otherwise run
  the scan and refuse non-conformal multi-fault input.
- `mesh/safs.geo` — switch back to `Mesh.Algorithm3D = 10` (HXT) for
  conformal multi-fault input. Keep Algorithm3D = 4 only as a manually-
  enabled fallback (`-setnumber Algorithm3D <n>`).

**Acceptance:**
- `bash run_two_crossing_2000m.sh` → `10/10 checks passed`. In
  particular `check_5_internal_interface`: `n_internal_with_2_tets`
  ≈ (#fault tris − #trace tris); `n_offending` = 0.
- `check_11_tube_uniformity` passes (R-004 acceptance).
- Run time on a laptop: < 60 s end-to-end.

## Phase 5 — Cascade to >2 faults (out of scope for this PR)

The conformalization in Phase 2 must already cascade (each new polyline
modifies the *current* triangulation, not the original) so that 3-, 4-,
… 8-fault inputs work without algorithm changes — only the cluster of
crossings grows. Validation: rerun `run_subsets_2000m.sh` with a new
"crossing subset" recipe (e.g., NW + Mill = 3 faults including a
crossing pair) and confirm 10/10 in that intermediate.

The all-8 build is the final goal but is gated on Pinto + Banning
CFM-source defects, which are independent (per README §Caveats 2) and
out of scope here.

---

## Why not the alternatives?

The README enumerates two alternative paths. Both have been prototyped
or considered and rejected for this code base:

- **Path B — `manifold` / CGAL boolean union:** requires closed manifolds;
  CFM faults are open surfaces with boundary. Closing each fault
  virtually before union and re-opening the seam is more work than the
  custom triangle splitter once the splitter exists. Also macOS install
  pain (CMake + C++17). Defer.
- **Path C — TetGen PLC:** has the same non-self-intersection requirement
  as HXT, so it does not solve the problem; the input still has to be
  conformalized first. The conformalization is the core work, regardless
  of mesher.

The custom Phase 1 + 2 pipeline (≈ Path A in README parlance) is the
shortest path to a working 2-fault test, and the primitives are
reusable for the all-8 goal once Pinto + Banning are addressed.

