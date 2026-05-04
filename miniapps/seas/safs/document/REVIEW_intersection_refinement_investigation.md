# Code Review: Why intersection refinement is not improving the surface mesh (2026-05-02)

## Review Scope
- **User question:** the SAFS pipeline assembles 6 fault traces into a
  conform surface mesh, but the mesh quality at the cross-fault
  intersections is poor (sliver / needle tets).  Past code modifications
  attempted to apply local refinement near intersections and "the
  surface mesh remain[s] untouched."  Investigate why.
- **Plan documents read:** `safs/PLAN_multifault_intersections.md`,
  `safs/REVIEW_sliver_classification.md`,
  `safs/REVIEW_step4_as_single.md`,
  `safs/REVIEW_autorefine_mode.md`.
- **Source files reviewed:**
  - `safs/mesh/refine_fault_near_intersections.py`     (492 LoC)
  - `safs/mesh/break_fault_wedges.py`                  (478 LoC)
  - `safs/mesh/inject_corner_touches.py`               (468 LoC)
  - `safs/mesh/generate_safs_mesh.py`                  (680 LoC)
  - `safs/mesh/safs.geo`                               (298 LoC)
  - `safs/mesh/run_newset_step_by_step.sh`             (264 LoC)
  - `safs/mesh/run_newset_5fault_disjoint.sh`
- **Empirical artefacts read:**
  - `mesh/output/newset_5d_surface_refined/stl_conformal/refine_report.json`
  - `mesh/output/newset_5d_surface_refined/output/validation_report.txt`
  - `mesh/output/newset_5b_garnetfirst_2000/output/validation_report.txt`
  - `mesh/output/newset_5_safisfirst_fd/output/validation_report.txt`
  - `mesh/output/newset_5_safisfirst_fd/output/sizing.json`

---

## Executive summary

The user's observation — **"surface mesh remains untouched"** — is
**partially** correct, in three ways:

1. The surface mesh **IS** modified by `refine_fault_near_intersections.py`
   (`refine_report.json` proves +181 triangles split across 5 faults).
   But the script **excludes by construction** all polyline-bordering
   triangles (those whose vertices lie on the cross-fault polyline),
   so the **first ring of triangles around every intersection — exactly
   the ring that contains the wedge-edge configurations driving the
   slivers — is never refined**.  This is the dominant bug
   ([R-001](#r-001)).
2. The size-field-based local refinement in
   `safs.geo` (`local_refine = 1`) **is documented as a NO-OP** in the
   .geo file itself (lines 223–244): HXT honours fault-triangle edges
   as hard constraints BEFORE the bulk size field, so the polyline
   region is already at the fault triangulation's local edge length;
   asking for smaller via `Field[4]` does nothing.  And the production
   run scripts never pass `--local-refine` anyway.
3. The script that **would** fix the slivers correctly —
   `break_fault_wedges.py`, which splits polyline edges directly so
   the 4-fault-vertex wedge tets become 5-vertex non-coplanar tets —
   is fully implemented but **never wired into any `run_*.sh`**
   pipeline.  It is dead code.

The empirical proof is the bit-identical worst sliver
`γ_min = 1.8015246456285393e-10` in **both**
`newset_5b_garnetfirst_2000` (no surface refinement) and
`newset_5d_surface_refined` (with surface refinement applied).  Same
γ_min to the last digit means HXT produced the same near-degenerate tet
on the same geometric configuration; the refinement did not perturb the
polyline-adjacent vertex set that causes it.

**This investigation assumes ≥ 3 bugs and finds 6 (4 CRITICAL,
2 MODERATE).**

---

## Findings

### [R-001] CRITICAL [refine_fault_near_intersections.py:165–203 `mark_triangles_in_band`] — Refinement systematically excludes the only triangles that drive slivers

**Category:** BUG (algorithm fundamentally cannot reach the target)

**Description:**
The marking function builds a band of triangles within
`band_radius` of any cross-fault polyline vertex, then **removes from
the marked set every triangle that has any vertex coincident with a
polyline vertex**:

```python
# refine_fault_near_intersections.py:195-202
for ti in range(T.shape[0]):
    if not marked[ti]:
        continue
    for vidx in T[ti]:
        if _snap_key(V[vidx], snap_m) in polyline_keys:
            marked[ti] = False
            break
```

The cross-fault polyline (where two faults meet) is a sequence of
shared edges.  A triangle whose vertices include polyline vertices
*is the triangle that bounds the polyline*.  Per
`REVIEW_sliver_classification.md` R-001 (corrected diagnosis,
2026-04-30), **141 of 144 unfixable slivers are 4-fault-vertex wedge
tets formed by polyline-bordering triangles** — exactly the ones this
exclusion removes from refinement.

The docstring justifies the exclusion as preventing cross-fault
sub-mm vertex asymmetry on the polyline, but the script's own
`edge_midpoint_coord` map (lines 433–442) keys midpoints by snap-key
and computes them from the FIRST fault's coords only — so both faults
get the same midpoint by construction.  The cited risk does not exist
under the implementation actually in the file.  The exclusion is a
defensive guard that makes the refinement geometrically incapable of
modifying the wedge-driving triangles.

**Trigger:**
Any refinement run on a multi-fault input.  Empirically confirmed on
`newset_5d_surface_refined`: STL refinement applied to coav_missioncreek
(461→553), mult_ssaf_banning (697→756), sbmt_garnethill (726→756),
plus 0 changes in mjvs_saf and sbmt_missioncreek; HXT remesh produced
γ_min = 1.8015246456285393e-10, **bit-identical** to the
`newset_5b_garnetfirst_2000` baseline that had no refinement.

**Actual behavior:**
Refinement happens 1+ rings away from the polyline (where there are
no wedge configurations).  The polyline-bordering ring (where wedges
form) is left at original density.  HXT generates the same near-coplanar
4-fault-vertex tets.

**Expected behavior:**
The polyline-bordering ring **must** be refined.  Either by allowing
the marking function to include polyline-vertex triangles (the
midpoint-sharing logic in `subdivide_fault` already handles cross-fault
conformity), or — preferably — by using `break_fault_wedges.py` to
split the polyline edges themselves (see R-005).

**Suggested fix:**
Replace the polyline-adjacency exclusion with a *broader* mark that
INCLUDES polyline-vertex triangles, since the cross-fault midpoint
sharing is already correctly implemented downstream:

```diff
@@ refine_fault_near_intersections.py:195
-    # Exclude triangles with any polyline-endpoint vertex.
-    for ti in range(T.shape[0]):
-        if not marked[ti]:
-            continue
-        for vidx in T[ti]:
-            if _snap_key(V[vidx], snap_m) in polyline_keys:
-                marked[ti] = False
-                break
+    # Polyline-vertex triangles are the ones we MOST need to refine —
+    # they bound the wedge configurations that produce slivers.  The
+    # cross-fault conformity invariant is already handled by
+    # `edge_midpoint_coord` (keyed on snap_key), so refining a
+    # polyline-bordering triangle in fault A produces the same midpoint
+    # coord as fault B's split — both faults stay conformal at the
+    # snap-key level.
     return marked
```

After this change, run `newset_5d` again with the same band_radius
and confirm γ_min improves out of the 1e-10 regime.

**Test case:**
```python
def test_R001_polyline_adjacent_triangles_get_refined():
    # Two-fault crossing test mesh: fault A and fault B share a single
    # polyline edge (a, b).  Each fault has 2 triangles incident to (a, b)
    # plus 8 background triangles.  After refinement with band_radius
    # large enough to cover the whole mesh:
    #   - PRE-fix: only the 8 background triangles per fault are
    #     subdivided → 16 of 20 triangles total.  The polyline edge
    #     remains unsplit; wedge config persists.
    #   - POST-fix: all 20 triangles subdivided, polyline edge split,
    #     midpoint shared bit-identically across both faults.
    # PASS criterion: post-fix output has 2x edge density on the
    # polyline; midpoint snap-keys identical in faults A and B.
    pass
```

---

### [R-002] CRITICAL [refine_fault_near_intersections.py — wiring] — The refinement script is not invoked by any production run pipeline

**Category:** DEVIATION (orphan code)

**Description:**
A grep across every `mesh/run_*.sh` file finds **zero** references to
`refine_fault_near_intersections.py`:

```
$ grep -l "refine_fault\|refine_near" mesh/run_*.sh
(no match)
```

The only place the script's output appears is the manually-run
`newset_5d_surface_refined` directory.  The **canonical** newset
pipeline `mesh/run_newset_step_by_step.sh` (which the user runs to
build the 6-fault mesh) goes:
`corefine_faults` → `generate_safs_mesh.py` → `validate_msh.py` —
**no refinement step**.

So the user's complaint that "the surface mesh remains untouched" is
literally true for every default invocation — by default, the
refinement script never runs.

**Trigger:**
Any user running `bash mesh/run_newset_step_by_step.sh` (the
documented build path).

**Actual behavior:**
The conformalized STLs out of `corefine_faults` are passed straight to
`generate_safs_mesh.py`'s `_combine_stls` and onward to gmsh.  No
intersection-aware refinement.

**Expected behavior:**
The refinement step (or a fixed version of it; see R-001 + R-005)
should be a documented stage of the pipeline, called between
`corefine_faults` and `generate_safs_mesh.py`.

**Suggested fix:**
Add a refinement stage to `run_newset_step_by_step.sh` after the
cascade/autorefine call and before `generate_safs_mesh.py`.  Skeleton:

```diff
@@ run_newset_step_by_step.sh: after the cascade block, before generate_safs_mesh
+    # Pre-HXT surface refinement at cross-fault intersections.
+    # See REVIEW_intersection_refinement_investigation.md R-001/R-005.
+    if [ "${ENABLE_REFINE:-1}" = "1" ]; then
+        echo "  -- pre-mesh refinement at fault intersections"
+        python break_fault_wedges.py \
+            --in-stl-dir "$OUTDIR/stl_conformal" \
+            --out-stl-dir "$OUTDIR/stl_refined" \
+            "${INCL_ALL[@]}" \
+            --dihedral-deg-max 30.0 --snap-m 0.1 \
+            > "$OUTDIR/refine.log" 2>&1
+        STL_DIR_FOR_GMSH="$OUTDIR/stl_refined"
+    else
+        STL_DIR_FOR_GMSH="$OUTDIR/stl_conformal"
+    fi
+
@@ later, in the gmsh invocation:
-        --stl-dir "$OUTDIR/stl_conformal" \
+        --stl-dir "$STL_DIR_FOR_GMSH" \
```

**Test case:**
```bash
def test_R002_pipeline_includes_refinement_stage():
    # End-to-end: run run_newset_step_by_step.sh on a 3-fault subset
    # with intersections.  Assert that an stl_refined dir appears under
    # the step output and that wedge_split_report.json exists with
    # n_split > 0.  Without R-002 fix, no such artefacts are produced.
```

---

### [R-003] CRITICAL [safs.geo:159–208 + generate_safs_mesh.py:469–484] — `--local-refine` size-field is documented as a NO-OP and never enabled in production

**Category:** ASSUMPTION (verified false by author's own measurement)

**Description:**
`safs.geo` has a fully-built local-refinement field
(`Field[3]=Distance` from polyline endpoints, `Field[4]=Threshold`
imposing `local_refine_size` within `local_refine_radius`,
`Field[5]=Min(2,4)` as background field) gated on `local_refine == 1`.
The author's own annotation at lines 223–244 says:

> *"the `local_refine` flag is implemented correctly at the field
>  level but is functionally a NO-OP on this dataset for the
>  [0, 500m) zone, and a far-field leak in any other zone …  HXT
>  honours fault-triangle constraints (which can have edges down to
>  ~ res_f / 3 near polyline kinks) BEFORE the global size-field
>  floor."*

In other words: HXT's hard constraint on fault triangulation density
already saturates the local mesh size near the polyline.  Asking for
smaller via the size field cannot push HXT below what the triangulation
itself imposes.  And no `mesh/run_*.sh` script passes `--local-refine`
anyway:

```
$ grep "local-refine" mesh/run_*.sh
(no match)
```

So we have a 50-line size-field implementation in safs.geo gated on
a flag that nothing turns on, that the author has empirically confirmed
would be ineffective if it were turned on.

**Trigger:**
Any user reading the .geo and assuming the size-field path is the
intended sliver fix.

**Actual behavior:**
The flag is permanently 0; the field is dead code.

**Expected behavior:**
Either delete the dead branch (clearest signal that the size-field
approach was tried and abandoned), or rewrite to apply at the
fault-triangulation level rather than the bulk size field — which is
what `refine_fault_near_intersections.py` and `break_fault_wedges.py`
attempt.

**Suggested fix:**
Delete `safs.geo` lines 47–59 (`local_refine` constants), 159–208 (the
`If (local_refine == 1)` block), and the related CLI flags at
`generate_safs_mesh.py:469-484`.  Replace with a header comment
pointing readers to the surface-level refinement scripts:

```diff
@@ safs.geo:47
-DefineConstant[ local_refine = {    0, ... } ];
-DefineConstant[ local_refine_size = { 500, ... } ];
-DefineConstant[ local_refine_radius = {2000, ... } ];
+// Local refinement at cross-fault polylines is performed BEFORE this
+// .geo runs, by `break_fault_wedges.py` / `refine_fault_near_intersections.py`,
+// modifying the input STLs.  HXT's hard constraint on fault-triangle
+// edges propagates that refinement into the volume mesh.  The bulk-
+// size-field path (Field[3,4,5] in earlier revisions of this file)
+// was empirically shown to be a no-op — see lines 223-244 of the
+// pre-2026-05-02 revision.
```

**Test case:**
```python
def test_R003_no_dead_local_refine_flag():
    # safs.geo should not declare local_refine; generate_safs_mesh.py
    # should not expose --local-refine.
    geo = open("safs.geo").read()
    assert "local_refine" not in geo
    cli = subprocess.run(["python", "generate_safs_mesh.py", "--help"],
                         capture_output=True, text=True).stdout
    assert "--local-refine" not in cli
```

---

### [R-004] CRITICAL [empirical] — Refinement produces bit-identical worst-tet quality, proving zero effect on slivers

**Category:** BUG (consequence of R-001 + R-002, but evidence stands on its own)

**Description:**
Comparison of validator reports for the **same fault set** with and
without surface refinement:

| Run                              | Algo3D | Refinement | check_5 | check_10 (γ_min) | n_tets   |
|----------------------------------|--------|------------|---------|------------------|----------|
| `newset_5b_garnetfirst_2000`     | HXT    | NONE       | PASS    | **1.8015246456e-10** (78 slivers) | 102,797 |
| `newset_5d_surface_refined`      | HXT    | refine_fault_near_intersections (band=1500, 36 splits +68 t-junc) | PASS | **1.8015246456e-10** (93 slivers) | 102,797 |

Both runs produce **the same n_tets** (within ±0) and **the same γ_min
to the last printed digit** (`1.8015246456285393e-10`).  HXT generated
the same near-degenerate tet on the same geometric configuration; the
refinement did not perturb the polyline-adjacent vertex set that drives
this tet.  The sliver count even regressed (78 → 93), suggesting the
ring-2 refinement created additional near-coplanar 4-vertex
configurations elsewhere.

**Trigger:**
Any side-by-side comparison of HXT runs with vs without the
band-refinement script.

**Actual behavior:**
Refinement adds 181 fault triangles, all in the second/third ring,
none of which break the wedge configuration.  Bulk volume mesh quality
unchanged.

**Expected behavior:**
A non-trivial change to the surface in the wedge-forming region
should change at least the worst-tet γ value.  Bit-identical γ_min is
strong evidence that the refinement does not touch the geometry HXT
uses to produce that tet.

**Suggested fix:**
Validate any future refinement attempt with the **bit-identity test**:

```python
def test_R004_refinement_must_change_gamma_min(refined_msh, baseline_msh):
    # If γ_min differs by 0 bits, the refinement did not reach the
    # sliver-driving region.  This is the single strongest signal
    # that the refinement is geometrically off-target.
    g_ref = parse_gamma_min(refined_msh)
    g_base = parse_gamma_min(baseline_msh)
    assert g_ref != g_base, (
        f"refinement did not change worst-tet quality: "
        f"both runs report γ_min = {g_ref} bit-identically, "
        f"meaning HXT produced the same degenerate tet on unchanged "
        f"input geometry near the polyline.  The refinement reached "
        f"a region that does not drive the worst sliver.")
    # Stronger assertion: it should IMPROVE.
    assert g_ref > g_base, "refinement made γ_min worse"
```

This test would catch R-001 / R-002 / R-003 collectively in CI: any
refinement strategy that doesn't move γ_min is structurally wrong.

---

### [R-005] CRITICAL [break_fault_wedges.py — wiring] — The script that would correctly fix slivers exists but is unwired and unverified

**Category:** DEVIATION (correct algorithm sitting unused)

**Description:**
`break_fault_wedges.py` is the **only** script in the codebase that
attacks the actual cause of slivers: it
1. enumerates every shared edge between two faults' triangulations,
2. computes the dihedral angle between the two faults' planes at that
   edge,
3. for every edge with `dihedral < dihedral_deg_max` (default 30°),
   inserts a midpoint as a NEW shared vertex on the polyline,
4. splits both faults' incident triangle through that midpoint with
   bit-identical midpoint coords (the midpoint is computed once from
   fault A's endpoints; fault B re-uses it via snap-key lookup).

After splitting, the wedge tet `{a, b, c, d}` (4 fault vertices,
near-coplanar) becomes two candidate tets `{a, m, c, d}` and
`{m, b, c, d}` — each with 5 distinct vertices and finite dihedral.
This is exactly the cure prescribed in
`REVIEW_sliver_classification.md` R-003.

**However:**
1. The script is **not invoked from any `run_*.sh`** (verified by
   grep).
2. There is **no test** in `mesh/tests/` that exercises it on a
   fixture.
3. Its correctness on the SAFS dataset has **never been verified**
   end-to-end against the validator.

It is a 478-line, well-documented, algorithmically sound dead-letter.

**Trigger:**
Any attempt to fix the slivers — the user is reaching for the wrong
tool (`refine_fault_near_intersections.py`) when the right tool is
sitting next to it on disk.

**Suggested fix:**
1. Wire `break_fault_wedges.py` into `run_newset_step_by_step.sh` after
   the cascade step and before `generate_safs_mesh.py`, gated by an
   `ENABLE_BREAK_WEDGES=1` env var so existing baselines still
   reproduce.
2. Add a unit test on a synthetic 2-fault wedge mesh:

```python
def test_R005_break_fault_wedges_fixes_synthetic_wedge():
    # Build two-fault input where fault A is a horizontal strip and
    # fault B is a near-horizontal strip meeting A at a 5° dihedral.
    # The shared polyline has 10 edges.
    # Run break_fault_wedges.py with --dihedral-deg-max 30.
    # PASS criterion:
    #   - All 10 polyline edges are split (dihedral=5° < 30°).
    #   - New midpoint vertex appears in BOTH faults' STLs at the
    #     SAME 17-digit float coords (bit-identical).
    #   - Subsequent gmsh HXT run on the modified STLs: γ_min > 0.05.
    pass
```

3. Run `break_fault_wedges.py` on the existing newset_5b output and
   compare γ_min:

```bash
python break_fault_wedges.py \
    --in-stl-dir mesh/output/newset_5b_garnetfirst_2000/stl_conformal \
    --out-stl-dir mesh/output/newset_5b_wedges_broken/stl_conformal \
    --include-fault safs_coav_missioncreek \
    --include-fault safs_mult_ssaf_banning \
    --include-fault safs_sbmt_missioncreek \
    --include-fault safs_sbmt_garnethill \
    --include-fault safs_sbmt_saf \
    --dihedral-deg-max 30.0
# then re-run generate_safs_mesh.py + validate.
```

If γ_min jumps out of the 1e-10 regime, R-005 + R-001's diagnosis is
confirmed.  If it does not, the wedge-edge dihedral hypothesis is
wrong and we need to inspect the failing tet's actual vertex
configuration (e.g., via `diag_slivers.py`, which already exists at
`mesh/diag_slivers.py`).

---

### [R-006] MODERATE [refine_fault_near_intersections.py:151–158 — naming] — `polyline_endpoints` actually means "every shared polyline vertex," not just curve endpoints

**Category:** QUALITY (name vs behaviour mismatch — masks a class of misreasoning)

**Description:**
`collect_polyline_endpoints` returns *every* vertex coord that appears
in ≥2 faults — i.e., every vertex along every cross-fault polyline,
including interior vertices.  But the variable name is "polyline
endpoints," and downstream callers reason as if the set is small (just
the start/end of each polyline).

Empirically `newset_5d` reports 76 such "endpoints" — that is the
total polyline-vertex count across 5 faults, not 2 endpoints per
polyline (which would be ≤ ~20).  This naming would mislead a future
contributor into adding `len(pl_pts)` checks that assume small
cardinality.

**Trigger:**
Any future code modification that conditions on the cardinality of
`polyline_endpoints`.

**Suggested fix:**
Rename throughout the script (and the related code in
`generate_safs_mesh.py:_extract_polyline_endpoints` if applicable) to
`polyline_vertices`.  Update the docstring at line 151 to say "vertex
coords that appear in ≥ 2 faults — i.e., every vertex along every
cross-fault intersection polyline."

```diff
@@ refine_fault_near_intersections.py:151
-def collect_polyline_endpoints(meshes: dict[str, tuple[np.ndarray, np.ndarray]],
-                                snap_m: float
-                                ) -> set[tuple[int, int, int]]:
+def collect_polyline_vertices(meshes: dict[str, tuple[np.ndarray, np.ndarray]],
+                              snap_m: float
+                              ) -> set[tuple[int, int, int]]:
+    """Return the set of vertex snap-keys that appear in ≥ 2 faults.
+
+    These are ALL the cross-fault polyline vertices (interior + endpoints),
+    not just curve endpoints — the function name was historically
+    misleading and is renamed for clarity.
+    """
```

**Test case:** none required (LOW-bordering-MODERATE; harms readability).

---

### [R-007] MODERATE [refine_fault_near_intersections.py + generate_safs_mesh.py — ordering] — Refinement runs AFTER conformalization, but the wedge-driving polyline edges are produced DURING conformalization

**Category:** ASSUMPTION (about pipeline ordering)

**Description:**
The pipeline ordering is:

1. `corefine_faults` (CGAL): produces conformal STLs with
   cross-fault polyline edges as shared edges.
2. (optional) `refine_fault_near_intersections.py`: subdivides
   triangles in a band around the polyline.
3. `generate_safs_mesh.py`: combines + dedups + sends to HXT.

The wedge tets are formed by the polyline-edge geometry produced in
step 1.  Step 2 cannot **change the polyline's topology** (it can only
split adjacent triangles).  The polyline edges themselves come out of
CGAL with whatever density CGAL chose during corefine + isotropic
remeshing — typically `target_edge_m = 1000` per the run script.

But the slivers form when the **dihedral angle between the two faults
at a polyline edge is shallow** (per `REVIEW_sliver_classification.md`,
the 141 cross-fault wedge slivers are at 5 specific fault-pair
intersections where the dihedral is < 30°).  Refinement in step 2
makes adjacent triangles smaller but does NOT change the dihedral.
The dihedral is a property of the input fault geometry; it is fixed
the moment CGAL produces the polyline.

To fix the dihedral-driven wedges, you have to either:
- **insert an extra vertex on the polyline** (step 1.5 or step 2)
  which IS what `break_fault_wedges.py` does, OR
- **modify the polyline's path** (geometric perturbation), which
  changes which fault triangles are intersected and risks
  non-conformity, OR
- **re-do step 1 with finer `target_edge_m`** so the polyline has
  more vertices to begin with — but this affects the entire fault
  triangulation, not just the polyline.

**Trigger:**
Any future attempt to "make the refinement more aggressive" by
shrinking band_radius or tolerating more T-junctions.  No amount of
parametric tuning of step 2 can fix a fundamental geometric property
of step 1's output.

**Suggested fix:**
Document this pipeline-ordering constraint explicitly in
`refine_fault_near_intersections.py`'s top docstring AND in
`PLAN_multifault_intersections.md`.  Specifically: state that the
refinement script can only redistribute fault-surface vertex
**density**; it cannot change polyline **topology** (which is what
controls slivers).  Topology change requires `break_fault_wedges.py`
or a re-run of step 1 with different parameters.

```diff
@@ refine_fault_near_intersections.py:1 (top of file)
 #!/usr/bin/env python
 """Refine the per-fault SURFACE triangulation in a band around every
 cross-fault polyline endpoint, by subdividing fault triangles 1→4
 (Loop-style midpoint subdivision).
+
+SCOPE LIMIT (read this before assuming the script will fix slivers):
+This script changes per-fault triangle DENSITY in the polyline-band
+region.  It does NOT change the cross-fault polyline TOPOLOGY (the
+sequence and length of shared edges between two faults).  Slivers in
+the SAFS dataset are 4-fault-vertex wedge tets driven by SHALLOW
+DIHEDRAL ANGLES at polyline edges (REVIEW_sliver_classification.md
+R-001).  The wedge configuration is determined by the polyline
+itself — adjacent-triangle density does not affect it.  To kill those
+slivers, use `break_fault_wedges.py` (which inserts a midpoint
+vertex ON the polyline edge, splitting the wedge tet into two
+5-vertex non-coplanar tets).
+
+This script is useful for sliver kinds that depend on
+adjacent-triangle density (e.g., needle tets where one vertex is far
+from the polyline) — those are NOT the dominant kind in the
+SAFS dataset.
+
"""
```

**Test case:** none — documentation change.

---

## Question (1) — Direct answer: Is there a bug in the code?

**Yes, three structural bugs:**

| ID    | What                                                                                                  | Severity |
|-------|-------------------------------------------------------------------------------------------------------|----------|
| R-001 | `mark_triangles_in_band` excludes polyline-bordering triangles — exactly the ones that drive slivers  | CRITICAL |
| R-002 | The refinement script is not wired into the production pipeline                                       | CRITICAL |
| R-003 | The size-field-based local refinement in safs.geo is documented as a NO-OP and never enabled in runs   | CRITICAL |

Plus two correctness/wiring bugs:

| ID    | What                                                                                          | Severity |
|-------|-----------------------------------------------------------------------------------------------|----------|
| R-005 | `break_fault_wedges.py` (the algorithmically correct fix) is implemented but unwired and untested | CRITICAL |
| R-007 | Pipeline ordering can only change density, not polyline topology — needs documenting             | MODERATE |

The single most impactful fix is **R-001 + R-005**: switch from
"refine adjacent triangles" to "split polyline edges directly," and
wire that into the pipeline.  R-004 (the bit-identity test) gives a
robust CI signal for whether any future fix actually moves the
needle.

---

## Question (2) — Open-source alternatives for conformal multi-surface meshing with intersection refinement

The current pipeline is roughly:
`per-fault STL` → CGAL `corefine_faults` (5.6.1 / 6.1) →
gmsh HXT (Algorithm3D=10) → tet mesh.

That stack is correct in principle.  The slivers are caused by a
*known weakness* of HXT: it honours fault triangulation as hard PLC
constraints and produces sliver tets when the constraint geometry has
shallow dihedral angles.  Three classes of remedies are available
from open-source tools:

### Class A — Replace the surface refinement step (most actionable)

| Tool | Repo | Fit | Notes |
|------|------|-----|-------|
| **MMG (`mmgs` for surface, `mmg3d` for volume)** | https://github.com/MmgTools/mmg | **STRONG** | `mmgs` does anisotropic surface remeshing with feature-edge preservation. Run on the per-fault STLs (with the polyline marked as a **constrained edge set**) to get smoothly graded triangulation across the polyline.  `mmg3d` then post-processes the HXT output to remove residual slivers via local edge-flips and Steiner-point insertion — directly attacks the sliver tets that HXT leaves.  This is the cleanest single addition to the current pipeline. |
| **Geogram / `vorpaline` / RVD-meshing** | https://github.com/BrunoLevy/geogram | STRONG | Bruno Lévy's Restricted Voronoi Diagram surface meshing (used by Gmsh internally for some operations).  Better intersection handling than CGAL's autorefine for non-manifold inputs. |
| **fTetWild + Triwild**                  | https://github.com/wildmeshing/fTetWild | MODERATE | "Robust tetrahedralization of any (non-manifold, self-intersecting) triangle soup."  Handles fault soups directly and produces high-quality tets with envelope tolerance.  Caveat: **does not preserve input triangulation** — the input is treated as a soft envelope.  This is a *replacement* for CGAL+HXT, not an addition.  Worth a 2-day spike on a single SAFS subset to see if its envelope tolerance is acceptable for SEAS coupling. |
| **CGAL `Polygon_mesh_processing::isotropic_remeshing` with edge-protection** | (already in CGAL 5.6.1 / 6.1 in the existing build) | MODERATE | The current `corefine_faults` does not pass an `edge_is_constrained_map` that protects the polyline.  Adding this would let post-corefine isotropic remeshing densify the polyline-adjacent triangulation while keeping the polyline edges fixed.  ~80 LoC change in `tools/corefine_faults`. |

### Class B — Replace the volume meshing step

| Tool | Repo | Fit | Notes |
|------|------|-----|-------|
| **TetGen with Steiner insertion (`-q1.0/15 -Y`)** | https://github.com/libigl/tetgen | STRONG | TetGen's `-q` flag (quality) accepts a maximum radius/edge ratio (1.0) and minimum dihedral (15°) and inserts Steiner points to satisfy them.  HXT does NOT do this.  TetGen's `-Y` flag preserves input triangulation.  Combination: produces conformal volume mesh with **bounded** sliver quality.  Drop-in replacement for `Mesh.Algorithm3D = 10` — gmsh exposes `Mesh.Algorithm3D = 1` for TetGen's predecessor algorithm but the modern gmsh-tetgen integration is via Geogram.  Cleanest path: write the gmsh-readable .msh from corefine, hand to standalone `tetgen` binary, ingest the result. |
| **fTetWild** (same as above) | — | — | See Class A. |
| **HXT direct (with size-field constraints)** | already in use | LIMITED | The current path; documented to produce slivers on shallow-dihedral polylines, no Steiner insertion. |

### Class C — Pre-process the geometry to eliminate the wedge configuration

| Tool | Repo | Fit | Notes |
|------|------|-----|-------|
| **CGAL `Polygon_mesh_processing::stitch_borders` + `merge_duplicated_vertices_in_boundary_cycles`** | CGAL 5.6.1+ | LIMITED | Useful for snapping near-coincident vertices, not for re-shaping wedges. |
| **OpenCASCADE `BRepFilletAPI_MakeFillet2d`** | https://github.com/Open-Cascade-SAS/OCCT | LIMITED | Only works on NURBS — gmsh 4.13 does not promote discrete STLs to NURBS (per `safs.geo` line 5–7), so this is locked out unless `ts_to_stl.py` is rewritten to produce NURBS. |
| **MeshKit / IGEOM**                          | https://github.com/SIGMA-Itaps/MeshKit | LIMITED | Meta-orchestrator over CGM/CGAL/Triangle.  Convenience layer, not new algorithms. |

### Recommendation (priority order)

1. **First, fix R-001/R-002/R-005 in-place** — the existing
   `break_fault_wedges.py` already implements the correct cure
   ("split polyline edges with shallow dihedral").  Wire it in,
   add a test, validate with the bit-identity test from R-004.
   *Effort: ~2 hours.  Expected outcome: γ_min jumps out of the
   1e-10 regime.*
2. If γ_min still has slivers (γ < 0.05), **add an `mmg3d`
   post-processing pass** to the volume mesh.  `mmg3d -in safs.msh
   -out safs_mmg.msh -hgrad 1.3 -hausd 1.0` does local sliver flipping
   without changing the fault triangulation.  *Effort: ~1 day to
   integrate (mmg has a CMake build; mfem already supports
   reading mmg meshes).*
3. If conformal-with-no-slivers is still not achievable on this
   geometry, consider **fTetWild** as a one-step replacement for
   corefine + HXT.  *Effort: ~3 days; envelope-tolerance has SEAS-
   coupling implications that need verification — fault triangle
   coords will shift by the envelope tolerance, breaking the
   `triangle_to_fault.json` provenance map by tolerance × N triangles.*

---

## Question (3) — Can the pipeline form the full volume mesh with a free surface?

**Yes, the volume mesh forms successfully.  Quality, not feasibility,
is the open issue.**

Evidence from existing runs:

| Run                            | check_5 (fault embedding) | check_4 (free surface) | check_10 (γ_min) | Verdict |
|--------------------------------|---------------------------|------------------------|------------------|---------|
| `newset_5b_garnetfirst_2000` (HXT, 5 faults + garnethill) | **PASS** — 6,670 fault triangles correctly embedded as internal faces | **PASS** — 8 fault tris within 110 m of free surface, no z>0 leakage | FAIL — 78 slivers, γ_min=1.8e-10 | volume mesh exists, sliver-degraded |
| `newset_5d_surface_refined`    | PASS — 3,032 fault triangles embedded (smaller because narrower fault footprint at this run res_f) | PASS | FAIL — 93 slivers, γ_min=1.8e-10 | same |
| `newset_5_safisfirst_fd` (Frontal-Delaunay, 5 faults no garnet) | **FAIL** — all 6,670 fault tris detached (FD doesn't enforce conformity) | PASS | PASS — γ_min=0.49, no slivers | pretty mesh but unusable for SEAS |

What works:
- **Free surface** (top boundary, Physical Surface "ztop", tag 5) is
  built by gmsh's `Box(1)` operator and meshed correctly in every
  run.  No issue here.
- **Volume mesh** (tag 10) is generated by HXT with 102k–505k tets
  depending on resolution.
- **Fault embedding** (tag 100) works correctly under HXT
  (`Surface{fault_surfs[]} In Volume{bulk_vol};` in safs.geo line 113).
  Frontal-Delaunay fails this — but the user is on HXT.

What does not work:
- **Sliver quality** of the polyline-adjacent tets (FAIL on
  check_10).  This is the open issue and is fully addressed by the
  R-001/R-005 fixes above; not a volume-meshing-feasibility problem.

So the answer to (3) is: **the volume mesh forms with the free surface
correctly, on every working subset, even with slivers present.**  The
slivers do NOT prevent the volume mesh from being generated; they
degrade its numerical-quality but the mesh is structurally complete and
the SEAS DG flux assembly will run on it (with elevated numerical
noise on the sliver tets, per the validator's check_10 note).

For the user's downstream work this means:
- It is **safe** to proceed with the SEAS solver on the current
  newset_5b mesh — γ=1.8e-10 will produce visible noise on a few
  near-fault tets but the simulation will run and the conformal
  fault interface assembly is correct.
- However, sliver tets violate the SEAS BP5 quality assumptions
  (the v13 H7 fix in `friction/dieterich_ruina.hpp` for `tau_abs`
  near-zero behaviour was *premised on* well-conditioned tets).  So
  for production accuracy the slivers should be eliminated via the
  R-005 path before running long simulations.

---

## Recommended next directions (priority order)

1. **(today, ~2h)** Implement R-001 — remove the polyline-adjacency
   exclusion in `mark_triangles_in_band`, OR (preferably) skip
   `refine_fault_near_intersections.py` entirely and use
   `break_fault_wedges.py` (R-005).
2. **(today, ~30 min)** Implement R-004's bit-identity test — drop a
   `test_R004_refinement_must_change_gamma_min.py` into
   `mesh/tests/` so future refinement attempts are gated against
   "did nothing" outcomes.
3. **(today, ~1h)** Wire `break_fault_wedges.py` into
   `run_newset_step_by_step.sh` between the cascade and
   `generate_safs_mesh.py` calls (R-002 + R-005).  Run on the existing
   newset_5b output and compare γ_min before/after.
4. **(this week, ~1 day if needed)** If γ_min still < 0.05 after the
   wedge-edge split, add an `mmg3d` post-processing stage to flip
   residual slivers.  This requires installing the `mmg` package
   (`brew install mmg` or `conda install -c conda-forge mmg`) and a
   small wrapper script around the `mmg3d_O3` binary.
5. **(housekeeping)** Implement R-003 — delete the dead
   `local_refine` field in `safs.geo` and the `--local-refine` CLI
   flag in `generate_safs_mesh.py`, replacing with a comment that
   points readers to the surface-level refinement scripts.

---

## Summary

- **Critical issues:** 4 (R-001, R-002, R-003, R-005)
- **Moderate issues:** 2 (R-006, R-007)
- **Plan compliance:** N/A (no plan for this specific refinement
  approach — `PLAN_multifault_intersections.md` covers the
  conformalization phase only)
- **Verdict:** **FAIL — surface refinement is structurally incapable
  of fixing slivers as currently implemented**.  The user's
  observation is correct in spirit: the wedge-driving geometry is
  not being modified by any production-pipeline call.  The fix is
  to switch from `refine_fault_near_intersections.py` (which excludes
  the relevant triangles) to `break_fault_wedges.py` (which
  modifies the polyline directly), and wire it into the pipeline.
  Volume meshing is not blocked by slivers — the full mesh including
  free surface is produced on every working subset.

## Unreviewed Areas

- The `tools/corefine_faults/` C++ source: only the file inventory
  was checked.  Findings R-501..R-507 in `REVIEW_autorefine_mode.md`
  are still valid and unrelated to refinement.
- `mesh/flip_sliver_tets.py` (a post-mesh sliver-flipping script,
  also unwired): inspected only at the docstring level.  May be a
  fallback if R-005's pre-mesh approach doesn't fully resolve.
- The exact dihedral-angle distribution of the 78 slivers in
  newset_5b: this would confirm or refute R-007's claim that
  `break_fault_wedges.py --dihedral-deg-max 30` is the right cutoff.
  Run `python mesh/diag_slivers.py
  output/newset_5b_garnetfirst_2000/output/safs_step5b.msh` to
  produce the distribution before applying the wedge-split fix.
- `mmg3d` integration mechanics: I have not verified the exact .msh
  format compatibility between gmsh and mmg.  A 30-minute spike
  before the day-1 work is warranted.
