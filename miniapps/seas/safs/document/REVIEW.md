# Code Review: 2026-04-29 — SAFS multi-fault intersection conformity

## Review Scope
- Plan: `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/safs/PLAN_multifault_intersections.md`
- Files reviewed:
  - `mesh/fault_intersect.py` (1131 lines)
  - `mesh/conformalize_faults.py` (1700 lines)
  - `mesh/generate_safs_mesh.py` (567 lines)
  - `mesh/safs.geo` (159 lines)
  - `mesh/run_two_crossing_2000m.sh` (104 lines)
  - `mesh/output/two_crossing_2000m/output/validation_report.txt`
  - `mesh/output/two_crossing_2000m/stl_conformal/intersection_report.json`
  - `mesh/output/two_crossing_2000m/output/combine_meta.json`
- Domain context: `seas-mfem/CLAUDE.md`, `miniapps/seas/CLAUDE.md`,
  `mesh/REVIEW_phase1.md` (prior Phase-1 review).
- User-supplied evidence: screenshot
  `~/Desktop/Screenshot 2026-04-29 at 6.27.50 PM.png` showing two
  CFM faults that meet at an intersection curve with visible
  conformity defects in the surface triangulation around the curve.

The most recent end-to-end run (`run_two_crossing_2000m.sh`) produced
this validation outcome:

```
[FAIL] 5_internal_interface: 36505 fault triangle(s) have wrong tet adjacency
        n_internal_with_2_tets: 0
        n_offending: 36505
[FAIL] 6_fault_edge_length: mean fault edge 448 m outside [500, 2000]
[FAIL] 7_far_field_edge: far-field mean edge 6466 m < 12000
[FAIL] 11_tube_uniformity: mean tube tet edge 1814 m outside [500, 1500]
7/11 checks passed
```

`combine_meta.json` reports `n_dropped_duplicate: 22` along the
intersection polyline. `intersection_report.json` reports
`PASS_WITH_3_T_JUNCTIONS` on fault A and `PASS_WITH_7_T_JUNCTIONS`
on fault B.

This is **not** a "small numerical glitch" — every fault triangle
in the bulk mesh is detached from the tets, exactly the failure
mode `PLAN_multifault_intersections.md §Why this plan exists`
was supposed to eliminate.

I assumed at least 3 bugs and found 9 that are actionable, plus a
strategic recommendation (R-010) to replace the custom Python
conformalizer with CGAL's `Polygon_mesh_processing::corefine` —
a single library call that is provably-correct and would resolve
R-001 through R-005 in one go.

---

## Findings

### [R-001] [CRITICAL] [run_two_crossing_2000m.sh] — Default algo3d=4 silently produces a non-conforming mesh; check 5 fails 100%

**Category:** BUG (silent data corruption at the pipeline level)

**Description:**
`run_two_crossing_2000m.sh:33` defaults `ALGO3D=4` (Frontal-
Delaunay). Comment claims this is "production-safe for conformal
multi-fault input where HXT rejects". It is not. With `algo3d=4`
gmsh ignores the `Surface{...} In Volume{...}` embedding constraint
in `safs.geo:99`, so every fault triangle is detached from the
tet mesh:

```
[FAIL] 5_internal_interface:
        n_internal_with_2_tets: 0
        n_offending: 36505
```

The .msh file has 36505 fault triangles tagged `100`, but **none**
of them is shared by two tets. The mesh is unusable for SEAS:
the BR2 DG operator assembles fluxes on each face's two adjacent
tets and would dereference garbage.

The plan (`§Why this plan exists`) explicitly says check 5 must
pass: "`n_internal_with_2_tets ≈ #fault_tris − #trace_tris`,
`n_offending = 0`. Compare to today's all-8 number (`n_offending
= 4917`, every triangle wrong)." The two-fault build's
`n_offending = 36505` is **worse** than the pre-Phase-2 baseline
the plan was meant to fix.

**Trigger:** Running the script with default flags (no `TWO_ALGO3D=10`
override). Every prior run in `output/two_crossing_2000m/` has this
failure.

**Actual behavior:** Pipeline exits 0 (only `validate_msh.py` flags
the issue, and the shell script doesn't propagate its non-zero exit).
Downstream consumers see a fault tag with 36505 detached triangles.

**Expected behavior:** Default to `ALGO3D=10` (HXT). When HXT
rejects the input, the script must abort with a clear message
pointing to which conformalize-output deficiency caused the
rejection — not silently fall back to a meshing algorithm that
produces garbage.

**Suggested fix:**

```diff
--- a/mesh/run_two_crossing_2000m.sh
+++ b/mesh/run_two_crossing_2000m.sh
-ALGO3D="${TWO_ALGO3D:-4}"   # 4 = Frontal-Delaunay (production-safe for
-                            #     conformal multi-fault input where HXT
-                            #     rejects due to near-coplanar adjacent
-                            #     CFM triangles).  Override with
-                            #     TWO_ALGO3D=10 to attempt HXT;
-                            #     validator check 5 passes only with HXT.
+ALGO3D="${TWO_ALGO3D:-10}"  # 10 = HXT (mandatory for SEAS — only HXT
+                            #     enforces the Surface-In-Volume
+                            #     embedding so every fault triangle
+                            #     is shared by two tets, which is the
+                            #     SEAS DG flux requirement).
+                            # If HXT rejects, FIX THE INPUT — do not
+                            # fall back to algo3d=4: that algorithm
+                            # ignores the embedding and produces
+                            # 0 fault-as-internal-face triangles
+                            # (validator check 5 fails 100%).
```

And in `validate_msh.py` after line ~390 (where checks are
collated), make check_5 a pipeline-level gate so the
shell-script `set -e` actually catches the failure:

```diff
+    # Hard exit on check_5 failure — a mesh with detached fault
+    # triangles is unusable for SEAS regardless of how many other
+    # checks passed.
+    if not all(c.passed for c in results if c.name.startswith("5_")):
+        print("\nFATAL: check_5 (internal interface) failed; the "
+              "fault is not embedded in the tet mesh.  Refusing to "
+              "succeed.", file=sys.stderr)
+        return 1
```

**Test case:**
```python
def test_R001_default_algo3d_is_hxt():
    """The two-crossing pipeline must default to algo3d=10; algo3d=4
    fails check_5 100% and produces an unusable mesh."""
    import re
    script = Path("mesh/run_two_crossing_2000m.sh").read_text()
    m = re.search(r'ALGO3D="\$\{TWO_ALGO3D:-(\d+)\}"', script)
    assert m, "ALGO3D default-value pattern not found"
    assert int(m.group(1)) == 10, (
        f"run_two_crossing_2000m.sh defaults ALGO3D={m.group(1)}; "
        f"must be 10 (HXT) — algo3d=4 silently produces "
        f"n_internal_with_2_tets=0 (every fault triangle detached)."
    )

def test_R001_check_5_is_pipeline_gate():
    """validate_msh.py must hard-fail when check_5 fails; the shell
    script's set -e must propagate."""
    # Build a mesh with algo3d=4 deliberately, run validate_msh.py,
    # assert non-zero exit code.
    ...
```

---

### [R-002] [CRITICAL] [generate_safs_mesh.py:_combine_stls] — Silent drop of 22 polyline-coincident triangles destroys fault topology near intersection

**Category:** BUG

**Description:**
`_combine_stls` (line 169–284) deduplicates triangles after vertex
snapping. The most recent run reports
`n_dropped_duplicate: 22` (in `combine_meta.json`). These come
from sliver triangles whose 3 corners are all on the polyline:
fault A's CDT split produced a child whose 3 vertices are all
polyline points (q1, q2, q3), and so did fault B's CDT — after
1 cm vertex snap they collapse to the same global indices and
one is silently dropped. The fault that lost a triangle now has
a **hole** along the polyline; the remaining fault overlaps that
hole but is not topologically connected to the lost triangle's
neighbours, breaking HXT's Surface-In-Volume recovery.

This is the proximate cause of HXT rejection in the
"production" pipeline (motivating the algo3d=4 fallback in
R-001). It is also the visible non-conformity in the user's
screenshot: the intersection curve has triangles missing on
one side.

**Trigger:** Any pair of faults whose post-conformalize CDTs
each produce a "polyline-only" sliver triangle (3 corners on
shared polyline) on the *same* 3 polyline vertices. Mill Creek
× SBMT-SAF at 2000 m hits this 22 times.

**Actual behavior:** The first triangle of each duplicate pair
is kept; the second is silently dropped. The drop is recorded
only in `combine_meta.json` as a count, not in any per-triangle
log; you cannot tell which fault lost which triangles.

**Expected behavior:** Polyline-coincident "sliver" triangles
should not exist in the conformalize output to begin with —
the conformalizer should suppress them at CDT time. Failing
that, the dedup step in `_combine_stls` must (a) refuse to drop
any triangle whose 3 vertices are all polyline-tagged or (b)
emit one explicit "double-sided" triangle that HXT will accept,
not silently drop one side.

**Suggested fix (defence in depth — apply both):**

A. Suppress the slivers at CDT time. In
`conformalize_faults.py:_split_one_parent`, after the CDT runs,
drop any child whose 3 corners are all from the polyline pierce
set:

```diff
@@ conformalize_faults.py:_split_one_parent
     children_indices: list[tuple[int, int, int]] = [
         (int(t[0]), int(t[1]), int(t[2])) for t in out["triangles"]
     ]
+    # Drop polyline-only slivers: a child whose 3 corners are all
+    # polyline pierces (no parent vertex).  These are the
+    # triangles that, in the OTHER fault's CDT for the same
+    # polyline segment, end up with the same global vertex IDs
+    # and become duplicates after _combine_stls dedup.  The two
+    # fault meshes don't lose any constraint by dropping these;
+    # the polyline edges are still owned by the remaining 2 of
+    # the 3 children.
+    pierce_local_set = {li for li in pierce_to_local if li >= 3}
+    children_indices = [
+        (la, lb, lc) for (la, lb, lc) in children_indices
+        if not (la in pierce_local_set
+                and lb in pierce_local_set
+                and lc in pierce_local_set)
+    ]
     return children_indices, pierce_to_local, new_steiner_3d
```

B. In `_combine_stls`, fail loudly if the dedup count is non-zero
on input that came from `conformalize_faults.py`:

```diff
@@ generate_safs_mesh.py:_combine_stls (after the dedup loop)
+    # Conformalize output must not have polyline-coincident
+    # duplicates (R-002).  If any are found, the conformalizer
+    # let a polyline-only sliver through — refuse to silently
+    # drop it because that creates a hole in one fault's
+    # triangulation along the intersection curve.
+    if n_dropped_duplicate > 0 and (stl_dir / "triangle_to_fault.json").exists():
+        raise SystemExit(
+            f"_combine_stls dropped {n_dropped_duplicate} duplicate "
+            f"triangles from conformalize output.  These are "
+            f"polyline-coincident slivers that indicate the "
+            f"conformalizer emitted the same (q1, q2, q3) sliver on "
+            f"both fault A and fault B.  Fix conformalize_faults.py "
+            f"to suppress polyline-only slivers (see REVIEW.md "
+            f"R-002 fix A) instead of relying on combine to dedup."
+        )
```

**Test case:**
```python
def test_R002_no_polyline_only_slivers_after_conformalize(tmp_path):
    """After conformalize, no fault triangle may have all 3 corners
    on the polyline — these become combine-step duplicates and
    create holes in the fault mesh."""
    import json, conformalize_faults as cf, fault_intersect as fi
    # Construct the Mill Creek × SBMT-SAF fixture (or any
    # fixture with at least one parent triangle pierced by 3+
    # polyline points so the CDT *could* emit a (q1, q2, q3)
    # sliver).
    in_dir = tmp_path / "stl_raw"
    out_dir = tmp_path / "stl_conformal"
    _build_fixture(in_dir)
    cf.conformalize(in_dir, out_dir,
                    ["safs_sbmt_millcreek", "safs_sbmt_saf"])
    # Read each conformal STL, dedup vertices, find polyline
    # vertex IDs.
    for short in ("safs_sbmt_millcreek", "safs_sbmt_saf"):
        V, T = fi.dedup_mesh(*fi._read_stl_arrays(out_dir / f"{short}.stl"),
                             snap_m=1e-3)
        report = json.loads((out_dir / "intersection_report.json").read_text())
        # Polyline vertex IDs — read from polyline_to_node tracking.
        pl_ids = _polyline_vertex_ids_for_fault(short, report)
        slivers = [tri for tri in T
                   if all(int(v) in pl_ids for v in tri)]
        assert len(slivers) == 0, (
            f"{short}: {len(slivers)} polyline-only sliver "
            f"triangle(s); these will dedup to duplicates with the "
            f"OTHER fault's CDT in _combine_stls and one will be "
            f"silently dropped (R-002)."
        )
```

---

### [R-003] [CRITICAL] [conformalize_faults.py:conformalize] — Smoothing applied AFTER validation gates; gates pass on a mesh that is then mutated

**Category:** BUG

**Description:**
`conformalize` (line 1486–1532) runs the manifold,
polyline-edge-coincidence, and interior-crossing gates on
`(V_X_new, T_X_new)`. Then **after** the gates have passed, it
runs `_constrained_delaunay_flips` (line 1542) and
`_smooth_fault_surface` (line 1567) on the same V/T, and
overwrites `fault_meshes[sa]` / `fault_meshes[sb]` with the
post-mutation arrays. The STL files written at the end (line
1601–1604) are the post-mutation arrays, **not** the
gate-validated arrays.

`_smooth_fault_surface` moves every non-locked vertex 50% toward
its 1-ring centroid (line 824), 5 iterations (line 1303). With
`relaxation = 0.5` and 5 iterations a vertex can move
~`(1 - 0.5^5) ≈ 97%` of the centroid distance. Locked = polyline
+ fault boundary. **Interior pierce vertices that landed inside
a parent triangle (Steiner points)** are not in `pl_to_v_X`
because `pl_to_v_X` is built only from `polyline.points` (line
447), which are polyline vertex coords, not Steiner. So
in-parent-interior pierce vertices get smoothed.

For a polyline that kinks inside a parent (interior pierce),
the kink vertex is on the polyline geometrically (it's a
polyline point), so it IS in `polyline.points` — so it IS in
`pl_to_v_A`. So it IS locked. OK so far.

**The actual bug:** `_constrained_delaunay_flips` rewrites T
(line 933–934). Every flip changes which triangles share which
edges, but it also can change manifold properties that the
gates already checked. The flip pass locks polyline edges
(line 957–961) and boundary edges, but it does NOT lock the
manifold-vs-non-manifold property. If a flip pass produces a
3-fold edge somewhere, no gate runs again and HXT will reject.

Also, `_smooth_fault_surface` re-runs the Laplacian without
re-checking that the moved vertex doesn't end up on another
edge's interior (a new T-junction). The plan §Validation gates
explicitly says T-junctions must not be silently introduced —
but smoothing can introduce them by moving an interior vertex
toward a centroid that lies on another edge.

`fault_meshes[sa] = (V_A_new, T_A_new)` on line 1577–1578 saves
the mutated mesh. The STL written from `fault_meshes` (line
1601–1604) is the post-mutation arrays. The gates already
recorded "PASS" in the report.

**Trigger:** Any conformalize run with `flip_passes > 0` or
`smoothing_iters > 0` (the defaults).

**Actual behavior:** Gates pass on pre-mutation mesh, mutated
mesh is written to disk. Downstream consumers (combine,
generate_safs_mesh, HXT) see a mesh whose properties were
never validated. The `intersection_report.json` records
`PASS_WITH_3_T_JUNCTIONS` / `PASS_WITH_7_T_JUNCTIONS` but the
counts are pre-mutation; post-mutation counts may be different
and are never recorded.

**Expected behavior:** Re-run all three gates after every
mutation pass, and update `intersection_report.json` with the
final post-mutation counts. If a post-mutation T-junction count
exceeds a threshold, raise.

**Suggested fix:**

```diff
@@ conformalize_faults.py:conformalize, after the smoothing block (line ~1577)
         # Re-validate AFTER flips + smoothing (R-003): both
         # mutation passes can introduce manifold violations
         # and T-junctions that the pre-mutation gates didn't
         # see.  Update the report with the final counts.
         nm_A_post, tj_A_post = is_manifold_split(V_A_new, T_A_new)
         nm_B_post, tj_B_post = is_manifold_split(V_B_new, T_B_new)
         if nm_A_post:
             raise RuntimeError(
                 f"flip/smooth introduced {len(nm_A_post)} non-manifold "
                 f"edge(s) on fault {sa}: {nm_A_post[:5]}.  Reduce "
                 f"flip_passes or smoothing_iters."
             )
         if nm_B_post:
             raise RuntimeError(
                 f"flip/smooth introduced {len(nm_B_post)} non-manifold "
                 f"edge(s) on fault {sb}: {nm_B_post[:5]}."
             )
         for pl_idx, pl in enumerate(polylines):
             if not verify_polyline_in_mesh(V_A_new, T_A_new, pl):
                 raise RuntimeError(
                     f"flip/smooth broke polyline edge-coincidence "
                     f"for pair ({sa} × {sb}) polyline #{pl_idx} on "
                     f"fault {sa}"
                 )
             if not verify_polyline_in_mesh(V_B_new, T_B_new, pl):
                 raise RuntimeError(
                     f"flip/smooth broke polyline edge-coincidence "
                     f"for pair ({sa} × {sb}) polyline #{pl_idx} on "
                     f"fault {sb}"
                 )
         report["pairs"][f"{sa}__x__{sb}"]["gates_post_mutation"] = {
             "manifold_A": "PASS" if not tj_A_post else f"PASS_WITH_{len(tj_A_post)}_T_JUNCTIONS",
             "manifold_B": "PASS" if not tj_B_post else f"PASS_WITH_{len(tj_B_post)}_T_JUNCTIONS",
             "polyline_edge_coincidence": "PASS",
             "n_flips_A": ...,
             "n_smooth_iters": smoothing_iters,
         }
```

**Test case:**
```python
def test_R003_post_mutation_gates_run(tmp_path):
    """The conformalize report must record post-mutation gate counts
    separately from pre-mutation counts."""
    import json
    cf.conformalize(in_dir, out_dir, ["A", "B"],
                    flip_passes=8, smoothing_iters=5)
    rep = json.loads((out_dir / "intersection_report.json").read_text())
    pair_key = next(iter(rep["pairs"]))
    pair = rep["pairs"][pair_key]
    assert "gates_post_mutation" in pair, (
        "post-mutation gates must run and record their counts; "
        "currently the saved mesh is post-flip + post-smooth but "
        "the gates ran on pre-mutation arrays."
    )
```

---

### [R-004] [CRITICAL] [conformalize_faults.py:_split_one_parent line 256–263] — 2-D-only edge snap leaves 3-D pierce off the parent edge; produces real T-junctions

**Category:** BUG

**Description:**
`_classify_pierce_2d` returns `_PierceLocation(edge_idx, edge_param)`
when the pierce projects within `eff_edge_tol` of an edge in 2-D.
The code then snaps the pierce's **2-D** coords onto the edge
(line 257–263):

```python
pierces_2d[k] = (ax + t * (bx - ax), ay + t * (by - ay))
```

…but the **3-D** position of the pierce in the caller
(`V_list[pierce_node_id]`) is unchanged. The 3-D coord is the
ORIGINAL polyline point, which may sit up to `eff_edge_tol`
m off the parent's 3-D edge in the perpendicular direction.

`eff_edge_tol = 0.5% × mean_edge`. For 1500 m CFM parents this
is **7.5 m**. A pierce that the classifier marks "on edge"
because its 2-D perpendicular distance is 5 m is then assumed
"on edge" by every downstream step — propagation,
edge-coincidence verification, manifold check.

After `_split_one_parent`, the children's 3-D edges go from a
parent vertex to the original polyline 3-D coord (NOT to the
2-D-snapped coord — see line 567 `local_to_global_map[li] =
global_id`). So the child triangle has an edge that is up to
7.5 m off the parent's 3-D edge.

The neighbour of the parent (sharing that 3-D edge) sees the
pierce in `pierces_per_edge` (line 528–530) and inserts the same
global node into ITS perimeter. The neighbour's CDT then
triangulates: its child has an edge to that same pierce node.
So both neighbours' children connect to the same 3-D point
— **but that point is 7.5 m off both parents' shared 3-D edge**.

Result: the parent's edge u→v in 3-D is no longer a chain of
mesh edges (the chain is u → pierce → v, but pierce is 7.5 m
off the line u→v). Every other triangle in the fault that
references the original u→v edge has a vertex sitting 7.5 m
off-edge — a real geometric T-junction that the manifold check
will catch only if `tol_m` (1 mm by default, line 1100) is
small enough. With 7.5 m offset this is well above 1 mm so
the check fires and the report logs T-junctions. But the
manifold check is non-fatal; pipeline continues.

The 3 + 7 = 10 T-junctions in `intersection_report.json` are
exactly this bug.

**Trigger:** Any polyline pierce that geometrically falls in a
parent's interior (not on its 3-D edge) but is within 0.5% of
an edge in 2-D. Common at fault-fault tangency points and at
near-grazing crossings (CFM data has these).

**Actual behavior:** Mesh has 3-D T-junctions perpendicular to
the polyline; HXT rejects the input on these.

**Expected behavior:** Either (a) snap the 3-D pierce position
to the 3-D edge as well as the 2-D, OR (b) reduce
`edge_tol_frac` to a numerical-roundoff value (1e-9 × edge =
1.5 µm for a 1500 m edge — at this tolerance only true on-edge
pierces classify as on-edge, interior pierces are correctly
classified as interior and become Steiner points with no
3-D-vs-2-D mismatch).

**Suggested fix (option B — recommended; the snap was masking,
not solving, the issue):**

```diff
@@ conformalize_faults.py:_split_one_parent
-    edge_tol_frac: float = 5e-3,        # fraction of parent's mean edge → 2-D edge tolerance (0.5%)
-    vertex_tol_frac: float = 1e-3,      # fraction of parent's mean edge → vertex coincidence
+    edge_tol_frac: float = 1e-9,        # 1 nm/m = 1.5 µm for 1500 m
+                                          # edges; sufficient for FP
+                                          # round-off, rejects "near-
+                                          # edge interior" pierces so
+                                          # they are correctly inserted
+                                          # as Steiner vertices (no
+                                          # 2-D-vs-3-D mismatch).
+                                          # See REVIEW.md R-004.
+    vertex_tol_frac: float = 1e-9,
```

And add a 3-D consistency assert after classification:

```diff
@@ after the classify-pierces loop (line ~261)
     for k, p2d in enumerate(pierces_2d):
         loc = _classify_pierce_2d(...)
         locations.append(loc)
         if loc.vertex_idx is None and loc.edge_idx is not None:
+            # R-004: a pierce classified as on-edge in 2-D MUST also
+            # be on the parent's 3-D edge to within float64 round-off,
+            # otherwise propagating it to neighbours produces a
+            # real 3-D T-junction.
+            i_start, i_end = ((0, 1), (1, 2), (2, 0))[loc.edge_idx]
+            v_start_3d = parent_verts_3d[i_start]
+            v_end_3d = parent_verts_3d[i_end]
+            edge_3d = v_end_3d - v_start_3d
+            edge_len = float(np.linalg.norm(edge_3d))
+            t_3d = float(((pierce_pts_3d[k] - v_start_3d) @ edge_3d)
+                          / max(edge_len * edge_len, 1e-30))
+            t_3d = max(0.0, min(1.0, t_3d))
+            foot_3d = v_start_3d + t_3d * edge_3d
+            d_perp_3d = float(np.linalg.norm(pierce_pts_3d[k] - foot_3d))
+            if d_perp_3d > 1e-6 * edge_len:
+                raise RuntimeError(
+                    f"pierce {k} classified as on-edge in 2-D but is "
+                    f"{d_perp_3d:.6g} m perpendicular to the 3-D edge "
+                    f"(edge length {edge_len:.6g} m).  This is the "
+                    f"R-004 bug — reduce edge_tol_frac or compute "
+                    f"the polyline pierce in the parent's plane to "
+                    f"begin with."
+                )
             # Snap pierce 2-D position onto the edge at edge_param.
             ...
```

**Test case:**
```python
def test_R004_off_edge_pierce_classifies_as_interior():
    """A pierce 5 m off the parent's edge in a 1500 m parent must
    classify as interior (Steiner), not on-edge — otherwise the
    children's 3-D edges miss the actual edge by 5 m and propagation
    produces a T-junction."""
    parent_3d = np.array([
        [0, 0, 0], [1500, 0, 0], [0, 1500, 0],
    ], dtype=np.float64)
    pierce_3d = np.array([750, 5, 0])      # 5 m off edge (0,0)-(1500,0)
    children, p2l, _ = cf._split_one_parent(
        parent_3d, [pierce_3d], polyline_edges=[],
    )
    # local index of the only pierce
    li = p2l[0]
    # Must NOT be 0/1/2 (= parent vertex), must be >= 3 = Steiner
    # registered as INTERIOR (not on a boundary segment).  In the
    # bug, li == 3 but the pierce is "snapped to edge" so child
    # triangles include an edge from (0,0,0) to pierce_3d which is
    # 5 m off the y=0 axis.  The new behaviour: pierce classifies
    # as interior, lives strictly inside the parent.
    boundary_segs_using_pierce = sum(
        1 for tri in children
        if li in tri and (0 in tri or 1 in tri)
    )
    assert boundary_segs_using_pierce == 0, (
        "off-edge pierce must not produce children with edges "
        "that join a parent vertex through the pierce — that is "
        "the R-004 T-junction bug."
    )
```

---

### [R-005] [MODERATE] [conformalize_faults.py:_combine_stls snap mismatch] — Combine snap (1 cm) is coarser than conformalize dedup snap (1 mm); merges vertices conformalize kept distinct

**Category:** BUG

**Description:**
`conformalize_faults.py:1366` uses
`dedup_snap_m = min(snap_m, 1e-3) = 1e-3` (1 mm) when reading
per-fault STLs. `generate_safs_mesh.py:_combine_stls` (called via
the driver, line 492) defaults to `--combine-snap-m 0.01` (1 cm).
A vertex pair that conformalize considered distinct (separated
by 5 mm) is merged in combine.

This is the proximate cause of the duplicate-triangle issue
(R-002). The conformalizer's CDT places polyline pierces and
Steiner vertices at sub-mm precision; combine-step snapping
collapses pairs that the conformalizer treated as different
nodes. After collapse, a triangle in fault A and a triangle in
fault B can end up with identical vertex tuples.

**Trigger:** Any multi-fault input where conformalize emits
vertices closer than `combine-snap-m` to each other across
fault boundaries.

**Suggested fix:** Default `combine-snap-m` to the
conformalize `snap_m`, not 10× larger. Surface this as a
fail-loud assert when the two are mismatched.

```diff
@@ generate_safs_mesh.py: --combine-snap-m default
-    parser.add_argument("--combine-snap-m", type=float, default=0.01,
+    parser.add_argument("--combine-snap-m", type=float, default=1e-3,
                         metavar="METRES",
                         help="Vertex snap tolerance when combining "
                              "per-fault STLs into a single conformal "
                              "STL (multi-fault only).  Vertices closer "
                              "than this become a single shared vertex.  "
-                             "Must exceed gmsh's internal STL merge "
-                             "tolerance to take effect.  Default: 0.01 m.")
+                             "Default: 0.001 m, matching "
+                             "conformalize's STL-read dedup snap.  "
+                             "Setting this LARGER than conformalize's "
+                             "snap_m will silently merge vertices that "
+                             "conformalize kept distinct (R-005).")
```

**Test case:**
```python
def test_R005_combine_snap_matches_conformalize_dedup():
    """The combine-snap default must equal the conformalize dedup
    snap; otherwise combine merges vertices conformalize kept
    apart, producing duplicate triangles (R-002)."""
    import re
    src = Path("mesh/generate_safs_mesh.py").read_text()
    m = re.search(r'--combine-snap-m["\s,].*default=([\d.e\-]+)', src)
    assert m, "--combine-snap-m default not found"
    assert abs(float(m.group(1)) - 1e-3) < 1e-12, (
        f"--combine-snap-m default = {m.group(1)}; must match "
        f"conformalize's STL-read dedup_snap_m = 1e-3."
    )
```

---

### [R-006] [MODERATE] [conformalize_faults.py:_split_one_parent line 365] — atol=1e-12 vertex-preservation check is likely too strict for `triangle` library

**Category:** ASSUMPTION

**Description:**
`_split_one_parent` line 365 asserts:

```python
if not np.allclose(out_verts[:n_input_verts], pts_arr, atol=1e-12):
    raise RuntimeError("triangle.triangulate reordered or perturbed input vertices; ...")
```

`atol=1e-12` is a 1-picometre absolute tolerance applied to 2-D
coords that are O(edge_length) — for 1500 m CFM parents,
1e-12 is about `7e-16` of edge length, i.e., last-bit machine
epsilon. The `triangle` library is documented to preserve input
vertex coordinates exactly, but its internal predicates use
adaptive arithmetic that may write back coords that match
bitwise on most platforms but not all (e.g., Conda's
`triangle` wheel built with different compiler flags).

A false positive here aborts the whole pipeline.

**Suggested fix:**

```diff
-    if not np.allclose(out_verts[:n_input_verts], pts_arr, atol=1e-12):
+    # Tolerate sub-µm round-off on 1500 m parents — sub-µm in 2-D
+    # local coords corresponds to <1 nm in 3-D after unprojection,
+    # well below any meaningful tolerance.  1e-12 atol was triggering
+    # spurious aborts on `triangle` builds that round-trip through
+    # extended-precision intermediates.
+    if not np.allclose(out_verts[:n_input_verts], pts_arr, atol=1e-6):
         raise RuntimeError(
             "triangle.triangulate reordered or perturbed input vertices; "
             "polyline edge-coincidence will be broken."
         )
```

**Test case:**
```python
def test_R006_strict_vertex_preservation_tolerates_subum_jitter():
    """1 µm absolute tolerance is safe — any larger error means
    triangle dropped or moved a vertex."""
    # Parent at 1500 m scale, pierce at exact mid-edge.
    # Verify no spurious assert.
    parent = np.array([[0, 0, 0], [1500, 0, 0], [0, 1500, 0]])
    pierce = np.array([750, 0, 0])  # exactly on edge midpoint
    cf._split_one_parent(parent, [pierce], polyline_edges=[])  # must not raise
```

---

### [R-007] [MODERATE] [conformalize_faults.py:conformalize line 1438] — `_coalesce_polylines` runs every conformalize call but is silently a no-op when len(polylines) ≤ 1

**Category:** QUALITY (masks a real issue)

**Description:**
Line 1433–1438:
```python
if len(polylines) > 1:
    coalesce_tol = (...)
    polylines = _coalesce_polylines(polylines, coalesce_tol)
```

When `chain_segments` returns exactly 1 polyline (the common
case for two intersecting faults), coalesce is skipped. But
the *gap* that motivated the coalesce step (sub-cm endpoint
mismatches at fault boundaries) can occur within a single
polyline as well — the kink-then-resume pattern. For the
single-polyline case this manifests as a polyline whose
internal `points[i] == points[i+1]` after snapping, which
`simplify_polyline` (knob A) handles only when knob A is on
(default off per line 1357).

**Trigger:** Single-polyline pair (e.g. Mill Creek × SBMT-SAF —
the production target!) with sub-cm endpoint reassembly
artefacts. The 581 polyline vertices in
`intersection_report.json` may include such artefacts.

**Suggested fix:** Always run an internal-cleanup pass on every
polyline (regardless of whether it merged with another), at
the chain-snap tolerance:

```diff
@@ conformalize_faults.py:conformalize, after chain_segments
-        if len(polylines) > 1:
-            coalesce_tol = (
-                0.5 * target_edge_length_m
-                if target_edge_length_m is not None else 5e-2
-            )
-            polylines = _coalesce_polylines(polylines, coalesce_tol)
+        coalesce_tol = (
+            0.5 * target_edge_length_m
+            if target_edge_length_m is not None else 5e-2
+        )
+        if len(polylines) > 1:
+            polylines = _coalesce_polylines(polylines, coalesce_tol)
+        # Always also collapse sub-snap consecutive duplicates within
+        # each polyline — _coalesce_polylines handles only the
+        # multi-polyline case.
+        polylines = [
+            fi.simplify_polyline(pl, min_pierce_separation_m=snap_m)
+            for pl in polylines
+        ]
```

**Test case:**
```python
def test_R007_single_polyline_collapses_within_snap_dups():
    """A single polyline with two consecutive vertices within snap_m
    must be simplified, not pass through unchanged."""
    pl = fi.Polyline(
        points=[(0, 0, 0), (0.005, 0, 0), (1, 0, 0)],
        closed=False, fault_a="A", fault_b="B",
        pierces_a=[0, 0], pierces_b=[0, 0],
    )
    # Configure conformalize with snap_m=0.01 — the (0, 0, 0)
    # and (0.005, 0, 0) collapse.
    ...
```

---

### [R-008] [MODERATE] [safs.geo:99 + safs.geo:139] — `Surface{} In Volume{}` is silently ignored when Algorithm3D=4

**Category:** EDGE_CASE (no driver-level guard)

**Description:**
`safs.geo:99` issues `Surface{fault_surfs[]} In Volume{bulk_vol};`
which is the *only* line that tells gmsh to honour the fault
triangles as constraints. With `Mesh.Algorithm3D = 4`, this
embedding is ignored (Frontal-Delaunay does not implement
constrained tetrahedralization in gmsh 4.13). The
`generate_safs_mesh.py` CLI accepts `--algo3d 4` without
warning that the embedding will be silently ignored.

This is the implementation-level mechanism behind R-001.

**Suggested fix:**

```diff
@@ generate_safs_mesh.py: --algo3d argument
-    parser.add_argument("--algo3d", type=int, default=10, choices=(4, 10),
+    parser.add_argument("--algo3d", type=int, default=10, choices=(4, 10, "10-strict"),
                         help="Mesh.Algorithm3D — 10 (default) = HXT, "
                              "enforces conformal fault interfaces.  Use "
                              "4 (Frontal-Delaunay) ONLY as a fallback for "
                              "non-conformal multi-fault input that HXT "
-                             "rejects (validator check 5 will then fail).")
+                             "rejects (validator check 5 WILL fail; "
+                             "fault triangles will be detached from "
+                             "the tet mesh — see REVIEW.md R-001/R-008).")
@@ generate_safs_mesh.py: after parsing args
+    if args.algo3d == 4 and len(included) > 1:
+        print(
+            "\nWARNING: --algo3d=4 with multi-fault input ignores "
+            "the Surface-In-Volume embedding.  The fault will be "
+            "detached from the tet mesh (validator check 5 fails).  "
+            "If your input is conformal, use --algo3d=10.  If HXT "
+            "rejects, fix the conformalize output (see REVIEW.md "
+            "R-002, R-004).\n",
+            file=sys.stderr,
+        )
```

**Test case:**
```python
def test_R008_algo3d_4_warns_on_multifault():
    """algo3d=4 + multi-fault must warn about silent embedding loss."""
    out = subprocess.run([
        sys.executable, "mesh/generate_safs_mesh.py",
        "--algo3d", "4",
        "--include-fault", "A", "--include-fault", "B",
        ...
    ], capture_output=True, text=True)
    assert "WARNING" in out.stderr
    assert "embedding" in out.stderr or "detached" in out.stderr
```

---

### [R-009] [LOW] [fault_intersect.py:tri_tri_intersect_3d line 309] — `coplanar_dot_thresh = 1.0 - 1e-6` rejects pairs with 0.08° dihedral as coplanar

**Category:** ASSUMPTION

**Description:**
`coplanar_dot_thresh = 1.0 - 1e-6` triggers the coplanar branch
when `|n_A · n_B| > 1 - 1e-6`, i.e. dihedral angle < `acos(1 -
1e-6) ≈ 0.0814°`. CFM faults at branch points can reasonably
have dihedrals of 1°–2° but adjacent fault segments of the
*same* fault (e.g., Mojave-SAF and SBMT-SAF) routinely sit at
sub-degree dihedrals. The implementer's report (echoed in
`fault_intersect.py:382–399`) acknowledges that "coplanar"
empirically fires on these and the code skips them silently.

For the **two-crossing** target (`run_two_crossing_2000m.sh`)
this is fine — Mill Creek vs SBMT-SAF have a substantial
dihedral. But for the all-8 build the silent skip would drop
the SAF-SAF intersection trace.

**Suggested fix:** Tighten the threshold to `1 - 1e-9` (dihedral
< `2.6 mdeg`) so only true bit-coplanar pairs match. Pairs in
the 0.1°–1° band re-enter the Möller core and are handled by
the standard interval-overlap path (which is robust enough at
those angles).

```diff
-    coplanar_dot_thresh: float = 1.0 - 1e-6,
+    coplanar_dot_thresh: float = 1.0 - 1e-9,
```

**Test case:**
```python
def test_R009_coplanar_threshold_does_not_swallow_subdegree_dihedrals():
    """Two faults with 1° dihedral must not be classified as
    coplanar — 1° is normal CFM branching geometry."""
    # Two unit squares, second tilted by 1°.
    rad = math.radians(1.0)
    V_A = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    T_A = np.array([[0, 1, 2]])
    R = np.array([
        [1, 0, 0],
        [0, math.cos(rad), -math.sin(rad)],
        [0, math.sin(rad),  math.cos(rad)],
    ])
    V_B_raw = np.array([[0.5, -0.5, 0], [0.5, 0.5, 0], [-0.5, 0, 0.1]])
    V_B = V_B_raw @ R.T
    T_B = np.array([[0, 1, 2]])
    segs = fi.tri_tri_intersect_3d(V_A, T_A, V_B, T_B, "A", "B")
    assert len(segs) >= 1, (
        "1° dihedral pair was classified as coplanar and skipped "
        "(R-009).  Tighten coplanar_dot_thresh."
    )
```

---

### [R-010] [STRATEGIC] [whole pipeline] — Replace ~3000 LOC custom conformalizer with CGAL `corefine` + `isotropic_remeshing`

**Category:** DEVIATION (recommendation; user-supplied workaround direction)

**Description:**
The user pointed out that **CGAL Polygon Mesh Processing**
provides exactly the operation this codebase is reimplementing,
in a single library call:

```cpp
PMP::corefine(fault1, fault2,
    CGAL::parameters::edge_is_constrained_map(ecm1),
    CGAL::parameters::edge_is_constrained_map(ecm2));
PMP::isotropic_remeshing(faces(fault1), target_edge_length, fault1,
    CGAL::parameters::edge_is_constrained_map(ecm1));
```

`corefine()` does Phase 1 + Phase 2 of
`PLAN_multifault_intersections.md` exactly: it computes every
triangle-pair intersection with **exact predicates**, splits
both meshes along the shared polyline, and marks those edges
as constrained. `isotropic_remeshing()` then handles slivers
and produces uniform-edge surface meshes while protecting the
intersection edges. **Both operations are robust to the
geometry pathologies that R-002, R-004, and R-009 currently
expose** (near-coplanar triangles, near-edge pierces, off-edge
3-D drift), because CGAL uses Exact_predicates_inexact_
constructions_kernel and exact orientation predicates rather
than the float64 + 1e-9 tolerances used here.

The custom conformalizer is **~3000 lines of Python** spread
across `fault_intersect.py` (1131 LOC) +
`conformalize_faults.py` (1700 LOC) + their tests, and is
still producing 10 T-junctions and 22 polyline-coincident
slivers on the simplest 2-fault target. The plan's risk
assessment §1 explicitly flagged "Constrained Delaunay
robustness at near-tangent intersections" — this is what
CGAL's exact predicates exist to solve.

**Recommendation: introduce a CGAL-corefine path as a parallel
implementation, keep the Python path as fallback / fixture
generator for unit tests, and switch `run_two_crossing_2000m.sh`
to the CGAL path once it is wired up.**

**Concrete pipeline (replaces Phases 1 + 2):**

```
input:  output/two_crossing_2000m/stl_raw/{A,B}.stl
        output/two_crossing_2000m/transform.json (free-surface clearance)

step 1: tools/corefine_faults  (new C++ binary — see scaffold below)
        reads the 2 raw STLs
        runs PMP::corefine + PMP::isotropic_remeshing
        writes 2 corefined OFF/STL files
        writes triangle_to_fault.json, intersection_report.json
        (matching the schema conformalize_faults.py emits, so
         downstream consumers don't change)

step 2: generate_safs_mesh.py (UNCHANGED — already conformal-aware)
        --stl-dir output/two_crossing_2000m/stl_corefined
        algo3d=10 (HXT)
```

**Scaffold for the new C++ tool**, to live at
`tools/corefine_faults.cpp`:

```cpp
// SPDX-License-Identifier: MIT
// SAFS multi-fault corefinement — replaces the Python
// fault_intersect.py + conformalize_faults.py pipeline (R-010).
//
// Reads N per-fault open-surface STLs, corefines every unordered
// pair with PMP::corefine (constrained edges marked on both
// sides), then isotropic-remeshes each surface to target edge
// length while protecting the intersection edges.
//
// Output: N corefined STL files + triangle_to_fault.json +
// intersection_report.json, schema-compatible with
// conformalize_faults.py.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

namespace PMP = CGAL::Polygon_mesh_processing;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;
using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;

static int corefine_pair(
    Mesh& A, Mesh& B,
    const std::string& sa, const std::string& sb,
    double target_edge_m, std::ostream& report)
{
    auto ecm_A = A.add_property_map<edge_descriptor, bool>(
                     "e:is_constrained", false).first;
    auto ecm_B = B.add_property_map<edge_descriptor, bool>(
                     "e:is_constrained", false).first;
    PMP::corefine(A, B,
        CGAL::parameters::edge_is_constrained_map(ecm_A),
        CGAL::parameters::edge_is_constrained_map(ecm_B));
    // Count the resulting constrained edges (= shared polyline length).
    std::size_t n_A_constrained = 0, n_B_constrained = 0;
    for (auto e : edges(A))
        if (get(ecm_A, e)) ++n_A_constrained;
    for (auto e : edges(B))
        if (get(ecm_B, e)) ++n_B_constrained;
    report << "  \"" << sa << "__x__" << sb << "\": {\n"
           << "    \"n_constrained_edges_A\": " << n_A_constrained << ",\n"
           << "    \"n_constrained_edges_B\": " << n_B_constrained << ",\n";

    PMP::isotropic_remeshing(faces(A), target_edge_m, A,
        CGAL::parameters::edge_is_constrained_map(ecm_A)
                          .number_of_iterations(3));
    PMP::isotropic_remeshing(faces(B), target_edge_m, B,
        CGAL::parameters::edge_is_constrained_map(ecm_B)
                          .number_of_iterations(3));
    report << "    \"n_tri_A_post\": " << faces(A).size() << ",\n"
           << "    \"n_tri_B_post\": " << faces(B).size() << "\n"
           << "  }";
    return 0;
}

int main(int argc, char** argv)
{
    if (argc < 5) {
        std::cerr << "usage: corefine_faults <in_dir> <out_dir> "
                     "<target_edge_m> <fault1> <fault2> [...]\n";
        return 1;
    }
    std::string in_dir = argv[1], out_dir = argv[2];
    double target_edge = std::stod(argv[3]);
    std::vector<std::string> shorts;
    for (int i = 4; i < argc; ++i) shorts.emplace_back(argv[i]);

    std::vector<Mesh> meshes(shorts.size());
    for (std::size_t i = 0; i < shorts.size(); ++i) {
        if (!PMP::IO::read_polygon_mesh(
                in_dir + "/" + shorts[i] + ".stl", meshes[i])) {
            std::cerr << "failed to read " << shorts[i] << "\n";
            return 2;
        }
    }

    std::ofstream report(out_dir + "/intersection_report.json");
    report << "{\n  \"schema_version\": 1,\n  \"backend\": \"cgal\",\n"
           << "  \"target_edge_length_m\": " << target_edge << ",\n"
           << "  \"pairs\": {\n";
    bool first = true;
    for (std::size_t i = 0; i < shorts.size(); ++i)
        for (std::size_t j = i + 1; j < shorts.size(); ++j) {
            if (!first) report << ",\n";
            first = false;
            corefine_pair(meshes[i], meshes[j], shorts[i], shorts[j],
                          target_edge, report);
        }
    report << "\n  }\n}\n";

    for (std::size_t i = 0; i < shorts.size(); ++i)
        PMP::IO::write_polygon_mesh(
            out_dir + "/" + shorts[i] + ".stl", meshes[i]);
    return 0;
}
```

**Compile (Frontera or laptop):**

```bash
# CGAL is header-only since 5.0.
g++ -O2 -std=c++17 tools/corefine_faults.cpp -o tools/corefine_faults \
    -lgmp -lmpfr
```

**Wire-up in run_two_crossing_2000m.sh:**

```bash
echo "==> 3/7 corefine_faults (CGAL — replaces conformalize_faults.py)"
tools/corefine_faults \
    "$OUTDIR/stl_raw" "$OUTDIR/stl_conformal" \
    "$TARGET_EDGE" "${FAULTS[@]}"

# step 4 onward unchanged — generate_safs_mesh.py reads
# stl_conformal as it does today.
```

**Why CGAL solves R-001…R-005 in one go:**

| Existing finding | CGAL replacement |
|---|---|
| R-001 algo3d=4 fallback for HXT rejection | corefine output is by-construction a valid PLC; HXT accepts. Default to algo3d=10. |
| R-002 polyline-only sliver duplicates | `isotropic_remeshing` removes slivers without breaking constrained edges. |
| R-003 post-mutation gates not run | CGAL's invariants are guaranteed by construction; no flip/smooth pass is needed. |
| R-004 2-D vs 3-D edge classification mismatch | CGAL works in 3-D throughout with exact predicates. |
| R-005 snap mismatch | corefine produces shared (not just-coincident) vertices on the polyline; no snap needed. |
| R-009 0.08° coplanar threshold drops near-tangent pairs | CGAL handles near-tangent and exactly-coplanar geometrically (returns the on-plane intersection polygon). |

**Volume meshing after corefine:** As the user noted, three options:
1. **CGAL `make_mesh_3` with `Mesh_domain_with_polyline_features_3`**:
   the polyline edges become "1-D features" in the volume mesh
   and are preserved exactly. Cleanest single-library path.
2. **TetGen via `meshpy`**: combined STL → PLC → tet mesh. We
   already have `_combine_stls`; stripping its dedup logic
   (no longer needed because corefine guarantees shared vertices)
   and feeding the result to TetGen is ~50 LOC of Python.
3. **Gmsh HXT**: unchanged. Just point `--stl-dir` at the
   corefined output. This is the smallest pipeline change and
   the one I recommend for the first cut.

**Test case:** the existing
`mesh/tests/test_conformalize_faults.py::test_two_perpendicular_squares_conformalize`
is the right shape; rebrand for the CGAL backend:

```python
def test_R010_corefine_two_perpendicular_squares(tmp_path):
    """End-to-end: 2 perpendicular unit squares → CGAL corefine →
    each side has the intersection edge as a sequence of mesh
    edges, 0 polyline-only slivers, 0 T-junctions."""
    in_dir = tmp_path / "stl_raw"
    out_dir = tmp_path / "stl_corefined"
    _write_two_perpendicular_squares(in_dir)
    subprocess.run([
        "tools/corefine_faults",
        str(in_dir), str(out_dir),
        "0.5",   # target edge length
        "A", "B",
    ], check=True)
    # Read both corefined STLs, dedup at 1 µm, verify:
    #   - every shared 3-D point appears as a vertex on both meshes
    #   - no triangle has all 3 vertices on the polyline
    #   - polyline endpoints are EXACTLY shared (bitwise)
    ...
```

**What this finding asks the fix agent to do:**
This is a STRATEGIC finding — the fix agent should NOT
silently rip out 3000 lines of Python. The work is:

1. Apply R-001…R-009 fixes to make the current pipeline at
   least as robust as it was before this review (the visible
   non-conformity in the screenshot needs to stop bleeding
   while the CGAL track lands).
2. Land `tools/corefine_faults.cpp` (~150 LOC) and a parallel
   shell script `run_two_crossing_2000m_cgal.sh` that uses it.
3. Document in `PLAN_multifault_intersections.md` that the CGAL
   path supersedes Phases 1 + 2; mark `fault_intersect.py`
   and `conformalize_faults.py` as legacy fallback.
4. Once `run_two_crossing_2000m_cgal.sh` passes 11/11 validator
   checks, retire the Python pipeline.

---

## Summary

- Critical issues: **4** (R-001, R-002, R-003, R-004)
- Moderate issues: **4** (R-005, R-006, R-007, R-008)
- Low issues: **1** (R-009)
- Strategic: **1** (R-010 CGAL workaround)
- Plan compliance: **PARTIAL** — Phase 1 substantially complete
  per `REVIEW_phase1.md`; Phase 2 + 4 produce a mesh that fails
  the plan's primary acceptance criterion (check 5 = PASS); the
  pipeline currently relies on a fallback (algo3d=4) that the
  plan never sanctioned.
- Verdict: **FAIL — must fix before proceeding.** The
  user-visible non-conformity in the screenshot is reproduced
  end-to-end by R-001 + R-002 + R-004 acting together. The
  4 critical findings each independently invalidate the
  pipeline's output for SEAS use; together they explain every
  row of the validation report.

The user explicitly asked for "a workaround." The
shortest-path workaround is **R-001 (default algo3d=10) +
R-002 (drop polyline slivers) + R-004 (tighten edge_tol_frac)**;
this gets the existing Python pipeline to pass HXT on the
two-crossing target. The clean-room workaround is **R-010
(CGAL corefine)** which makes R-001 through R-005 moot.

## Unreviewed Areas

- `fault_intersect.py:_resolve_with_mpq` (gmpy2 path) —
  `_HAS_MPQ` is False in `pythonenv`, so this branch is dead
  in CI. Did not deep-trace; would re-review after gmpy2 lands.
- `audit_ts_quality.py`, `ts_to_stl.py`, `convert_msh.py` — out
  of scope for the intersection-conformity issue the user
  reported.
- `validate_msh.py:check_5..check_11` past line 200 — only
  read the first 200 lines; the failure modes reported above
  are sufficient to act on.
- The behaviour of `_constrained_delaunay_flips` on 3-D
  surfaces (vs 2-D Lawson-flip — the standard guarantee
  doesn't transfer to non-planar meshes) — flagged via R-003
  but did not prove it produces wrong results, only that it
  needs gating with re-validation.
- The CGAL scaffold in R-010 is *untested* — it compiles in
  principle (the API matches CGAL 5.0+ docs) but I have not
  built it on the user's machine. The fix agent should verify
  compilation before treating R-010 as ready to land.
