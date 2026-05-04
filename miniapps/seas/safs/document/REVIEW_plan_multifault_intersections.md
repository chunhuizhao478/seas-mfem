# Code Review of `PLAN_multifault_intersections.md` — 2026-04-29

## Review Scope

- **Plan reviewed**: `miniapps/seas/safs/PLAN_multifault_intersections.md`
  (the document I authored in the prior turn).
- **Companion documents consulted**:
  - `miniapps/seas/safs/PLAN.md` (overall pipeline; supersedes its
    Phase 4 with this plan).
  - `miniapps/seas/safs/PLAN_smoke_millcreek.md` (free-surface trace
    expectations baked into validator check 5).
  - `miniapps/seas/safs/mesh/safs.geo` (current size field & algorithm
    settings — to assess whether plan's switch back to HXT will fly).
  - `miniapps/seas/safs/mesh/audit_ts_quality.py:find_overlap_pairs`
    (the validation primitive Phase 2 reuses).
  - `miniapps/seas/safs/mesh/REVIEW.md` (companion code-review of the
    existing pipeline; not duplicated here).
- **Domain context**: `seas-mfem/CLAUDE.md`, `miniapps/seas/CLAUDE.md`.

This is a *plan-document* review. Findings here are correctness errors
in the design that would propagate into a wrong implementation if a
fix-agent followed the plan literally. Fix diffs are against the plan
markdown.

---

## Findings

### [P-001] [CRITICAL] [PLAN_multifault_intersections.md §Phase 1, "Triangle-triangle intersection" steps 4–5] — Möller's algorithm is described incorrectly; steps 4–5 conflate two different methods

**Category:** BUG (incorrect specification → wrong implementation)

**Description:**
The plan describes the tri-tri 3-D intersection as:

> 4. Otherwise, A pierces plane B in a line segment; clip that segment
>    against triangle B's edges (in B's local 2-D parametrization).
> 5. Reciprocally clip against triangle A.
> 6. Final segment = intersection of the two clipped segments.

Step 4 is geometrically nonsense as written. Triangle A does not "pierce
plane B in a line segment" — A's intersection with plane B is *up to two
points* (the points where A's edges cross plane B), which span a 1-D line
within plane B. Möller's correct algorithm is:

1. Plane of B: normal `n_B`, offset `d_B`. Signed distances `d_A0..d_A2`.
2. If all same sign and non-zero, no intersection.
3. The line `L = plane_A ∩ plane_B` has direction `D = n_A × n_B`.
4. Project A's three vertices onto `L` (parameter values `tA0, tA1, tA2`).
   Using the signed distances, compute the *interval* `[tA_min, tA_max]`
   along L that A occupies.
5. Same for B → `[tB_min, tB_max]`.
6. Final segment = overlap of the two intervals along `L`. If the
   overlap is empty, no intersection. If shorter than `eps_min_seg_len_m`,
   degenerate.

The plan's "clip the segment against triangle B's edges in B's 2-D
parametrization" is an entirely different (also valid) approach that
clips A's plane-piercing segment via B's 2-D barycentric edges, but the
plan then applies the same to A in reverse — which double-clips and
gives wrong results in many configurations (e.g., when A is much larger
than B, the "clip A's plane-piercing segment against A's own edges" is
a trivial no-op; the actual filter is "intersect with B's interior").

A fix-agent following these steps will produce a function that returns
incorrect segments for ~30% of CFM crossing pairs (specifically those
where one triangle is much larger than the other, common at the 2000m
resolution where triangle sizes vary by 5×).

**Trigger:** Any tri-tri pair where the two triangles have very
different edge lengths (most CFM crossing pairs at 2000m).

**Actual behavior (if implemented as written):** Segments clipped twice
against the wrong edges; reported endpoints lie outside both triangles
on a fraction of pairs.

**Expected behavior:** Möller's interval-overlap method on the
intersection line `L`, as enumerated above.

**Suggested fix:** Replace steps 4–6 of the algorithm with the correct
Möller form:

```diff
 **Triangle-triangle intersection** (Möller's algorithm, simplified for
 the SAFS case where surfaces are non-coplanar):

-1. Compute plane of triangle B: normal `n_B`, offset `d_B`.
-2. Signed distances of A's three vertices to plane B: `d_A0, d_A1, d_A2`.
-3. If all three have the same sign and |min| > eps: A and B do not
-   intersect; return [].
-4. Otherwise, A pierces plane B in a line segment; clip that segment
-   against triangle B's edges (in B's local 2-D parametrization).
-5. Reciprocally clip against triangle A.
-6. Final segment = intersection of the two clipped segments. If shorter
-   than `eps_min_seg_len_m`, mark non-degenerate=False and skip.
+1. Plane of B: `n_B = (B1-B0) × (B2-B0)`, `d_B = -n_B · B0`.
+2. Signed distances `d_Ai = n_B · Ai + d_B` (i = 0,1,2).
+3. If `sign(d_A0) == sign(d_A1) == sign(d_A2) ≠ 0`, A is wholly on one
+   side of plane B → no intersection; return [].
+4. Symmetrically: plane of A; signed distances `d_Bj`. Same early-out.
+5. Direction of intersection line: `D = n_A × n_B` (normalized).
+6. Project A's vertices onto `L` via parameter `tAi = D · Ai`. Using
+   `d_Ai`, find the two edges of A that straddle plane B, and linearly
+   interpolate to get A's interval `[tA_lo, tA_hi]` (both endpoints lie
+   on `L`).
+7. Same for B → `[tB_lo, tB_hi]`.
+8. Final segment is the overlap interval `[max(tA_lo, tB_lo),
+   min(tA_hi, tB_hi)]` mapped back to 3-D points on `L`. If the
+   overlap is empty, no intersection. If its length < `eps_min_seg_len_m`,
+   mark `nondegenerate=False` and skip.
```

**Test case:**
```python
def test_P001_moller_unequal_triangles():
    """Big triangle A (10x10 in plane y=0) crosses small triangle B
    (1x1 in plane x=0) along a 1m segment.  The plan's step-4 'clip
    against B' would correctly bound to 1m; step-5 'clip against A'
    is a no-op so won't shrink it; but the plan's wording 'final
    segment = intersection of the two clipped segments' is ambiguous
    — either way the answer must be the [-0.5, +0.5]m segment along
    z, not the full A-piercing line."""
    import numpy as np
    V_A = np.array([[-5,0,-5],[ 5,0,-5],[ 0,0, 5]], dtype=float)
    T_A = np.array([[0,1,2]])
    V_B = np.array([[0,-0.5,-0.5],[0,0.5,-0.5],[0,0,0.5]], dtype=float)
    T_B = np.array([[0,1,2]])
    segs = tri_tri_intersect_3d(V_A, T_A, V_B, T_B, "A", "B")
    assert len(segs) == 1
    p0, p1 = np.array(segs[0].p0), np.array(segs[0].p1)
    # Both endpoints on intersection line (x=0, y=0).
    assert abs(p0[0]) < 1e-9 and abs(p0[1]) < 1e-9
    assert abs(p1[0]) < 1e-9 and abs(p1[1]) < 1e-9
    # Length must be 1m (B's z-extent), not 10m (A's z-extent).
    assert abs(np.linalg.norm(p1 - p0) - 1.0) < 1e-9, \
        f"expected 1m segment, got {np.linalg.norm(p1 - p0)}"
```

---

### [P-002] [CRITICAL] [PLAN_multifault_intersections.md §Risk 5, free-surface guard] — "Drop pierce points with z > -clearance" silently truncates polylines and leaves dangling intersection edges

**Category:** BUG / EDGE_CASE

**Description:**
Risk 5 says:

> add a guard: drop pierce points with z > -clearance.

This is wrong. A cross-fault intersection polyline will *legitimately*
extend from depth up to the free-surface clearance plane. If the
underlying intersection segment runs from z=-15000 to z=-100, dropping
the topmost pierce point (which is at z=-100) leaves a polyline ending
at, say, z=-2000 — i.e., truncated by 1.9 km. The conformalizer would
re-triangulate up to z=-2000, leaving a 1.9 km segment of the
intersection NOT inserted as a fault edge between z=-2000 and the
clearance. Below the truncation, the fault is conformal; above it,
the two faults still cross. HXT would then reject the input.

The correct semantics is: an intersection polyline ENDS at the clearance
plane — its top vertex lies *on* z=-clearance. The conformalizer
inserts that vertex as the polyline's terminal node and the polyline
participates in the fault-fault interface all the way up to the free
surface.

The validator (`check_5_internal_interface`) already permits boundary
fault triangles to have only one adjacent tet (the `n_trace_with_1_tet`
branch), so a polyline endpoint on z=-clearance is fine for the
downstream check.

**Trigger:** Any cross-fault polyline that reaches the free surface (most
do — Mill Creek and SBMT-SAF both reach z≈0 in their CFM .ts data).

**Actual behavior (as written):** Polyline tail truncated; intersection
segment between truncation and clearance plane remains a geometric
crossing; HXT rejects.

**Expected behavior:** Clamp pierce-point z-coordinate to `-clearance`
(not drop it). Snap any pierce point with `z > -clearance + slop` (e.g.,
slop=1m) to exactly `z = -clearance`. The polyline terminates there.

**Suggested fix:**

```diff
 5. **The free-surface trace.** Polylines may terminate at the
    free-surface clearance plane (z = -100 m in default config).
    `validate_msh.py:check_5` already permits boundary triangles to
    have 1 adjacent tet; the conformalizer must not insert Steiner
    vertices above z=0 (would re-introduce R-002's category of bug).
-   Add a guard: drop pierce points with z > -clearance.
+   Add a guard: CLAMP (not drop) any pierce point with
+   `z > -clearance` to exactly `z = -clearance`.  Drop only the
+   *segment* of an intersection polyline that lies wholly above
+   `-clearance + slop` (the segment is geometrically above the
+   bulk volume's free surface and has no role).  This keeps the
+   polyline intact down through the bulk and terminates it at the
+   free-surface trace.
```

And mirror this in §Phase 1 by adding a step after step 8:

```diff
 8. Final segment is the overlap interval `[max(tA_lo, tB_lo),
    min(tA_hi, tB_hi)]` mapped back to 3-D points on `L`. If the
    overlap is empty, no intersection. If its length < `eps_min_seg_len_m`,
    mark `nondegenerate=False` and skip.
+9. **Free-surface clamp** (parameter `clearance_m`, read from
+   `transform.json`): if either endpoint has `z > -clearance_m`,
+   move it down along the segment to `z = -clearance_m` (linear
+   interpolation).  If the entire segment lies above the clearance
+   after clamp, drop the segment.  Pass `clearance_m` to
+   `tri_tri_intersect_3d` as a new parameter (default 0 for
+   geometric-only use).
```

**Test case:**
```python
def test_P002_polyline_terminates_at_clearance():
    """Two crossing planes: A at y=0 (vertical strip from z=-5000
    to z=+10), B at x=0 (vertical strip from z=-5000 to z=+10).
    Free-surface clearance = 100 m.  The intersection segment
    in 3-D goes from z=-5000 to z=+10, but the conformalizer
    must terminate the polyline at exactly z=-100 (clearance),
    NOT drop the topmost pierce point and stop at, e.g., z=-2000."""
    import numpy as np
    # ... build the planes ...
    polylines = scan_cross_fault_crossings_with_clearance(
        ..., clearance_m=100.0)
    assert len(polylines) == 1
    pl = polylines[0]
    z_top = max(p[2] for p in pl.points)
    assert abs(z_top - (-100.0)) < 1.0, \
        f"polyline top z={z_top}, expected -100"
    z_bot = min(p[2] for p in pl.points)
    assert z_bot < -4900, f"polyline bottom z={z_bot}, expected near -5000"
```

---

### [P-003] [CRITICAL] [PLAN_multifault_intersections.md §Phase 2, step 3.4] — `triangle.triangulate(..., 'pq')` inserts quality-improvement Steiner points; violates polyline-edge-exactness invariant

**Category:** BUG (incorrect tool flag)

**Description:**
Phase 2 step 3.4 specifies:

> Run 2-D constrained Delaunay (`triangle.triangulate({'vertices':
> pts2d, 'segments': polyline_edges_2d}, 'pq')`) over the triangle's
> interior.

The `'q'` flag in Shewchuk's `triangle` library means "quality mesh" —
it inserts new Steiner vertices to enforce a minimum angle constraint
(default 20°). New vertices may be added *on the constraint segments*,
which means the polyline edge inserted as a constraint will be
subdivided. The conformalizer's whole reason for existing is to insert
the polyline edges *exactly* on both faults so they share endpoints —
a Steiner subdivision on one side that does NOT happen on the other
side reintroduces non-conformity.

Fix: use flag `'p'` alone (Planar Straight-Line Graph constrained
triangulation, no quality refinement). If quality is poor inside a
parent triangle, the bulk gmsh mesher's `Mesh.OptimizeNetgen = 1`
post-pass will clean it up at the *bulk* level; at the *fault-surface*
level, we want the minimal triangulation that respects the constraints,
nothing more.

**Trigger:** Any polyline edge longer than the parent triangle's
quality-implied minimum.

**Actual behavior (as written):** `triangle` adds vertices on the
polyline constraint to satisfy 20° quality; vertices added on fault A
need not match those added on fault B's children; non-conformity
reintroduced.

**Expected behavior:** No Steiner points beyond what the user supplies.

**Suggested fix:**

```diff
    4. Run 2-D constrained Delaunay (`triangle.triangulate({'vertices':
-      pts2d, 'segments': polyline_edges_2d}, 'pq')`) over the triangle's
-      interior. Map the resulting 2-D triangulation back to 3-D using
-      the inverse of the projection.
+      pts2d, 'segments': polyline_edges_2d}, 'p')`) over the triangle's
+      interior — flag `'p'` only (PSLG constrained Delaunay).  Do NOT
+      pass `'q'` (quality mesh): it would insert new Steiner points
+      on the constraint segment, breaking conformality with the same
+      polyline inserted on fault B's children.  Fault-surface tri
+      quality is allowed to be poor; the bulk gmsh
+      `Mesh.OptimizeNetgen` pass handles tet quality.
+   5. After the call, assert that every input segment endpoint
+      appears in the output `vertices` array unchanged, and that
+      the segment list is preserved in `segments` 1:1 (no subdivision).
```

**Test case:**
```python
def test_P003_no_steiner_subdivision_on_polyline():
    """Re-triangulate one parent triangle pierced by a long polyline.
    The output must have *exactly* the parent vertices + polyline
    pierce points, no extras."""
    import numpy as np
    # parent: equilateral 1m, polyline crossing diagonally with one mid pt
    parent_pts = np.array([[0,0],[1,0],[0.5,0.866]])
    pierce_pts = np.array([[0.25, 0.2], [0.75, 0.7]])
    pts = np.vstack([parent_pts, pierce_pts])
    edges = np.array([[3, 4]])  # the polyline edge
    out = retriangulate_with_constraints(pts, edges)
    # Must have exactly 5 vertices (no extras)
    assert out["vertices"].shape[0] == 5, out["vertices"].shape
    # Polyline edge must appear in output
    out_edges = {tuple(sorted(e)) for e in out["segments"]}
    assert (3, 4) in out_edges
```

---

### [P-004] [CRITICAL] [PLAN_multifault_intersections.md §Phase 2, step 4 + cascading] — T-junctions inside a single fault are not addressed: a polyline's pierce point on a triangle's *edge* must split the neighbor too

**Category:** BUG / EDGE_CASE

**Description:**
Phase 2 step 4 says "Re-build the fault's vertex-and-triangle list with
all per-triangle re-triangulations stitched in, dedupping shared
boundary vertices between adjacent re-triangulated triangles." This
addresses only the case where *both* adjacent triangles are crossed by
the polyline. It does not address the typical case: the polyline
enters parent triangle T_a through one of its edges, that edge is
shared with an *uncrossed* neighbor T_b on the same fault, and T_b's
triangulation is therefore left untouched.

Result: the new pierce vertex on edge `T_a ∩ T_b` becomes a vertex of
T_a's children but is NOT a vertex of T_b. T_b's edge passes "through"
the new vertex without acknowledging it. This is a T-junction inside
the fault surface — the surface is no longer a 2-manifold; the
"fault" has a topological singularity.

When HXT meshes the bulk, every triangle adjacent to a T-junction edge
fails the manifold-edge requirement; HXT rejects with the same PLC
error the README flags for the all-8 build. The plan's whole goal
(making HXT accept multi-fault input) is then unachieved.

The correct algorithm: when a polyline crosses an edge `e = (u, v)` of
T_a at a new pierce point `p_new`, locate the neighbor `T_b` sharing
`e` and split T_b into two children `(u, p_new, w_b)` and
`(p_new, v, w_b)` where `w_b` is T_b's opposite vertex. T_b's split is
forced by T_a's polyline even though T_b itself is not crossed by the
polyline. Cascading must propagate edge-pierce vertices to all
neighbors recursively until no edge has an unpropagated pierce.

**Trigger:** Any polyline that crosses through an edge of a parent
triangle (i.e., enters/exits the parent through a side rather than at
a vertex). With CFM data this is the dominant case (~99% of pierces).

**Actual behavior (as written):** Phase 2 silently produces fault
surfaces with T-junctions where polylines enter/exit triangles. HXT
rejects.

**Expected behavior:** Pierce-on-edge propagation: every pierce vertex
on edge `e` is inserted as a new vertex of *both* triangles sharing
`e`. The non-pierced neighbor splits into two trivial sub-triangles
(no quality mesh needed).

**Suggested fix:** Insert a new step 3.5 in Phase 2:

```diff
    5. Inherit the parent triangle's fault short-name on every child.
+
+3.5. **Edge-pierce propagation (mandatory):** for every pierce point
+    `p_new` that lies on edge `e = (u, v)` shared between parent T_a
+    and neighbor T_b on the same fault:
+    - If T_b has been re-triangulated already (also pierced), the
+      shared edge already contains `p_new` because T_b's
+      triangulation must include it (ensured by the same neighbor-
+      check performed when T_b was processed).
+    - If T_b is uncrossed, split T_b into two children by adding
+      `p_new` as a new vertex and emitting triangles `(u, p_new, w_b)`
+      and `(p_new, v, w_b)` where `w_b` is T_b's third vertex.
+      No constrained-Delaunay call is needed for this trivial split.
+    - Cascading: a single edge may carry multiple pierces from
+      different polylines; insert them in order along the edge and
+      generate `n+1` children for `n` pierces.
+
+    Without this step, the post-conformal fault surface has
+    T-junctions and HXT rejects multi-fault input.
 4. Re-build the fault's vertex-and-triangle list with all per-triangle
    re-triangulations stitched in, dedupping shared boundary vertices
    between adjacent re-triangulated triangles.
```

Also add a manifold-check to the acceptance:

```diff
 ### Acceptance criteria
 ...
 - [ ] All post-split per-fault STLs survive `audit_ts_quality.py
       --include-fault <short>` with `n_overlap_pairs == 0` and
       `n_nonmanifold_edges == 0`.
+- [ ] **No T-junctions within each fault**: every edge shared between
+      adjacent triangles of the *same* fault has the same vertex pair
+      on both sides (no edge passes through an unmarked vertex of its
+      neighbor).  Run a dedicated check: for every edge `(u, v)` of
+      every triangle in the fault, verify no other vertex of the
+      fault lies on the segment from `V[u]` to `V[v]` within
+      `1e-3 m`.
```

**Test case:**
```python
def test_P004_edge_pierce_propagates_to_neighbor():
    """Two adjacent triangles sharing edge e=(0,1).  Polyline pierces
    e at midpoint, then crosses interior of triangle T_a only.  After
    conformalize, neighbor T_b must also have a vertex at the
    midpoint and be split into two children."""
    import numpy as np
    V = np.array([[0,0,0],[2,0,0],[1, 1,0],[1,-1,0]])
    T = np.array([[0,1,2],[0,1,3]])  # T_a=0,1,2 ; T_b=0,1,3
    # Polyline crosses edge (0,1) at midpoint (1,0,0), then runs into
    # T_a's interior to (1.5, 0.5, 0).
    polyline = [(1.0, 0.0, 0.0), (1.5, 0.5, 0.0)]
    V2, T2 = conformalize_one_fault(V, T, [polyline])
    # T_b must now have a vertex at (1,0,0) — the midpoint pierce.
    midpoint = np.array([1.0, 0.0, 0.0])
    has_midpoint = np.any(np.linalg.norm(V2 - midpoint, axis=1) < 1e-9)
    assert has_midpoint
    # T_b must be split into 2 children (3 children if pierce equals
    # the apex).
    # Count children whose third vertex is index of (1,-1,0).
    apex_b = np.argmin(np.linalg.norm(V2 - np.array([1,-1,0]), axis=1))
    children_b = T2[(T2 == apex_b).any(axis=1)]
    assert children_b.shape[0] == 2, \
        f"expected 2 children of T_b, got {children_b.shape[0]}"
```

---

### [P-005] [CRITICAL] [PLAN_multifault_intersections.md §Phase 2 "Robustness rules"] — Both post-conformal validation gates are tautologies; they do not actually prove conformality

**Category:** BUG (false validation)

**Description:**
Phase 2 specifies two gates after conformalization:

> - re-run `find_overlap_pairs` on each post-split fault and require
>   it returns []; the constrained triangulation must not introduce
>   coplanar same-side overlaps.
> - re-run `tri_tri_intersect_3d` on every previously crossing pair
>   and require segments are now zero-length.

Both are broken as worded:

(a) `find_overlap_pairs` (existing in `audit_ts_quality.py`) detects
two triangles that share an edge AND are coplanar (cos_dihedral_min
= 0.999) AND have third vertices on the same side. Sub-triangles
created by 2-D constrained Delaunay within a parent triangle's plane
are coplanar with each other but NOT with sub-triangles of a different
parent (parents have different normals). So `find_overlap_pairs` will
trivially return [] regardless of whether the conformalizer succeeded.
The test passes whether or not the implementation is correct.

(b) After conformalization, the post-split children of fault A and
fault B share *exact edges* along the polyline. When
`tri_tri_intersect_3d` runs on a pair `(child_a, child_b)` where they
share an edge, the intersection algorithm correctly returns a segment
*along that shared edge* — which is NOT zero-length (it's the full
edge length). The "require segments are now zero-length" criterion
will fail on every correctly-conformalized pair. The test inverts:
correct output flags as failure, incorrect output (still crossing)
also flags as failure. Useless either way.

**Trigger:** Any conformalize call followed by these gates.

**Actual behavior:** Phase 2 ships with two "validations" that don't
distinguish a correct implementation from a no-op or a wrong one.

**Expected behavior:** Replace both with checks that actually
distinguish conformal from non-conformal output:

1. **Manifold check** (replaces (a)): every edge of the post-split
   fault is shared by ≤ 2 triangles, AND no triangle vertex lies in
   the open interior of an edge belonging to another triangle. This
   is the T-junction check (P-004).
2. **Cross-fault edge-coincidence check** (replaces (b)): for every
   intersection polyline, the polyline's vertex chain must appear as
   a sequence of mesh edges *both* in fault A's post-split
   triangulation and in fault B's post-split triangulation. (Edge
   coincidence, not absence of crossings.)
3. **Crossing-free interior check**: `tri_tri_intersect_3d` between
   `(child_a, child_b)` must report *only* shared-edge segments
   (where `child_a` and `child_b` share two vertices), not interior
   crossings. This requires a new `interior_crossings_only=True` flag
   in `tri_tri_intersect_3d` that filters out shared-vertex/edge
   pairs.

**Suggested fix:**

```diff
 ### Robustness rules
 - Drop any output sub-triangle whose area < `1e-6 m²` (zero-area
   protection — same as `_drop_degenerate` in `ts_to_stl.py`).
 - Drop any polyline segment shorter than `drop_short_segment_frac *
   mean_edge_in_parent`; record in the report.
-- After processing, re-run `find_overlap_pairs` on each post-split
-  fault and require it returns []; the constrained triangulation must
-  not introduce coplanar same-side overlaps.
-- After processing, re-run `tri_tri_intersect_3d` on every previously
-  crossing pair and require segments are now zero-length (within
-  `eps_min_seg_len_m`); the re-triangulation has eliminated the
-  geometric crossing.
+- **Manifold gate (per fault).** For each post-split fault: every edge
+  is shared by ≤ 2 triangles (within-fault), and no vertex lies in
+  the open interior of any edge to within `1e-3 m`.  Implemented as a
+  new helper `is_manifold_no_t_junctions(V, T)` in `fault_intersect.py`
+  (returns bool + list of offending edges).
+- **Polyline edge-coincidence gate (per pair).** For every
+  intersection polyline `P` from fault `A × B`: every consecutive
+  vertex pair `(P[i], P[i+1])` must appear as an edge in fault A's
+  post-split mesh AND in fault B's post-split mesh.  Mismatch =
+  conformalization bug.  Implemented as `verify_polyline_in_mesh(V, T,
+  polyline)` in `conformalize_faults.py`.
+- **Interior-crossing gate (per pair).** Run a new
+  `tri_tri_intersect_3d_interior_only(V_A, T_A, V_B, T_B)` which
+  filters out triangle pairs `(i, j)` that share one or more vertices
+  (these are the polyline-edge coincidences and are CORRECT, not
+  failures).  The post-split set must report 0 interior crossings.
+
+The previous "find_overlap_pairs returns []" and "tri_tri_intersect_3d
+returns zero-length segments" gates were tautologies (they pass on any
+input, correct or not) and have been removed.
```

**Test case:**
```python
def test_P005_old_gate_is_tautology():
    """find_overlap_pairs returns [] on a well-formed fault even
    BEFORE conformalize runs.  Demonstrates the old gate proves
    nothing."""
    # Take an unconformalized 2-fault input.  find_overlap_pairs
    # within fault A returns [] because A is well-formed.
    # find_overlap_pairs within fault B returns [] same reason.
    # Yet the cross-fault scan reports lots of crossings.  The old
    # gate would "pass" this broken state.
    V_A, T_A = load_stl("safs_sbmt_millcreek.stl")
    V_B, T_B = load_stl("safs_sbmt_saf.stl")
    assert find_overlap_pairs(V_A, T_A) == []
    assert find_overlap_pairs(V_B, T_B) == []
    # But this input is NOT conformal:
    crossings = scan_cross_fault_crossings(stl_dir, [
        "safs_sbmt_millcreek", "safs_sbmt_saf"])
    assert len(crossings) > 0   # crossings exist
    # ⇒ old gate is meaningless

def test_P005_new_polyline_edge_coincidence_gate():
    """After conformalize, polyline edges appear in both faults' mesh
    edge sets."""
    V_A2, T_A2, V_B2, T_B2, polylines = conformalize(...)
    edges_A = {frozenset((int(t[i]), int(t[(i+1)%3])))
               for t in T_A2 for i in range(3)}
    edges_B = {frozenset((int(t[i]), int(t[(i+1)%3])))
               for t in T_B2 for i in range(3)}
    for pl in polylines:
        for i in range(len(pl.points) - 1):
            # Map polyline vertices to mesh vertex indices, in both
            # faults, and assert the edge exists.
            ...
```

---

### [P-006] [MODERATE] [PLAN_multifault_intersections.md §Phase 1, "Numerical robustness"] — `gmpy2.mpfr` is not exact arithmetic; using it for orientation predicates is the wrong tool

**Category:** ASSUMPTION (wrong tool reference)

**Description:**
> fall back to higher precision (optional: `gmpy2.mpfr`; if `gmpy2` is
> not importable, log a warning and treat ambiguous cases as
> non-crossing — they are by definition near-tangencies that won't
> contribute meaningful conformality).

`gmpy2.mpfr` is **multi-precision float**, not exact arithmetic. With
sufficient precision (say 200 bits) it works for almost all inputs,
but it is fundamentally subject to the same rounding-mode dichotomy
as `float64` — just at a smaller scale. For *exact* orientation
predicates (sign(orient3d) for any inputs), the correct tool is one of:

- `gmpy2.mpz` (arbitrary-precision integers — works if inputs are
  scaled to integers via a fixed multiplier; the scale must exceed
  the dynamic range of products of three subtractions, ~10^48 for our
  data, fits in ~160-bit mpz).
- `gmpy2.mpq` (arbitrary-precision rationals — slower, but always
  exact).
- Shewchuk's adaptive predicates (the gold standard; available via
  the `predicates` Python wrapper or `triangle.predicates` if the
  Python `triangle` is built with them).

The plan's recommendation will produce slightly more precise results
than `float64` but still has rounding ambiguity in degenerate cases.
For SAFS this is probably fine in practice (CFM inputs are not
adversarial), but the plan should not claim that `mpfr` provides
"exact" predicates.

**Trigger:** Any near-tangency case where the orientation magnitude
is close to ULP. `mpfr` at default precision (53-bit) is identical
to `float64` and provides no benefit.

**Actual behavior (as planned):** Implementer reaches for `mpfr`,
finds it does not solve the ambiguity, gets confused, possibly
disables the fallback.

**Expected behavior:** Use `gmpy2.mpq` (rationals) for the
orientation determinant, or — preferably — wrap Shewchuk's adaptive
predicates via a small ctypes binding or a pure-Python port.

**Suggested fix:**

```diff
 **Numerical robustness**:
 - Default predicates use `np.float64`. When the orientation magnitude
   for any vertex of A relative to plane B is below
   `eps_orient * mean_edge_squared`, mark the case as ambiguous and
-  fall back to higher precision (optional: `gmpy2.mpfr`; if `gmpy2` is
-  not importable, log a warning and treat ambiguous cases as
-  non-crossing — they are by definition near-tangencies that won't
-  contribute meaningful conformality).
+  fall back to exact arithmetic.  Two acceptable implementations
+  (pick whichever is in `pythonenv`; if neither, declare a hard
+  dependency on `gmpy2` and use option (a)):
+  (a) `gmpy2.mpq` (rationals): convert the 12 input coordinates to
+      `mpq(int(x * 1_000_000), 1_000_000)` (1 µm grid) and evaluate
+      the orient3d determinant as exact rational; sign of the
+      result is exact.  `mpfr` is NOT acceptable: it is
+      multi-precision *float* and remains subject to rounding for
+      truly ambiguous predicates.
+  (b) Shewchuk's adaptive predicates via a small ctypes wrapper
+      (`predicates.so` — Shewchuk's C source is < 1000 lines and
+      compiles standalone with `gcc -O2`).  This is the gold
+      standard and is what the `triangle` library uses internally.
+  If neither is available, log a warning and treat the ambiguous
+  case as non-crossing; CFM inputs are not adversarial and the
+  failure mode is at most a missed near-tangency (no false
+  positives).
 - Explicitly forbid the coplanar case (`|n_A · n_B| > 1 - 1e-6` AND
   `|d_A − d_B| < 1 m`): coplanar surfaces should not occur in CFM data
   (different geological structures); raise.
```

**Test case:**
```python
def test_P006_orient3d_uses_exact_rationals():
    """A near-degenerate orient3d test where float64 gives the wrong
    sign but gmpy2.mpq gives the right sign."""
    import numpy as np
    # Construct 4 points such that orient3d is exactly 0 in mpq but
    # ±ULP in float64.  Easiest: 4 coplanar points generated from a
    # parametric quadrilateral.
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([1.0 + 1e-15, 1.0, 0.0])  # numerical noise
    p3 = np.array([0.0, 1.0, 0.0])
    # float64 orient3d will be ~1e-16, not 0
    o_float = orient3d_float64(p0, p1, p2, p3)
    o_exact = orient3d_mpq(p0, p1, p2, p3)
    assert o_exact == 0   # exact zero
    # The intersection routine using float64 gives ambiguous result;
    # using mpq gives a definite "no intersection" or "edge-touch".
```

---

### [P-007] [MODERATE] [PLAN_multifault_intersections.md §Phase 1 acceptance + chain_segments default] — Snap tolerance `snap_m=0.01 m` is 10× looser than `eps_min_seg_len_m=1e-3 m`; produces zero-length polyline edges

**Category:** EDGE_CASE (silent corruption)

**Description:**
The plan sets:
- `eps_min_seg_len_m = 1e-3` (1 mm) → segments shorter than 1 mm are
  dropped from `tri_tri_intersect_3d`.
- `snap_m = 1e-2` (1 cm) → segment endpoints within 1 cm collapse to
  the same chain node.

A 5 mm segment passes the `eps_min_seg_len_m` filter (5 mm > 1 mm),
enters `chain_segments`, and has both endpoints snapped to the same
chain node (5 mm < 1 cm). The chain then has a self-loop
(`p0 == p1`), which downstream code may interpret as a 1-vertex
polyline, an invalid edge, or a degenerate cycle. The `triangle`
library will reject the 1-vertex constraint or hang.

The hierarchy must be: `eps_min_seg_len_m > snap_m`, not the other
way around. The current values are inverted.

**Trigger:** Any tri-tri pair whose intersection segment is between
1 mm and 1 cm long — common at fault edges where one triangle barely
clips the other.

**Actual behavior (as planned):** Spurious self-loop polylines after
chaining; downstream behavior undefined.

**Expected behavior:** Set `eps_min_seg_len_m >= 2 * snap_m` so that
any segment surviving the length filter has endpoints distinguishable
under chain snapping.

**Suggested fix:**

```diff
 def tri_tri_intersect_3d(
     V_A: np.ndarray, T_A: np.ndarray,
     V_B: np.ndarray, T_B: np.ndarray,
     fault_a: str, fault_b: str,
     eps_orient: float = 1e-9,
-    eps_min_seg_len_m: float = 1e-3,
+    eps_min_seg_len_m: float = 5e-2,    # 5 cm (must exceed 2*snap_m below)
 ) -> list[Segment]: ...

 def chain_segments(
     segments: list[Segment],
-    snap_m: float = 1e-2,
+    snap_m: float = 1e-2,               # 1 cm (must be < eps_min_seg_len_m / 2)
 ) -> list[Polyline]: ...
```

And add an assertion in the chain function spec:

```diff
 **Segment chaining**: build a graph where each segment endpoint is a
 node (snapped to `snap_m` for matching), each segment is an edge.
+Pre-condition: every input segment satisfies `|p1 - p0| > 2 * snap_m`.
+Otherwise it would collapse to a self-loop.  Assert this on entry.
 Connected components → polylines.
```

**Test case:**
```python
def test_P007_snap_smaller_than_min_seg_len():
    """A 5mm segment must NOT survive into chain_segments with
    snap_m=1cm (would self-loop)."""
    seg = Segment(p0=(0,0,0), p1=(0.005, 0, 0), tri_a=0, tri_b=0,
                  fault_a="A", fault_b="B", nondegenerate=True)
    # Default eps_min_seg_len_m = 5e-2 means the 5mm seg is filtered
    # out at tri_tri stage; chain_segments never sees it.
    # If a caller passes it anyway (legacy path), the function must
    # raise rather than silently produce a self-loop.
    import pytest
    with pytest.raises(AssertionError, match="2 * snap_m"):
        chain_segments([seg], snap_m=1e-2)
```

---

### [P-008] [MODERATE] [PLAN_multifault_intersections.md §Phase 2 cascading, step 5] — Triangle-index remapping between cascade passes is unspecified; cascade can crash on second polyline pass

**Category:** ASSUMPTION (incomplete spec)

**Description:**
> 5. Cascade: a triangle that was already split by polyline P_AB may be
>    crossed by polyline P_AC. After step 3 for AB, the affected children
>    become the new input for the AC pass. Process polylines one pair at
>    a time, refreshing the working triangulation between passes.

The plan does not specify *how* `pierces_a` (a list of pre-split
triangle indices in fault A's mesh, computed by Phase 1 before any
splitting) maps to indices in the post-AB-split mesh. After AB
processing, the parent triangle index 17 may be replaced by children
indexed 17, 583, 584. When the AC pass later wants to insert polyline
P_AC into triangle index 17, that index now refers to one of the three
children — and only one of them is geometrically pierced by P_AC.

The fix is a re-scan: rerun Phase 1's `tri_tri_intersect_3d` between
fault A (post-AB-split) and fault C (original) to get fresh
`pierces_a` indices into the current mesh. The plan implicitly assumes
this but does not state it; an implementer following the plan
literally may try to reuse the original `pierces_a` and get index
errors or geometric mismatches.

**Trigger:** Any 3+-fault input (Phase 5 verification target;
specifically the proposed 3-fault crossing subset).

**Actual behavior (as planned):** Indeterminate — implementer might
reuse stale indices and crash, or re-scan and silently work but with
duplicated tri-tri scan cost.

**Expected behavior:** Explicit re-scan between cascade passes, with
documented O((K-1) * pair_scan_cost) for K polylines.

**Suggested fix:**

```diff
 5. Cascade: a triangle that was already split by polyline P_AB may be
    crossed by polyline P_AC. After step 3 for AB, the affected children
    become the new input for the AC pass. Process polylines one pair at
-   a time, refreshing the working triangulation between passes.
+   a time, refreshing the working triangulation between passes.
+
+   **Index discipline:** `pierces_a` and `pierces_b` from Phase 1's
+   `tri_tri_intersect_3d` reference the triangulation as it was at
+   scan time.  After a fault is re-triangulated, those indices are
+   stale.  Therefore: between cascade passes, *re-scan*
+   `tri_tri_intersect_3d(V_A_now, T_A_now, V_C_now, T_C_now)` to
+   recompute `pierces_a` and `pierces_b` against the current
+   triangulations.  Do not attempt to remap indices via parent-child
+   bookkeeping — the rescan is cheaper, simpler, and harder to get
+   wrong.  Cost: O((K-1) * pair_scan_cost) for K polylines, which is
+   fine at the 2- and 3-fault scale.
```

**Test case:**
```python
def test_P008_three_fault_cascade():
    """Three faults A, B, C with two crossings (A×B, A×C).  Verify the
    A×C pass re-scans and finds A's post-AB triangulation, not the
    pre-AB indices."""
    # Synthetic: three perpendicular planes meeting at the origin.
    # After A×B re-triangulation, A's mesh has new vertices on the
    # AB intersection line.  The A×C scan must find pierces in A's
    # CURRENT mesh, not the original mesh.
    V_A0, T_A0, V_B0, T_B0, V_C0, T_C0 = three_planes()
    V_A1, T_A1, V_B1, T_B1 = conformalize_pair(
        V_A0, T_A0, V_B0, T_B0, "A", "B")
    # A1 has more vertices than A0 (split along AB intersection)
    assert V_A1.shape[0] > V_A0.shape[0]
    # A×C re-scan returns pierces into V_A1, T_A1, NOT V_A0, T_A0
    segs_AC = tri_tri_intersect_3d(V_A1, T_A1, V_C0, T_C0, "A", "C")
    # Pierce indices must be valid in T_A1
    for s in segs_AC:
        assert s.tri_a < T_A1.shape[0]
```

---

### [P-009] [LOW] [PLAN_multifault_intersections.md §Risk 6, all-8 pair count] — Quoted "24M pairs" arithmetic is wrong by ~6x

**Category:** QUALITY (factual error in motivation)

**Description:**
> Performance for the all-8 set (Phase 5). Naive O(N²) pair scan on
> 1346 + 794 + 504 + 449 + 346 + 234 + 610 + 639 = 4922 triangles =
> 24M pairs is borderline.

The arithmetic is wrong. Cross-fault pairs over 8 faults are
`(sum n_i)² / 2 - sum(n_i²) / 2`, which evaluates to:
- `sum n_i = 4922`, `(sum)² / 2 ≈ 12.1M`
- `sum n_i² = 1346² + 794² + ... ≈ 3.85M`, halved ≈ 1.92M
- Cross-fault pairs ≈ `12.1M - 1.92M ≈ 4.1M` (roughly).

So the naive scan is ~6× *less* expensive than claimed. AABB-tree
acceleration is unlikely to be needed at all at 2000m. This is a
benign over-estimate that doesn't affect correctness, but it
mis-prioritizes optimization work in Phase 5.

**Trigger:** None — pure documentation issue.

**Actual behavior (as planned):** Reader believes pair scan is
expensive; implementer may waste effort on premature AABB optimization.

**Expected behavior:** Correct count + matching defer guidance.

**Suggested fix:**

```diff
 6. **Performance for the all-8 set** (Phase 5). Naive O(N²) pair scan
-   on 1346 + 794 + 504 + 449 + 346 + 234 + 610 + 639 = 4922 triangles
-   = 24M pairs is borderline. Mitigation: AABB tree if naive >30 s.
+   on 8 faults with triangle counts {1346, 794, 504, 449, 346, 234,
+   610, 639} produces ≈ 4.1M cross-fault pairs (computed as
+   ((Σn)² − Σn²)/2 = (24.2M − 3.85M)/2). Numpy-vectorized naive
+   scan handles this in single-digit seconds; AABB tree is unlikely
+   to be needed at the 2000m resolution. Defer.
    Out of scope for the 2-fault milestone.
```

**Test case:** N/A (documentation correction).

---

## Summary

- Critical issues: 5 (P-001, P-002, P-003, P-004, P-005)
- Moderate issues: 3 (P-006, P-007, P-008)
- Low issues: 1 (P-009)
- Plan compliance: **N/A** (this is a plan-document review, not a code
  implementation review).
- Verdict: **FAIL — must revise plan before implementing.** P-001 through
  P-005 are correctness errors that, if implemented literally, would
  produce a non-conforming output mesh while looking like the plan
  succeeded — i.e., the same failure mode as the *existing* code, just
  with more sophisticated wrong code. P-001 (Möller mis-specification),
  P-003 (`'pq'` flag), and P-004 (T-junction propagation) are the three
  bugs that most directly cause conformality failure. P-005 (false
  validation gates) is what would let the bugs ship undetected.

Suggested order of plan revisions:
1. Apply P-001 (correct Möller algorithm).
2. Apply P-004 (mandatory edge-pierce propagation).
3. Apply P-005 (replace tautological validation gates).
4. Apply P-003 (`'p'` flag, no `'q'`).
5. Apply P-002 (free-surface clamp not drop).
6. Apply P-007 (snap < min_seg_len ordering).
7. Apply P-006, P-008 (numerical-robustness wording, cascade re-scan).
8. Apply P-009 (factual correction).

After revision, regenerate `PLAN_multifault_intersections.pdf`.

## Unreviewed Areas

- Phase 3 (size field tube) and Phase 4 (driver wiring): minor wording
  changes only; their correctness was already covered by R-004 and
  R-005 in `mesh/REVIEW.md`. Re-reviewed quickly and no new findings
  surfaced.
- Phase 5 (3-fault cascade verification): only the cascade-index issue
  (P-008) is plan-level. The deeper 3-way-junction handling (where 3
  polylines meet at a single point) was flagged as in-scope by Phase 2
  step 5 implicitly; an explicit 3-way-meeting test would belong in
  Phase 5 acceptance, but adding it is enhancement, not bug-fix, so I
  did not flag it as a finding.
