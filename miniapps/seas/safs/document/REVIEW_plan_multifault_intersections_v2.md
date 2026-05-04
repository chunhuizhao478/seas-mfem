# Code Review of `PLAN_multifault_intersections.md` — Round 2 (2026-04-29)

## Review Scope

- **Plan reviewed (current state):**
  `miniapps/seas/safs/PLAN_multifault_intersections.md` (after the
  Round-1 fixes from `REVIEW_plan_multifault_intersections.md` and
  the user-driven readability rewrite of Phase 1's algorithm and
  Phase 2's algorithm).
- **Companion reviews (informational, not duplicated):**
  `miniapps/seas/safs/REVIEW_plan_multifault_intersections.md` (Round 1
  — applied), `miniapps/seas/safs/mesh/REVIEW.md` (existing-code review,
  separate scope).
- **Domain context**: `seas-mfem/CLAUDE.md`,
  `miniapps/seas/CLAUDE.md`, `miniapps/seas/safs/PLAN.md`,
  `miniapps/seas/safs/mesh/safs.geo`,
  `miniapps/seas/safs/mesh/audit_ts_quality.py:find_overlap_pairs`.
- **Audit method:** re-executed all three review passes from scratch.
  This is **not** a "verify the previous fixes" review — Round 1's
  edits introduced a worked example, ASCII pictures, and a complete
  Phase 2 rewrite, all of which need fresh adversarial reading.

---

## Findings

### [P2-001] [CRITICAL] [PLAN_multifault_intersections.md §Testing strategy table] — Test row for Phase 2 still references the OLD tautological gates that Round 1 replaced

**Category:** DEVIATION

**Description:**
The "Testing strategy" table (line ~826) Row 2 reads:

> | 2 | 90°-cross squares + cascade with a 3rd plane | Mill × SBMT-SAF
> post-split STLs are crossing-free | `find_overlap_pairs` empty;
> cross-scan empty |

The "Validator gate" column for Phase 2 is `` `find_overlap_pairs`
empty; cross-scan empty ``. But Round 1's P-005 replaced exactly those
two gates with the manifold + polyline-edge-coincidence +
interior-crossing gates because the originals were tautologies.

The §Validation gates subsection now correctly says:
> The previous "find_overlap_pairs returns []" and "tri_tri_intersect_3d
> returns zero-length segments" gates were tautologies …

But the summary table in §Testing strategy still tells the implementer
to use them. A fix-agent reading this table mechanically will
implement the wrong gates and the conformality check will be
meaningless — the same failure mode that motivated P-005 in the
first place.

**Trigger:** Anyone reading the Testing strategy summary table to
figure out what to test for.

**Actual behavior (as written):** Plan contradicts itself —
prescribes the new gates in one section and the old tautologies in
another.

**Expected behavior:** Table updated to reference the new gates.

**Suggested fix:**

```diff
 | Phase | Synthetic test | Real fixture | Validator gate |
 |---|---|---|---|
 | 1 | 90°-cross squares; 100 m parallel; near-tangent | Mill × SBMT-SAF crossing-pair count vs README | scan returns expected pairs |
-| 2 | 90°-cross squares + cascade with a 3rd plane | Mill × SBMT-SAF post-split STLs are crossing-free | `find_overlap_pairs` empty; cross-scan empty |
+| 2 | 90°-cross squares + cascade with a 3rd plane | Mill × SBMT-SAF post-split STLs pass all three Phase-2 gates | (a) manifold gate per fault, (b) polyline-edge-coincidence gate per pair, (c) interior-crossing-only gate per pair |
 | 3 | None (size-field config only) | Single-fault smoke retains 10/10; 2-fault tube uniformity | check_11 passes |
```

**Test case:**
```python
def test_P2_001_testing_strategy_matches_validation_gates():
    """The testing-strategy summary table must list the gates that
    the §Validation gates subsection actually defines, not the
    deprecated ones."""
    plan = open("PLAN_multifault_intersections.md").read()
    # The deprecated gate names must NOT appear inside the Phase 2
    # row of the testing strategy table.
    table_pos  = plan.find("## Testing strategy")
    risk_pos   = plan.find("## Risk assessment")
    table_text = plan[table_pos:risk_pos]
    phase2_row = [r for r in table_text.splitlines()
                  if r.startswith("| 2 |")][0]
    assert "find_overlap_pairs" not in phase2_row, (
        "Phase 2 test row still references the deprecated tautological "
        "gate; should mention the manifold/edge-coincidence/interior-"
        "crossing gates instead.")
    assert "manifold" in phase2_row.lower() or \
           "edge-coincidence" in phase2_row.lower() or \
           "interior-crossing" in phase2_row.lower()
```

---

### [P2-002] [CRITICAL] [§Phase 2 step (3) cascade rule vs §Phase 2 Edge cases item "A triangle is crossed by 2+ polylines"] — The two prescriptions for handling multiple polylines per parent contradict each other

**Category:** BUG (incompatible specs)

**Description:**
Step (3) "Cascade across 3+ faults" (line ~587) prescribes:

> Process polylines one pair at a time, refreshing the working
> triangulation between passes. … Rule: between cascade passes,
> re-scan.

Edge case "A triangle is crossed by 2+ polylines" (line ~705) prescribes:

> re-triangulate with all relevant polyline segments as constraints
> in the same Delaunay call. The `triangle` library handles multiple
> constraint segments natively.

These are **two different mechanisms for the same case**: one parent
triangle pierced by polylines from two distinct fault pairs (e.g.,
`P_AB` and `P_AC`). Step (3) says do `P_AB` first, then re-scan
to find which (post-split) child of the parent is now pierced by
`P_AC`, then process. The Edge-cases item says do both polylines'
constraint segments in a single CDT call on the original parent.

A fix-agent will get confused about which to implement. Either
mechanism alone is correct, but they cannot coexist:

- "All at once" requires knowing all polylines that touch a parent
  *before* re-triangulating. That breaks if cascade pass 2 introduces
  a new pierce that wasn't visible in pass 1.
- "One at a time + re-scan" requires re-triangulating after each
  polyline. That works but the Edge-cases item explicitly says "the
  same Delaunay call" — directly contradicting it.

For 3-fault cascading both polylines through one parent are visible
*after* the second tri-tri scan (Step 3's re-scan), which means the
"all at once" approach is fine *if and only if* you re-scan first
and gather all polylines per current parent before each CDT call.
The plan does not say this.

**Trigger:** Any 3-fault input where one parent triangle of fault A
is pierced by both `P_AB` and `P_AC`. Likely in the Phase 5 3-fault
target (Mill Creek + SBMT-SAF + Mission Creek SBMT, San Gorgonio
knot region).

**Actual behavior (as written):** Implementer picks one of two
mechanisms; the other branch's wording becomes wrong; plan is
internally inconsistent.

**Expected behavior:** A single, unambiguous prescription:

> Process polyline-pairs one pair at a time. Between pairs, re-scan
> with `tri_tri_intersect_3d` to refresh pierce indices against the
> current triangulation. **Within a pair**, if a single parent
> triangle is pierced by multiple segments of the same polyline (or
> by two polylines whose pierces co-locate on the same parent), pass
> all the constraint segments to one CDT call.

**Suggested fix:**

```diff
 - **A triangle is crossed by 2+ polylines** (cascading): re-triangulate
-  with all relevant polyline segments as constraints in the same
-  Delaunay call. The `triangle` library handles multiple constraint
-  segments natively.
+  with all relevant polyline segments as constraints in the same
+  Delaunay call **within a single cascade pass**.  The `triangle`
+  library handles multiple constraint segments natively.  This is
+  the case where, after Step (3)'s re-scan, the *current* re-scan
+  shows a parent triangle pierced by multiple segments of the same
+  polyline (e.g. a polyline with a kink that re-enters the parent),
+  not the case of two polylines from different pairs landing on the
+  same parent — those two polylines belong to different cascade
+  passes by Step (3), so they are processed sequentially with a
+  re-scan in between.
```

And cross-reference Step (3):

```diff
 Rule: **between cascade passes, re-scan**. Don't try to remap indices
 via parent-child bookkeeping; the rescan is `O(pair-scan-cost)` per
 pass and is the simpler, safer choice.
+
+Within a single cascade pass, *all* constraint segments touching one
+parent triangle (multiple segments from the same polyline, plus any
+already-recorded pierces from preceding within-pass propagation) are
+collected and passed in one `triangle.triangulate({segments: …},
+'p')` call.  Mixing the "one CDT call" rule (within pass) with the
+"re-scan between passes" rule is correct; only the mixing must be
+made explicit.
```

**Test case:**
```python
def test_P2_002_three_fault_cascade_same_parent():
    """Three faults A, B, C; A's triangle 5 is pierced by both P_AB
    and P_AC.  Cascade processes P_AB first; after re-scan, A's
    triangle 5 is gone (split into children).  P_AC's pierces must
    land on those children, NOT triangle 5."""
    V_A, T_A, V_B, T_B, V_C, T_C = three_planes_meeting_at_origin()
    # P_AB pierces A's triangle 5 along a known segment.
    # P_AC also pierces A's triangle 5 along a different known segment.
    out_A1, out_B = conformalize_pair(V_A, T_A, V_B, T_B, "A", "B")
    # After AB, original triangle 5 of A no longer exists.
    assert len(out_A1.T) > len(T_A)
    # P_AC re-scan must find pierces in out_A1.T, not the original T_A.
    segs_AC = tri_tri_intersect_3d(out_A1.V, out_A1.T,
                                    V_C, T_C, "A", "C")
    # Pierces on A side reference the post-AB indices.
    for s in segs_AC:
        assert s.tri_a < len(out_A1.T)
    # And the conformalize_pair("A", "C") with out_A1 as input must
    # not crash on stale indices.
    out_A2, out_C = conformalize_pair(out_A1.V, out_A1.T,
                                       V_C, T_C, "A", "C")
    assert len(out_A2.T) >= len(out_A1.T)
```

---

### [P2-003] [CRITICAL] [§Phase 1 step (6) "Recover 3-D endpoints by mapping back through the interpolation"] — The 3-D endpoint recovery formula is not stated; "linearly blend between consecutive p_hits" is ambiguous

**Category:** ASSUMPTION (incomplete spec)

**Description:**
Step (6) is supposed to map the overlap interval `[t_start, t_end]`
back to 3-D points. The current wording:

> Recover 3-D endpoints by mapping back through the interpolation in
> step (4): keep, alongside each `t`, the `p_hit` it came from, and
> **linearly blend** between consecutive p_hits on each triangle to
> recover the 3-D endpoint at any new `t` in `[tA_lo, tA_hi]`. (Or
> equivalently: `p = (project A onto L)` evaluated at `t = t_start`,
> which is just `p_hit` of whichever boundary defined `t_start`.)

This is muddled in three ways:

1. "Linearly blend between consecutive p_hits" — there are exactly
   two p_hits per triangle (the two interpolated edge crossings).
   "Consecutive" implies an ordering; "blend" doesn't say what
   parameter governs the blend.
2. "p_hit of whichever boundary defined `t_start`" — this is fine
   when `t_start = tA_lo` or `t_start = tA_hi` (boundary of A's
   interval), but if `t_start = tB_lo` and `tB_lo` is strictly
   inside `[tA_lo, tA_hi]`, the boundary that defined `t_start` is
   on triangle B, not A. We need a 3-D point on `L` at `t = tB_lo`
   that's interpolated *along A's chord*.
3. The linear-interpolation formula in 3-D is straightforward:
   ```
   α = (t - tA_lo) / (tA_hi - tA_lo)
   p_at_t = (1 - α) · p_hit_A_lo + α · p_hit_A_hi
   ```
   But the plan never states this. An implementer who reads only
   the second clause ("p = p_hit of whichever boundary defined
   t_start") will copy the wrong endpoint when the overlap is
   strictly inside both intervals.

**Trigger:** Whenever the overlap interval is strictly inside one
triangle's chord, e.g.:
- A's interval = [-5, +5]
- B's interval = [-0.5, +0.5]
- Overlap = [-0.5, +0.5]
- `t_start = tB_lo = -0.5`, defined by B, not A.

The worked example *is* exactly this case (A goes [-5, 5]; B goes
[-0.5, 0.5]; overlap = [-0.5, 0.5] is fully inside A's interval).
Under the plan's "p_hit of whichever boundary defined t_start"
clause, the implementer should pick B's `p_hit_B_lo = (0, 0, -0.5)`
— which happens to also be on A's chord because the answer is on
both. But this is geometric coincidence, not a derivation from the
formula. For non-axis-aligned cases, picking the "boundary that
defined" `t_start` from B and using it directly is fine *only because
both triangles' chords are subsets of `L` and meet at the same
3-D point at the same `t`*. The plan should make this principle
explicit, not leave it to "linearly blend".

**Actual behavior (as written):** Implementer guesses what "blend"
means, possibly returning the wrong p_hit when the overlap doesn't
align with one of A's two computed crossings.

**Expected behavior:** Explicit formula:

```
For overlap endpoint at parameter t:
    α = (t - tA_lo) / (tA_hi - tA_lo)         in [0, 1]
    p = (1 - α) · p_hit_A_lo + α · p_hit_A_hi  3-D position on L

Either A's chord or B's chord can be used; both pass through L at
the same 3-D point at any given t (this is the definition of L).
Pick the one whose p_hits are already in hand to avoid an extra
interpolation.
```

**Suggested fix:**

```diff
-Recover 3-D endpoints by mapping back through the interpolation in
-step (4): keep, alongside each `t`, the `p_hit` it came from, and
-**linearly blend** between consecutive p_hits on each triangle to
-recover the 3-D endpoint at any new `t` in `[tA_lo, tA_hi]`. (Or
-equivalently: `p = (project A onto L)` evaluated at `t = t_start`,
-which is just `p_hit` of whichever boundary defined `t_start`.)
+Recover the 3-D position at any parameter `t` along `L` by linearly
+interpolating along **A's chord** (or B's — both give the same
+answer because both chords are subsets of `L`):
+
+```
+α = (t - tA_lo) / (tA_hi - tA_lo)            # in [0, 1]
+p_at_t = (1 - α) · p_hit_A_lo + α · p_hit_A_hi
+```
+
+Apply the formula at `t = t_start` and `t = t_end` to get the two
+3-D endpoints of the overlap segment.  When `t_start` happens to
+equal one of the four boundary t-values (`tA_lo`, `tA_hi`, `tB_lo`,
+`tB_hi`), the corresponding `p_hit` is the answer directly; the
+formula reduces to that p_hit (with α = 0 or α = 1), so a fast
+path that returns the matching p_hit is correct but optional.
```

**Test case:**
```python
def test_P2_003_overlap_interior_to_one_triangle():
    """A's chord on L spans t = [-5, +5].  B's chord spans
    [-0.5, +0.5].  Overlap is fully inside A's interval, so both
    overlap endpoints come from B's pierce points — but the implementer
    must compute them via the formula, not via 'pick A's p_hit'."""
    import numpy as np
    # Triangles from the plan's worked example.
    V_A = np.array([[-5,0,-5], [5,0,-5], [0,0,5]])
    T_A = np.array([[0,1,2]])
    V_B = np.array([[0,-0.5,-0.5], [0,0.5,-0.5], [0,0,0.5]])
    T_B = np.array([[0,1,2]])
    segs = tri_tri_intersect_3d(V_A, T_A, V_B, T_B, "A", "B")
    assert len(segs) == 1
    p0, p1 = np.array(segs[0].p0), np.array(segs[0].p1)
    expected = np.array([(0, 0, +0.5), (0, 0, -0.5)])
    pts = np.array([p0, p1])
    pts_sorted = pts[np.argsort(pts[:, 2])][::-1]
    np.testing.assert_allclose(pts_sorted, expected, atol=1e-9)
    # If implementer returned A's p_hit_A_lo or p_hit_A_hi instead,
    # they'd get (0, 0, +5) or (0, 0, -5) — fail.
```

---

### [P2-004] [MODERATE] [§Phase 2 step (1) "before/after" picture] — ASCII picture is geometrically incoherent; vertex labels (q1, q2, q1') don't match the listed children

**Category:** QUALITY (broken documentation; will mislead the implementer)

**Description:**
The picture (line ~528) reads:

```
        before                          after step (1)
                                         q1, q2 = polyline pierce points
            v2                              v2
           /  \                            /| \
          /    \                          / |  \
         /      \                        /  q2  \         ← 4 children:
        /   ●━━━━━━━━━━ polyline      ─/───●    \         ─────────────
       /        \                      /    \    \         (v0, q1, v2)
      /          \                    / q1   \    \        (q1, q2, v2)
     /            \                  /  ●─────\    \       (q1, q1', q2)
    v0────────────v1                v0─────────────v1      … etc.
                                          (q1' = entry on edge v0v1)
```

Multiple problems:
1. The legend says **q1, q2 = polyline pierce points** — only two — but
   the children list mentions a third vertex `q1'` ("entry on edge
   v0v1"), which the implementer can't reconcile.
2. The "before" panel shows the polyline as `●━━━━━━━━━━` — a
   dot followed by a horizontal line. There is no second dot. So
   it's unclear what the polyline's geometry actually is.
3. The "after" panel shows `q1` *both* on edge `v0v1` (the `●─────`
   on line 540) *and* somewhere in the interior (the `q1` label on
   line 539). The same vertex appears in two places.
4. The children list ends with `… etc.` — the implementer doesn't
   know how many children to expect.
5. The "4 children" header conflicts with `… etc.` — either it's
   exactly 4, or it's "etc."

**Trigger:** An implementer following the picture as a worked
specification.

**Actual behavior:** Confusion; risk of building a wrong stitching.

**Expected behavior:** A simple, unambiguous picture for **the
canonical case**: polyline enters one edge, exits another, with **two
pierce points and three children**. State the geometry explicitly:

**Suggested fix:**

```diff
 Picture (one parent triangle, one polyline crossing it):
 
 ```
-        before                          after step (1)
-                                         q1, q2 = polyline pierce points
-            v2                              v2
-           /  \                            /| \
-          /    \                          / |  \
-         /      \                        /  q2  \         ← 4 children:
-        /   ●━━━━━━━━━━ polyline      ─/───●    \         ─────────────
-       /        \                      /    \    \         (v0, q1, v2)
-      /          \                    / q1   \    \        (q1, q2, v2)
-     /            \                  /  ●─────\    \       (q1, q1', q2)
-    v0────────────v1                v0─────────────v1      … etc.
-                                          (q1' = entry on edge v0v1)
+    BEFORE                            AFTER step (1)
+    ──────                            ──────────────
+
+         v2                                v2
+         /\                                /\
+        /  \                              /  \
+       /    \                            / q2 \      q1, q2 = the two
+      /      \                          /───●  \     pierce points
+     /   ●←━━━━━━━━━●                  /   /│   \    where the polyline
+    /  q1          q2 \               /   / │    \   enters & leaves
+   /                   \             /q1 /  │     \  the parent triangle
+  v0────────────────────v1          v0──●───┴──────v1
+                                        |
+                                  3 children:
+                                    (v0, q1, v2)
+                                    (q1, q2, v2)
+                                    (q1, v1, q2)
+                                  — q1 lies on edge v0-v1
+                                  — q2 lies on edge v1-v2
+                                  — the polyline edge q1-q2
+                                    is shared by children 2 and 3
 ```
+
+The number of children is `2 + (# pierce points strictly interior
+to the parent)`.  In the canonical "enters one edge, exits another"
+case shown here, both pierces are on edges, so 0 interior, hence
+exactly 3 children.  A polyline that enters through an edge and
+terminates inside the parent (e.g., at a fault boundary) gives
+1 interior pierce ⇒ 4 children.  Always: count is determined by
+the geometry, not a fixed number.
```

**Test case:**
```python
def test_P2_004_canonical_split_produces_three_children():
    """Parent triangle pierced by a polyline that enters edge v0v1 and
    exits edge v1v2.  Result: exactly 3 children, polyline edge q1-q2
    shared between children 2 and 3."""
    import numpy as np
    V = np.array([[0,0,0], [2,0,0], [1,1,0]])     # parent v0,v1,v2
    T = np.array([[0,1,2]])
    polyline = [(1.0, 0.0, 0.0), (1.5, 0.5, 0.0)]  # q1 on v0v1, q2 on v1v2
    V_out, T_out = conformalize_one_parent(V, T, polyline)
    assert T_out.shape[0] == 3, T_out.shape[0]
    # Polyline edge appears in two children.
    q1_idx = np.where(np.linalg.norm(V_out - polyline[0], axis=1) < 1e-9)[0][0]
    q2_idx = np.where(np.linalg.norm(V_out - polyline[1], axis=1) < 1e-9)[0][0]
    edge_count = sum(
        1 for tri in T_out
        if {q1_idx, q2_idx}.issubset({int(x) for x in tri}))
    assert edge_count == 2, f"polyline edge in {edge_count} children, expected 2"
```

---

### [P2-005] [MODERATE] [§Phase 2 step (2) "BAD vs GOOD" picture] — Picture has duplicate `p_new` symbols and conflicting geometry

**Category:** QUALITY (broken documentation)

**Description:**
The propagation picture (line ~555):

```
    BAD (no propagation)             GOOD (propagation)

         T_a (split)                      T_a (split)
         /|\                              /|\
        / | \                            / | \
       /  |  \                          /  |  \
      /   ●   \   ← p_new on edge      /   ●   \
     /    │    \                      /    │\   \
    u─────────v                      u─────●─────v
     \         /                      \    │\   /     ← p_new is now a
      \       /                        \   │ \ /        vertex of T_b too;
       \ T_b /                          \  │  ●         T_b is split into
        \   /                            \ │ /          (u, p_new, w_b) and
         \ /                              \│/           (p_new, v, w_b).
          w_b                              w_b
```

Problems:
1. The "GOOD" panel shows **three** `●` symbols but the legend
   identifies only one `p_new`. Reader has to guess which one is
   p_new.
2. The `●` on line 562 is described as "p_new on edge" (correct),
   but the additional `●` symbols on the diagonal lines down to
   `w_b` are unlabeled — they look like additional vertices, but
   they're meant to indicate the *split lines* (T_b → two children).
3. The vertical bars (`│` and `\`) are out of alignment, making it
   look like the diagonals don't actually go to the same vertex.
4. `w_b` is supposed to be the vertex of T_b opposite edge `(u, v)`.
   In the picture, `w_b` is at the *bottom* of T_b. That's
   geometrically reasonable, but the location of the new `●`
   symbols on the diagonal connecting to `w_b` suggests they sit on
   the diagonals, not the edge — contradicting "p_new is on the
   shared edge".

**Trigger:** Anyone trying to follow the picture.

**Actual behavior:** Reader cannot tell which point is p_new.

**Expected behavior:** Single `●` symbol labeled `p_new`, and the
split lines drawn without additional dots.

**Suggested fix:**

```diff
 ```
-    BAD (no propagation)             GOOD (propagation)
-
-         T_a (split)                      T_a (split)
-         /|\                              /|\
-        / | \                            / | \
-       /  |  \                          /  |  \
-      /   ●   \   ← p_new on edge      /   ●   \
-     /    │    \                      /    │\   \
-    u─────────v                      u─────●─────v
-     \         /                      \    │\   /     ← p_new is now a
-      \       /                        \   │ \ /        vertex of T_b too;
-       \ T_b /                          \  │  ●         T_b is split into
-        \   /                            \ │ /          (u, p_new, w_b) and
-         \ /                              \│/           (p_new, v, w_b).
-          w_b                              w_b
+    BAD  (no propagation)             GOOD  (propagation applied)
+
+        T_a (split)                       T_a (split)
+        /  \                              /  \
+       /    \                            / |  \
+      /      \                          /  |   \
+     /        \                        /   |    \
+    u─────●────v   p_new is here      u────●─────v   p_new is here
+     \   T_b   /     (only on T_a)     \  /│\    /    (on BOTH T_a and T_b)
+      \       /                         \/ │ \  /
+       \     /                          /  │  \/
+        \   /                          /   │   \
+         \ /                          /    │    \
+          w_b                              w_b
+
+    T_b is one undivided triangle    T_b is split into two children:
+    even though p_new lies on its    (u, p_new, w_b) and (p_new, v, w_b)
+    edge — that's the T-junction.    — the T-junction is gone.
 ```
```

**Test case:** Visual / documentation check; no code test (the
picture is prose-equivalent).

---

### [P2-006] [MODERATE] [§Phase 2 step (1) sub-step b. "near_vertex_snap_frac × mean_edge"] — Snap-to-corner can drop a polyline pierce point silently, splitting the polyline into disconnected pieces and breaking conformality

**Category:** EDGE_CASE / BUG

**Description:**
Sub-step b. snaps a pierce point `q_i` to a parent vertex `v_j` when
`|q_i - v_j| < near_vertex_snap_frac * mean_edge`. With the default
`near_vertex_snap_frac = 1e-2` and CFM 2000m triangles
(mean_edge ≈ 2 km), the snap radius is **20 m**.

A polyline whose neighbouring pierce points `q_{i-1}` and `q_{i+1}`
both land within 20 m of triangle vertex `v_j` will have `q_i` and
its neighbours all snapped to `v_j` — collapsing a polyline edge to
zero length.

Even worse: when the snap is applied independently on fault A's
parent and fault B's parent, the same polyline vertex `q_i` may snap
to vertex `v_j_A` on fault A but **not** snap on fault B (because
fault B's parent triangle is geometrically different and `v_j_A`
isn't a vertex of fault B). Now fault A's mesh has a `v_j_A` vertex
on the polyline; fault B's mesh has a different `q_i` vertex at a
different 3-D location. The polyline edge-coincidence gate (P-005)
will fail.

The `drop_short_segment_frac = 5e-2` rule (≈ 100 m at 2000m mean
edge) catches the *short-segment* case but does not coordinate with
sub-step b. Specifically, sub-step b. is a **per-fault, per-parent**
local decision that doesn't see the segment length on the OTHER
fault. The two faults can disagree on which pierces to snap.

**Trigger:** Cross-fault polyline that grazes a CFM vertex on one
fault but not the other (common because CFM 2000m triangles on
fault A and fault B have different vertex placements).

**Actual behavior (as planned):** Sub-step b. silently snaps,
producing per-fault triangulations that disagree on the polyline
geometry. Polyline edge-coincidence gate fails post-conformalize.

**Expected behavior:** Snap decisions for polyline pierces must be
**globally consistent across both faults**. Either:

(a) Don't snap pierces. Always insert as Steiner. Accept the small
slivers (sub-step a. project + sub-step c. CDT will keep them
non-degenerate as long as `eps_min_seg_len_m` is enforced upstream).

(b) Snap centrally before fault-by-fault re-triangulation: pre-
process the polyline, decide globally which pierces to snap (and to
which fault's vertex), and modify the polyline to use that vertex
on **both** faults' parent triangles. This requires the snapped-to
vertex to be a vertex of both — generally only possible at chain
endpoints lying on the free-surface trace or fault boundary.

**Suggested fix (recommend option (a) — simpler and matches the
"trust constrained Delaunay" philosophy already adopted for `'p'`):**

```diff
 | sub-step | what | why |
 |---|---|---|
 | **a. project** | parent's 3 vertices + P's pierce points → 2-D | constrained Delaunay is a 2-D operation |
-| **b. snap to corner** | if a pierce point is within `near_vertex_snap_frac × mean_edge` of a parent vertex, replace it with that vertex | avoids zero-area sliver children |
+| **b. (no per-parent snap)** | Insert every pierce as a Steiner vertex unconditionally.  Per-parent corner-snapping is dropped: the same polyline pierce can be snap-distance from a vertex of fault A's parent without being snap-distance from any vertex of fault B's parent, breaking cross-fault polyline-edge coincidence.  Slivers from short edges are caught upstream by `eps_min_seg_len_m` (Phase 1 step 7) — anything that survives is guaranteed long enough to constrain. |
 | **c. CDT** | call `triangle.triangulate({verts, segments}, 'p')` | flag `'p'` only — `'q'` would insert Steiner points on the constraint and re-break conformality |
```

If retaining per-parent snapping is preferred, then sub-step b.
must be replaced with a global pre-processing step that records the
post-snap polyline once and uses that on both faults; flag this as
follow-up work and add a per-pair sanity check.

**Test case:**
```python
def test_P2_006_pierce_near_vertex_on_one_fault_only():
    """Cross-fault polyline whose middle pierce point lies 5m from a
    vertex of fault A's parent (would be snapped under
    near_vertex_snap_frac = 0.01 with mean_edge = 2km), but >100m
    from any vertex of fault B's parent.  After conformalize, fault
    A's mesh has the snapped vertex on the polyline; fault B's mesh
    has the original pierce.  The polyline-edge-coincidence gate
    must catch this — and the recommended fix (no snap) makes it not
    happen in the first place."""
    V_A, T_A = make_2km_triangle_with_vertex_at(np.array([1000., 0., -50.]))
    V_B, T_B = make_2km_triangle_offset_by(50.)
    polyline = [(0, 0, -50), (1005, 0, -50), (2000, 0, -50)]
    out_A, out_B = conformalize_pair(V_A, T_A, V_B, T_B, "A", "B",
                                      near_vertex_snap_frac=0.01)
    # Polyline-edge-coincidence: middle vertex (1005, 0, -50) must
    # exist in BOTH meshes and the consecutive polyline edges must
    # appear as mesh edges in both.
    p_mid = np.array([1005., 0., -50.])
    in_A = np.any(np.linalg.norm(out_A.V - p_mid, axis=1) < 1e-2)
    in_B = np.any(np.linalg.norm(out_B.V - p_mid, axis=1) < 1e-2)
    assert in_A == in_B, (
        f"Polyline pierce {p_mid} present in A:{in_A}, B:{in_B} — "
        f"per-parent snap drift; faults disagree on polyline geometry.")
```

---

### [P2-007] [MODERATE] [§Phase 1 "Numerical robustness" — `eps_orient * mean_edge_squared`] — Dimensional analysis is wrong: orient3d has units of length³ but the threshold uses length²

**Category:** BUG (wrong units)

**Description:**
The plan says (line ~411):

> When the orientation magnitude for any vertex of A relative to
> plane B is below `eps_orient * mean_edge_squared`, mark the case
> as ambiguous

The *signed distance* to a plane has units of [length]. The
*orient3d determinant* is `(B1-B0) × (B2-B0) · (Ai-B0)` which has
units of [length]³ (= 6 × tetrahedron volume).

If the implementer interprets "orientation magnitude" as signed
distance, the threshold `eps_orient × mean_edge²` (units length²)
is dimensionally wrong by one factor of length.

If the implementer interprets it as orient3d, `eps_orient ×
mean_edge²` is wrong by one factor of length the *other* way (orient3d
is length³, not length²).

For SAFS data (mean_edge ≈ 2000 m), `eps_orient × mean_edge² = 1e-9
× 4e6 = 4e-3` — and we don't know what units that is. If it's
meant for signed distance, 4 mm — extremely tight, almost everything
is "ambiguous". If for orient3d, 4 mm³ on triangles of 2000 m scale,
massively too lax (orient3d on real CFM triangles is on order
1e9 m³).

This will silently mis-flag near-tangent cases and either spam the
"ambiguous → exact arithmetic" branch (wasting time) or never enter
it (defeating the safety net).

**Trigger:** Every call to `tri_tri_intersect_3d` with default
`eps_orient`.

**Actual behavior:** Wrong threshold; either too strict or too lax
depending on what "orientation magnitude" means in the implementer's
head.

**Expected behavior:** Pick one quantity and use the correct units.
For signed distance threshold:
- `signed_distance_threshold = eps_orient × mean_edge` (length × dimensionless = length).

For orient3d threshold:
- `orient3d_threshold = eps_orient × mean_edge³` (length³ × dimensionless = length³).

**Suggested fix:**

```diff
 **Numerical robustness**:
 - Default predicates use `np.float64`. When the orientation magnitude
   for any vertex of A relative to plane B is below
-  `eps_orient * mean_edge_squared`, mark the case as ambiguous and
+  `eps_orient * mean_edge` (for signed-distance predicates, units of
+  length) or `eps_orient * mean_edge**3` (for orient3d-determinant
+  predicates, units of length³) — pick the threshold that matches
+  the quantity you actually compute.  Mark the case as ambiguous and
   fall back to **exact arithmetic**.
```

**Test case:**
```python
def test_P2_007_ambiguity_threshold_dimensional():
    """Construct a triangle pair where vertex A_i is at signed
    distance d from plane B with d = 0.5 * eps_orient * mean_edge
    (i.e. just below the dimensionally-correct threshold).  The
    ambiguity branch must fire."""
    eps = 1e-9
    edge = 2000.0
    threshold = eps * edge  # = 2e-6 m, the correct length-units threshold
    V_A = np.array([[-1000, 0, threshold * 0.5],
                    [+1000, 0, threshold * 0.5],
                    [0,    0, +1000]])
    T_A = np.array([[0, 1, 2]])
    V_B = np.array([[0, -1, 0], [0, 1, 0], [0, 0, 1000]])
    T_B = np.array([[0, 1, 2]])
    n_ambiguous = 0
    def hook(): nonlocal n_ambiguous; n_ambiguous += 1
    tri_tri_intersect_3d(V_A, T_A, V_B, T_B, "A", "B",
                          eps_orient=eps, _ambiguity_hook=hook)
    assert n_ambiguous >= 1
```

---

### [P2-008] [MODERATE] [§Phase 2 §Validation gates §2 "Polyline edge-coincidence gate"] — Vertex-match tolerance is unspecified; gate is under-defined

**Category:** ASSUMPTION (incomplete spec)

**Description:**
The gate (line ~645) says:

> every consecutive vertex pair `(P[i], P[i+1])` must appear as an
> edge in fault A's post-split mesh AND in fault B's post-split
> mesh

But "appear as an edge" requires identifying a polyline 3-D point
with a mesh vertex. The mesh has float64 coordinates; the polyline
has float64 coordinates from `tri_tri_intersect_3d`'s linear
interpolation of input vertices. The 2-D project / 3-D unproject
roundtrip in Phase 2 introduces rounding on the order of
`1e-9 × max_coord` (≈ 1e-4 m at SAFS scale). After conformalize,
fault A's vertex at the pierce point may be 1e-7 m offset from the
original polyline coordinate.

The plan does not specify a matching tolerance. If a fix-agent
defaults to `==`, the gate will report failures for every pierce due
to round-off. If they pick a too-loose tolerance (1 m), unrelated
mesh vertices may match.

**Trigger:** Gate verification on any non-trivial conformalize output.

**Actual behavior (as written):** Implementer guesses tolerance; gate
either always-fails or has false positives.

**Expected behavior:** Specify the tolerance explicitly; align with
the chain `snap_m = 1e-2`. A 2× margin (`2 × snap_m = 2e-2 m =
2 cm`) is generous and cannot collide with unrelated CFM vertices
(min spacing in 2000m CFM is ~100s of m).

**Suggested fix:**

```diff
 2. **Polyline edge-coincidence gate (per pair).** For every
    intersection polyline `P` from fault `A × B`: every consecutive
    vertex pair `(P[i], P[i+1])` must appear as an edge in fault A's
    post-split mesh AND in fault B's post-split mesh. A mismatch
    indicates a conformalizer bug (a pierce point made it onto only
    one side). Implemented as `verify_polyline_in_mesh(V, T,
-   polyline)` in `conformalize_faults.py`.
+   polyline, match_tol_m=2e-2)` in `conformalize_faults.py`.  The
+   tolerance is `2 × snap_m` (2 cm by default) — large enough to
+   absorb 2-D projection / 3-D unprojection round-off, far smaller
+   than the minimum CFM vertex spacing (~100 m), so no false
+   positives on real data.  Any mesh vertex within `match_tol_m` of
+   a polyline point is considered "the same vertex" for the purpose
+   of the gate.
```

**Test case:**
```python
def test_P2_008_polyline_gate_tolerance():
    """Conformalize a 2-fault crossing.  Verify that the gate's
    match tolerance accepts mesh-vertex round-off (~1e-7 m) as
    matching the polyline coordinate, and rejects a deliberate 1-m
    offset of one mesh vertex."""
    out = conformalize(...)   # produces V_A2, T_A2, polylines
    polyline = out.polylines[0]
    assert verify_polyline_in_mesh(out.V_A, out.T_A, polyline,
                                     match_tol_m=2e-2), \
        "round-off-level mismatch should be tolerated"
    # Inject a 1-m offset on one polyline vertex's mesh location.
    bad = out.V_A.copy()
    pl_v = polyline.points[1]
    idx = np.argmin(np.linalg.norm(bad - pl_v, axis=1))
    bad[idx] += np.array([1.0, 0, 0])
    assert not verify_polyline_in_mesh(bad, out.T_A, polyline,
                                         match_tol_m=2e-2), \
        "1-m offset should be detected"
```

---

### [P2-009] [MODERATE] [§Phase 1 "Worked numeric example" Step (4)] — Example uses a degenerate triangle where vertex A2 lies exactly on plane B; this case is precisely what the algorithm's edge-case handling is supposed to cover, but the example skips over it

**Category:** ASSUMPTION (misleading example)

**Description:**
In the worked example:
```
A = (-5, 0, -5), (5, 0, -5), (0, 0, 5)
B = (0, -0.5, -0.5), (0, 0.5, -0.5), (0, 0, 0.5)
```

Plane B is `x = 0`. Vertex A2 = `(0, 0, 5)` has x = 0, so
`d_A2 = 0` exactly. By the `tri_tri_intersect_3d`'s sign-test
semantics, `sign(d_A2) = 0`, *not* `+` or `−`. The all-same-sign
check `sign(d_A0) == sign(d_A1) == sign(d_A2)` is "−, +, 0" — three
distinct values. Step (1) does not reject. Good.

But Step (4) says "exactly two of A's three edges cross plane B".
With `d_A0 = -5, d_A1 = +5, d_A2 = 0`:
- Edge (A0, A1): signs `-, +` — crosses.
- Edge (A1, A2): signs `+, 0` — endpoint on plane; depending on
  implementation, this may or may not count as "crossing".
- Edge (A2, A0): signs `0, -` — same ambiguity.

The example then says "A's two crossing edges are (A0,A1) and
(A1,A2)" — picking arbitrarily between (A1,A2) and (A2,A0). Both
edges have an endpoint exactly on plane B (at A2 itself), so both
yield the same `p_hit = A2`. The implementer can pick either; the
example happens to pick (A1,A2) and gets `p_hit = (0, 0, 5), t = -5`.

This is **a degenerate case the algorithm prose explicitly does
not handle**: Step (4) reads "for each of those two edges (Ai, Aj)
whose signed distances differ in sign". With `d_A2 = 0`, edge
(A1, A2) has `d_A1 = +5` and `d_A2 = 0` — does `+5 vs 0` "differ
in sign"? Strictly, `sign(0)` is conventionally `0`, distinct from
`+1` or `−1`, so technically yes; but a `> 0`-vs-`< 0` test will
miss the edge entirely.

A reader who follows the example numerically will either:
- Get the right answer by accident (because A2 happens to be one of
  the answer endpoints).
- Get a wrong answer because their sign convention treats `d_A2 = 0`
  differently than the example's prose implies.

The example should either:
1. Use non-degenerate data (move A2 to e.g. `(0, 0.01, 5)` so
   `d_A2 = 0.01 ≠ 0`).
2. Explicitly say "this example is degenerate; A2 is on plane B,
   and the algorithm picks edges (A0,A1) and one of (A1,A2)/(A2,A0)
   arbitrarily — both give the same p_hit."

**Trigger:** Reader follows the example to verify their
implementation; their implementation handles `d_Ai = 0` differently.

**Actual behavior:** Reader is confused or copies the wrong
convention.

**Expected behavior:** Non-degenerate example, OR explicit
acknowledgement of the degenerate case.

**Suggested fix:** Replace A's third vertex to break the degeneracy:

```diff
-A = (-5, 0, -5), (5, 0, -5), (0, 0, 5)         (a triangle in the xz-plane, y=0)
+A = (-5, 0, -5), (5, 0, -5), (0, 0.01, 5)      (very nearly the xz-plane; y=0.01 at A2 keeps A2 off plane B at x=0… wait, A2's x is 0, not y; the degenerate condition is d_A2 = n_B · A2 + d_B; need A2.x ≠ 0 to break it)
+A = (-5, 0, -5), (5, 0, -5), (0.01, 0, 5)      (a triangle approximately in the xz-plane; A2.x = 0.01 keeps A2 strictly off plane B at x=0)

 B = (0, -0.5, -0.5), (0, 0.5, -0.5), (0, 0, 0.5)  (a triangle in the yz-plane, x=0)

 Step (1):  n_B = (1, 0, 0), d_B = 0
-           d_A0 = -5,  d_A1 = +5,  d_A2 = 0   → mixed signs ✓
+           d_A0 = -5,  d_A1 = +5,  d_A2 = +0.01   → mixed signs ✓
 Step (2):  n_A = (0, -50, 0), d_A = 0
-           d_B0 = +25, d_B1 = -25, d_B2 = 0   → mixed signs ✓
+           d_B0 ≈ +25, d_B1 ≈ -25, d_B2 ≈ 0   → strictly speaking still degenerate; pick a less-pathological B too
```

Or — and I recommend this — pick a clean fully non-degenerate
example.  Suggested replacement:

```
A = (-5, -1, -5), ( 5, -1, -5), (0,  3,  5)
B = ( 0, -2, -2), ( 0,  2, -2), (0,  0,  2)
```
A is a triangle tilted in y; B is a triangle in plane x = 0. None
of A's vertices are on x = 0. The intersection is along the line
{x = 0, y from interpolation, z from interpolation}, the math is
slightly less round, but the example is non-degenerate and
demonstrably exercises every algorithm step without ambiguity.

**Test case:**
```python
def test_P2_009_worked_example_is_non_degenerate():
    """The worked example in the plan must not require any
    'd_Ai = 0' branch — every signed distance must be strictly
    non-zero."""
    A = np.array([[-5, -1, -5], [5, -1, -5], [0, 3, 5]])
    B = np.array([[0, -2, -2], [0, 2, -2], [0, 0, 2]])
    n_B = np.cross(B[1] - B[0], B[2] - B[0])
    d_B = -n_B.dot(B[0])
    d_A = n_B.dot(A.T) + d_B
    n_A = np.cross(A[1] - A[0], A[2] - A[0])
    d_A_off = -n_A.dot(A[0])
    d_B_local = n_A.dot(B.T) + d_A_off
    assert all(abs(d) > 1e-6 for d in d_A), d_A
    assert all(abs(d) > 1e-6 for d in d_B_local), d_B_local
```

---

### [P2-010] [LOW] [§Phase 1 §Edge cases] — "Triangle B's edge lies in A's plane (1-D intersection). Treat as zero-area intersection; do not emit a segment." — wrong; this is a positive-length-segment case

**Category:** BUG (wrong specification)

**Description:**
The edge-case rule (line ~457) says:

> **Triangle B's edge lies in A's plane (1-D intersection).** Treat
> as zero-area intersection; do not emit a segment.

This is geometrically wrong. If B has an *edge* lying in A's plane,
that edge is a 1-D segment in 3-D, not a 2-D area. Specifically:
- The two endpoints of B's edge have `d_Bi = d_Bj = 0` (both on
  plane A).
- B's third vertex has `d_Bk ≠ 0`.
- Step (2)'s all-same-sign check sees `0, 0, ≠0` — not all same
  sign, so does not reject.

The intersection is the segment from `B_i` to `B_j` along that
edge — a real, length-positive segment that should be emitted (after
clipping against triangle A's interior).

The plan's "Treat as zero-area; do not emit" is the wrong response.
For SAFS data this case is vanishingly rare (CFM faults aren't
constructed to share planes), but the specification is incorrect as
written. Worse, if it occurs in the cascade phase (Phase 2's edge-
pierce propagation creates new edges that COULD coincide with another
fault's plane), the algorithm silently drops a real intersection.

**Trigger:** Triangle B has an edge with both endpoints exactly on
plane A. Rare in raw CFM, possible after Phase 2 cascade
re-triangulation.

**Actual behavior (as written):** Drops the intersection segment.

**Expected behavior:** Either (a) detect this case and emit B's
edge as the intersection segment, or (b) in `tri_tri_intersect_3d`'s
default float64 path, treat `d_Bi == 0 && d_Bj == 0` as
`d_Bi = +eps && d_Bj = +eps` (perturbation) and let the standard
algorithm handle it.

**Suggested fix:**

```diff
 ### Edge cases to handle
-- **Triangle B's edge lies in A's plane (1-D intersection).** Treat as
-  zero-area intersection; do not emit a segment.
+- **Triangle B's edge lies in A's plane (B has two consecutive
+  vertices with d_Bi = d_Bj = 0).** This is a positive-length 1-D
+  intersection along that edge, not a zero-area case.  Detect it
+  before running the main interpolation: if any pair of B's vertices
+  has |d| < eps_orient × mean_edge, clip the connecting edge against
+  triangle A's interior (using A's barycentric coordinates) and
+  emit the resulting sub-segment.  If the case occurs in cascade
+  phase post-Step-(3), it indicates two faults have started sharing
+  a plane — log loudly, do not silently drop.
+- **Both faults nearly coplanar but not exactly (dihedral < 1°).**
+  Already rejected by the `|n_A · n_B| > 1 − 1e-6` guard in
+  Numerical robustness §.
```

**Test case:** Synthetic — construct a case and verify the segment
is emitted (omitted; LOW severity).

---

### [P2-011] [LOW] [§Phase 5 acceptance criteria] — "polyline #2's input triangulation differs from polyline #1's input" is a tautology

**Category:** QUALITY (weak acceptance criterion)

**Description:**
The acceptance metric (line ~812) reads:

> Cascading metric: `intersection_report.json` shows that polyline
> #2's input triangulation differs from polyline #1's input (i.e.,
> the conformalizer is operating on the already-split surface, not
> the original).

This is true *whenever cascading runs at all*. Any non-zero polyline
#1 changes the triangulation before polyline #2 is processed. The
metric does not distinguish a *correct* cascade from one that
silently uses stale indices anyway (e.g., re-uses original `pierces_a`
into a triangulation where they no longer make sense).

**Trigger:** N/A (documentation issue).

**Actual behavior:** Acceptance passes on incorrect implementations.

**Expected behavior:** Replace with: at least one of polyline #2's
pierces lands on a triangle that did NOT exist in the pre-cascade
triangulation (i.e., a child created by polyline #1's split). This
is the property cascading is supposed to enable.

**Suggested fix:**

```diff
 - [ ] 3-fault crossing subset passes 11/11 validator checks.
-- [ ] Cascading metric: `intersection_report.json` shows that
-      polyline #2's input triangulation differs from polyline #1's
-      input (i.e., the conformalizer is operating on the
-      already-split surface, not the original).
+- [ ] **Cascade actually cascades**: at least one of polyline #2's
+      pierces lands on a child triangle created by polyline #1's
+      split (verifies that the re-scan produced fresh indices into
+      the post-#1 triangulation, not stale indices into the
+      pre-#1 original).  Implemented by recording in
+      `intersection_report.json` the parent-id (or "post-#1 child")
+      of each pierce, and asserting at least one cascade-derived
+      pierce in the 3-fault target.
```

---

## Summary

- Critical issues: 3 (P2-001, P2-002, P2-003)
- Moderate issues: 6 (P2-004, P2-005, P2-006, P2-007, P2-008, P2-009)
- Low issues: 2 (P2-010, P2-011)
- Plan compliance: **N/A** (this is a plan-document review).
- Verdict: **FAIL — must revise plan before implementing**, but the
  scale of failure is smaller than Round 1.  The Round-1 fixes did
  resolve the algorithmic correctness issues (Möller, free-surface
  clamp, T-junction propagation, validation gates).  Round-2's
  findings concentrate on:
  - Internal inconsistencies (P2-001, P2-002) — the plan now
    contradicts itself in two places.
  - Worked-example and ASCII-picture quality (P2-003, P2-004,
    P2-005, P2-009) — the readability rewrite introduced new
    documentation bugs.
  - Tolerance / units gaps (P2-006, P2-007, P2-008) — silent-failure
    risks in the implementation, all addressable with one-line
    spec edits.

Suggested revision order:
1. P2-001 (testing-strategy table) — pure search-and-replace.
2. P2-002 (cascade vs edge-case conflict) — clarifying paragraph.
3. P2-003 (Step (6) interpolation formula) — paste explicit formula.
4. P2-004, P2-005 (pictures) — redraw with 1 unambiguous case.
5. P2-006 (drop snap-to-corner) — one bullet, big behavioural impact.
6. P2-007 (units) — one-word fix per call site.
7. P2-008 (gate tolerance) — add `match_tol_m=2e-2` to one signature.
8. P2-009 (worked example) — replace input data.
9. P2-010 (1-D edge case) — rewrite one bullet.
10. P2-011 (cascade acceptance) — rewrite one bullet.

After revision, regenerate the PDF.

## Unreviewed Areas

- **Phase 3 §Edge cases "Two faults closer than 2 × tube_radius"**:
  the plan says "tubes overlap — fine; size field takes the
  minimum". This is correct but does not address whether the bulk
  mesher correctly resolves the *between-fault* region. If two
  faults are 2km apart and `tube_radius = 4km`, every tet between
  them is in both tubes; `Field[2].SizeMin = res_f` gives uniform
  res_f throughout. That's fine. No new finding.
- **Phase 4 §Files to modify "`Mesh.Algorithm3D = 10` (HXT) for
  conformal input"**: the README documents that HXT can fail even
  on apparently-conformal input due to PLC-recovery quirks. The
  plan does not specify a fallback ladder. This was Round-1's
  Recommendation §Path C; deferred there as out-of-scope. Still
  out-of-scope here.
- **Phase 4 §Files to modify "single-fault smoke pass-through"**:
  works correctly when `len(included) == 1` (no pairs → no work).
  No bug.
- **gmpy2.mpq vs Shewchuk predicates choice**: both are correct;
  the plan recommends `mpq` if `gmpy2` is in `pythonenv`. No bug.
