# Code Review of Phase 1 Implementation — 2026-04-29

## Review Scope

- **Plan**: `miniapps/seas/safs/PLAN_multifault_intersections.md` §Phase 1
- **Files reviewed**:
  - `miniapps/seas/safs/mesh/fault_intersect.py` (~440 lines)
  - `miniapps/seas/safs/mesh/tests/test_fault_intersect.py` (~280 lines)
- **Implementation report consulted**: the `/code-implement` summary in
  the conversation (notes a coplanar-skip deviation from plan and that
  the README ±10% real-fixture criterion is borderline).
- **Domain context**: `seas-mfem/CLAUDE.md`,
  `miniapps/seas/CLAUDE.md`, `PLAN_multifault_intersections.md`,
  `mesh/REVIEW.md` (existing-code review, separate scope).

This is a **fresh adversarial review** of the Phase 1 implementation.
I traced every public function, every helper, and every test through
explicit numerical scenarios. I assumed at least 3 bugs; I found 9.

---

## Findings

### [F-001] [CRITICAL] [fault_intersect.py:chain_segments] — Branched components silently orphan all but one branch

**Category:** BUG

**Description:**
For a connected component with a degree-3+ junction node (e.g., a Y, T,
or star — possible in Phase 5 cascade or whenever an intersection
polyline has a true T-junction), `chain_segments` walks one branch
from a leaf through the junction to another leaf, then **breaks**
because the next adjacency at the junction is to an already-`used`
segment in the *current* walk's set. After the walk, line 446 marks
**every** segment in the component as visited:

```python
visited_segs |= component
```

This drops the orphaned branches from the output entirely — they are
neither emitted as a separate polyline nor flagged as a warning. The
plan's edge case "Polyline branches at a T-junction" claims chaining
"produces 1 polyline per connected component" — but the implementation
produces 1 polyline that contains *only some* of the component's
segments, and silently discards the rest.

**Trigger:** Three or more segments meeting at a single snapped node
(degree ≥ 3 junction). Constructible synthetically; rare in 2-fault
data but expected in Phase 2 cascade and 3-fault tests (Phase 5).

**Actual behavior:**
For a Y of 3 segments meeting at the origin, `chain_segments` returns
**1 polyline containing 2 of the 3 segments**; the third is silently
gone.

**Expected behavior:**
Either (a) emit one polyline per branch (so a Y produces 3 polylines
sharing a common endpoint), or (b) emit a single polyline that walks
all branches with documented multi-arm semantics. Option (a) matches
the plan's "T-junction in the 3-fault case manifests as 3 polylines
meeting at a point and is fine".

**Suggested fix:**
Change the post-walk bookkeeping to mark only what was actually
walked, and let the outer loop pick up the remaining branches as
new polylines:

```diff
-        visited_segs |= component
+        visited_segs |= used
```

Then the within-walk filter must also exclude already-visited segments
across outer iterations (currently it only checks the local `used`):

```diff
         while True:
             next_sidx: int | None = None
             next_node: int | None = None
             for sidx, other in adj[cur]:
-                if sidx in component and sidx not in used:
+                if (sidx in component
+                        and sidx not in used
+                        and sidx not in visited_segs):
                     next_sidx = sidx
                     next_node = other
                     break
```

Endpoint-finding (the degree-1 search) must also exclude already-
visited segments from the degree count:

```diff
         endpoint_node: int | None = None
         for sidx in component:
+            if sidx in visited_segs:
+                continue
             for node_id in seg_endpoints[sidx]:
                 deg_in_component = sum(
-                    1 for (sx, _) in adj[node_id] if sx in component
+                    1 for (sx, _) in adj[node_id]
+                    if sx in component and sx not in visited_segs
                 )
                 if deg_in_component == 1:
                     endpoint_node = node_id
                     break
```

**Test case:**
```python
def test_F001_branched_component_emits_all_segments():
    """A T-junction (3 segments meeting at the origin) must yield
    polylines covering all 3 segments — currently the implementation
    drops one branch silently."""
    p_origin = (0.0, 0.0, 0.0)
    p1 = (1.0, 0.0, 0.0)
    p2 = (-1.0, 0.0, 0.0)
    p3 = (0.0, 1.0, 0.0)
    segs = [
        fi.Segment(p_origin, p1, 0, 0, "A", "B", True),
        fi.Segment(p_origin, p2, 1, 1, "A", "B", True),
        fi.Segment(p_origin, p3, 2, 2, "A", "B", True),
    ]
    plines = fi.chain_segments(segs, snap_m=0.01)
    n_segs_total = sum(len(p.points) - 1 for p in plines)
    assert n_segs_total == 3, (
        f"chain_segments dropped a branch: "
        f"3 input segments yielded only {n_segs_total} polyline edges"
    )
```

---

### [F-002] [MODERATE] [fault_intersect.py:tri_tri_intersect_3d, lines ~424–445] — Plan's Phase-1 step ordering reversed: free-surface clamp executes before the length filter

**Category:** DEVIATION

**Description:**
PLAN_multifault_intersections.md Phase 1 specifies the order:

> **(7) Length filter.** if `|p_end − p_start| < eps_min_seg_len_m`:
>   mark `nondegenerate=False`; skip
> **(8) Free-surface clamp** …

The implementation does step (8) first, then step (7):

```python
if clearance_m > 0.0:
    clamped = _free_surface_clamp(p_start, p_end, clearance_m)
    if clamped is None:
        continue
    p_start, p_end = clamped

seg_len = float(np.linalg.norm(p_end - p_start))
if seg_len < eps_min_seg_len_m:
    continue
```

The implementer's `/code-implement` report acknowledged this deviation
("clamp first, then length-filter") and argued it is "more sensible".
That may be defensible, but the deviation is silent — neither the
docstring nor a comment flags it, so an adversarial reviewer like
me had to spot it from the source.

**Trigger:** Any segment longer than `eps_min_seg_len_m` *before* the
free-surface clamp but shorter *after* (e.g., a 60 mm segment with
57 mm above the clearance plane → 3 mm post-clamp).

**Actual behavior:** Post-clamp 3 mm segment is dropped (`< 5 cm`).
Plan would have kept the pre-clamp 60 mm segment, then clamped to 3
mm, then emitted the 3 mm segment.

**Expected behavior:** Either match the plan order, or document the
deviation in the function docstring and the implementation report.
Both behaviors are arguably correct; the inconsistency is the bug.

**Suggested fix (option A — match plan order):**

```diff
-        if clearance_m > 0.0:
-            clamped = _free_surface_clamp(p_start, p_end, clearance_m)
-            if clamped is None:
-                continue
-            p_start, p_end = clamped
-
         seg_len = float(np.linalg.norm(p_end - p_start))
         if seg_len < eps_min_seg_len_m:
             continue
+
+        if clearance_m > 0.0:
+            clamped = _free_surface_clamp(p_start, p_end, clearance_m)
+            if clamped is None:
+                continue
+            p_start, p_end = clamped
```

**Suggested fix (option B — keep current order, document):**

```diff
+    # NOTE: Plan §Phase-1 specifies length-filter (step 7) BEFORE
+    # free-surface clamp (step 8); this implementation does clamp
+    # first because a long pre-clamp segment can become a sub-cm
+    # post-clamp sliver that no caller wants.  Both orderings drop
+    # spurious sub-eps segments; only segments whose length straddles
+    # eps across the clamp see different behaviour.
     if clearance_m > 0.0:
         clamped = _free_surface_clamp(p_start, p_end, clearance_m)
```

Pick option B if the user wants the safer behaviour; pick option A
if literal plan compliance matters more.

**Test case:**
```python
def test_F002_length_filter_runs_before_clamp_per_plan():
    """A pre-clamp 60 mm segment that becomes a 3 mm segment post-
    clamp.  Plan says: emit the 3 mm clamped segment.  Implementation
    currently drops it (clamp-then-length).  Whichever order is
    chosen, this test pins the contract."""
    # Construct triangles whose intersection runs from z = -0.060 to
    # z = +0.057.  Clearance = 0.060 → top endpoint clamps to z = -0.060.
    # Pre-clamp length = 0.117 m (> 0.05 = eps_min_seg_len_m → keep).
    # Post-clamp length = 0.000 m → ambiguous behaviour.
    # Adjust the geometry to land exactly in the contested zone.
    eps_min = 0.05
    pre_len = 0.10
    post_len = 0.005   # below eps after clamp
    # ... (construct two crossing squares with intersection of pre_len
    # spanning z = -clearance + post_len down to z = -clearance + post_len - pre_len)
    # Run with clearance = (top z + 0.001) so we clamp.
    # Assert the chosen contract.
```

---

### [F-003] [MODERATE] [test_fault_intersect.py:test_two_perpendicular_unit_squares_intersection_endpoints] — Acceptance assertion is too loose; allows a regression that drops half the segments to pass

**Category:** QUALITY (test gap)

**Description:**
The plan acceptance criterion says "in ≤ 2 of the 4 triangles".
Geometric truth on the 90°-cross test: **exactly 2 of the 4
(i, j) pairs** produce 1 m segments. The other 2 pairs produce
zero-length overlaps that are filtered out.

The test asserts:

```python
assert 1 <= len(segs) <= 2, f"expected 1 or 2 segments, got {len(segs)}"
```

A regression that emits only 1 segment would still pass this
assertion. The acceptance criterion is "≤ 2", but "exactly 2" is
required for full conformality coverage (each triangle of A must
register a crossing with at least one triangle of B).

**Trigger:** Any regression that drops segments while keeping the
unit-square fixture identical.

**Expected behavior:** Tighten to `len(segs) == 2`. There is no
geometric ambiguity in the fixture — the answer is exact.

**Suggested fix:**

```diff
-    assert 1 <= len(segs) <= 2, f"expected 1 or 2 segments, got {len(segs)}"
+    assert len(segs) == 2, (
+        f"expected exactly 2 segments (the 4 triangle pairs decompose "
+        f"into 2 zero-overlap and 2 one-metre segments); got {len(segs)}"
+    )
```

**Test case:** the test itself is the test. After tightening, run it.

---

### [F-004] [MODERATE] [test_fault_intersect.py:test_chain_segments_one_polyline_from_perpendicular_squares] — Polyline-length assertion too loose

**Category:** QUALITY (test gap)

**Description:**
The 90°-cross intersection runs from z = -1 to z = +1, total length
**exactly 2.0 m** along the z-axis. The test asserts:

```python
total = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
assert total >= 1.0, f"polyline total length {total} < 1.0 m"
```

`>= 1.0` permits half the diagonal to be lost. If F-001 (branched
component) regresses to drop half a polyline, this test still passes
because half a 2 m polyline is 1 m, which is `>= 1.0`.

**Suggested fix:**

```diff
-    assert total >= 1.0, f"polyline total length {total} < 1.0 m"
+    assert abs(total - 2.0) < 1e-9, (
+        f"polyline total length {total} ≠ 2.0 m (the geometric "
+        f"truth for the 90°-crossing unit-square fixture)"
+    )
```

---

### [F-005] [MODERATE] [test_fault_intersect.py — missing] — Real-fixture acceptance criterion is not asserted by any test

**Category:** DEVIATION

**Description:**
The plan's Phase-1 acceptance criteria include:

> `scan_cross_fault_crossings("output/all8_2000m/stl", [<all 8>])`
> returns the README-documented crossing-pair counts to within ±10%
> (Mill × SBMT-SAF ≈ 163, Mission × SAF ≈ 135, Pinto × Mill ≈ 42, etc.)

> Run time of `scan_cross_fault_crossings` on the 8-fault 2000m
> set: ≤ 30 s on a laptop.

Neither is asserted in `test_fault_intersect.py`. The implementer ran
the scan as a manual smoke check, reported 0.51 s and the per-pair
segment counts in the conversation, but did not pin any of those
numbers as a test.

Consequence: a future regression that, e.g., breaks the AABB filter
(making the scan 100× slower) or skips an entire pair (making
Mill × SBMT-SAF return 0 segments) would pass the entire CI suite.

**Trigger:** Any regression in scan correctness or performance.

**Expected behavior:** A test that runs the all-8 scan when the
fixture is present, asserts the runtime budget and pair-count
budget, and is `pytest.mark.skipif`'d when the fixture is absent.

**Suggested fix:** add to `tests/test_fault_intersect.py`:

```python
import time

@pytest.mark.skipif(
    not (Path(__file__).resolve().parent.parent
         / "output" / "all8_2000m" / "stl"
         / "safs_sbmt_saf.stl").exists(),
    reason="all-8 fixture not generated; run run_all8_2000m.sh first"
)
def test_F005_real_fixture_scan_runtime_and_counts():
    """PLAN Phase-1 real-fixture acceptance: 8-fault scan completes in
    ≤ 30 s and reports approximately README-documented pair counts."""
    stl_dir = (Path(__file__).resolve().parent.parent
               / "output" / "all8_2000m" / "stl")
    included = [
        "safs_coav_missioncreek", "safs_mjvs_saf",
        "safs_mult_banning", "safs_mult_ssaf_banning",
        "safs_pmfz_pinto", "safs_sbmt_millcreek",
        "safs_sbmt_missioncreek", "safs_sbmt_saf",
    ]
    t0 = time.perf_counter()
    out = fi.scan_cross_fault_crossings(stl_dir, included, clearance_m=100.0)
    elapsed = time.perf_counter() - t0
    assert elapsed <= 30.0, f"scan took {elapsed:.1f} s (budget 30 s)"

    # PLAN cites these README-documented counts; allow ±25% to absorb
    # the implementation's coplanar-skip + clearance-clamp filters
    # (documented deviation in /code-implement report).  Tighten to
    # ±10% once gmpy2 is available and the coplanar-skip is replaced.
    expected_min = {
        ("safs_sbmt_millcreek",     "safs_sbmt_saf"):       100,    # README ~163
        ("safs_sbmt_missioncreek",  "safs_sbmt_saf"):        80,    # README ~135
        ("safs_pmfz_pinto",         "safs_sbmt_millcreek"): 25,    # README ~ 42
    }
    counts = {pair: sum(len(p.points) - 1 for p in plines)
              for pair, plines in out.items()}
    for pair, lower in expected_min.items():
        assert pair in counts, f"expected pair {pair} not detected"
        assert counts[pair] >= lower, (
            f"{pair} reported only {counts[pair]} segments; expected >= "
            f"{lower} (README-documented)"
        )
```

---

### [F-006] [MODERATE] [test_fault_intersect.py:test_aabb_filter_skips_far_apart_pairs] — Test does not isolate AABB-filter behavior; coplanar guard would also produce empty output

**Category:** QUALITY (test gap)

**Description:**
The test fixture uses two triangles, both with z = 0 (in the same
plane). If the AABB filter were broken (so the pair entered the Möller
core), the coplanar guard would still warn-and-skip. The test cannot
distinguish "AABB filter rejected" from "coplanar guard rejected".

**Trigger:** A regression that makes the AABB filter always-true
(returns full overlap matrix). The test still passes because the
coplanar guard catches the otherwise-broken case.

**Suggested fix:** Move B's triangle out of the z=0 plane so it is
non-coplanar with A AND has a non-overlapping AABB.

```diff
     V_A = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float64)
     T_A = np.array([[0, 1, 2]], dtype=np.int64)
-    V_B = np.array([[100, 0, 0], [101, 0, 0], [100, 1, 0]], dtype=np.float64)
+    V_B = np.array([[100, 0, 5], [101, 0, 5], [100, 1, 5]], dtype=np.float64)
     T_B = np.array([[0, 1, 2]], dtype=np.int64)
     segs = fi.tri_tri_intersect_3d(V_A, T_A, V_B, T_B, "A", "B")
     assert segs == []
```

Now A is in z=0, B is in z=5 — non-coplanar, non-overlapping AABBs.
A regression in the AABB filter would actually cause the per-pair
loop to run (and Möller to correctly find no crossing); regression in
Möller would be caught separately. The test now isolates AABB.

---

### [F-007] [LOW] [fault_intersect.py:_categorise vs ambiguity test] — Inconsistent strict-vs-inclusive boundary semantics

**Category:** ASSUMPTION

**Description:**
`_categorise` (line 95) treats `d == thresh` as "in band" (s = 0):

```python
s[d > thresh] = 1     # strict >
s[d < -thresh] = -1   # strict <
```

The ambiguity test (line 392) uses strict less-than:

```python
ambiguous = (
    (np.abs(d_A_to_B) < thresh_A).any()
    or (np.abs(d_B_to_A) < thresh_B).any()
)
```

A vertex with `|d| == thresh` is categorized as 0 (band) but **not**
flagged as ambiguous (`|d| < thresh` is False). Such a vertex bypasses
the gmpy2 fallback and proceeds to interval calculation, which handles
it correctly in practice (the vertex becomes a p_hit on the plane).

So this is a correctness-preserving inconsistency, but the asymmetry
will trip up future maintainers and may interact subtly with FP
edge cases.

**Suggested fix:** Use `<=` consistently:

```diff
     ambiguous = (
-        (np.abs(d_A_to_B) < thresh_A).any()
-        or (np.abs(d_B_to_A) < thresh_B).any()
+        (np.abs(d_A_to_B) <= thresh_A).any()
+        or (np.abs(d_B_to_A) <= thresh_B).any()
     )
```

(Or change `_categorise` to strict `>=`/`<=`. Either choice is fine
as long as both match.)

---

### [F-008] [LOW] [fault_intersect.py:Segment.nondegenerate] — Field is dead; always True for emitted segments

**Category:** QUALITY (dead code)

**Description:**
The plan's wording "mark `nondegenerate=False`; skip" was ambiguous.
The implementation chose "skip emit", so every emitted Segment has
`nondegenerate=True`. The field is unused information.

```python
out.append(Segment(
    ...
    nondegenerate=True,   # always; never False
))
```

**Suggested fix:** Remove the field, or repurpose it (e.g., emit
`nondegenerate=False` segments alongside True ones for caller
inspection). Removing is simpler.

---

### [F-009] [LOW] [fault_intersect.py:_interval_on_L line 158] — `edge_on_plane` branch is unreachable in the gmpy2-absent code path

**Category:** QUALITY (dead branch in current configuration)

**Description:**
`_interval_on_L` raises `edge_on_plane = True` when both endpoints of
an edge are within `thresh_zero` of the other plane. But the caller
`tri_tri_intersect_3d` runs the ambiguity check **before** calling
`_interval_on_L`. With even one vertex within `thresh`, ambiguity
fires; without `gmpy2`, we `continue` past `_interval_on_L` entirely.

So the `edge_on_plane` branch (and its accompanying warning) is dead
in the gmpy2-absent default configuration. With gmpy2 it would be
reachable but the test suite has no `gmpy2`-installed run.

**Suggested fix:** Either install gmpy2 in `pythonenv` and add a
gmpy2-on test for the edge case, or document the unreachability in
the comment block above the branch.

---

## Summary

- Critical issues: **1** (F-001)
- Moderate issues: **5** (F-002, F-003, F-004, F-005, F-006)
- Low issues: **3** (F-007, F-008, F-009)
- Plan compliance: **PARTIAL**
  - Plan's "Phase 1 acceptance: scan returns README-documented counts
    ±10%" is reported in the implementation summary but not asserted
    by any test (F-005).
  - Plan's "step ordering — length-filter then clamp" is reversed
    (F-002).
  - Coplanar-pair handling deviates from the plan's `raise` to a
    `warn + skip` (documented in the implementer's report; not a new
    finding).
- Verdict: **PASS WITH FIXES** — no incorrect-result bug at the
  Phase-1 acceptance fixture level (16/16 tests pass, all-8 scan
  completes), but F-001 will cause silent data loss as soon as the
  pipeline ever sees a branching polyline (Phase 2 cascade or hand-
  written test). Fix F-001 + F-005 before considering Phase 1 done.

Suggested fix order:
1. **F-001** (chain_segments branch orphaning) — apply the diff,
   add the test, verify the test now passes and pre-fix it failed.
2. **F-005** (real-fixture acceptance test) — add the test; pin the
   current counts as a lower bound so future regressions are caught.
3. **F-002** (step ordering) — pick option A or option B; document.
4. **F-003**, **F-004** (test tightening) — one-line changes.
5. **F-006** (AABB test isolation) — one-line change.
6. **F-007** (band semantics) — cosmetic; one-line.
7. **F-008** (dead field), **F-009** (dead branch) — cleanup; defer
   to a tidy-up pass if convenient.

## Unreviewed Areas

- **`_resolve_with_mpq`**: gmpy2 fallback path. Implemented but no
  test exercises it (gmpy2 absent in `pythonenv`). Did not deeply
  trace the predicate logic — the function is unreachable in the
  current CI configuration. Re-review when gmpy2 lands.
- **The 1-D-edge-in-plane intersection case** (`edge_on_plane`
  branch + plan §Phase 1 §Edge cases): not implemented; logged-
  warn-and-skip. The implementer noted this as a known limitation.
  Out of scope for this review round.
- **Threshold scope choice (`mean_edge_A` vs `mean_edge_B`)**:
  implementation uses `mean_edge` of the *same* triangle as the
  vertex being categorised. Plan does not specify; defensible.
  Worth revisiting if real-fixture results show systematic
  under-detection on small-triangle pairs.
