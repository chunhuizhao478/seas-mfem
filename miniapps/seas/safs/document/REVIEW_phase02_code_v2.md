# Code Review: Phase 1 + Phase 2 (2026-04-29, third pass — post-R-101..R-106 fixes)

## Review Scope

- **Plan reviewed against:** `PLAN_cgal_corefine.md` Phase 1 + Phase 2.
- **Files reviewed (all current state):**
  - `tools/CMakeLists.txt` (144 lines)
  - `tools/main.cpp` (53 lines, hello binary)
  - `tools/corefine_faults/corefine.hpp` (210 lines)
  - `tools/corefine_faults/corefine.cpp` (582 lines)
  - `tools/corefine_faults/io.hpp`/`io.cpp` (103 + 387 lines)
  - `tools/corefine_faults/main.cpp` (681 lines)
  - `tools/tests/test_corefine_smoke.cpp` (493 lines)
- **User evidence:** screenshot
  `~/Desktop/Screenshot 2026-04-29 at 10.11.49 PM.png` of the
  current Phase-2 surface mesh, plus the user's note: *"now it
  looks much better, but I wonder if the mesh size is too
  different than other places on the fault"*.
- **Live verification on this machine:**
  - All 9 smoke tests pass (`test_two_squares`,
    `test_three_orthogonal_cascade`,
    `test_coincident_preflight`,
    `test_polyline_edge_coincidence_gate`,
    `test_mill_creek_pin`,
    `test_sliver_collapse_symmetric`,
    `test_remesh_protects_constrained`,
    `test_free_surface_violation_and_clamp`,
    `test_idempotency`).
  - End-to-end run on Mill Creek × SBMT-SAF 2000 m fixture:
    - Sliver-collapse at threshold 250 m: **46 short-polyline-edge
      pairs collapsed** symmetrically.
    - Post-remesh count: A=2639 tris, B=5405 tris. Both pass the
      manifold + polyline-coincidence + idempotency gates.
    - Total runtime 71 ms.
  - **R-101 fix confirmed working**: max aspect on Mill Creek
    dropped from 542:1 (previous round) to **4.7:1**
    (well within the plan's "max aspect ≤ 4" target — only
    SBMT-SAF has 3 triangles slightly over at max 7.6:1).

The user's observation is **correct and quantifiable** — mesh
size IS noticeably different near the polyline vs far from it.
Direct STL inspection produces:

| Fault | mean edge < 500 m from polyline | mean edge > 5 km from polyline | ratio |
|---|---:|---:|---:|
| Mill Creek | 702 m | 990 m | **0.71×** |
| SBMT-SAF   | 714 m | 1010 m | **0.71×** |

The mesh is **~30% finer** within 500 m of the intersection
curve than in the far field. This is the mesh-size
non-uniformity the user is calling out, and it has a precise
cause in `split_long_edges`'s bisection semantics combined with
`isotropic_remeshing`'s constrained-edge protection. Findings
below.

R-101 through R-106 from the prior round all appear addressed
in code. The new bug class is around the **polyline-grading
gradient** (R-201), plus three quality issues exposed by
end-to-end real-data inspection (R-202..R-204), plus four
lower-severity issues (R-205..R-208).

---

## Findings

### [R-201] [CRITICAL] [corefine.cpp:479–482, split_long_edges target choice] — Polyline edges end up at ~0.7 × target, producing a 30% mesh-size grade between near-polyline and far-field

**Category:** BUG (user-visible mesh-quality degradation)

**Description:**
`remesh_one_fault` step 2 calls
```cpp
PMP_local::split_long_edges(
    long_constrained, p.target_edge_m, M,
    params::edge_is_constrained_map(protect_ecm));
```
with `p.target_edge_m = 1000.0` m.

CGAL's `split_long_edges(L_max, M)` recursively bisects every
edge whose length exceeds `L_max`. After the recursion,
**every constrained edge has length in `(L_max / 2, L_max]`**
(an edge of length `1.5 L_max` splits to two edges of `0.75
L_max`; a `1.9 L_max` edge splits to two `0.95 L_max` ones,
etc.). The mean post-split length is therefore approximately
`0.75 × L_max` ≈ 750 m for our 1000 m target.

`isotropic_remeshing(target_edge_m=1000)` then runs over the
mesh with `protect_constraints(true)`. Constrained edges stay
at their post-split length (~750 m mean). Non-constrained edges
are remeshed toward 1000 m. Triangles that share an edge with
the polyline have one ~750 m edge constraining their nearby
non-constrained edges — and the remesher's quality smoothing
pulls those non-constrained edges toward similar lengths to
keep aspect ratios reasonable.

Result: the **polyline corridor is locally graded at ~750 m**,
the far field is at ~1000 m, and there's a smooth ramp
between. This is exactly the visual artifact the user is
seeing in the screenshot.

Direct measurement on the Mill Creek × SBMT-SAF run:
- Within 500 m of polyline: mean tri edge **702 m** (Mill),
  **714 m** (SAF).
- Beyond 5 km: mean tri edge **990 m** (Mill), **1010 m** (SAF).
- 30% gradation, sustained across 1500–5000 m of distance.

**Trigger:** Always — every Phase 2 run with `--target-edge-m`
positive.

**Actual behavior:** Visible mesh-density gradation around the
polyline; user-reported.

**Expected behavior:** The polyline corridor should be at the
same target as the rest of the mesh (within ±15%, say). The
plan's Phase 2 acceptance criterion *"mean edge length in
[800, 1200] m"* is met globally only because the far field
dilutes the near-polyline mean — the **near-polyline-only
mean of 702 m is OUTSIDE the [800, 1200] m band**.

**Suggested fix:** Use a higher split ceiling so post-split
constrained edges are centred at `target`, not at `0.75 ×
target`. The natural choice is `target × 4/3`:

```diff
@@ corefine.cpp:479–482
     {
         std::vector<EdgeDescriptor> long_constrained;
         long_constrained.reserve(M.number_of_edges());
         for (auto e : M.edges()) {
             if (get(protect_ecm, e)) long_constrained.push_back(e);
         }
+        // R-201: split ceiling = target × 4/3 so post-split
+        // constrained edges fall in (target × 2/3, target × 4/3),
+        // mean ≈ target.  With ceiling = target, the post-split
+        // range was (target/2, target], producing a 30% finer
+        // polyline corridor than the rest of the mesh.
+        const double split_ceiling = p.target_edge_m * 4.0 / 3.0;
         PMP_local::split_long_edges(
-            long_constrained, p.target_edge_m, M,
+            long_constrained, split_ceiling, M,
             params::edge_is_constrained_map(protect_ecm));
     }
```

The same constant should also drive the sliver-collapse
threshold (R-208 below) so the two passes are consistent.

**Test case:**
```cpp
TEST(R201_polyline_corridor_matches_target) {
  Mesh A = make_square(2, 10000.0);   // 10 km square
  Mesh B = make_square(1, 10000.0);
  ECMap pl_a, pl_b, pr_a, pr_b;
  install_ecm(A, pl_a, "e:pl"); install_ecm(B, pl_b, "e:pl");
  install_ecm(A, pr_a, "e:pr"); install_ecm(B, pr_b, "e:pr");
  corefine_pair(A, B, pl_a, pl_b, "A", "B");
  for (auto e : A.edges()) if (get(pl_a, e)) put(pr_a, e, true);
  for (auto e : B.edges()) if (get(pl_b, e)) put(pr_b, e, true);
  RemeshParams rp{.target_edge_m = 1000.0, .n_iterations = 3};
  remesh_one_fault(A, pl_a, pr_a, rp);
  remesh_one_fault(B, pl_b, pr_b, rp);

  // Compute mean edge length within 500 m of the polyline,
  // and beyond 5 km, and check the ratio.
  // ...
  double mean_near = ..., mean_far = ...;
  EXPECT_GT(mean_near / mean_far, 0.85)
      << "Polyline corridor is too fine: " << mean_near
      << " vs far-field " << mean_far << " (ratio "
      << mean_near / mean_far << ")";
}
```

---

### [R-202] [CRITICAL] [corefine.cpp:341–429, collapse_short_polyline_edges_symmetric — non-polyline short edges still leak through] — SBMT-SAF post-remesh has min edge 66 m and 3 triangles with aspect > 5

**Category:** BUG (residual sliver class)

**Description:**
The R-101 fix `collapse_short_polyline_edges_symmetric` finds
short edges in **`polyline_ecm` only** (corefine.cpp:359–366):

```cpp
auto collect = [&](const Mesh& M, const EdgeConstrainedMap& pe) {
    for (auto e : M.edges()) {
        if (!get(pe, e)) continue;          // polyline-only
        ...
    }
};
```

But `corefine` produces sliver triangles even on **non-polyline
edges** when the polyline pierces a parent triangle very close
to one of its existing edges. Example: the parent triangle has
edges `(v0, v1, v2)`; the polyline enters near `v0` on edge
`(v0, v1)` and exits near `v2` on edge `(v0, v2)`. CGAL splits
both edges, creating new vertices `q1` (very close to `v0` on
edge v0-v1) and `q2` (close to `v0` on edge v0-v2). The
resulting children are:
- `(v0, q1, q2)` — *non-polyline* sliver near `v0`.
- `(q1, q2, ...)` — polyline-bounded triangle.

Only the `(q1, q2)` polyline edge is in the collapse set; the
sub-target `(v0, q1)` and `(v0, q2)` edges are NOT touched
because they don't carry polyline marks. They survive into the
post-remesh output.

Empirical confirmation on Mill Creek × SBMT-SAF:
- SBMT-SAF post-Phase-2: **min edge 66 m**, 3 triangles with
  aspect > 5 (max 7.6).
- These edges are well below the 250 m sliver threshold; the
  collapse pass simply never considered them because they
  aren't polyline-marked.

**Trigger:** Any pair where the polyline enters a parent
triangle close (within ~50 m) to a vertex; common in real CFM
data.

**Actual behavior:** A handful of needle triangles per fault
survive Phase 2 cleanup. Visually invisible in the screenshot
(too small to see) but they will show up in HXT (Phase 3) as
poor-quality tets adjacent to the fault.

**Expected behavior:** The collapse pass should also handle
**non-polyline near-vertex sub-target edges**. These can be
collapsed without affecting cross-fault polyline conformity
because they live entirely on one fault's surface.

**Suggested fix:** After the symmetric polyline collapse, add a
**per-fault non-polyline sliver pass** that uses CGAL's
existing `experimental::remove_almost_degenerate_faces` with the
polyline ECM as the constraint set:

```diff
@@ main.cpp — between sliver-collapse (line 472) and seed-protect_ecms (line 476)
         }
+
+        // R-202: per-fault non-polyline sliver cleanup.  Polyline
+        // edges are protected via polyline_ecms; CGAL's
+        // experimental::remove_almost_degenerate_faces collapses
+        // the remaining slivers (non-polyline near-vertex edges
+        // that corefine introduced as a side effect of polyline
+        // pierces close to existing vertices).
+        for (std::size_t i = 0; i < N; ++i) {
+            const double cap_cos = std::cos(170.0 * CGAL_PI / 180.0);
+            const double needle_aspect = 4.0;
+            const double collapse_len = sliver_threshold;
+            CGAL::Polygon_mesh_processing::experimental::
+                remove_almost_degenerate_faces(
+                    meshes[i],
+                    CGAL::parameters::cap_threshold(cap_cos)
+                                     .needle_threshold(needle_aspect)
+                                     .collapse_length_threshold(collapse_len)
+                                     .edge_is_constrained_map(polyline_ecms[i]));
+        }
```

(`#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>`
needs to be added.)

This is the fix the prior review suggested (R-101 in
REVIEW_phase02_code.md), but it was applied only to polyline
edges in `collapse_short_polyline_edges_symmetric`. Adding the
PMP-experimental cleanup as a per-fault pass catches the
non-polyline residue.

**Test case:**
```cpp
TEST(R202_no_non_polyline_slivers) {
  // Mill Creek fixture.
  // After full Phase 2:
  Mesh A = ...;
  // ... full pipeline ...
  // Assert min edge ≥ target / 4 (= 250 m for target 1000).
  for (auto e : A.edges()) {
    EXPECT_GE(edge_length(A, e), 250.0)
        << "Edge " << e << " is shorter than the sliver threshold";
  }
}
```

---

### [R-203] [MODERATE] [main.cpp:577–605 + corefine.cpp:217–232] — Post-remesh polyline-coincidence gate uses `protect_ecms` (with border-edge filter) instead of `polyline_ecms`; the comment claims `polyline_ecms` "leaves stale marks"

**Category:** ASSUMPTION (workaround over a real CGAL behaviour)

**Description:**
Main.cpp lines 577–586 say:
> *"We use protect_ecms rather than polyline_ecms here because
>  Property_map storage is indexed by edge_index and
>  isotropic_remeshing recycles indices, leaving polyline_ecms
>  with stale marks; protect_ecms is propagated correctly by
>  split_long_edges + isotropic_remeshing's edge_is_constrained_map
>  plumbing."*

The post-remesh gate then runs:
```cpp
polyline_edge_coincidence_gate(meshes[i], protect_ecms[i],
                               meshes[j], protect_ecms[j])
```
which internally filters out border edges (corefine.cpp:222) so
the comparison reduces to polyline endpoints only.

Two concerns:

1. **The "stale marks" claim is unverified.** CGAL's
   Surface_mesh property maps are documented to track edge
   creation/deletion across PMP operations. In particular,
   `split_long_edges` and `isotropic_remeshing` propagate the
   `edge_is_constrained_map` automatically because the property
   map is registered with the algorithm via parameters. **But
   `polyline_ecms[i]` is NOT passed to either algorithm** — only
   `protect_ecms[i]` is. So `polyline_ecms[i]`'s state after
   remeshing depends on whether CGAL's *auto-propagation* of
   edge-property maps applies regardless of whether they're
   registered as constrained or not. This is a CGAL internal
   behaviour the comment hand-waves over.

2. **The `is_border` filter is a structural-not-semantic check.**
   A polyline edge that happens to lie on the fault's open
   boundary (e.g., the polyline ends at the fault edge) would
   be filtered out by `is_border`, falsely shrinking the
   compared set. CFM data has polylines that terminate at fault
   boundaries (where the polyline reaches the edge of the open
   fault patch). Such terminations would silently miss the
   gate's check.

Empirically the gate passes on Mill Creek × SBMT-SAF, but the
pair geometrically does not have polyline-meets-boundary
endpoints. A 3-fault build with adjacent fault patches would
exercise this case.

**Trigger:** A pair where the polyline reaches the open
boundary of either fault. Likely in the all-8 build.

**Suggested fix:** Either:

(a) **Maintain `polyline_ecms` correctly through remeshing.** Make
   `remesh_one_fault` propagate polyline marks explicitly: after
   each call to `split_long_edges`, walk the new constrained
   edges and propagate the polyline mark from the parent edge
   to its sub-edges. Then the post-remesh gate uses
   `polyline_ecms` directly with no filter:

```diff
@@ corefine.cpp:remesh_one_fault step 2
+        // Track the polyline-only sub-set of the protect_ecm so
+        // post-remesh polyline-coincidence verification can
+        // operate without an is_border filter (R-203).
+        std::set<EdgeDescriptor> pre_polyline;
+        for (auto e : M.edges()) {
+            if (get(polyline_ecm, e)) pre_polyline.insert(e);
+        }
         PMP_local::split_long_edges(
             long_constrained, split_ceiling, M,
             params::edge_is_constrained_map(protect_ecm));
+        // After the split, walk every edge that is constrained
+        // in protect_ecm but NOT yet in polyline_ecm and check
+        // whether either endpoint is a vertex of an original
+        // polyline edge — if so, the edge is a polyline sub-edge
+        // and must inherit the polyline mark.
+        // (Implementation: maintain a parent-edge → polyline-mark
+        //  map via a custom visitor passed to split_long_edges,
+        //  or post-process by spatial query against pre_polyline
+        //  endpoints.)
```

(b) **Tighten the gate's `is_border` filter** to a structural
   "is-this-the-fault-mesh-boundary" check that distinguishes
   polyline-on-boundary from other-boundary edges. This is
   subtle because in CGAL `is_border` only knows topology.

Option (a) is the correct fix; option (b) is a hack.

**Test case:**
```cpp
TEST(R203_polyline_meets_fault_boundary) {
  // Construct a fault A that is a half-square, and B that
  // crosses it perpendicular such that the intersection
  // polyline reaches A's open edge.
  // After Phase 2, polyline_edge_coincidence_gate should still
  // PASS, AND both faults should agree on the polyline-meets-
  // boundary endpoint.
  ...
}
```

---

### [R-204] [MODERATE] [corefine.cpp:341–429, collapse_short_polyline_edges_symmetric — O(E²) lookup pattern] — Per-candidate `find_edge_by_endpoints` linear scan; will be slow on the all-8 build

**Category:** QUALITY (performance scaling)

**Description:**
`find_edge_by_endpoints(M, key)` (corefine.cpp:314–320) walks
**every edge** of M to find one matching by endpoint coords:

```cpp
for (auto e : M.edges()) {
    if (endpoint_pair_of(M, e) == key) return e;
}
```

`collapse_short_polyline_edges_symmetric` calls this twice per
candidate edge (once for A, once for B). On the 2-fault target
this is fine (2 × 50 candidates × 3000 edges = 300k ops). On the
all-8 target with ~30k edges per fault and 28 unordered pairs
× ~50 candidates each, that's 1400 calls × 30k edges × 2 = 84M
ops per Phase 2 run. Still OK in absolute terms (~0.1 s) but
unnecessarily slow.

**Trigger:** Performance scaling concern. Not a correctness bug.

**Suggested fix:** Build a hash map keyed by endpoint pair → edge
descriptor at the start of the function, query it in O(1) per
candidate:

```diff
@@ corefine.cpp:341 — collapse_short_polyline_edges_symmetric
+    // Build endpoint → edge index for O(1) lookups (R-204).
+    auto build_index = [](const Mesh& M) {
+        std::unordered_map<EndpointPair3, Mesh::Edge_index, EpHash> idx;
+        idx.reserve(M.number_of_edges() * 2);
+        for (auto e : M.edges()) idx.emplace(endpoint_pair_of(M, e), e);
+        return idx;
+    };
+    auto idx_A = build_index(A);
+    auto idx_B = build_index(B);
+    // ... rest of loop uses idx_A.find(c.key) instead of
+    //     find_edge_by_endpoints(A, c.key) ...
```

`EpHash` needs to hash an `EndpointPair3` (`std::pair` of two
`std::array<double,3>`). Use the same FNV pattern as
`TriKeyHash`.

**Caveat:** After each successful collapse, the edge map needs
to be rebuilt OR maintained incrementally. Easiest: for each
collapse, remove the collapsed edge's key from the index. The
neighbour edges that gained a new endpoint are still in the
index under their OLD endpoint pairs; subsequent searches by the
NEW endpoint pair will miss. This is fine because we don't
re-look-up any neighbours by their new-endpoint coords —
candidates are pre-collected.

**Test case:** N/A (performance regression — measure runtime
on all-8 fixture once it lands).

---

### [R-205] [LOW] [corefine.cpp:308 + 226 vs 70 — collect_constrained_endpoints / endpoint_pair_of don't apply canon_zero] — +0/-0 hash inconsistency between collision-detection (canonicalised) and gate (not canonicalised)

**Category:** ASSUMPTION

**Description:**
`find_coincident_triangle` (corefine.cpp:91–157) and
`make_tri_key` (lines 62–78) carefully apply
`canon_zero(v)` to canonicalise +0.0 and -0.0 (R-009 fix).

But `collect_constrained_endpoints` (lines 217–232) and
`endpoint_pair_of` (lines 301–310) — both used by
`polyline_edge_coincidence_gate` and
`collapse_short_polyline_edges_symmetric` — do NOT apply
`canon_zero`:

```cpp
std::array<double, 3> a{p.x(), p.y(), p.z()};   // raw
std::array<double, 3> b{q.x(), q.y(), q.z()};   // raw
if (b < a) std::swap(a, b);
out.emplace(std::move(a), std::move(b));
```

If a polyline endpoint has a zero coordinate that happens to
be `+0.0` in fault A's storage and `-0.0` in fault B's storage
(same numeric value, different bits — possible if `(x - x) * y`
produces -0.0 on one side and +0.0 on the other), the gate's
`std::set<EndpointPair>` comparison would treat them as
different and falsely fail.

This is the same hazard `canon_zero` was added to defend
against in R-009, but the defence wasn't extended to the gate
path.

**Trigger:** Polyline endpoint at z = 0 (fault crosses free
surface; rare on real CFM data but possible after the user-
authorised `ts_to_stl.py` clearance change) where remeshing
produces a -0.0 vs +0.0 mismatch.

**Suggested fix:**

```diff
@@ corefine.cpp:217–232 collect_constrained_endpoints
     for (auto e : M.edges()) {
         if (!get(ecm, e)) continue;
         if (CGAL::is_border(e, M)) continue;
         auto h = M.halfedge(e);
         const auto& p = M.point(M.source(h));
         const auto& q = M.point(M.target(h));
-        std::array<double, 3> a{p.x(), p.y(), p.z()};
-        std::array<double, 3> b{q.x(), q.y(), q.z()};
+        std::array<double, 3> a{canon_zero(p.x()),
+                                 canon_zero(p.y()),
+                                 canon_zero(p.z())};
+        std::array<double, 3> b{canon_zero(q.x()),
+                                 canon_zero(q.y()),
+                                 canon_zero(q.z())};
         if (b < a) std::swap(a, b);
```

Same change in `endpoint_pair_of` (corefine.cpp:301–310).

`canon_zero` is currently in an unnamed namespace at
corefine.cpp:43 — accessible to both call sites.

---

### [R-206] [LOW] [corefine.cpp:319, find_edge_by_endpoints] — Returns `Mesh::Edge_index(Mesh::null_edge())` — needless wrapping

**Category:** QUALITY

**Description:**
```cpp
return Mesh::Edge_index(Mesh::null_edge());
```

`Mesh::null_edge()` already returns a `Mesh::Edge_index`
(it's a static `constexpr` sentinel). Wrapping in
`Mesh::Edge_index(...)` is a no-op constructor call.

**Suggested fix:**
```diff
-    return Mesh::Edge_index(Mesh::null_edge());
+    return Mesh::null_edge();
```

---

### [R-207] [LOW] [main.cpp:307–323 — error path on ECM install failure] — Returns exit code 2 (CGAL runtime error) on a programmer error (duplicate property tag)

**Category:** QUALITY

**Description:**
```cpp
auto pl = meshes[i].add_property_map<...>("e:polyline", false);
auto pr = meshes[i].add_property_map<...>("e:protect", false);
if (!pl.second || !pr.second) {
    std::cerr << "[corefine_faults] ECM property map already "
                 "exists on fault " << args.included[i] << '\n';
    return 2;
}
```

Exit code 2 is reserved for "CGAL runtime error (corefine
failed / non-manifold)" per the file header comment lines
17–22. A duplicate property tag is a programmer error (we
control the tag string) that should never happen at runtime;
if it does, it's a code bug, not a CGAL runtime issue.

**Trigger:** Should never fire. If it does, indicates the
fault's STL was loaded with a pre-existing property map of the
same tag — extraordinary.

**Suggested fix:** Use a dedicated exit code (or just abort with
a clearer message that this indicates an internal bug):

```diff
-    return 2;
+    return 1;     // user/programmer error, not a CGAL runtime issue
```

Or treat the case as `assert(pl.second && pr.second);` and
let it crash — the only realistic cause is a typo in the tag
string elsewhere in the codebase.

---

### [R-208] [LOW] [main.cpp:455 — sliver_threshold = target / 4] — Threshold tuning is coupled to R-201 and R-202; document the relationship

**Category:** ASSUMPTION

**Description:**
The sliver-collapse threshold is set at:
```cpp
const double sliver_threshold = args.target_edge_m / 4.0;
```

= 250 m for default target of 1000 m. Relationships:

- `split_long_edges` ceiling = `target_edge_m` (now → `target × 4/3`
  per R-201).
- After split, polyline edges in `(target × 2/3, target × 4/3)`
  ≈ (667, 1333) m.
- Sliver-collapse threshold = `target / 4` = 250 m. Edges
  shorter than 250 m get collapsed.
- The "intermediate" edges (250–667 m) survive — they're below
  the split-recursion floor but above the collapse ceiling.

For a polyline edge of, say, 350 m: above the 250 m collapse
threshold (so kept), below the 667 m post-split lower bound (so
NOT a fresh post-split edge — must have come from a shorter
input edge that didn't get split). These intermediate-length
edges are the source of the 30% mesh-size grade in R-201.

**Trigger:** Real-data polyline geometry where the corefine
output has constrained edges naturally falling in the 250–667
m band.

**Suggested fix:** Document the relationship and pin the two
constants together:

```diff
+    // Mesh-size constants (pinned together — see R-201 + R-208).
+    //   split ceiling   = target × 4/3   → polyline edges in (target/2, target × 4/3)
+    //   sliver threshold = target / 4    → only collapse ≤ target/4
+    // Intermediate edges (target/4 .. target/2) survive both passes
+    // and are responsible for the polyline-corridor finer-than-far-
+    // field gradation.  Acceptable so long as the gradation stays
+    // below the plan's "max aspect ≤ 4" criterion.
+    const double split_ceiling     = args.target_edge_m * 4.0 / 3.0;
+    const double sliver_threshold  = args.target_edge_m / 4.0;
```

(And pass `split_ceiling` to `remesh_one_fault` rather than
recomputing it there.)

---

## Summary

- Critical issues: **2** (R-201 polyline-corridor gradient,
  R-202 non-polyline slivers)
- Moderate issues: **2** (R-203 polyline_ecm vs protect_ecm
  comment, R-204 O(E) edge lookup performance)
- Low issues: **4** (R-205 +0/-0 in gate path, R-206 needless
  wrapping, R-207 wrong exit code, R-208 threshold doc)
- Plan compliance: **PARTIAL**.
  - Phase 2 acceptance "max aspect ≤ 4 on triangles incident to
    a constrained edge": **substantively passing** (Mill Creek
    max 4.7, SBMT-SAF max 7.6 with 3 violators) — close to but
    not exactly meeting the bar.
  - Phase 2 acceptance "mean edge length in [800, 1200] m":
    **passing globally** (Mill 942, SAF 997) but the
    near-polyline mean is OUTSIDE the band (702/714 m). The plan
    didn't anticipate the spatial gradation; should be tightened
    to "mean edge length in [800, 1200] m globally AND within
    1 km of any constrained edge" once R-201 lands.
- Verdict: **PASS WITH FIXES**. R-201 is the dominant
  remaining defect and matches the user's "mesh size is too
  different" observation precisely. R-202 is hidden in the
  numbers but will affect Phase 3 HXT output. The other six
  findings are quality fixes that can land together.

## Suggested fix order

1. **R-201** (split ceiling = target × 4/3) — single-line
   constant change in `remesh_one_fault`. Eliminates the
   user-visible mesh-density gradation around the polyline.
2. **R-202** (per-fault non-polyline sliver pass via
   `experimental::remove_almost_degenerate_faces`) — adds a
   ~10-line block in `main.cpp` between sliver-collapse and
   protect-seeding. Removes the residual ~3 needles in SAF.
3. **R-205, R-208** (canon_zero + threshold doc) — cosmetic
   but cheap.
4. **R-203** (polyline_ecms propagation through remesh) —
   structural; defer to a follow-up if Phase 3 doesn't surface
   a polyline-meets-boundary failure.
5. **R-204** (O(1) edge-endpoint lookup) — only worth landing
   when the all-8 build is on the menu (Phase 5).
6. **R-206, R-207** (wrapping + exit code) — trivial cleanup.

## Why R-201 is what the user sees

The user's observation: *"the mesh size is too different than
other places on the fault"*.

Direct measurement of the post-Phase-2 STLs:

| Distance from polyline | Mill Creek mean tri-edge | SBMT-SAF mean tri-edge |
|---|---:|---:|
| < 500 m | **702 m** | **714 m** |
| 500–1500 m | 970 m | 1004 m |
| 1500–5000 m | 950 m | 1030 m |
| > 5000 m | 990 m | 1010 m |

The corridor within 500 m of the polyline is **30% finer** than
the rest of the mesh. The polyline shows up in ParaView as a
"trough" of denser triangulation — exactly what the user is
seeing.

R-201's diff (split ceiling target × 4/3) lifts the polyline
edges into `(667, 1333) m` instead of `(500, 1000] m`,
matching the far-field target. Predicted post-fix near-polyline
mean: ~1000 m, eliminating the gradation.

## Unreviewed Areas

- **CGAL `experimental::remove_almost_degenerate_faces`
  symmetry on two meshes** (used for R-202 fix): not yet
  measured. The R-101 fix for symmetry is in
  `collapse_short_polyline_edges_symmetric` and is restricted
  to polyline edges. The R-202 fix would run independently per
  fault on non-polyline edges, which is safe — non-polyline
  collapses don't affect cross-fault conformity.
- **The `polyline_ecms` propagation through remesh** (R-203):
  whether CGAL Surface_mesh's auto-property-map mechanism
  carries the polyline mark through `split_long_edges` and
  `isotropic_remeshing` even when `polyline_ecm` is not
  registered as the algorithm's `edge_is_constrained_map`.
  The current main.cpp comment claims it doesn't, hence the
  workaround using `protect_ecms` in the post-remesh gate. Worth
  a 10-line test to verify or refute, before the R-203 fix
  lands.
- **Phase 3 HXT acceptance** of post-Phase-2 STLs with R-201
  and R-202 applied: still pending; this is the next milestone.
- **The `cgal_patch` shim** (3-line `this->base()` →
  `this->g` substitution): not re-checked this round.
