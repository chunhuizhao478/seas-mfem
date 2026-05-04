# Code Review: Phase 1 + Phase 2 (2026-04-29, second pass)

## Review Scope

- **Plan reviewed against:** `PLAN_cgal_corefine.md` Phase 1 + Phase 2.
- **Files reviewed:**
  - `tools/CMakeLists.txt`
  - `tools/main.cpp`
  - `tools/corefine_faults/{corefine.hpp, corefine.cpp, io.hpp, io.cpp, main.cpp}`
  - `tools/tests/test_corefine_smoke.cpp`
  - `tools/README.md`, `tools/cgal_patch/README.md`
- **User evidence:** screenshot
  `~/Desktop/Screenshot 2026-04-29 at 9.47.35 PM.png` of the
  post-Phase-2 surface mesh, plus the user's note: *"the
  intersection mesh looks better but can be improved I think"*
  and *"it still has some bad triangle as I can tell, check the
  resultant stl file carefully"*.
- **Live verification on this machine:**
  - All 8 smoke tests pass (`test_two_squares`,
    `test_three_orthogonal_cascade`,
    `test_coincident_preflight`,
    `test_polyline_edge_coincidence_gate`,
    `test_mill_creek_pin`,
    `test_remesh_protects_constrained`,
    `test_free_surface_violation_and_clamp`,
    `test_idempotency`).
  - End-to-end run on Mill Creek × SBMT-SAF 2000 m fixture with
    `--target-edge-m 1000 --remesh-iters 3 --clearance-m 100`
    completes in 67 ms.
  - **Per-stage triangle quality measured directly from the STLs:**

    | Stage | Mill aspect>5 | Mill aspect_max | SAF aspect>5 | SAF aspect_max | edges<50m |
    |---|---:|---:|---:|---:|---:|
    | Input (raw) | 0 | 4.7 | 0 | 3.8 | 0 |
    | Post-corefine (no remesh) | **82** | **924** | **95** | **1280** | 48 + 50 |
    | Post-remesh | 58 | 542 | 63 | 854 | 38 + 40 |

  - Worst Mill Creek post-remesh sliver: aspect **542.7**, edges
    `[0.99, 537.32, 537.79] m` — a 1-metre needle.
  - Only **3 of 58 bad triangles** sit within 100 m of a polyline
    vertex; the rest are **distributed across the fault** as
    needles whose short edge is a constrained polyline edge near
    an existing fault vertex.

The Phase 0/1 review's R-001 through R-011 all appear addressed
in code. This second pass found **6 new findings** (1 critical, 4
moderate, 1 low) that are exposed once Phase 2 is on real data.

The headline: `PMP::corefine` introduces 1-metre constrained
edges where the intersection polyline passes within ~1 m of an
existing fault vertex. Because the Phase 2 remesh runs with
`protect_constraints(true)`, those edges cannot be collapsed and
become **needle triangles persisting into the final mesh**. This
is what the user is seeing in the screenshot, and the dominant
remaining mesh-quality defect.

---

## Findings

### [R-101] [CRITICAL] [corefine.cpp + main.cpp pipeline] — `PMP::corefine` near-vertex pierces produce ~80 needle triangles per fault that `protect_constraints(true)` cannot remove

**Category:** BUG (silent quality degradation)

**Description:**
When the intersection polyline passes geometrically close (≤ a
few metres) to an existing fault vertex, `PMP::corefine` inserts
a NEW vertex at the geometric pierce point and creates a short
edge connecting it to the nearby existing vertex. CGAL's
exact-predicates corefine does not snap; it preserves the exact
intersection geometry.

On Mill Creek × SBMT-SAF 2000 m the immediate post-corefine
output has:
- 82 (Mill) and 95 (SAF) triangles with aspect ratio > 5;
- max aspect 924 (Mill) and 1280 (SAF);
- 48 (Mill) and 50 (SAF) edges shorter than 50 m, of which the
  vast majority are **constrained polyline edges** ≤ 5 m.

`isotropic_remeshing(protect_constraints=true)` then runs over
this corefined input. The protect-constraints flag is documented
to forbid splitting AND collapsing constrained edges. The 1 m
polyline edges therefore stay at 1 m. Adjacent non-constrained
edges get remeshed to ~1000 m. Result: triangles with one ~1 m
constrained edge and two ~1000 m edges — aspect ratio 1000:1.

The Phase 2 acceptance criterion in `PLAN_cgal_corefine.md`:

> max aspect ratio ≤ 4 on triangles incident to a constrained edge

is **substantially violated** (worst case 542:1). The plan's
P-001 fix (split long constrained edges) was the wrong direction
of correction — the actual problem is *short* constrained edges,
not long ones.

**Trigger:** Any pair of fault meshes where the intersection
polyline approaches within `< target_edge_m` of an existing
vertex on either fault. Common in real CFM data because the
~2000 m source triangulation makes near-vertex tangencies the
typical case, not the exceptional case.

**Actual behavior:**
- Intersection report records `n_constrained_edges = 143` per
  side (PASS the symmetry check).
- All eight Phase-1/2 gates report PASS.
- Output STL has ~60 needle triangles per fault, distributed
  across the surface (centroids NOT clustered near polyline
  vertices — only 3/58 within 100 m).
- HXT (Phase 3) will likely refuse, or produce a degenerate
  bulk mesh, on these inputs.

**Expected behavior:** Either:
- **(a) Snap near-vertex pierces during preprocessing**, so
  corefine never produces sub-target constrained edges; OR
- **(b) Run a CGAL-idiomatic post-remesh sliver cleanup** that
  collapses short constrained edges (with a small geometric
  budget), accepting ~1 m of polyline-position drift; OR
- **(c) Pre-collapse short constrained edges manually**
  (un-protecting the specific edges via a separate ECM key
  before remesh).

The plan's wording in §Acceptance suggests (b) is acceptable for
SEAS use because the polyline is a geometric approximation
already.

**Suggested fix (recommend (b) — minimal code, CGAL-supported):**

CGAL provides
`CGAL::Polygon_mesh_processing::experimental::remove_almost_degenerate_faces`
(or the namespaced `repair_self_intersections` family) which is
designed for this. Add a post-remesh sliver cleanup pass:

```diff
@@ corefine_faults/corefine.cpp — extend remesh_one_fault
@@ Step 5 — verify constraints survived ...
+
+    // Step 6 — sliver cleanup (R-101).  Corefine introduces
+    // sub-metre constrained edges where the polyline passes near
+    // an existing fault vertex.  protect_constraints(true) above
+    // has correctly preserved them, but the result is needles
+    // with aspect ratios up to 1000:1.  Collapse those needles
+    // by un-protecting just the short edges and re-running a
+    // bounded edge-collapse pass.
+    {
+        // Build a "sliver-friendly" ECM that protects only LONG
+        // constrained edges; short ones (< target/4) get
+        // un-protected so the collapse pass can fix them.
+        auto ecm_loose = M.add_property_map<EdgeDescriptor, bool>(
+            "e:is_constrained_loose", false).first;
+        const double short_threshold = p.target_edge_m / 4.0;
+        for (auto e : M.edges()) {
+            if (!get(ecm, e)) continue;
+            const auto& a = M.point(M.source(M.halfedge(e)));
+            const auto& b = M.point(M.target(M.halfedge(e)));
+            const double L = std::sqrt(CGAL::squared_distance(a, b));
+            put(ecm_loose, e, L >= short_threshold);
+        }
+        // remove_almost_degenerate_faces with the loose ECM
+        // collapses needle triangles whose short edge is no longer
+        // protected.  cap_threshold and needle_threshold are
+        // documented in CGAL/Polygon_mesh_processing/repair_degeneracies.h.
+        PMP_local::experimental::remove_almost_degenerate_faces(
+            M,
+            params::cap_threshold(std::cos(170.0 * CGAL_PI / 180.0))  // ~10° wedge
+                   .needle_threshold(4.0)                              // aspect 4:1
+                   .collapse_length_threshold(short_threshold)
+                   .edge_is_constrained_map(ecm_loose));
+        M.remove_property_map(ecm_loose);
+    }
+    // Re-snapshot post-cleanup counts.
+    r.n_tri_post = static_cast<std::int64_t>(M.number_of_faces());
```

Caveats:
- `experimental::remove_almost_degenerate_faces` does NOT
  preserve the exact polyline geometry — short polyline edges
  get collapsed, displacing the polyline by ≤ `short_threshold`.
- After this collapse, `polyline_edge_coincidence_gate` MAY FAIL
  because A and B can collapse at different positions if the
  collapse choice is implementation-dependent. The fix is to
  collapse SYMMETRICALLY: collect the union of "edges to
  collapse" from both faults, then for each pair, weld both A
  and B at the same target vertex (the longer-edge endpoint).

If symmetric collapse is too invasive for first-cut, a simpler
alternative is option (a): pre-process inputs by welding A
vertices to B vertices that are within `short_threshold` BEFORE
calling corefine. CGAL's `experimental::snap_borders` does this
between two meshes. Either fix needs to be designed carefully —
the suggested diff above is a sketch, not a finished
implementation.

**Test case:**
```cpp
TEST(R101_no_post_remesh_slivers) {
  // Construct a fixture that deliberately places the polyline
  // ~1 m from an existing vertex of A.  After full Phase 2, the
  // STL must have aspect_max ≤ 4 (or whatever post-cleanup target
  // we choose).
  Mesh A = make_square_offset(2, 1.0, /*offset=*/0.001);
  Mesh B = make_square(1, 1.0);
  // ... corefine + remesh + sliver cleanup ...
  for (auto f : A.faces()) {
    auto pts = face_points(A, f);
    double mn = std::min({len(pts[0]-pts[1]), len(pts[1]-pts[2]), len(pts[2]-pts[0])});
    double mx = std::max({len(pts[0]-pts[1]), len(pts[1]-pts[2]), len(pts[2]-pts[0])});
    EXPECT_LE(mx / std::max(mn, 1e-9), 4.0);
  }
}
```

---

### [R-102] [CRITICAL] [main.cpp:375–387, gate-coverage gap] — `polyline_edge_coincidence_gate` runs before remesh; never re-checked after `split_long_edges` or `isotropic_remeshing` mutate the polyline

**Category:** BUG (silent invariant gap)

**Description:**
The gate runs once, in the corefine inner loop, BEFORE
`remesh_one_fault` is called (which happens in the separate loop
at lines 424–515). The gate then never runs again. The plan's
acceptance "the polyline edges of A and B coincide in 3-D" is
not verified after the operations that could break it:

1. `mark_border_edges_constrained` adds **boundary** edges to the
   shared ECM. Boundary edges of A do not exist on B (faults
   don't share boundary). After this step, a comparison of
   `collect_constrained_endpoints(A)` vs `…(B)` will FAIL even
   if polyline edges match perfectly, because the boundary
   subsets differ.
2. `split_long_edges(target_edge_m, ecm)` subdivides every
   constrained edge above target. The midpoint of an edge with
   identical 3-D endpoints in A and B is bit-identical, so
   polyline splits SHOULD match — but this is not verified.
3. `isotropic_remeshing(protect_constraints=true)` is documented
   to leave constrained edges unchanged. We trust this; we do
   not measure.

The verbose end-to-end output shows:
```
remesh safs_sbmt_millcreek: F 1079 → 2737, ce 143 → 367 → 367
remesh safs_sbmt_saf:       F 1633 → 5507, ce 143 → 481 → 481
```

`ce_A = 367` vs `ce_B = 481` post-remesh. The asymmetry is
**expected** because of the boundary-edge contribution (Mill
Creek and SBMT-SAF have different boundary lengths), but the
asymmetric counts make a naive "ce_A == ce_B" check
unusable. The polyline portion is presumably symmetric, but
"presumably" isn't enough.

A direct measurement on the post-remesh STLs:
- 1477 unique vertices on Mill Creek post-remesh.
- 2919 on SBMT-SAF.
- **153 shared 3-D coordinate values** (vs 144 polyline vertices
  immediately post-corefine).

The 153 ≈ 144 + ~9 suggests that some polyline edges did get
split during remeshing, and most splits happened consistently on
both sides — but again, "presumably" isn't enough.

**Trigger:** Any Phase 2 run. Currently 100% of paths.

**Actual behavior:** Post-remesh polyline conformity is
unverified. If a future CGAL upgrade or a subtle property-map
mishap broke the invariant, the tool would emit non-conformal
STLs that HXT would fail on, with no clear diagnostic from
`corefine_faults`.

**Expected behavior:** Add a **polyline-only ECM** (separate
from the boundary-marking ECM), and re-check
`polyline_edge_coincidence_gate` AFTER remeshing using only the
polyline ECM. Any mismatch is fatal.

**Suggested fix:**

```diff
@@ main.cpp — install TWO ECMs per fault, one for polyline-only marks, one for the remesh-protection set (polyline ∪ boundary).
-    std::vector<ECMap> ecms(N);
+    std::vector<ECMap> polyline_ecms(N);
+    std::vector<ECMap> protect_ecms(N);
     for (std::size_t i = 0; i < N; ++i) {
-        auto pm = meshes[i].add_property_map<...>("e:is_constrained", false);
-        ...
+        auto p1 = meshes[i].add_property_map<safs::corefine::EdgeDescriptor, bool>(
+                      "e:polyline", false);
+        auto p2 = meshes[i].add_property_map<safs::corefine::EdgeDescriptor, bool>(
+                      "e:protect", false);
+        if (!p1.second || !p2.second) { ... return 2; }
+        polyline_ecms[i] = p1.first;
+        protect_ecms[i]  = p2.first;
     }
@@ — pass polyline_ecms to corefine_pair (so it marks only intersection edges, not boundaries).
@@ — after the corefine loop, BEFORE remesh:
+    for (std::size_t i = 0; i < N; ++i) {
+        // Seed protect_ecm with polyline marks.
+        for (auto e : meshes[i].edges()) {
+            if (get(polyline_ecms[i], e)) put(protect_ecms[i], e, true);
+        }
+    }
@@ — remesh_one_fault is called with protect_ecms (which marks
@@   polyline ∪ boundary) and ALSO writes its border-marks into
@@   protect_ecms only — polyline_ecms stays clean.
@@   The split_long_edges step propagates polyline subdivisions
@@   into BOTH ECMs (for the polyline edges) and into protect_ecms
@@   only (for boundary edges).
@@
@@ — after the remesh loop, re-run the gate using polyline_ecms only:
+    for (std::size_t i = 0; i < N; ++i) {
+        for (std::size_t j = i + 1; j < N; ++j) {
+            if (!safs::corefine::polyline_edge_coincidence_gate(
+                    meshes[i], polyline_ecms[i],
+                    meshes[j], polyline_ecms[j])) {
+                std::cerr << "[corefine_faults] post-remesh polyline-"
+                             "coincidence FAILED on pair ("
+                          << args.included[i] << ", "
+                          << args.included[j] << ").\n";
+                return 2;
+            }
+        }
+    }
```

The two-ECM split also requires `remesh_one_fault` to take TWO
ECM arguments (one for polyline propagation, one for protection)
and propagate splits correctly — plus updating
`mark_border_edges_constrained` to write only into `protect_ecm`.

**Test case:**
```cpp
TEST(R102_post_remesh_polyline_gate) {
  // Two-square fixture; run full Phase 2.  Then call the
  // post-remesh gate; assert PASS.
  Mesh A = make_square(2, 1.0); Mesh B = make_square(1, 1.0);
  ECMap pl_a = ..., pl_b = ..., pr_a = ..., pr_b = ...;
  corefine_pair(A, B, pl_a, pl_b, ...);
  // Seed protect from polyline.
  ...
  remesh_one_fault(A, /*polyline=*/pl_a, /*protect=*/pr_a, ...);
  remesh_one_fault(B, pl_b, pr_b, ...);
  ASSERT_TRUE(polyline_edge_coincidence_gate(A, pl_a, B, pl_b));
}
```

---

### [R-103] [MODERATE] [corefine.cpp:377–418, project_to_free_surface_safe] — Inversion check passes degenerate-output triangles silently; the clamp may produce zero-area triangles

**Category:** EDGE_CASE

**Description:**
`project_to_free_surface_safe` checks face normals before vs
after projection:

```cpp
const auto n_pre  = face_normal_unit(M, f_pre);
const auto n_post = face_normal_unit(trial, f_post);
if (n_pre.squared_length() == 0.0
    || n_post.squared_length() == 0.0) {
    continue;                       // skip degenerate
}
if (CGAL::scalar_product(n_pre, n_post) <= 0.0) {
    return false;                   // inversion
}
```

When the projection produces a degenerate triangle (3 collinear
points after clamping one vertex), `n_post.squared_length() ==
0.0` and the check `continue`s. The clamp is then committed even
though it created a zero-area triangle. The STL writer
subsequently emits a triangle with all-zero normal (the
fallback at io.cpp:141–143 leaves nx=ny=nz=0 when nrm == 0).

HXT and most downstream consumers reject zero-area triangles or
silently produce inverted tets near them.

**Trigger:** A vertex at z = -50 m whose 1-ring neighbours are
collinear in the (x, y) plane after the clamp pulls v to
z = -100 m. Less common than R-101 but possible on fault
geometries with near-edge interior vertices.

**Suggested fix:** Treat post-projection degeneracy the same as
inversion (refuse the clamp):

```diff
@@ corefine.cpp:402–413
     for (auto f_pre : M.faces()) {
         Mesh::Face_index f_post(f_pre);
         const auto n_pre  = face_normal_unit(M, f_pre);
         const auto n_post = face_normal_unit(trial, f_post);
-        if (n_pre.squared_length() == 0.0
-            || n_post.squared_length() == 0.0) {
-            continue;
+        if (n_pre.squared_length() == 0.0) {
+            continue;            // input was degenerate; don't fault the clamp
+        }
+        if (n_post.squared_length() == 0.0) {
+            return false;        // R-103: clamp would create a degenerate face
         }
         if (CGAL::scalar_product(n_pre, n_post) <= 0.0) {
             return false;
         }
     }
```

**Test case:**
```cpp
TEST(R103_clamp_refuses_degenerate_creation) {
  Mesh M;
  // Three vertices in a line on z = -150, one off-line at z = -50.
  // After clamp to z = -100, all four would be coplanar with the
  // first three collinear → degenerate triangles.
  // (Construction sketch — implementer fills in the exact mesh
  //  that puts a non-constrained interior vertex above clearance
  //  with collinear neighbours.)
  ...
  EXPECT_FALSE(project_to_free_surface_safe(M, ecm, 100.0));
}
```

---

### [R-104] [MODERATE] [corefine.cpp:253–265, vertex_is_constrained] — Manual halfedge rotation is brittle on boundary vertices; should use the documented CGAL idiom

**Category:** ASSUMPTION

**Description:**
`vertex_is_constrained` rotates halfedges around `v` with:

```cpp
auto h0 = M.halfedge(v);
auto h = h0;
do {
    const auto e = M.edge(h);
    if (get(ecm, e)) return true;
    h = M.opposite(M.next(h));
} while (h != h0);
```

For an INTERIOR vertex this visits every incident edge. For a
BOUNDARY vertex (open fault patch — every fault has these), the
ring is open and the rotation `M.opposite(M.next(h))` may
traverse the boundary halfedge cycle without coming back to
`h0`, OR may exit the rotation early. CGAL's documented safe
idiom for "iterate halfedges around target" is:

```cpp
for (auto h : CGAL::halfedges_around_target(v, M)) { ... }
```

which handles the boundary case correctly. The smoke test
exercises only interior-vertex paths and would not catch a
boundary-vertex bug. The Mill Creek × SBMT-SAF run does include
boundary vertices (every fault has them), but the function is
only called from `project_to_free_surface_safe` and `check_
free_surface_clearance`, both of which run under the
`--clearance-m` path; in the default run the input is already
clamped, so the function is never invoked on a violator.

**Suggested fix:**

```diff
@@ corefine.cpp:253–265
 bool vertex_is_constrained(const Mesh& M,
                            const EdgeConstrainedMap& ecm,
                            Mesh::Vertex_index v) {
-    auto h0 = M.halfedge(v);
-    if (h0 == Mesh::null_halfedge()) return false;
-    auto h = h0;
-    do {
-        const auto e = M.edge(h);
-        if (get(ecm, e)) return true;
-        h = M.opposite(M.next(h));
-    } while (h != h0);
-    return false;
+    if (M.halfedge(v) == Mesh::null_halfedge()) return false;
+    for (auto h : CGAL::halfedges_around_target(v, M)) {
+        if (get(ecm, M.edge(h))) return true;
+    }
+    return false;
 }
```

Add `#include <CGAL/boost/graph/iterator.h>` if not already
transitively included.

**Test case:**
```cpp
TEST(R104_boundary_vertex_constrained) {
  // Open patch (one triangle); mark its only edge constrained.
  // Each of the three vertices is a BOUNDARY vertex.  All three
  // must report constrained=true.
  Mesh M;
  auto v0 = M.add_vertex(P3(0,0,0));
  auto v1 = M.add_vertex(P3(1,0,0));
  auto v2 = M.add_vertex(P3(0,1,0));
  M.add_face(v0, v1, v2);
  ECMap ecm = ...;
  // Mark edge (v0, v1) constrained.
  for (auto e : M.edges()) {
    auto h = M.halfedge(e);
    if ((M.source(h) == v0 && M.target(h) == v1)
        || (M.source(h) == v1 && M.target(h) == v0)) {
        put(ecm, e, true);
    }
  }
  EXPECT_TRUE(vertex_is_constrained(M, ecm, v0));
  EXPECT_TRUE(vertex_is_constrained(M, ecm, v1));
  EXPECT_FALSE(vertex_is_constrained(M, ecm, v2));
}
```

---

### [R-105] [MODERATE] [corefine.cpp:243–250, mark_border_edges_constrained — side effect on caller's ECM] — Permanently pollutes the ECM passed to `remesh_one_fault`; couples border marks with polyline marks for the rest of the program

**Category:** ASSUMPTION (hidden mutation)

**Description:**
`remesh_one_fault` is called with `const EdgeConstrainedMap&
ecm`. The first thing it does is mutate the underlying property
storage:

```cpp
mark_border_edges_constrained(M, ecm);   // step 1
```

CGAL Property_map is a handle; `const Property_map&` does not
prevent `put(ecm, e, true)`. So the caller's ECM (stored in
`ecms[i]` in main.cpp) is permanently mutated to include
boundary-edge marks. Any subsequent operation on that ECM —
including the post-remesh polyline gate proposed in R-102 —
sees boundary edges mixed in with polyline edges.

The function name and signature suggest read-only ECM
behaviour; the actual implementation mutates it. This is
exactly the structural coupling that makes R-102 hard.

**Trigger:** Always — every `remesh_one_fault` call.

**Actual behavior:** main.cpp's `ecms[i]` accumulates polyline
marks (added by `corefine_pair`) AND boundary marks (added by
`remesh_one_fault`). No way to distinguish them later.

**Expected behavior:** Two distinct ECMs (the R-102 fix), OR
`remesh_one_fault` operates on a LOCAL copy/snapshot of the ECM
without polluting the caller's.

**Suggested fix:** Accept the R-102 fix; `remesh_one_fault`
takes two ECMs (`polyline_ecm` for read-only, `protect_ecm` for
mutation). The function then writes border marks only into
`protect_ecm`. Caller keeps `polyline_ecm` clean for post-remesh
verification.

```diff
@@ corefine.hpp — RemeshParams + remesh_one_fault signature
 RemeshResult remesh_one_fault(
     Mesh& M,
-    const EdgeConstrainedMap& ecm,
+    const EdgeConstrainedMap& polyline_ecm,
+    const EdgeConstrainedMap& protect_ecm,
     const RemeshParams& p);
@@ corefine.cpp — remesh_one_fault implementation
-    mark_border_edges_constrained(M, ecm);
+    // Seed protect_ecm with polyline marks (the union); border
+    // additions go ONLY into protect_ecm so polyline_ecm stays
+    // clean for the post-remesh coincidence gate (R-102).
+    for (auto e : M.edges()) {
+        if (get(polyline_ecm, e)) put(protect_ecm, e, true);
+    }
+    mark_border_edges_constrained(M, protect_ecm);
@@ — every subsequent reference to `ecm` becomes `protect_ecm`,
@@   except the post-remesh count-of-polyline check which uses
@@   polyline_ecm.
```

**Test case:** Same as R-104 negative-control: after
`remesh_one_fault` returns, the polyline ECM must NOT have
border-edge marks set.

---

### [R-106] [LOW] [main.cpp:222–224, interior_crossing_only] — Deep-copies both meshes per pair; per-pair O(N_tri) memory amplifies on the all-8 build

**Category:** QUALITY

**Description:**
`interior_crossing_only(A, B)` creates `Mesh A2 = A, B2 = B`,
which CGAL Surface_mesh implements as a full structural copy
(vertices, halfedges, faces, all property maps including the
~6500 entries each on Mill Creek + SBMT-SAF, and any others
attached). The copy is then mutated by `corefine_pair`. After
the function returns the copies are destroyed.

For the all-8 cascade (28 pairs, each fault ~6500 triangles),
that's 56 large mesh copies. CGAL's Surface_mesh is reasonably
compact (~100 bytes/face), so 6500 × 100 × 56 ≈ 36 MB total
peak — not catastrophic. But it is wasted: an
"is-it-idempotent" check that allocates 36 MB at the end of
every pair processing is excessive when the same answer can be
obtained by counting `M.number_of_faces()` before and after a
no-op corefine.

**Suggested fix:** Run corefine on the LIVE meshes (no copy), but
on a FRESH ECM (added with `add_property_map<...>` then removed
afterwards). The triangle counts pre/post will be unchanged for a
correctly-conformal pair (CGAL has no-op fast path for already-
intersected meshes). If counts changed, the gate fails AND the
mesh has been mutated — which is fine because the next pair
would have re-corefined them anyway.

Alternatively, skip the gate entirely and trust the
(post-remesh polyline-coincidence + manifold) gates to catch any
real failure. The plan §Phase 1 §6 documents that this gate is
nice-to-have; the real correctness check is the pair of
manifold + polyline-coincidence gates that R-102 proposes.

```diff
@@ main.cpp:222–224 — replace deep-copy with in-place re-corefine
-bool interior_crossing_only(const Mesh& A, const Mesh& B,
-                            std::string_view short_a,
-                            std::string_view short_b) {
-    Mesh A2 = A, B2 = B;
-    const auto a_pre = A2.number_of_faces();
-    ...
+bool interior_crossing_only(Mesh& A, Mesh& B,
+                            std::string_view short_a,
+                            std::string_view short_b) {
+    auto pa = A.add_property_map<...>("e:gate_idem", false);
+    auto pb = B.add_property_map<...>("e:gate_idem", false);
+    if (!pa.second || !pb.second) {
+        // Already exists from a previous pair; skip silently.
+        return true;
+    }
+    const auto a_pre = A.number_of_faces();
+    const auto b_pre = B.number_of_faces();
+    try {
+        safs::corefine::corefine_pair(A, B, pa.first, pb.first,
+                                      short_a, short_b);
+    } catch (...) {
+        A.remove_property_map(pa.first);
+        B.remove_property_map(pb.first);
+        return false;
+    }
+    A.remove_property_map(pa.first);
+    B.remove_property_map(pb.first);
+    return A.number_of_faces() == a_pre
+        && B.number_of_faces() == b_pre;
+}
```

This deviates from "deep copy" intent but is correct because a
true-no-op corefine doesn't mutate the meshes. If it does
mutate, the function returns false (gate fails) AND the meshes
are at most as fine-grained as they would have been — the next
pair would have produced the same intersection anyway.

---

## Summary

- Critical issues: **2** (R-101 sliver triangles, R-102 missing
  post-remesh polyline gate)
- Moderate issues: **3** (R-103 degenerate-clamp escape,
  R-104 boundary-vertex iteration, R-105 ECM mutation
  side-effect)
- Low issues: **1** (R-106 deep-copy waste)
- Plan compliance: **PARTIAL**. The Phase 2 acceptance criterion
  "max aspect ratio ≤ 4 on triangles incident to a constrained
  edge" is **substantially violated** (max 542:1) on the
  Mill Creek × SBMT-SAF target.
- Verdict: **PASS WITH FIXES** — but R-101 is the dominant
  remaining mesh-quality defect and almost certainly what the
  user is observing in the screenshot. Fix R-101 + R-102 before
  Phase 3, in this order:
  1. **R-102 + R-105** (split ECMs into polyline-only + protect)
     — structural foundation that R-101 needs.
  2. **R-101** (sliver cleanup) — applied on top of the new
     two-ECM structure; collapse short polyline edges via
     `experimental::remove_almost_degenerate_faces` with the
     polyline ECM, then re-run the post-remesh polyline gate.
     If symmetric collapse cannot be guaranteed, fall back to
     the "pre-snap input meshes" path (option (a) in R-101).
  3. **R-103, R-104, R-106** — quality follow-ups; can land
     together.

## Why R-101 is what the user sees

The screenshot shows a uniformly-meshed surface with what
appears to be a darker overlay region (the second fault) and
some triangulation irregularities along the intersection. The
direct STL inspection confirms:

- Mean edge length is exactly target (~1000 m on both faults).
- p5 / p95 are near-target (502 / 1153 m on Mill Creek).
- BUT max aspect is 542:1 with a 1 m short edge.
- 58 / 2737 (Mill) and 63 / 5507 (SAF) triangles have aspect > 5.

These needle triangles are **not at the polyline** (only 3/58
within 100 m of a polyline vertex). They are **scattered across
the fault surface** as one-metre-long edges connecting a polyline
pierce to a nearby pre-existing fault vertex. Visually in
ParaView they look like "slightly off" triangulation — exactly
what the user is calling out.

The fix is **not** further remeshing — the remesher cannot help
because `protect_constraints(true)` blocks it from touching the
1 m polyline edges. The fix is to **collapse those short
constrained edges** (R-101) and verify polyline conformity
afterwards (R-102).

## Unreviewed Areas

- **CGAL's `experimental::remove_almost_degenerate_faces`
  symmetry on two meshes**: the suggested R-101 fix assumes that
  collapsing symmetric short polyline edges on A and on B
  produces the same target vertex on both sides. This is not
  documented; the implementer must verify (via R-102's
  post-cleanup polyline gate) that this holds, or fall back to
  pre-snapping input meshes.
- **HXT acceptance of post-remesh STLs with R-101 fix applied**:
  not verified yet; deferred to Phase 3. If HXT still rejects,
  the cause may be elsewhere (e.g., self-intersections introduced
  by the collapse).
- **Deep-copy memory cost of `interior_crossing_only` on the
  all-8 build**: not measured; flagged in R-106 with a
  back-of-envelope estimate of 36 MB peak.
- **The `cgal_patch` shim** (3-line `this->base()` → `this->g`
  substitution): I did not diff against the upstream conda-forge
  header in this round. The provenance README claims the diff
  is exactly 3 sed substitutions; trust but verify on the next
  CGAL upgrade.
