# Code Review: Phase 1 + Phase 2 (2026-04-29, fourth pass — post-R-201..R-208)

## Review Scope

- **Plan:** `PLAN_cgal_corefine.md` Phase 1 + Phase 2.
- **Files reviewed (all current state):**
  - `tools/CMakeLists.txt`, `tools/main.cpp`, `tools/README.md`
  - `tools/corefine_faults/{corefine.hpp, corefine.cpp, io.hpp, io.cpp, main.cpp}`
  - `tools/tests/test_corefine_smoke.cpp`
- **User context:** *"the improvements are small as suggested by
  agent, can we improve it further"* — referring to the R-201
  fix (split ceiling target × 4/3) yielding only marginal
  reduction in the polyline-corridor mesh-density gradient.
- **Live verification on this machine:**
  - All 11 smoke tests pass.
  - End-to-end on Mill Creek × SBMT-SAF 2000 m (target 1000 m,
    remesh-iters 3, clearance 100 m): 69 ms total. 46 short
    polyline edges collapsed at 250 m threshold.
  - **Direct STL measurement vs prior round (v2):**

    | Metric | v2 (target ceiling) | v3 (4/3 × target ceiling) |
    |---|---:|---:|
    | Mill near-polyline (< 500 m) mean tri-edge | 702 m | 723 m |
    | Mill far-field (> 5 km) mean tri-edge | 990 m | 998 m |
    | Mill near/far ratio | 0.71 | **0.72** |
    | SAF near-polyline mean tri-edge | 714 m | 755 m |
    | SAF far-field mean tri-edge | 1010 m | 1014 m |
    | SAF near/far ratio | 0.71 | **0.74** |
    | Mill aspect_max | 4.7 | 4.93 |
    | SAF aspect_max | **7.6** | **15.19** |
    | SAF triangles aspect > 5 | 3 | 3 |
    | SAF min edge | 66 m | 66 m |

  - **The R-201 fix gave a 1–3% improvement in the gradient AND
    doubled the SAF aspect_max from 7.6 to 15.19**.
  - **R-202 (`cleanup_non_polyline_slivers`) was implemented as
    code but DISABLED in main.cpp** because empirically the
    per-fault cleanup breaks cross-fault polyline coincidence on
    real CFM data (10 vs 9 triangles removed asymmetrically per
    side).

The user's intuition is correct: the improvement was marginal.
The deeper diagnosis (next two findings) explains *why* and
proposes the structural fix.

I assumed the code contains at least 3 bugs. I found 4
(2 critical, 1 moderate, 1 low). The headline:

- The R-201 fix attacks only the upper tail of the polyline-edge
  distribution; the actual mean is 615 m (target 1000 m), driven
  by the LOWER tail.
- R-201 made max aspect ratio worse on SBMT-SAF (7.6 → 15.19)
  because the larger split ceiling lets surrounding edges grow,
  magnifying the aspect ratio of any short polyline edge that
  escapes collapse.
- Real fix: replace `split_long_edges + collapse_short` with
  **uniform polyline resampling** (place new vertices at exactly
  `target_edge_m` arc-length spacing along each polyline chain).

---

## Findings

### [R-301] [CRITICAL] [corefine.cpp:594–604, remesh_one_fault step 2 + main.cpp:465 sliver_threshold] — `split_long_edges + collapse_short_polyline_edges` only addresses the long and very-short tails; the bulk of polyline edges (250–667 m) survives untouched, producing a mean of 615 m for a 1000 m target

**Category:** BUG (root cause of "improvements are small")

**Description:**
Direct measurement of the polyline-edge length distribution on
the Mill Creek × SBMT-SAF post-Phase-2 STL output, target = 1000 m:

```
polyline edges: 99 (matched both sides)
length stats: min=176 p5=272 p25=447 p50=590 mean=615 p75=725 p95=1158 max=1315 m

distribution:
  [0,    250) m:  3 ( 3.0%)   — should have been collapsed; link-condition saved them
  [250,  500) m: 31 (31.3%)   — DEAD ZONE: above collapse, below split-recursion floor
  [500,  667) m: 30 (30.3%)   — DEAD ZONE
  [667, 1000) m: 25 (25.3%)   — split-recursion lower half
  [1000,1333) m: 10 (10.1%)   — split-recursion upper half
  [1333, ...) m:  0 ( 0.0%)   — split-recursion upper bound
```

The R-201 fix (split ceiling = target × 4/3) lifted the
post-recursion bound from 1000 m to 1333 m. **It did nothing for
the 60% of polyline edges that fall in the 250–667 m "dead
zone"** — these edges are below the split threshold (1333 m)
and above the collapse threshold (250 m), so neither pass
touches them.

The actual mean of 615 m → corridor mean tri-edge of ~720 m →
30% finer than far-field. Almost identical to the pre-fix value
of 702 m. The marginal improvement the user observed.

Root cause: `PMP::split_long_edges(L_max)` is a **one-sided**
operation. It only addresses edges ABOVE `L_max`. Edges below
the target stay below. To get a polyline distribution centred at
target, you need a **resampling** operation that handles BOTH
tails: short-edge merging AND long-edge subdivision, with the
same target.

The closest CGAL primitive is `PMP::isotropic_remeshing` itself,
but that requires un-protecting the polyline (which would let
the remesher move polyline vertices off the original geometry,
breaking cross-fault coincidence).

**Trigger:** Real CFM data, where the input polyline-edge length
distribution is wide (175 m to 1500 m) and not centred at the
target.

**Actual behavior:** Mean polyline edge = 0.62 × target. Mean
tri-edge in the 500-m corridor around the polyline = 0.72 ×
target. Visible in screenshots as a "darker band" of denser
mesh along the intersection.

**Expected behavior:** Mean polyline edge should be within ±10%
of target. Mean tri-edge in the corridor should match far-field
within ±5%.

**Suggested fix:** Replace the
`split_long_edges + collapse_short_polyline_edges_symmetric` pair
with a custom **uniform polyline resampling** that:

1. Identifies connected polyline chains (graph-walk over polyline
   edges using the polyline ECM).
2. For each chain, computes the cumulative arc length.
3. Places new vertices at evenly-spaced arc-length values
   `target_edge_m × (1, 2, 3, ...)`, by linear interpolation
   between original polyline vertices.
4. Applies the same set of new vertex positions on BOTH faults
   simultaneously (driving cross-fault coincidence by
   construction — the new vertex coordinates are computed from
   the polyline geometry, identical between A and B).
5. Replaces the old polyline edges with the new chain.

Sketch — add to `corefine.cpp`:

```cpp
namespace safs::corefine {

// Replace the chain of constrained edges between two polyline
// "junction" vertices with a uniform resampling at target_edge_m
// arc-length spacing.  Applied symmetrically to both faults.
//
// Returns the number of polyline edges in the input; the output
// chain has approx round(total_length / target_edge_m) edges.
struct ResampleResult {
    int n_chains_processed = 0;
    int n_edges_in  = 0;
    int n_edges_out = 0;
};

ResampleResult resample_polyline_uniform_symmetric(
    Mesh& A, const EdgeConstrainedMap& polyline_ecm_A,
    Mesh& B, const EdgeConstrainedMap& polyline_ecm_B,
    double target_edge_m);

} // namespace
```

Implementation sketch (~80 LOC):

```cpp
// 1. Walk the polyline graph on A; find connected components.
//    Each component is a chain (or a closed loop) of polyline
//    edges.  Use a DFS with bookkeeping by (vertex, edge) → next.
std::vector<std::vector<Mesh::Vertex_index>> chains_a = walk_polyline(A, polyline_ecm_A);
std::vector<std::vector<Mesh::Vertex_index>> chains_b = walk_polyline(B, polyline_ecm_B);

// 2. Match A's chains to B's by canonical endpoint coords.
//    Both have bit-identical polyline endpoints post-corefine.
auto matched = match_chains_by_endpoints(A, chains_a, B, chains_b);

// 3. For each matched chain pair (cA, cB):
for (auto& [cA, cB] : matched) {
    // Compute arc length of cA (== arc length of cB by construction).
    double L = total_arc_length(A, cA);
    int n_new = std::max(1, (int)std::round(L / target_edge_m));
    double seg = L / n_new;

    // Compute new vertex positions by linear interp on cA.
    std::vector<Kernel::Point_3> new_pts;
    for (int k = 1; k < n_new; ++k) {
        double t = k * seg;
        new_pts.push_back(arc_length_to_point(A, cA, t));
    }

    // Replace chain in A and B with the new vertices/edges.
    // CGAL provides Euler::add_vertex_to_face / split_edge for this;
    // alternatively, mark the entire old chain as "to-delete"
    // and add new vertices + new edges in one pass.
    rebuild_chain(A, cA, new_pts, polyline_ecm_A);
    rebuild_chain(B, cB, new_pts, polyline_ecm_B);
}
```

The crucial property: `new_pts` is computed once per chain
(from cA's arc length), then applied identically to both A and
B. Polyline vertex coordinates are bit-equal across faults by
construction, so cross-fault coincidence is preserved without
relying on the link-condition check.

This replaces both `split_long_edges` (step 2 of `remesh_one_
fault`) AND `collapse_short_polyline_edges_symmetric` for the
polyline portion. Border-edge handling stays as-is (boundaries
don't need to match across faults).

**Test case:**
```cpp
TEST(R301_polyline_uniformly_distributed) {
  // Pre-load Mill Creek × SBMT-SAF; corefine; resample.
  auto stats = stat_polyline_edge_lengths(A, polyline_ecm_A);
  // After uniform resampling at target=1000:
  EXPECT_NEAR(stats.mean, 1000.0, 100.0);    // ±10% of target
  EXPECT_GT(stats.min, 700.0);               // no edges < 0.7 × target
  EXPECT_LT(stats.max, 1300.0);              // no edges > 1.3 × target
  EXPECT_LT(stats.cv, 0.10);                 // low coefficient of variation
}
```

**Why the simpler split-ceiling tweak (R-201) was insufficient:**
`split_long_edges` is one-sided. To get bilateral uniformity
you either (a) collapse short edges aggressively (which moves
polyline vertices and is risky for cross-fault conformity), or
(b) resample uniformly (the suggested fix). The plan's strategy
of "split + collapse" with conservative thresholds left a 60%-
of-edges dead zone that no pass touched.

---

### [R-302] [CRITICAL] [corefine.cpp:594 split_ceiling and end-to-end SAF aspect] — Raising split ceiling from target to target × 4/3 doubled the worst-case aspect on SBMT-SAF (7.6 → 15.19)

**Category:** BUG (regression introduced by the R-201 fix)

**Description:**
Direct measurement, target_edge_m = 1000:

| Round | Split ceiling | Mill aspect_max | **SAF aspect_max** |
|---|---|---:|---:|
| v2 | target (1000 m) | 4.7 | **7.6** |
| v3 | target × 4/3 (1333 m) | 4.93 | **15.19** |

The 175 m sliver polyline edge that escapes collapse (link-
condition fails) now sits adjacent to non-constrained
neighbour edges that the remesher pulled up to ~1333 m (the
new split ceiling, also tied to the remesher's upper edge
limit). Aspect ratio = 1333 / 175 ≈ 7.6 (theory), but with
multiple short-polyline-edges adjacent to the same vertex, the
worst-case longest non-constrained edge is even longer, giving
the observed 15:1.

The R-201 fix's analysis assumed the post-split distribution
would be tight around `target`. Empirically (R-301) it is
spread from 175 to 1315 m. R-201 made the LARGE end larger
(now 1333) while doing nothing for the SMALL end (still 175).
That maximised the spread → worse aspect.

This is a **net regression** vs v2 on the worst-case quality
metric. The user-visible "improvements are small" comment
applies to the gradient; on aspect ratio, it's strictly worse.

**Trigger:** Same input as the prior round (Mill Creek × SBMT-
SAF 2000 m).

**Actual behavior:** SAF post-Phase-2 has 3 triangles with
aspect > 5, max 15.19; min edge 66 m. Visible only at extreme
zoom in ParaView, but HXT in Phase 3 will object.

**Expected behavior:** Aspect ratios should be IMPROVED by the
R-201 fix, not regressed.

**Suggested fix:** R-301 (uniform resampling) eliminates the
short-edge tail entirely, so the aspect-ratio regression
disappears. As a STOP-GAP if R-301 lands later, revert R-201's
split ceiling to `target_edge_m`:

```diff
@@ corefine.cpp:594
-    const double split_ceiling = p.target_edge_m * 4.0 / 3.0;
+    // R-201/R-302: split ceiling reverted to `target_edge_m`
+    // pending R-301 (uniform resampling).  The 4/3 multiplier
+    // empirically widened the polyline-edge spread, doubling
+    // the worst-case aspect on SBMT-SAF (7.6 → 15.19).  Until
+    // resampling lands, the conservative `target` ceiling
+    // produces fewer slivers at the cost of a still-30%
+    // corridor gradient (which R-301 is the right fix for).
+    const double split_ceiling = p.target_edge_m;
```

Note the test `test_polyline_corridor_matches_target` will then
need its bounds widened back to `[0.65, 1.30] × target` (or
better, replaced by a real-data fixture once R-301 lands).

**Test case:**
```cpp
TEST(R302_no_aspect_regression_vs_v2) {
  // Run end-to-end on the Mill Creek × SBMT-SAF fixture.
  // After Phase 2, no triangle on either fault may have
  // aspect > 8 (catches the v3 regression of 15.19; v2 had 7.6).
  auto stats = full_pipeline_aspect_stats();
  EXPECT_LT(stats.max_aspect_mill, 8.0);
  EXPECT_LT(stats.max_aspect_saf,  8.0);
}
```

---

### [R-303] [MODERATE] [main.cpp:489–507, R-202 cleanup is implemented but disabled] — Disabled cleanup means non-polyline slivers persist; the comment correctly identifies the cross-fault asymmetry but doesn't propose the right alternative

**Category:** ASSUMPTION (gives up on a class of fix instead of solving it)

**Description:**
The implementer added `cleanup_non_polyline_slivers` (corefine.cpp:
487–541) but disabled the call in main.cpp:489–507 with the
comment:

> *"R-202 cleanup pass — DEFERRED. ... A removes 10 triangles, B
>   removes 9 — different decisions per fault — and the resulting
>   `is_border`-vs-interior status of one polyline edge differs
>   between A and B (R-203's concern about the border filter
>   materialising in practice)."*

The diagnosis is correct, but the response (defer) leaves the
non-polyline sliver class untreated. SBMT-SAF still has min
edge 66 m and aspect 15.19 because of these slivers.

The right alternative is either:

(a) **Lock the cleanup's vertex set symmetrically across A and
    B** before running. For each fault pair, compute the
    intersection of "candidates for cleanup" on A and on B (by
    triangle endpoint coords), then run cleanup only on the
    intersection on each side. CGAL doesn't expose
    "constraint-by-triangle-set" directly, but you can constrain
    every vertex incident to a non-symmetric candidate.

(b) **Run cleanup once on a "merged" representation**: combine
    A's and B's non-polyline edges into a single set, decide
    which to collapse based on cross-fault metrics, then apply
    the decisions symmetrically.

(c) **Pre-process the input STLs** with `experimental::snap_borders`
    before corefine, so the polyline doesn't pierce within ε of
    a vertex to begin with. This is the R-101 alternative
    described in the v2 review and is the most CGAL-idiomatic.

The current code abandons all three. The slivers persist. R-301
(uniform resampling) is the fix for the **polyline** edges; R-303
needs a separate fix for **non-polyline** near-polyline-vertex
slivers.

**Trigger:** Same as R-302; SAF residual 3 needles.

**Suggested fix:** Implement option (c) in a new pre-corefine
helper. Add to `corefine.cpp`:

```cpp
// Pre-corefine vertex weld: for each pair (A, B), find vertices
// of A within `tol` of any vertex of B (or vice versa) and snap
// them to a common position.  Eliminates the near-vertex pierce
// case at its source.
struct WeldResult {
    int n_welded_pairs = 0;
};
WeldResult preweld_near_vertices(Mesh& A, Mesh& B, double tol);
```

Implementation sketch (~40 LOC):

```cpp
// Build a KD-tree of B's vertex coordinates; query each A vertex.
// For each A vertex within tol of a B vertex: set both to the
// midpoint coord.  Cross-fault conformity is preserved by
// construction (same midpoint on both sides).
```

Add to main.cpp BEFORE the corefine loop:

```cpp
const double weld_tol = args.target_edge_m / 4.0;  // ~250 m for target 1000
for (std::size_t i = 0; i < N; ++i)
    for (std::size_t j = i + 1; j < N; ++j)
        safs::corefine::preweld_near_vertices(meshes[i], meshes[j], weld_tol);
```

This eliminates the source of the slivers (input vertex
geometry), so corefine produces only polyline edges of length
≥ `weld_tol`. R-301 then uniform-resamples those to
`[0.7, 1.3] × target_edge_m`. Combined effect: no slivers,
polyline corridor matches far-field within 10%.

**Test case:**
```cpp
TEST(R303_pre_weld_eliminates_short_polyline_edges) {
  // Construct two faults where A has a vertex 5 m from one of
  // B's vertices.  Without pre-weld, corefine produces a 5 m
  // sliver edge.  With pre-weld at tol=10 m, the two vertices
  // are merged BEFORE corefine — no sliver.
  ...
  for (auto e : A.edges()) {
    if (get(polyline_ecm_A, e)) {
      EXPECT_GT(edge_length(A, e), weld_tol * 0.5);
    }
  }
}
```

---

### [R-304] [LOW] [test_polyline_corridor_matches_target — synthetic fixture only] — The smoke test passes on a 10 km uniform square fixture but does not catch the real-data ratio of 0.62 × target

**Category:** QUALITY (test-coverage gap)

**Description:**
`test_polyline_corridor_matches_target` (lines 303–369) uses a
`make_square(10000.0)` fixture where the polyline is the x-axis,
length 20 km. Input has ONE polyline edge of 20 km, which
bisects to 10 km, then 5 km, then 2.5 km, then 1.25 km
(stops because 1.25 km ≤ 1.333 km). All resulting polyline
edges are exactly 1.25 km. Mean = 1.25 km ≈ 1.25 × target.

Test bounds: `0.70 ≤ mean ≤ 1.30 × target`. **Passes
trivially** because synthetic input is uniform.

Real CFM data has a wide input polyline-edge distribution
(175 m to 1500 m). Most input edges are SHORTER than `target`,
so `split_long_edges` doesn't touch them. The post-pass mean is
0.62 × target, **outside the test's 0.70 lower bound**, but
the smoke test never sees this because it uses the synthetic
fixture.

The Mill Creek × SBMT-SAF pin test (`test_mill_creek_pin`,
lines 191–234) measures `n_constrained_edges` but does NOT
measure mean edge length. So the regression is invisible to
CI.

**Trigger:** Real CFM data only. CI passes; humans see the
bad mesh.

**Suggested fix:** Add a real-data acceptance check to
`test_mill_creek_pin`:

```diff
@@ test_corefine_smoke.cpp test_mill_creek_pin (after line ~233)
     CHECK(pin >= 133 && pin <= 153,
           "Mill Creek 2000 m pin: ce should be 143 ± 10 ...");
+
+    // R-304: also verify the real-data polyline-edge mean is
+    // within ±15% of target after full Phase 2.  The squares
+    // smoke test cannot catch this regression because synthetic
+    // input has uniform edge distribution.
+    sc::EdgeConstrainedMap pr_a, pr_b;
+    install_ecm(A, pr_a, "e:protect_real");
+    install_ecm(B, pr_b, "e:protect_real");
+    for (auto e : A.edges()) if (get(ea, e)) put(pr_a, e, true);
+    for (auto e : B.edges()) if (get(eb, e)) put(pr_b, e, true);
+    sc::RemeshParams rp;
+    rp.target_edge_m = 1000.0;
+    rp.n_iterations  = 3;
+    rp.clearance_m   = 100.0;
+    sc::remesh_one_fault(A, ea, pr_a, rp);
+    sc::remesh_one_fault(B, eb, pr_b, rp);
+    auto polyline_mean = [&](const Mesh& M,
+                             const sc::EdgeConstrainedMap& ecm) {
+        double s = 0.0; int n = 0;
+        for (auto e : M.edges()) {
+            if (!get(ecm, e)) continue;
+            if (CGAL::is_border(e, M)) continue;
+            auto h = M.halfedge(e);
+            s += std::sqrt((M.point(M.target(h)) -
+                            M.point(M.source(h))).squared_length());
+            ++n;
+        }
+        return n > 0 ? s/n : 0.0;
+    };
+    const double mean_a = polyline_mean(A, pr_a);
+    const double mean_b = polyline_mean(B, pr_b);
+    std::cerr << "[smoke]   Mill×SAF mean polyline edge "
+              << "A=" << mean_a << " B=" << mean_b
+              << " (target=" << rp.target_edge_m << ")\n";
+    CHECK(std::abs(mean_a - rp.target_edge_m) <= 0.15 * rp.target_edge_m,
+          "Mill polyline mean edge length is outside ±15% of target — "
+          "R-301 uniform resampling not landed yet");
+    CHECK(std::abs(mean_b - rp.target_edge_m) <= 0.15 * rp.target_edge_m,
+          "SAF polyline mean edge length is outside ±15% of target");
```

This will FAIL with the current implementation (mean ≈ 615 m
vs target 1000 m), correctly identifying the R-301 regression.
Once R-301 lands, the test passes.

---

## Summary

- Critical issues: **2** (R-301 polyline-edge distribution dead
  zone, R-302 R-201 made aspect_max worse)
- Moderate issues: **1** (R-303 R-202 cleanup disabled —
  alternate fix needed)
- Low issues: **1** (R-304 smoke-test coverage gap)
- Plan compliance: **PARTIAL**.
  - Phase 2 acceptance "max aspect ≤ 4 on triangles incident to
    a constrained edge": **REGRESSED** (v2 → v3: 7.6 → 15.19 on
    SBMT-SAF).
  - Phase 2 acceptance "mean edge length in [800, 1200] m":
    Globally still passing; near-polyline mean (720 m) is
    outside the band.
- Verdict: **FAIL — must fix before proceeding**. The R-201
  fix's net effect is roughly zero on the gradient AND a
  regression on aspect ratio. R-301 is the structural fix the
  user's instinct points at; without it the Phase 2 acceptance
  is unmet on real data.

## Why R-201 alone wasn't enough

R-201's analysis: *"split ceiling = target → post-split edges in
(target/2, target], mean = 0.75 × target → corridor 30% finer."*

The analysis assumed input polyline edges are LONG (everything
> target) so split_long_edges actually runs on them. On real
CFM data, input edges are MOSTLY shorter than target:

```
input polyline edge length (Mill Creek × SBMT-SAF, post-corefine):
  median ≈ 700 m, std ≈ 350 m
  60% of edges are in [250, 700) m — never touched by split_long_edges
```

`split_long_edges(L_max)` is a one-sided primitive. When most
edges start below `L_max`, increasing `L_max` does almost
nothing. The fix needs the OTHER direction (subdivide-and-merge)
or a complete replacement (uniform resampling = R-301).

## Suggested fix order

1. **R-302** (revert split ceiling to `target`) — one-line
   safety stop-gap. Restores v2's better aspect_max while R-301
   is being implemented. Land first.
2. **R-301** (uniform polyline resampling) — the structural
   fix. ~80 LOC of new helper plus rewiring `remesh_one_fault`
   step 2. Eliminates the polyline-edge dead zone AND the
   sliver class that R-101's collapse left behind. Most of the
   visible improvement the user wants.
3. **R-303** (pre-corefine vertex weld) — kills the residual
   non-polyline sliver class at the source. ~40 LOC. Land
   together with R-301 if possible; the two are
   complementary (R-301 fixes polyline edges, R-303 fixes
   non-polyline near-polyline-vertex edges).
4. **R-304** (smoke test extension) — must land WITH R-301 so
   future regressions are caught. The current smoke test would
   have passed even at 60% of target.

## Rough effort estimate

| Finding | LOC | Risk |
|---|---:|---|
| R-302 (revert ceiling) | 1 line | Low; reverts to verified v2 behaviour |
| R-301 (resampling) | ~80–100 | Medium; new geometric algorithm |
| R-303 (pre-weld) | ~40 | Low–medium; standard KD-tree query |
| R-304 (test) | ~30 | Low |

Total: ~150 LOC for a Phase-2 output that would actually meet
the plan's acceptance criteria.

## Unreviewed Areas

- **CGAL `Euler::collapse_edge`'s exact link-condition behaviour on the
  3 sub-250 m polyline edges that survive in the SAF run.** I did
  not trace which specific edges fail the link condition. If
  R-303's pre-weld eliminates them at source the question is
  moot.
- **Whether `experimental::snap_borders` (CGAL ≥ 5.4) does the
  pre-weld natively.** R-303 sketches a hand-rolled solution; if
  CGAL has a one-liner, prefer it.
- **The CGAL `cgal_patch` shim** — not re-checked this round.
- **HXT (Phase 3) acceptance.** Phase 3 still pending; once it's
  on the bench, the aspect-15.19 triangles will probably get
  rejected by HXT and the Phase-3 driver will need to surface
  that as a diagnostic.
