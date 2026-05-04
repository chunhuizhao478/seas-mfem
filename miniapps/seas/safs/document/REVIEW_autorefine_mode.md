# Code Review: autorefine_mode.cpp + main.cpp wiring (2026-04-30)

## Review Scope

- **Plan:** `PLAN_cgal_corefine.md` (Phases 1–5).  The autorefine path
  is documented in the file header (`autorefine_mode.hpp`) as the
  6-faults-or-more replacement for the pairwise cascade — it is NOT
  described in the plan, so this is **scope creep** and must satisfy
  the plan's invariants on its own.
- **Files reviewed:**
  - `tools/corefine_faults/autorefine_mode.hpp`  (58 lines)
  - `tools/corefine_faults/autorefine_mode.cpp`  (497 lines, NEW)
  - `tools/corefine_faults/main.cpp`             (autorefine
                                                  branch lines
                                                  364–436)
  - `tools/CMakeLists.txt`                       (target wiring)
- **Live verification:**
  - The new code builds cleanly.
  - I did NOT exercise the autorefine path end-to-end on any fixture;
    findings are from static analysis + cross-reference with CGAL
    5.6.1 documentation.

The autorefine path is a fundamentally different algorithm from the
pairwise-cascade path Phases 1–4 specify.  It is the production path
for `--mode autorefine`.  I assumed at least 3 bugs and found 7
(4 CRITICAL, 3 MODERATE).

---

## Findings

### [R-501] [CRITICAL] [autorefine_mode.cpp:105 + 119–136] — `orient_polygon_soup` may DUPLICATE polygons; their fault tags are never assigned and the demux silently loses them

**Category:** BUG (silent face loss)

**Description:**
Per CGAL 5.6.1 docs for
`Polygon_mesh_processing::orient_polygon_soup`:

> When it is not possible to produce a consistent orientation, some
> polygons are duplicated.  The duplicates are appended at the end of
> the polygon list.

The autorefine driver:

```cpp
PMP::orient_polygon_soup(all_points, all_polygons);   // line 105
// ...
PMP::polygon_soup_to_polygon_mesh(all_points, all_polygons, combined);
// ...
std::int64_t fi = 0;
for (auto f : combined.faces()) {
    if (fi < static_cast<std::int64_t>(face_to_fault.size())) {  // line 131
        put(fault_pm, f, face_to_fault[fi]);
    }
    ++fi;
}
```

`face_to_fault` was built BEFORE `orient_polygon_soup` ran.  Any
polygon `orient_polygon_soup` appended is at index ≥
`face_to_fault.size()`, so the `if` condition is false → those faces
are LEFT WITH `fault_pm = -1` (the default).

The demux loop at line 354 then drops them:
```cpp
const std::int64_t fid = get(fault_pm, f);
if (fid < 0) { ++n_skipped_no_tag; continue; }
```

Result: faces that should be in fault A's output STL are silently
omitted.  No warning, no error — just `n_skipped_no_tag++` in the
demux diagnostic line.  Per `safs_mult_banning` and
`safs_mult_ssaf_banning` are typical CFM cases where consistent
orientation across faults is impossible (each fault has its own
normal direction), so duplicates ARE expected on real input.

**Trigger:** Any all-8 build where two faults disagree on normal
direction (the documented all-8 case).

**Actual behavior:** Silent face loss, asymmetric per-fault output
sizes that look "almost right" in the diagnostic.

**Expected behavior:** Either (a) extend `face_to_fault` after
`orient_polygon_soup` returns by mapping appended duplicates back to
their pre-orient parent, or (b) fail-fast if duplicates are
appended (the orient step is documented to indicate when this
happened).

**Suggested fix:**

```diff
@@ autorefine_mode.cpp:105
-    PMP::orient_polygon_soup(all_points, all_polygons);
+    const std::size_t n_polys_pre_orient = all_polygons.size();
+    PMP::orient_polygon_soup(all_points, all_polygons);
+    if (all_polygons.size() != n_polys_pre_orient) {
+        // R-501: orient_polygon_soup duplicated polygons to make the
+        // soup orientable (e.g., two faults disagree on normal
+        // direction).  The duplicates are at indices [n_pre, end);
+        // tag them with the SAME fault as their pre-orient parent.
+        // This preserves face provenance through the post-autorefine
+        // demux.
+        if (n_polys_pre_orient != face_to_fault.size()) {
+            throw std::runtime_error(
+                "[autorefine_faults] internal: face_to_fault is out "
+                "of sync with all_polygons before orient_polygon_soup");
+        }
+        face_to_fault.reserve(all_polygons.size());
+        for (std::size_t i = n_polys_pre_orient;
+             i < all_polygons.size(); ++i) {
+            // Each appended duplicate is a copy of polygon (i - n_pre)
+            // — see CGAL::Polygon_mesh_processing::orient_polygon_soup
+            // implementation.  Pin to the same fault.
+            face_to_fault.push_back(face_to_fault[i - n_polys_pre_orient]);
+        }
+    }
```

If the "appended duplicates copy index (i - n_pre) of the original"
mapping is not actually documented, the safer alternative is option
(b) fail-fast:

```diff
+    if (all_polygons.size() != n_polys_pre_orient) {
+        throw std::runtime_error(
+            "[autorefine_faults] orient_polygon_soup duplicated "
+            + std::to_string(all_polygons.size() - n_polys_pre_orient)
+            + " polygon(s) to make the combined soup orientable.  "
+              "Demux cannot map duplicates back to their fault; aborting. "
+              "Investigate which fault(s) disagree on normal direction.");
+    }
```

**Test case:**
```cpp
TEST(R501_orient_appends_duplicates_no_silent_loss) {
    // Build two faults whose triangle pair shares all three vertex
    // coordinates but with opposite winding.  orient_polygon_soup
    // will duplicate the polygon to make the soup orientable.
    // After autorefine_faults, the count of skipped-no-tag faces
    // in the demux must be 0 (duplicates either correctly tagged
    // or the function aborts).
    Mesh A = make_triangle({0,0,0}, {1,0,0}, {0,1,0}, /*ccw*/true);
    Mesh B = make_triangle({0,0,0}, {1,0,0}, {0,1,0}, /*ccw*/false);
    write_stl("/tmp/A.stl", A); write_stl("/tmp/B.stl", B);
    auto rr = autorefine_faults({"/tmp/A.stl","/tmp/B.stl"},
                                 {"/tmp/Aout.stl","/tmp/Bout.stl"},
                                 {"A","B"});
    // Either the call throws (option b), or the per-fault outputs
    // each contain at least one face (option a).
    EXPECT_TRUE(rr.per_fault_n_output[0] > 0
             && rr.per_fault_n_output[1] > 0);
}
```

---

### [R-502] [CRITICAL] [autorefine_mode.cpp:73–117] — Cross-fault coincident vertices are not deduped in the combined soup; produce duplicate vertices in the Surface_mesh that autorefine cannot reduce

**Category:** BUG

**Description:**
Per-fault `stl_to_soup` runs `repair_polygon_soup(erase_all_duplicates=true)`,
which dedups vertices WITHIN each fault.  But the combined soup is
built by APPENDING each fault's `pts` into `all_points`:

```cpp
const std::size_t off = all_points.size();
for (auto& p : pts) all_points.push_back(p);
for (auto& poly : polys) {
    std::vector<std::size_t> shifted;
    shifted.reserve(poly.size());
    for (auto idx : poly) shifted.push_back(idx + off);
    all_polygons.push_back(std::move(shifted));
    // ...
}
```

If fault A and fault B both contain a vertex at coordinate
(123.0, 456.0, -789.0), fault A's index_in_all_points is, say, 17,
and fault B's is, say, 42.  These are two DISTINCT vertex_indices in
`all_points` despite being the same 3-D point.

`PMP::polygon_soup_to_polygon_mesh` uses the `all_polygons` indices
verbatim — it does NOT de-duplicate vertices by coordinate.  The
resulting Surface_mesh has TWO vertex_descriptors at the same XYZ.

When `autorefine` then runs:
- triangles in fault A and fault B that share the (123, 456, -789)
  vertex by COORDINATE do not share it by HALFEDGE — autorefine
  treats them as a near-zero-length edge between two coincident
  vertices and may produce zero-area faces or fail
  `does_self_intersect` checks unpredictably.
- The resulting per-fault demux produces inconsistent polyline
  geometry between A's STL and B's STL: both have the (123, 456,
  -789) point, but different vertex_indices on either side, with
  potentially different connectivity post-autorefine.

The comment at line 98–104 explicitly says "do NOT erase duplicate
POLYGONS across faults" — but the per-VERTEX dedup is what's
actually missing here, and it is necessary.

**Trigger:** Any input where two faults have a vertex at the same
3-D coordinate.  Common when faults share a CFM tile boundary (the
SAF-Banning system has this).

**Actual behavior:** Coincident vertices duplicated in the combined
mesh; autorefine sees a degenerate-edge configuration; the demux
output has one of the two faults with subtly broken polyline
endpoints.

**Expected behavior:** The combined soup should call
`repair_polygon_soup` with vertex-dedup but NOT polygon-dedup.

**Suggested fix:**

```diff
@@ autorefine_mode.cpp:103–105
     // CRITICAL: do NOT erase duplicate POLYGONS across faults.  Two
     // faults can have triangles whose three vertex coords coincide
     // (e.g., adjacent SAF segments share boundary triangles); the
     // legacy `erase_all_duplicates(true)` flag dedups them and the
     // per-fault demux loses one of the two faults' provenance.
     // Instead, only ensure the soup is orientable; autorefine itself
     // will then process the genuine self-intersections.
+    // R-502: but DO dedup VERTICES across faults.  Two faults that
+    // share a CFM tile boundary have bit-coincident vertex
+    // coordinates; without coordinate-based vertex dedup, the
+    // combined Surface_mesh has two distinct vertex_indices at the
+    // same XYZ and autorefine cannot produce a coherent shared
+    // polyline.
+    PMP::repair_polygon_soup(
+        all_points, all_polygons,
+        CGAL::parameters::erase_all_duplicates(false)
+                         .require_same_orientation(false));
     PMP::orient_polygon_soup(all_points, all_polygons);
```

`erase_all_duplicates(false)` keeps polygon duplicates (per the
comment's intent) but `repair_polygon_soup` STILL deduplicates
vertices by coordinate (its primary documented behavior).
Verify against CGAL 5.6.1 docs that this combination is supported.

**Test case:**
```cpp
TEST(R502_cross_fault_vertex_dedup) {
    // Two faults that share vertex (0,0,0).
    write_stl("/tmp/A.stl", triangle({0,0,0}, {1,0,0}, {0,1,0}));
    write_stl("/tmp/B.stl", triangle({0,0,0}, {0,0,1}, {-1,0,0}));
    auto rr = autorefine_faults(...);
    // Read both output STLs; the (0,0,0) vertex should appear exactly
    // once across the two combined sets (i.e., the polyline endpoint
    // is bit-identical between A's STL and B's STL).
    Mesh Aout = read("/tmp/Aout.stl"), Bout = read("/tmp/Bout.stl");
    auto a_zero_count = count_vertices_at(Aout, {0,0,0});
    auto b_zero_count = count_vertices_at(Bout, {0,0,0});
    EXPECT_EQ(a_zero_count, 1);
    EXPECT_EQ(b_zero_count, 1);
    // Distance between the two zero-vertices must be 0 (bit-equal),
    // not float-noise away.
    EXPECT_EQ(get_zero_vertex(Aout), get_zero_vertex(Bout));
}
```

---

### [R-503] [CRITICAL] [autorefine_mode.cpp:111–117] — `polygon_soup_to_polygon_mesh` is called with a soup that the precondition explicitly rejects; behaviour is undefined per CGAL docs

**Category:** BUG (precondition violation)

**Description:**

```cpp
if (!PMP::is_polygon_soup_a_polygon_mesh(all_polygons)) {
    std::cerr << "[autorefine_faults] WARN: combined soup is NOT a "
                 "polygon mesh after orient/repair — autorefine may "
                 "still proceed, but cross-fault conformity is at "
                 "risk.\n";
}
PMP::polygon_soup_to_polygon_mesh(all_points, all_polygons, combined);
```

CGAL's `polygon_soup_to_polygon_mesh` documents the
**precondition** that the soup must already be a polygon mesh.
Calling it on a non-mesh soup is undefined behaviour — it may
produce an invalid Surface_mesh, an assertion failure, or a
silently-broken mesh.

The current code prints a warning and proceeds anyway.  This is
exactly the path the all-8 build with cross-fault coincident
vertices (R-502) would hit.

**Trigger:** Any input that fails the
`is_polygon_soup_a_polygon_mesh` check.  R-502 can produce this on
real CFM data; so can any input with non-orientable cross-fault
geometry.

**Actual behavior:** Silently constructs a possibly-invalid
Surface_mesh and proceeds.  Downstream symptoms (autorefine failures,
demux face loss, HXT rejection) are difficult to attribute back to
this root cause.

**Expected behavior:** Fail-fast with a diagnostic that names the
specific culprit (which polygon/vertex pair makes the soup
non-mesh).

**Suggested fix:**

```diff
@@ autorefine_mode.cpp:111–117
-    if (!PMP::is_polygon_soup_a_polygon_mesh(all_polygons)) {
-        std::cerr << "[autorefine_faults] WARN: combined soup is NOT a "
-                     "polygon mesh after orient/repair — autorefine may "
-                     "still proceed, but cross-fault conformity is at "
-                     "risk.\n";
-    }
-    PMP::polygon_soup_to_polygon_mesh(all_points, all_polygons, combined);
+    if (!PMP::is_polygon_soup_a_polygon_mesh(all_polygons)) {
+        // R-503: polygon_soup_to_polygon_mesh's precondition is that
+        // the soup IS a polygon mesh; calling it otherwise is UB.
+        // Fail with a diagnostic listing the per-fault triangle
+        // counts so the user can attribute the failure to the right
+        // input fault.
+        std::ostringstream msg;
+        msg << "[autorefine_faults] FATAL: combined polygon soup is "
+               "NOT a polygon mesh after orient/repair.  Per-fault "
+               "input triangle counts:\n";
+        for (std::size_t k = 0; k < N; ++k) {
+            msg << "  " << shorts[k] << ": "
+                << res.per_fault_n_input[k] << " triangles\n";
+        }
+        msg << "Likely causes: cross-fault duplicate vertices "
+               "(R-502) or genuinely non-orientable input (some "
+               "CFM faults agree on neither orientation).  Aborting.";
+        throw std::runtime_error(msg.str());
+    }
+    PMP::polygon_soup_to_polygon_mesh(all_points, all_polygons, combined);
```

**Test case:**
```cpp
TEST(R503_non_polygon_mesh_aborts) {
    // Build a soup that fails is_polygon_soup_a_polygon_mesh
    // (e.g., three triangles sharing one edge).  Call autorefine_faults;
    // assert it throws with "FATAL: combined polygon soup".
    EXPECT_THROW(autorefine_faults(...), std::runtime_error);
}
```

---

### [R-504] [CRITICAL] [autorefine_mode.cpp:419–469] — `target_edge_m` and `remesh_iters` parameters are SILENTLY IGNORED in autorefine mode (the entire post-autorefine cleanup block is wrapped in `if (false)`)

**Category:** BUG (silent feature drop)

**Description:**

```cpp
// Step 6.5: per-fault sliver cleanup is DELIBERATELY OMITTED.
// ...
(void)target_edge_m;
(void)remesh_iters;
if (false) {
    for (std::size_t k = 0; k < N; ++k) {
        // ... remove_almost_degenerate_faces + isotropic_remeshing ...
    }
}
```

The autorefine driver's API takes `target_edge_m = 1000.0` and
`remesh_iters = 3` as parameters (with the same defaults the cascade
mode uses).  `main.cpp:381` passes them through:
```cpp
auto rr = safs::corefine::autorefine_faults(
    in_paths, out_paths, shorts,
    args.no_remesh ? 0.0 : args.target_edge_m,
    args.no_remesh ? 0   : args.remesh_iters);
```

But `autorefine_mode.cpp` accepts the values, casts them to void to
suppress unused-variable warnings, and never uses them.  The
sliver-cleanup block that WOULD use them is wrapped in `if (false)`
— dead code.

User effect:
- `corefine_faults --mode autorefine --target-edge-m 1000 --remesh-iters 3`
  ignores BOTH flags entirely; output is the raw autorefine result
  with no remeshing.  This produces a mesh with whatever edge-length
  distribution `autorefine` happens to leave behind, which is
  typically dominated by the input triangle sizes — NOT the user's
  requested target edge length.
- This silently breaks Phase 3's acceptance criterion (mean edge
  ∈ [500, 2000] m at `res_f`=1000 m) when autorefine mode is used.

The comments at lines 405–418 + 456–465 explain the rationale ("any
per-fault edit breaks HXT PLC recovery on real CFM data"), which IS
a real problem — but silently dropping the user-requested behaviour
is the wrong response.  The fix is either to (a) document the
limitation in `--help` AND `--mode autorefine` rejects non-default
target_edge_m / remesh_iters, or (b) print a clear runtime warning
when those flags are used in autorefine mode.

**Trigger:** Any `--mode autorefine` invocation with a positive
`--target-edge-m` (the default) and the user is unaware that the
flag is ignored.

**Actual behavior:** Flags accepted, values ignored, no warning.

**Expected behavior:** Either honour the flags (re-enable the
sliver-cleanup block — the comments correctly note this is unsafe
at the multi-fault HXT-PLC stage, so this is hard) OR explicitly
warn the user at runtime.

**Suggested fix (minimal — runtime warning):**

```diff
@@ autorefine_mode.cpp:419
-    (void)target_edge_m;
-    (void)remesh_iters;
-    if (false) {
+    if (target_edge_m > 0.0 || remesh_iters > 0) {
+        // R-504: autorefine mode does NOT honour --target-edge-m
+        // or --remesh-iters because per-fault remesh corrupts the
+        // bit-identical cross-fault polyline endpoints autorefine
+        // produces (HXT then rejects the combined PLC).  Warn loudly
+        // so the user knows their flags are not honoured.
+        std::cerr << "[autorefine_faults] WARN: --mode autorefine "
+                     "ignores --target-edge-m (" << target_edge_m
+                  << ") and --remesh-iters (" << remesh_iters
+                  << ").  The output STL retains the input fault "
+                     "triangulation's edge sizes; downstream gmsh "
+                     "size-field tuning is the only knob that "
+                     "controls the bulk-mesh edge size near faults "
+                     "in this mode.  Use `--mode cascade` if "
+                     "per-fault remeshing is required.\n";
+    }
+    if (false) {
         for (std::size_t k = 0; k < N; ++k) {
```

(Keeps the dead `if (false)` block as is — the comments around it
are valuable and the dead code communicates "this WAS tried, here is
why it's disabled".)

**Test case:**
```bash
# R-504 acceptance: --mode autorefine emits a WARN line about
# --target-edge-m being ignored.
out=$(corefine_faults --mode autorefine --target-edge-m 500 \
                       --remesh-iters 3 ... 2>&1)
echo "$out" | grep -q "WARN.*--mode autorefine.*ignores --target-edge-m"
```

---

### [R-505] [MODERATE] [autorefine_mode.cpp:228–270] — Point-in-triangle tie-breaking by perpendicular distance is noise-driven near cross-fault polylines

**Category:** ASSUMPTION (numerical fragility)

**Description:**
The post-autorefine demux tags each face by point-in-triangle search:

```cpp
const bool inside = (sa >= -kBaryEps)
                 && (sb >= -kBaryEps)
                 && (sc >= -kBaryEps);
if (!inside) continue;
const double perp = CGAL::scalar_product(pt.normal, dvec);
const double perp2 = perp * perp / pt.area2;
if (perp2 < best_in_perp2) {
    best_in_perp2 = perp2;
    best_in = pt.fault;
}
```

A post-autorefine face whose centroid lies on a cross-fault polyline
projects-into BOTH parents (one from fault A, one from fault B).
Both pass the `inside` check.  The tiebreaker selects the parent
with smaller perpendicular distance to the centroid's plane.

But for cross-fault polylines, the centroid lies on BOTH planes
(the polyline IS the plane intersection); both `perp2` values are
~0 ± float64 noise.  Whichever floats first wins, and the choice
is essentially random.

The downstream effect: a face that should belong to fault A may
get tagged as fault B.  The `add_face` call on B's mesh either
succeeds (B's mesh now has a triangle from A's actual surface,
producing a non-manifold mismatch) or fails-with-null (face lost
entirely — `n_skipped_add_face_null++`).  Either way, the per-fault
output STLs are not faithful to the input fault decomposition.

**Trigger:** Any face whose centroid lies on a cross-fault polyline
(geometrically: any face directly adjacent to the polyline within
~ε of either fault's plane).  Common.

**Suggested fix:** Break ties by the SIGN of the perpendicular
distance.  After autorefine, each post-face's normal is consistent
with ONE pre-tri's normal (autorefine preserves normal direction
through subdivision).  Pick the pre-tri with the same normal sign
on the post-face's vertices:

```diff
@@ autorefine_mode.cpp:262–270
-            const auto dvec = cent - pt.a;
-            const double perp = CGAL::scalar_product(pt.normal, dvec);
-            const double perp2 = perp * perp / pt.area2;
-            if (perp2 < best_in_perp2) {
-                best_in_perp2 = perp2;
-                best_in = pt.fault;
-            }
+            // R-505: tie-break inside-multiple parents by NORMAL
+            // ALIGNMENT, not perpendicular-distance noise.  The
+            // post-face inherits its normal direction from exactly
+            // one parent; the parent whose normal points the SAME
+            // way as the post-face's normal is the correct one.
+            const auto dvec = cent - pt.a;
+            const double perp = CGAL::scalar_product(pt.normal, dvec);
+            const double perp2 = perp * perp / pt.area2;
+            // Compute post-face normal locally.
+            const auto post_normal = CGAL::cross_product(b - a, c - a);
+            const double align = CGAL::scalar_product(pt.normal, post_normal);
+            // Skip parents whose normal direction is opposite the
+            // post-face's; only consider co-oriented parents.
+            if (align <= 0.0) continue;
+            if (perp2 < best_in_perp2) {
+                best_in_perp2 = perp2;
+                best_in = pt.fault;
+            }
```

**Test case:**
```cpp
TEST(R505_polyline_tiebreak_uses_normal_alignment) {
    // Two faults intersecting along a polyline.  Post-autorefine,
    // a face with centroid on the polyline must be assigned to the
    // fault whose normal it shares, NOT the other one.
    Mesh A = make_square_xy(...);
    Mesh B = make_square_xz(...);
    auto rr = autorefine_faults(...);
    // Sum of per-fault output triangle counts must equal the
    // post-autorefine combined count (no faces lost to the
    // wrong-fault add_face-null path).
    EXPECT_EQ(rr.per_fault_n_output[0] + rr.per_fault_n_output[1],
              rr.n_output_faces);
}
```

---

### [R-506] [MODERATE] [autorefine_mode.cpp:257] — `kBaryEps = 1e-7` is in barycentric units, giving a triangle-size-dependent world tolerance

**Category:** ASSUMPTION

**Description:**
```cpp
constexpr double kBaryEps = 1e-7;
const bool inside = (sa >= -kBaryEps)
                 && (sb >= -kBaryEps)
                 && (sc >= -kBaryEps);
```

Barycentric coordinates are in [0, 1] regardless of triangle size.
A negative-by-1e-7 barycentric coordinate corresponds to a world-space
distance of `1e-7 × edge_length` (for centroid offsets along an
edge).  For a typical CFM fault triangle (edge ~1500 m) this is
~0.15 mm — fine for absorbing float64 round-off.  But for a tiny
post-autorefine sliver (edge ~10 m) the same tolerance is ~1 µm —
TOO TIGHT, may falsely reject genuine inside-the-triangle centroids
near the polyline.  Conversely for large pre-autorefine fault tiles
(edge ~10 km), the tolerance corresponds to ~1 mm — looser than
strictly needed.

Symptom: scattered `n_outside_fallback` increments where the
inside-check should have succeeded; the fallback then uses the noisy
nearest-centroid tiebreaker.

**Suggested fix:** Convert to a world-space tolerance and divide by
the local edge length to get the equivalent barycentric:

```diff
@@ autorefine_mode.cpp:257
-            constexpr double kBaryEps = 1e-7;
+            // R-506: keep tolerance in WORLD units (1 mm), then
+            // convert to local barycentric units by dividing by a
+            // characteristic length of the triangle (sqrt(area2)).
+            constexpr double kWorldEpsM = 1e-3;   // 1 mm
+            const double bary_eps = kWorldEpsM / std::sqrt(pt.area2);
             const bool inside = (sa >= -kBaryEps)
                              && (sb >= -kBaryEps)
                              && (sc >= -kBaryEps);
```
(replace all three `kBaryEps` with `bary_eps`.)

**Test case:**
```cpp
TEST(R506_inside_check_robust_to_tri_size) {
    // Build two parents: one with edge 10 m, one with edge 10 km.
    // A post-face centroid that lies ~ 0.5 mm INSIDE the triangle
    // boundary (numerically) must be classified as inside in BOTH
    // size regimes, not just the medium-size case.
    EXPECT_TRUE(barycentric_inside(small_tri, p_05mm_in));
    EXPECT_TRUE(barycentric_inside(large_tri, p_05mm_in));
}
```

---

### [R-507] [MODERATE] [autorefine_mode.cpp:396–402] — Demux skip counts are logged but never fatal; silent face loss

**Category:** ASSUMPTION

**Description:**

```cpp
std::cerr << "[autorefine_faults] demux: "
          << combined.number_of_faces() << " in combined, "
          << n_added << " added, "
          << n_collapsed << " degenerate, "
          << n_skipped_no_tag << " no-tag, "
          << n_skipped_bad_tag << " bad-tag, "
          << n_skipped_add_face_null << " add_face-null\n";
```

If the demux loses 50 % of faces (e.g., R-501 fires on a
non-orientable input), `n_skipped_no_tag` is large and the per-fault
STL outputs are silently incomplete.  No threshold gates this; the
function returns success.  The downstream pipeline writes the
incomplete STLs as if they were correct, and the user sees mesh
artefacts in ParaView with no error.

**Suggested fix:** Compute the skip fraction and fail if it exceeds
a sanity threshold (say 1 %):

```diff
@@ autorefine_mode.cpp:402
+        const std::int64_t n_skipped_total =
+            n_skipped_no_tag + n_skipped_bad_tag + n_skipped_add_face_null;
+        const std::int64_t n_total_in = static_cast<std::int64_t>(
+            combined.number_of_faces());
+        if (n_total_in > 0
+            && n_skipped_total > std::max<std::int64_t>(10, n_total_in / 100)) {
+            // R-507: more than 1 % (or > 10 absolute) of faces lost
+            // in the demux is almost certainly an algorithmic
+            // failure (R-501, wrong tagging, etc.), not legitimate
+            // edge-case dropping.  Fail rather than ship an
+            // incomplete output.
+            throw std::runtime_error(
+                "[autorefine_faults] demux dropped "
+                + std::to_string(n_skipped_total) + " of "
+                + std::to_string(n_total_in)
+                + " faces; refusing to write incomplete output. "
+                  "Check the demux diagnostic line above for the "
+                  "specific skip category.");
+        }
```

**Test case:**
```cpp
TEST(R507_high_demux_skip_rate_aborts) {
    // Build a fixture that triggers R-501 (orient appends
    // duplicates) without the R-501 fix in place; expect the demux
    // to drop > 1% of faces; expect the autorefine call to throw.
    EXPECT_THROW(autorefine_faults(...), std::runtime_error);
}
```

---

## Summary

- Critical issues: **4** (R-501 silent face loss from orient
  append; R-502 cross-fault vertex dedup; R-503 precondition
  violation on non-mesh soup; R-504 `--target-edge-m` /
  `--remesh-iters` silently ignored)
- Moderate issues: **3** (R-505 polyline tiebreak; R-506
  triangle-size-dependent tolerance; R-507 silent demux skip)
- Low issues: 0
- Plan compliance: **NONE** — autorefine mode is not described
  in `PLAN_cgal_corefine.md`.  This entire path is undocumented
  scope creep; even if every finding above is fixed, the path
  needs a plan amendment before it can be considered production.
- Verdict: **FAIL — must fix before proceeding**.  R-501 + R-502 +
  R-503 are root causes that would silently corrupt the all-8
  build (the documented motivation for autorefine mode); R-504
  silently breaks the user's CLI contract.

## Direct answer to the user's concern about mesh generation

The new `autorefine_mode.cpp` is the production path on `--mode
autorefine` and is wired into `main.cpp:364`.  It has four CRITICAL
bugs that will silently corrupt the all-8 mesh:

1. **R-501**: `orient_polygon_soup` appends duplicates to make
   non-orientable cross-fault soups orientable; their fault tags
   are never assigned; the demux drops them.
2. **R-502**: Two faults sharing a 3-D vertex coordinate end up
   with TWO distinct vertex_indices in the combined mesh; autorefine
   sees a near-zero-edge configuration; demux output has subtly
   broken polyline endpoints between the two faults' STLs.
3. **R-503**: `polygon_soup_to_polygon_mesh` is called on a non-mesh
   soup (warning printed, proceeds anyway); per CGAL docs this is
   undefined behaviour.
4. **R-504**: User-supplied `--target-edge-m` and `--remesh-iters`
   are silently ignored — the post-autorefine cleanup block is
   dead code (`if (false)`).  The user's expectation of
   res_f-controlled fault mesh size is violated.

The pairwise-cascade path (`--mode cascade`, the default) is
unaffected by these findings — it is the path Phases 1–4 of the
plan describe and the path the 2-fault Phase-3 acceptance ran
through.

## Suggested fix order

1. **R-503** (one-line throw on non-mesh soup) — fail-fast surfaces
   R-501 + R-502 as their actual root cause rather than letting
   them produce silent corruption.
2. **R-502** (vertex dedup in combined soup) — fixes the
   "non-mesh soup" trigger from one of its sources.
3. **R-501** (face_to_fault extension or fail-fast on appended
   duplicates) — fixes the OTHER source of "non-mesh soup".
4. **R-504** (runtime warning when target_edge_m/remesh_iters > 0
   in autorefine mode) — restores user-CLI contract honesty.
5. **R-505 + R-506 + R-507** — quality / robustness; can land
   together.

## Unreviewed Areas

- **`PMP::experimental::autorefine` semantics in CGAL 5.6.1.**  The
  `experimental` namespace is documented as preliminary and
  subject-to-change.  The plan's pinned ABI window
  [5.4, 6.0) does not guarantee `experimental` API stability across
  minor bumps.  Any CGAL upgrade should re-validate this path's
  output.
- **End-to-end run of `--mode autorefine` on the all-8 fixture.**
  The implementer's last reported run (Phase 5) used `--mode
  cascade` (default) and exited 2 on `safs_mult_banning`'s
  self-intersection.  I did not verify that `--mode autorefine`
  produces a valid all-8 STL set on the same input — likely it
  doesn't, given R-501–R-503.
- **`combined_pts = combined.points();` (line 302).**  This copies
  all vertex coords into a local variable that is never used
  afterwards.  Dead code; not a bug, just clutter.
- **The `Surface_mesh::add_face` failure modes captured by
  `n_skipped_add_face_null`.**  Not enumerated in the diagnostic;
  could be (a) duplicate edge with opposite orientation, (b)
  non-manifold edge, (c) other.  If R-507's threshold is hit,
  knowing which sub-cause dominates would help diagnosis.
