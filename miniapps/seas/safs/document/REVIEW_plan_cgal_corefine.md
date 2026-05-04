# Plan Review: PLAN_cgal_corefine.md (2026-04-29)

## Review Scope

- **Plan reviewed:** `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/safs/PLAN_cgal_corefine.md`
- **Cross-referenced:**
  - `mesh/conformalize_faults.py` (the legacy implementation whose schema the plan claims to inherit)
  - `mesh/generate_safs_mesh.py` (downstream consumer of `triangle_to_fault.json`)
  - `mesh/write_fault_provenance.py` (claimed consumer; actually KDTree-based, does NOT read `triangle_to_fault.json`)
  - `mesh/validate_msh.py` (the 11/11 acceptance suite the plan promises)
  - `mesh/output/two_crossing_2000m/output/validation_report.txt` (current 7/11 baseline)
  - `mesh/output/two_crossing_2000m/stl_conformal/{triangle_to_fault.json,intersection_report.json}` (existing schema files)
  - `mesh/run_two_crossing_2000m.sh` (script being parallel-cloned and later modified)
  - CGAL Polygon Mesh Processing 5.x documentation (corefine, isotropic_remeshing, polygon-soup repair)
- **Domain context:** `seas-mfem/CLAUDE.md`, `miniapps/seas/CLAUDE.md`, project memories on Frontera approval and module-list discipline.
- **Output written to:** *this file* —
  `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/safs/REVIEW_plan_cgal_corefine.md`. The pre-existing
  `REVIEW.md` (the Python-pipeline review) is left intact.

This is an **adversarial plan review**, not a code review (no
binary exists yet). The findings below identify places where
the plan would mislead the implementer, omit a critical step,
overstate a guarantee, or set an acceptance bar that is not
achievable from the work the plan describes. I assumed the
plan contains at least 3 unstated bugs; I found 12 actionable
issues (4 critical, 6 moderate, 2 low).

---

## Findings

### [P-001] [CRITICAL] [Phase 2 §1–2] — `protect_constraints(true)` blocks splitting long polyline edges; output will have permanent slivers along the intersection

**Category:** ASSUMPTION (silent technical-debt; will fail acceptance)

**Description:**
The plan specifies (Phase 2 §2):
```cpp
PMP::isotropic_remeshing(
    faces(M), p.target_edge_m, M,
    params::edge_is_constrained_map(ecm)
           .number_of_iterations(p.n_iterations)
           .protect_constraints(true)
           .relax_constraints(false));
```

CGAL docs state: *"If `protect_constraints` is `true`, constraints
are protected from being collapsed AND split. Their length
remains unchanged."* The Mill Creek × SBMT-SAF intersection has
**581 polyline vertices over the polyline length** (per the
existing `intersection_report.json`); average constrained-edge
length is therefore on the order of the **CFM raw triangulation
edge** (~1500 m) — not the `target_edge_m = 1000 m` the plan
asks for.

After remeshing with `protect_constraints(true)`:
- Non-constrained edges → ~1000 m (target).
- Constrained polyline edges → still ~1500 m (untouched).
- Triangles bordering the polyline have a **1.5:1 aspect ratio**
  (best case) and aspect ratios degrade further when the polyline
  has any geometric kink because the kink edges are even shorter
  in 3-D after corefine.

This is exactly the sliver problem the plan claims to solve.

**Trigger:** Any input where the polyline-as-emitted-by-corefine
has edges that differ from `target_edge_m` by more than ~10%.
Mill Creek × SBMT-SAF hits this on every polyline edge.

**Actual behavior:** Phase 2 acceptance "mean edge in [800,
1200] m" passes (the *mean* is fine because 99% of edges are
non-constrained), but local mesh quality near the polyline is
exactly what corefine left behind — slivers. Validator check
6 may pass on mean but check 11 (tube uniformity, p95 metric)
may fail.

**Expected behavior:** Pre-split long constrained edges to
`target_edge_m` BEFORE the isotropic_remeshing call, so all
constrained edges are within the target band.

**Suggested fix:**
Insert a `PMP::split_long_edges` step in the remesh pipeline:

```diff
@@ Phase 2 §2 implementation
+    // P-001: split any constrained edge longer than target so
+    // protect_constraints(true) doesn't lock in slivers.  CGAL
+    // marks new sub-edges as constrained automatically when the
+    // original was constrained — see CGAL/PMP/repair.h docs.
+    {
+      // Build a filtered edge range containing only constrained edges.
+      std::vector<edge_descriptor> long_constrained;
+      for (auto e : edges(M))
+          if (get(ecm, e)) long_constrained.push_back(e);
+      PMP::split_long_edges(long_constrained, p.target_edge_m, M,
+          params::edge_is_constrained_map(ecm));
+    }
+
     PMP::isotropic_remeshing(
         faces(M), p.target_edge_m, M,
         params::edge_is_constrained_map(ecm)
                .number_of_iterations(p.n_iterations)
                .protect_constraints(true)
                .relax_constraints(false));
```

And update Phase 2 acceptance criteria with **min-edge** and
**aspect-ratio** bounds that catch slivers:

```diff
- [ ] After running `corefine_faults` with default flags
   (`--target-edge-m 1000`) on Mill Creek × SBMT-SAF 2000 m, every
   output STL has mean edge length in [800, 1200] m.
+ [ ] After running with default flags, every output STL has:
+      - mean edge length in [800, 1200] m,
+      - **min edge length ≥ 250 m** (= target_edge / 4 — catches
+        polyline-induced slivers),
+      - **max aspect ratio ≤ 4** on triangles incident to a
+        constrained edge (= median_edge / shortest_edge per tri).
```

**Test case:**
```cpp
// tools/tests/test_split_long_constrained.cpp
TEST(P001) {
  // Two 1 m squares that intersect in a single 1 m edge.
  // The intersection edge is 1 m long; target_edge = 0.2 m.
  Mesh A, B; build_two_squares(A, B);
  EdgeConstrainedMap ecmA = ..., ecmB = ...;
  PMP::corefine(A, B, params::edge_is_constrained_map(ecmA),
                       params::edge_is_constrained_map(ecmB));
  remesh_one_fault(A, ecmA, {.target_edge_m = 0.2});
  // After fix: the polyline has been split into ~5 sub-edges of
  // ~0.2 m each.  Without fix: it remains 1 edge of 1 m.
  std::int64_t n_constrained_post = 0;
  double max_constrained_len = 0;
  for (auto e : edges(A))
      if (get(ecmA, e)) {
          ++n_constrained_post;
          max_constrained_len = std::max(max_constrained_len,
              CGAL::sqrt(squared_length(A, e)));
      }
  ASSERT_GE(n_constrained_post, 4);    // split happened
  ASSERT_LE(max_constrained_len, 0.3); // ≤ 1.5 × target
}
```

---

### [P-002] [CRITICAL] [Phase 2 §4] — Per-pair ECM accumulation pattern is wrong; CGAL already propagates marks through splits

**Category:** ASSUMPTION (the plan's recipe will lose constraint marks)

**Description:**
Phase 2 §4 says:
> "Per-pair `EdgeConstrainedMap` for (A, C) only marks the edges
>  from corefining with C. The edges from (A, B) are no longer
>  marked in (A, C)'s ECM. Phase 2 says 'accumulate per-pair
>  ECMs into a per-fault ECM' — but the per-pair maps are CGAL
>  Property_maps tied to the mesh's lifetime."

The proposed `accumulate_ecm` function tries to union per-pair
ECMs into a per-fault ECM. This is the **wrong recipe**.

The correct CGAL pattern (from `corefinement_difference_remeshed.cpp`,
the example REVIEW.md R-010 was modeled on):
1. Install **one** ECM on each fault BEFORE the first corefine.
2. Pass the SAME ECM to every `PMP::corefine` call that fault
   participates in. CGAL accumulates: when a previously-marked
   edge is split, **both sub-edges are marked**; new corefine
   intersections add new marks; pre-existing marks are preserved.
3. After all corefines, pass that single accumulated ECM to
   `PMP::isotropic_remeshing`.

The plan's per-pair-then-accumulate pattern fails because:
- A per-pair ECM created AFTER (A, B) is run won't have marks
  for edges that already existed; only the ones (A, B) added.
- An ECM passed to (A, B) is mutated by (A, B). When (A, C)
  later subdivides edges of A that came from (A, B), the (A, B)
  ECM IS updated (marks propagate through splits) — but only if
  the same ECM is passed to (A, C). The plan's separate-ECMs
  recipe doesn't do that.

**Trigger:** Any 3+-fault input where one fault participates in
more than one pair. Phase 5's all-8 build hits this on every
fault.

**Actual behavior under the plan's recipe:** Faults that
participate in multiple pairs lose the constrained marks from
all pairs except the most recent. Remeshing relaxes the older
polylines back into ordinary edges; downstream fault topology
is destroyed.

**Expected behavior:** One ECM per fault, shared across all
corefine calls.

**Suggested fix:**

```diff
@@ Phase 2 §4 — replace the "accumulate per-pair ECMs" recipe
- 4. **Pairwise vs all-at-once**: corefine each pair (Phase 1
-    loop), then **after all pairs are corefined**, run
-    `remesh_one_fault` once per fault, accumulating constrained
-    edges from all pairs the fault participated in. The Phase 1
-    `EdgeConstrainedMap` per-pair gets merged into a per-fault
-    `EdgeConstrainedMap`:
-    ```cpp
-    // Per fault, union all per-pair ecms into one master ecm.
-    EdgeConstrainedMap accumulate_ecm(
-        Mesh& M, const std::vector<EdgeConstrainedMap>& per_pair);
-    ```
-    This is necessary because remeshing fault A protects only
-    the edges that came from corefining A with whichever pair
-    it was last in — without the union, an edge between A and C
-    could be relaxed when remeshing because A's last corefine was
-    with B.
+ 4. **One ECM per fault, shared across all corefine calls**:
+    Install a single `EdgeConstrainedMap` on each fault BEFORE
+    any corefine runs:
+    ```cpp
+    auto ecm = M.add_property_map<edge_descriptor, bool>(
+                   "e:is_constrained", false).first;
+    ```
+    Pass the SAME `ecm` to every `PMP::corefine` call that fault
+    participates in.  CGAL guarantees:
+      (i) Marks already set on an edge are preserved when that
+          edge is split during a subsequent corefine — both
+          sub-edges inherit the mark.
+      (ii) New intersection edges from later corefines are
+           marked True alongside the older marks.
+    No `accumulate_ecm` step is required — and writing one would
+    be wrong, because per-pair-local ECMs lose marks on edges
+    that the per-pair call didn't see.
+
+    Concrete loop in main.cpp:
+    ```cpp
+    std::vector<EdgeConstrainedMap> ecms(N);
+    for (std::size_t i = 0; i < N; ++i)
+        ecms[i] = meshes[i].add_property_map<edge_descriptor, bool>(
+                       "e:is_constrained", false).first;
+    for (std::size_t i = 0; i < N; ++i)
+        for (std::size_t j = i + 1; j < N; ++j)
+            corefine_pair(meshes[i], meshes[j], ecms[i], ecms[j], ...);
+    for (std::size_t i = 0; i < N; ++i)
+        remesh_one_fault(meshes[i], ecms[i], remesh_params);
+    ```
```

**Test case:**
```cpp
// tools/tests/test_three_fault_cascade.cpp
TEST(P002) {
  // Three 1 m squares: A in xy-plane, B in xz-plane, C in yz-plane.
  // All pairwise intersections are along coordinate axes.
  Mesh A, B, C;  build_three_orthogonal_squares(A, B, C);
  auto ecmA = A.add_property_map<edge_descriptor, bool>("e:c", false).first;
  auto ecmB = B.add_property_map<edge_descriptor, bool>("e:c", false).first;
  auto ecmC = C.add_property_map<edge_descriptor, bool>("e:c", false).first;

  PMP::corefine(A, B, params::edge_is_constrained_map(ecmA),
                       params::edge_is_constrained_map(ecmB));
  PMP::corefine(A, C, params::edge_is_constrained_map(ecmA),
                       params::edge_is_constrained_map(ecmC));

  // A now participates in TWO polylines: A∩B (x-axis) and A∩C (y-axis).
  // The shared ECM ecmA must mark BOTH polylines.
  std::int64_t n_marked_A = 0;
  for (auto e : edges(A)) if (get(ecmA, e)) ++n_marked_A;
  // x-axis polyline: ~1 edge after corefine before P-001 split.
  // y-axis polyline: ~1 edge.
  // So at least 2 marked edges.
  ASSERT_GE(n_marked_A, 2);
}
```

---

### [P-003] [CRITICAL] [Phase 3 acceptance criteria] — "11/11 validator checks" is unachievable from corefine alone; checks 7 and 11 depend on gmsh size-field tuning that the plan does not modify

**Category:** DEVIATION (plan over-promises)

**Description:**
Phase 3 acceptance says:
> `validate_msh.py` reports **11/11 checks passed** on the CGAL
> pipeline output.

But of the 4 checks currently failing (5, 6, 7, 11), only checks
**5** (`internal_interface`) and **6** (`fault_edge_length`)
are influenced by the corefine output. Checks **7**
(`far_field_edge`) and **11** (`tube_uniformity`) are functions
of the **gmsh `safs.geo` size field** (Distance + Threshold
fields), not the input STL geometry.

Concretely from the current `validation_report.txt`:
- check 7: mean far-field tet edge **6466 m** vs threshold
  **12000 m** (= 0.6 × `res_ff` = 0.6 × 20000). The far field
  is being meshed too fine because the `algo3d=4` fallback
  silently disables the embedding constraint and the size
  field doesn't ramp up enough away from the (effectively
  ignored) fault.
- check 11: mean in-tube tet edge **1814 m** vs target band
  **[500, 1500] m** (band centred on `res_f = 1000 m`). The
  tube isn't being honoured because the fault isn't actually
  embedded.

With `algo3d=10` (HXT) and the corefine output, **both checks
might flip to PASS for free** — because HXT actually embeds
the fault, the size field becomes effective, far tets reach
`res_ff`, and tube tets stay near `res_f`. That's the
**hopeful** outcome.

But it is not **guaranteed** by anything the plan does. The
plan must either (a) demote the acceptance bar to **9/11**
(checks 5, 6 are guaranteed; 7, 11 are stretch) and document
size-field tuning as a follow-up, OR (b) commit to actually
verifying this on the laptop and including a falsifiable
size-field acceptance.

**Trigger:** Phase 3 implementation completes and the developer
runs `validate_msh.py` on the CGAL output. Even with all of
Phase 1 + 2 perfect, checks 7 / 11 may still report FAIL if the
tube_radius+ramp_dist parameters interact poorly with the new
post-remesh fault edge distribution.

**Actual behavior:** Phase 3 is reported "complete" with 9/11
checks passed; the user thinks the plan failed.

**Expected behavior:** The plan distinguishes "must-pass"
checks (5, 6, 8 → corefine guarantees these) from "stretch"
checks (7, 11 → require size-field validation post-CGAL).

**Suggested fix:**

```diff
@@ Phase 3 acceptance criteria
- [ ] `bash mesh/run_two_crossing_2000m_cgal.sh` exits 0.
- [ ] `validate_msh.py` reports **11/11 checks passed** on the
-   CGAL pipeline output. In particular:
+ [ ] `bash mesh/run_two_crossing_2000m_cgal.sh` exits 0.
+ [ ] `validate_msh.py` reports **at least 9/11 checks passed**
+   (the 4 checks that must pass for SEAS usability, plus the
+   5 that already pass with the Python pipeline):
   - check 5 (`internal_interface`): `n_offending == 0`,
     `n_internal_with_2_tets ≈ n_fault_tris - n_trace_tris`.
+    [MUST PASS — this is the proof the corefine fix worked]
   - check 6 (`fault_edge_length`): mean edge in
     [500, 2000] m.
+    [MUST PASS — this is the proof remeshing worked]
+ - [ ] If checks 7 and 11 still fail, the failure is attributed to
+       the gmsh size-field (`safs.geo`'s Threshold field with
+       `tube_radius` + `ramp_dist`), NOT the corefine output.
+       Phase 3 is conditionally complete; size-field tuning is
+       a follow-up ticket.  The diagnostic: with
+       `validate_msh.py` reporting `mean_far_edge_m << 0.6 *
+       res_ff`, increase `--ramp-dist` or reduce `--tube-radius`
+       to extend the size-ramp into the far field.
+ - [ ] **Stretch goal — 11/11**: with one round of size-field
+       tuning (no code changes; CLI flag changes only), all 11
+       checks pass.  Exact `tube_radius` / `ramp_dist` /
+       `res_ff` values that achieve 11/11 are recorded in
+       `output/two_crossing_2000m_cgal/sizing.json` and pinned
+       in the run script as the new defaults.
```

**Test case:**
```python
def test_P003_corefine_owns_check5_and_6_only():
    """The CGAL plan's hard guarantees: check 5 and check 6
    flip from FAIL to PASS.  Checks 7 and 11 are stretch and
    depend on size-field tuning, not corefine output."""
    report = parse_validation_report(
        "output/two_crossing_2000m_cgal/output/validation_report.txt")
    assert report.check_5.passed, (
        "CGAL plan must deliver check 5 PASS; this is the "
        "primary acceptance bar.")
    assert report.check_6.passed, (
        "CGAL plan must deliver check 6 PASS via remeshing.")
    # Checks 7 and 11 are not asserted to pass at the corefine
    # acceptance gate.
    if not report.check_7.passed or not report.check_11.passed:
        sizing = json.loads(Path(
            "output/two_crossing_2000m_cgal/output/sizing.json"
        ).read_text())
        # Diagnostic must be informative for the size-field
        # follow-up ticket.
        print(f"Stretch checks not yet passing; size-field "
              f"settings: {sizing}")
```

---

### [P-004] [CRITICAL] [Phase 1 §2] — Polygon-soup → mesh recipe is incomplete; will leave duplicate vertices that make corefine non-deterministic

**Category:** BUG (in the plan's spec)

**Description:**
Phase 1 §2 says:
> `read_ascii_stl` uses CGAL's
> `CGAL::IO::read_polygon_soup` + `orient_polygon_soup_with_*`
> to handle non-manifold input from `ts_to_stl.py`.

This is **incomplete**. The standard CGAL recipe for STL → mesh
that produces a manifold `Surface_mesh` is **four steps**:

```cpp
std::vector<Point_3> points;
std::vector<std::vector<std::size_t>> polygons;
CGAL::IO::read_polygon_soup(path, points, polygons);
PMP::repair_polygon_soup(points, polygons);   // dedup, drop degenerate
PMP::orient_polygon_soup(points, polygons);
PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);
```

The plan misses both `repair_polygon_soup` and
`polygon_soup_to_polygon_mesh`, and gets the orient name wrong
(`orient_polygon_soup_with_*` is not a CGAL function — possibly
a confusion with `orient_polygon_soup_extension`, but the
plain `orient_polygon_soup` is what's needed).

**Without `repair_polygon_soup`:** ASCII STL stores three vertex
slots per triangle. After `read_polygon_soup`, every shared
vertex appears K times in `points` (once per incident
triangle). `polygon_soup_to_polygon_mesh` will reject the
input as non-manifold; even if it accepts, `corefine` on the
resulting mesh has undefined behaviour because half-edges
don't connect properly across "shared" vertices that are
actually distinct point indices.

**Trigger:** Reading any STL written by `ts_to_stl.py` (which
is every input the tool will see).

**Actual behavior under the plan's recipe:** `polygon_soup_to_polygon_mesh`
fails on the first STL; tool exits with a CGAL-internal error
message that doesn't point to the actual cause (missing
`repair_polygon_soup` step).

**Expected behavior:** The four-step recipe above. The
`repair_polygon_soup` step also handles the few degenerate
triangles `ts_to_stl.py` may emit (those are the same ones
the existing pipeline drops at line 230–242 of
`generate_safs_mesh.py`).

**Suggested fix:**

```diff
@@ Phase 1 §2 — STL I/O
2. **STL I/O** (`io.cpp`):
   ```cpp
   namespace safs::io {
   bool read_ascii_stl(const std::filesystem::path& p,
                       CGAL::Surface_mesh<K::Point_3>& out_mesh,
                       std::string* out_solid_name = nullptr);
   bool write_ascii_stl(const std::filesystem::path& p,
                        const CGAL::Surface_mesh<K::Point_3>& mesh,
                        std::string_view solid_name);
   } // namespace safs::io
   ```
-   - `read_ascii_stl` uses CGAL's
-     `CGAL::IO::read_polygon_soup` + `orient_polygon_soup_with_*`
-     to handle non-manifold input from `ts_to_stl.py`. Reject if
-     post-orient mesh has > 0 non-manifold edges; report which
-     fault failed.
+   - `read_ascii_stl` uses the standard CGAL polygon-soup recipe
+     to recover manifold connectivity from STL's per-triangle
+     vertex slots:
+     ```cpp
+     std::vector<K::Point_3> points;
+     std::vector<std::vector<std::size_t>> polygons;
+     if (!CGAL::IO::read_polygon_soup(p.string(), points, polygons))
+         return false;
+     PMP::repair_polygon_soup(points, polygons,
+         params::erase_all_duplicates(true)
+                .require_same_orientation(false));
+     PMP::orient_polygon_soup(points, polygons);
+     if (!PMP::is_polygon_soup_a_polygon_mesh(polygons)) {
+         std::cerr << "[corefine_faults] " << p
+                   << ": polygon soup is non-manifold even after "
+                   << "repair; ts_to_stl.py output is corrupt.\n";
+         return false;
+     }
+     PMP::polygon_soup_to_polygon_mesh(points, polygons, out_mesh);
+     return true;
+     ```
+     Required CGAL headers:
+     `<CGAL/Polygon_mesh_processing/repair_polygon_soup.h>`,
+     `<CGAL/Polygon_mesh_processing/orient_polygon_soup.h>`,
+     `<CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>`.
+     Reject if post-recipe mesh has > 0 non-manifold edges
+     (`PMP::is_polygon_mesh(out_mesh) == false`).
```

**Test case:**
```cpp
TEST(P004) {
  // ASCII STL with two adjacent triangles sharing an edge.
  // STL stores 6 vertices total (3 per triangle); after recipe
  // the Surface_mesh must report exactly 4 vertices (the shared
  // edge's two vertices appear once each).
  Mesh M;
  ASSERT_TRUE(safs::io::read_ascii_stl(fixture_path("two_tris.stl"), M));
  ASSERT_EQ(num_vertices(M), 4);
  ASSERT_EQ(num_faces(M), 2);
}
```

---

### [P-005] [MODERATE] [Constraints / Interface constraints + Phase 3] — Plan over-states which downstream tools consume `triangle_to_fault.json`

**Category:** ASSUMPTION (plan factual error)

**Description:**
The plan says:
> `triangle_to_fault.json` schema (consumed by
> `generate_safs_mesh.py:_combine_stls` line 445 and by
> `write_fault_provenance.py`)

Verified by reading both consumers:
- `generate_safs_mesh.py:445` only checks
  `(args.stl_dir / "triangle_to_fault.json").exists()`.  It does
  **not** read any keys from the JSON.  The file is a presence
  flag, nothing more.
- `write_fault_provenance.py` does **not** open
  `triangle_to_fault.json` at all.  It uses **KDTree centroid
  matching** against the per-fault STLs (read via
  `_read_stl_centroids`, line 27) — independent of any JSON
  produced by the conformalize step.

**Why this matters:** The plan elevates the schema to a contract
that constrains the implementation. It's not a contract — it's
just a marker file. The implementer who treats it as a
byte-equivalence requirement will spend time matching keys that
nobody reads.

**Trigger:** Any reader of the plan trying to estimate how
flexible the JSON output can be.

**Suggested fix:**

```diff
@@ Constraints / Interface constraints
- - **JSON schemas** consumed by downstream tools must match exactly:
-   - `triangle_to_fault.json` schema (consumed by
-     `generate_safs_mesh.py:_combine_stls` line 445 and by
-     `write_fault_provenance.py`):
+ - **JSON schemas** the new tool emits:
+   - `triangle_to_fault.json` is a **presence-only marker**
+     (`generate_safs_mesh.py:445` only calls `.exists()` on it;
+     `write_fault_provenance.py` ignores it and uses KDTree
+     centroid matching on the STLs).  We still emit it for
+     diagnostic and human-readability reasons, with the same
+     schema the Python conformalizer used:
     ```json
     {
       "schema_version": 1,
       "faults": {
         "<short_name>": {"n_triangles": <int>, "range": [<int>, <int>]}
       },
       "n_total_triangles": <int>
     }
     ```
+    No downstream tool reads these values; correctness is
+    diagnostic only.
@@
-  - `intersection_report.json` schema (consumed by humans + the
-    R-001 fix's hard-fail logic):
+  - `intersection_report.json` is a **diagnostic-only** report
+    consumed by humans.  It is NOT consumed by any of
+    `generate_safs_mesh.py`, `write_fault_provenance.py`,
+    `validate_msh.py`, or `convert_msh.py`.  The schema below
+    is for documentation continuity with the Python output —
+    the implementer is free to add or omit keys with the
+    constraint that the **per-fault** and **pairs** top-level
+    keys, plus per-pair `gates.{manifold_A,manifold_B,
+    polyline_edge_coincidence,interior_crossing_only}` strings,
+    are preserved (so existing dashboards / tail-of-log greps
+    continue to work).
```

**Test case:**
```python
def test_P005_no_consumer_reads_triangle_to_fault_keys():
    """No downstream tool reads triangle_to_fault.json keys —
    only its existence.  Confirm by static grep."""
    import re
    consumers = [
        "mesh/generate_safs_mesh.py",
        "mesh/write_fault_provenance.py",
        "mesh/validate_msh.py",
        "mesh/convert_msh.py",
    ]
    for c in consumers:
        src = Path(c).read_text()
        if "triangle_to_fault" not in src:
            continue
        # Allow .exists(); forbid .read_text() / json.load(...t2f).
        assert ".read_text()" not in src or "triangle_to_fault" not in (
            src[:src.rfind(".read_text()")].splitlines()[-1]
        ), f"{c} reads triangle_to_fault.json content"
```

---

### [P-006] [MODERATE] [Phase 1 acceptance + Phase 2 acceptance] — Acceptance criteria use mean-only metrics that hide local sliver damage

**Category:** EDGE_CASE

**Description:**
Phase 2 acceptance: "every output STL has mean edge length in
[800, 1200] m." This is a global mean; a fault with 99% of
edges at exactly 1000 m and 1% of edges at 10 m (slivers along
the polyline) passes this criterion.

Validator check 6 in `validate_msh.py:323` reports `mean_m`,
`max_m`, `p50_m`, `p95_m`. The plan should propagate these to
the C++ tool's acceptance.

**Trigger:** Any sliver pattern that survives Phase 2 — for
instance, a polyline pierce that lands at sub-target distance
from a parent vertex (CGAL doesn't unilaterally collapse to the
vertex; it inserts a new node).

**Suggested fix:**

```diff
@@ Phase 2 §Acceptance criteria
- [ ] After running `corefine_faults` with default flags
   (`--target-edge-m 1000`) on Mill Creek × SBMT-SAF 2000 m, every
   output STL has mean edge length in [800, 1200] m.
+ [ ] After running with default flags on Mill Creek × SBMT-SAF
+   2000 m, every output STL satisfies:
+      - mean edge length in [800, 1200] m
+      - p5 edge length ≥ 250 m   (= target / 4 — sliver gate)
+      - p95 edge length ≤ 2500 m (= 2.5 × target — long-tail gate)
+      - max aspect ratio ≤ 4 on any triangle incident to a
+        constrained edge (max-edge / min-edge per triangle).
```

**Test case:** the C++ smoke test in Phase 2 should compute
these statistics and assert them, not just the mean.

---

### [P-007] [MODERATE] [Phase 1 §3] — `intersection_report.json` schema in plan deviates from existing Python schema without migration notes

**Category:** DEVIATION (silently changes the on-disk format)

**Description:**
The plan emits a schema with these keys that the **Python output
does not have**: `backend`, `clearance_m`,
`n_constrained_edges_A`, `n_constrained_edges_B`, `remesh_iters`.
And **omits** these keys that the Python output has:
`drop_short_segment_frac`, `min_pierce_separation_m`,
`cdt_min_angle_deg`, `n_polylines`, `n_polyline_vertices`,
`n_dropped_short_A`, `n_dropped_short_B`, `smoothing_iters`.

The plan claims "byte-equivalent schemas to Python" (Phase 1
§Interfaces). That is false. The schema is intentionally
different to reflect the different backend.

**Trigger:** Any human or script that compares the two reports
side-by-side will see a key-set diff and may interpret it as
data loss.

**Suggested fix:** Add an explicit migration table and call out
the schema delta as a *deliberate* (not accidental) change.

```diff
@@ Constraints / Interface constraints
+   - **Schema delta vs Python output (deliberate):**
+
+     | Key                       | Python | CGAL | Notes |
+     |---|---|---|---|
+     | `schema_version`          | 1      | 1    | Bumped to 2 in a future ticket if needed |
+     | `backend`                 | absent | "cgal" | NEW — disambiguates output |
+     | `snap_m`                  | yes    | yes  | unchanged semantics |
+     | `clearance_m`             | absent | yes  | NEW — record clamp value |
+     | `target_edge_length_m`    | yes    | yes  | unchanged |
+     | `drop_short_segment_frac` | yes    | absent | n/a — CGAL handles via remeshing |
+     | `min_pierce_separation_m` | yes    | absent | n/a — corefine produces canonical pierces |
+     | `cdt_min_angle_deg`       | yes    | absent | n/a — no CDT step |
+     | `pairs.<x>.n_polylines`   | yes    | absent | replaced by `n_constrained_edges_*` |
+     | `pairs.<x>.n_polyline_vertices` | yes | absent | replaced by `n_constrained_edges_*` |
+     | `pairs.<x>.n_constrained_edges_{A,B}` | absent | yes | NEW — direct CGAL count |
+     | `pairs.<x>.smoothing_iters` | yes  | absent | replaced by `remesh_iters` |
+     | `pairs.<x>.remesh_iters`  | absent | yes  | NEW — CGAL remesh iteration count |
+     | `pairs.<x>.n_dropped_short_{A,B}` | yes | absent | n/a — corefine doesn't drop |
+     | `gates`                   | yes    | yes  | values: `"PASS"` always (CGAL) vs `"PASS_WITH_<n>_T_JUNCTIONS"` (Python) |
+     | `per_fault.<short>.{n_vertices,n_triangles,stl_path}` | yes | yes | unchanged |
+
+     The `intersection_report.json` is diagnostic-only (P-005),
+     so this delta does not break any downstream tool.  It IS
+     visible to humans diffing the two outputs; document the
+     delta in the README so operators don't waste time
+     reconciling it.
```

**Test case:**
```python
def test_P007_intersection_report_keysets_documented():
    """Ensure the documented schema delta in PLAN_cgal_corefine.md
    matches what the C++ tool actually emits."""
    py_out = json.loads((PY_OUTPUT_DIR / "intersection_report.json").read_text())
    cgal_out = json.loads((CGAL_OUTPUT_DIR / "intersection_report.json").read_text())
    py_keys = set(py_out.keys()) | {"top.<no-pair-keys>"}
    cgal_keys = set(cgal_out.keys()) | {"top.<no-pair-keys>"}
    expected_only_cgal = {"backend", "clearance_m"}
    expected_only_python = {
        "drop_short_segment_frac", "min_pierce_separation_m",
        "cdt_min_angle_deg",
    }
    assert (cgal_keys - py_keys) == expected_only_cgal, (
        f"unexpected new keys: {(cgal_keys - py_keys) - expected_only_cgal}")
    assert (py_keys - cgal_keys) == expected_only_python, (
        f"unexpected removed keys: {(py_keys - cgal_keys) - expected_only_python}")
```

---

### [P-008] [MODERATE] [Phase 1 §6 + Phase 2 §3] — Free-surface z-clamp can invert triangles or pull constrained vertices off the polyline

**Category:** EDGE_CASE

**Description:**
Phase 2 §3 says:
> After `isotropic_remeshing`, walk every vertex; if `v.z() >
> -p.clearance_m`, project it down to `z = -p.clearance_m` in
> place.

Two failure modes the plan doesn't address:
1. **Constrained-vertex motion.** A polyline endpoint that lies
   at `z = -50 m` (above the `z = -100 m` clearance plane)
   would be pulled to `z = -100 m`, which:
   - Breaks polyline edge-coincidence with the OTHER fault that
     shares the same polyline endpoint (CGAL guarantees both
     faults' polyline vertices are bit-identical *immediately
     after corefine* — a unilateral z-clamp on fault A breaks
     this).
   - Moves a constrained vertex despite `relax_constraints(false)`,
     defeating its purpose.
2. **Triangle inversion.** Pulling vertex `v` from `z = -50` to
   `z = -100` while its 1-ring neighbours are at `z = -110, z =
   -90` produces a triangle whose normal flips. Validator check
   1 (tag inventory) won't catch this; downstream (HXT, the SEAS
   solver) may produce wrong fluxes silently.

The input STL from `ts_to_stl.py` is **already** clamped to
`z ≤ -clearance_m`, and `isotropic_remeshing` with
`protect_constraints(true)` does not move boundary vertices. So
the post-remesh clamp is, in practice, a no-op — but it's a
silent corruption hazard if any of the assumptions above breaks.

**Trigger:** An input where (a) the user passes `--clearance-m`
larger than what `ts_to_stl.py` used, OR (b) CGAL's remesher
does in fact pull a non-locked interior vertex up by a few
metres (possible in edge cases of the smoothing iterations).

**Suggested fix:** Replace the z-clamp with an **assert** that
no vertex violates the invariant. If the assert fires, that's
a real bug that should be diagnosed, not silently fixed.

```diff
@@ Phase 2 §3 — Free-surface-clearance ceiling
- 3. **Free-surface-clearance ceiling** (R2-7):
-    - After `isotropic_remeshing`, walk every vertex; if `v.z() >
-      -p.clearance_m`, project it down to `z = -p.clearance_m` in
-      place. (Conservative; preserves the invariant `ts_to_stl.py`
-      already enforces.)
-    - Document: "this is a hard clamp; vertices above the
-      clearance plane are pulled down to it. The remesher's
-      interior-iteration smoothing can drift vertices ~m at most;
-      a clamp avoids R-001-class failures where a fault triangle
-      re-emerges above the free surface."
+ 3. **Free-surface-clearance assertion** (R2-7):
+    - After `isotropic_remeshing`, walk every vertex; if any
+      vertex has `v.z() > -p.clearance_m + 1e-3` (1 mm
+      tolerance for FP round-off), exit 2 with the offending
+      vertex's coordinates and which fault it belongs to.
+    - Rationale: the input STL is already clamped by
+      `ts_to_stl.py`.  `protect_constraints(true)` plus
+      `relax_constraints(false)` keeps boundary and constrained
+      vertices fixed.  Interior-vertex Laplacian smoothing in
+      the remesher cannot pull a vertex above its 1-ring's
+      maximum z, which is bounded by neighbours that are
+      themselves below the clearance plane.  If the assertion
+      fires, that's a CGAL bug or a pathological input — do
+      NOT silently clamp; surface it for diagnosis.
+    - The `--allow-z-clamp` CLI flag is added for the
+      explicit-override case (e.g., a user runs with
+      `--clearance-m` larger than what `ts_to_stl.py` produced).
+      In that mode, the clamp is applied but only after a
+      preflight check that no constrained vertex would be
+      moved.  If a constrained vertex would be moved, exit 1
+      with a "rerun ts_to_stl.py with the larger clearance" message.
```

**Test case:**
```cpp
TEST(P008) {
  // Mesh with one interior vertex at z = -50 m, neighbours at
  // z = -150 m.  clearance = 100.  Without the assertion, the
  // clamp pulls the vertex to z = -100, leaving its neighbours
  // at -150 — flipped triangles.
  Mesh M = build_one_high_vertex();
  RemeshParams p{.target_edge_m = 100, .clearance_m = 100};
  // Under the new assertion semantics: this exits 2 with a clear
  // message naming the offending vertex.
  ASSERT_DEATH(remesh_one_fault(M, ecm, p), "vertex above clearance");
}
```

---

### [P-009] [MODERATE] [Phase 4] — Deleting `run_two_crossing_2000m_cgal.sh` removes the A/B comparison harness; better to invert the relationship

**Category:** QUALITY (rollback-resistance)

**Description:**
Phase 4 says:
> `miniapps/seas/safs/mesh/run_two_crossing_2000m_cgal.sh` —
> delete (now redundant).

This destroys the ability to A/B compare future regressions
against the Python pipeline. If a future code change breaks the
CGAL output, you can no longer easily run the canonical CGAL
recipe; you'd have to read the git history to reconstruct it.

A better pattern: keep `run_two_crossing_2000m_cgal.sh` as the
canonical CGAL recipe, and make `run_two_crossing_2000m.sh`
either (a) a thin wrapper that calls the `_cgal` script, or (b)
a side-by-side runner that produces both outputs.

**Suggested fix:**

```diff
@@ Phase 4 §Files to modify
- `miniapps/seas/safs/mesh/run_two_crossing_2000m_cgal.sh` —
-   delete (now redundant).
+ `miniapps/seas/safs/mesh/run_two_crossing_2000m_cgal.sh` —
+   keep as the canonical CGAL recipe; document in its header
+   that it is the production path.
+ `miniapps/seas/safs/mesh/run_two_crossing_2000m.sh` —
+   either (a) thin wrapper that invokes the `_cgal` script
+   with the canonical OUTPUT_SUFFIX, or (b) renamed to
+   `_legacy.sh` and replaced by a 1-line wrapper that points
+   to the `_cgal` script.  Recommendation: option (b).  The
+   1-line wrapper:
+   ```bash
+   #!/usr/bin/env bash
+   exec "$(dirname "$0")/run_two_crossing_2000m_cgal.sh" "$@"
+   ```
```

**Test case:** N/A (file-existence check is the test).

---

### [P-010] [MODERATE] [Phase 0 + Phase 5] — Frontera CGAL build risk understated; module list and CGAL version not pinned

**Category:** ASSUMPTION

**Description:**
Phase 0 says CGAL is header-only since 5.0 and lists Boost 1.71
as compatible. The Frontera recipe is "TBD: requires
user-approved Frontera login session." Per project memory
`feedback_sbatch_modules.md`, when Frontera build fails the fix
is to copy the working module list from existing sbatch files
verbatim — including `fftw3` and the `LD_LIBRARY_PATH` lines
that drivers commonly omit. The plan does not mention this.

CGAL 5.6.1 (the version cited once in the plan's Frontera
section) requires C++17. Frontera's default `intel/19.1.1`
supports `-std=c++17` but not all C++17 library features (e.g.,
`<filesystem>` requires `-lstdc++fs` linkage on older
toolchains; the plan uses `std::filesystem` extensively).

**Trigger:** Phase 5 hits Frontera, build fails on either C++17
filesystem linkage or a missing module dep, and the developer
has no recipe to fall back on.

**Suggested fix:**

```diff
@@ Phase 0 README §3 (Frontera recipe)
- 3. Frontera recipe placeholder (per project memory, do not run
-    anything on Frontera in this phase; document the **expected**
-    module list — typically `intel/19.1.1`, `boost/1.71`, `gcc/9.1`,
-    plus a manual CGAL header fetch via `git clone --depth 1 -b
-    v5.6.1 https://github.com/CGAL/cgal.git $WORK/cgal-5.6.1` and
-    `-DCGAL_DIR=$WORK/cgal-5.6.1`). Mark with **TBD: requires
-    user-approved Frontera login session per project memory
-    `feedback_frontera_approval.md`**.
+ 3. Frontera recipe (DO NOT EXECUTE in Phase 0 — per project
+    memory `feedback_frontera_approval.md`, Frontera commands
+    require explicit user approval).  Document the **starting
+    point** that the user can run when ready:
+
+    a. **Module list** — copy from a known-working sbatch
+       (per project memory `feedback_sbatch_modules.md`).  The
+       `seas_tpv102` driver currently uses:
+       ```bash
+       module reset
+       module load intel/19.1.1 impi/19.0.9 boost/1.71 \
+                   gcc/9.1.0 cmake/3.24.2 fftw3/3.3.10 \
+                   gmp mpfr
+       export LD_LIBRARY_PATH=$TACC_GMP_LIB:$TACC_MPFR_LIB:$LD_LIBRARY_PATH
+       ```
+       (Verify the exact module versions by reading
+       `safs_origin.py`'s sibling `seas_tpv102` sbatch when the
+       user authorises a login session.)
+
+    b. **CGAL fetch** (CGAL is header-only since 5.0):
+       ```bash
+       cd $WORK
+       git clone --depth 1 -b v5.6.1 https://github.com/CGAL/cgal.git
+       export CGAL_DIR=$WORK/cgal/lib/cmake/CGAL
+       ```
+       Pin **5.6.1**: tested with the `corefinement_difference_remeshed`
+       example referenced in REVIEW.md R-010.  Versions 5.4 and
+       5.5 also work; 6.0 introduces ABI changes and is NOT
+       on the tested list.
+
+    c. **CMake build** with explicit C++17-filesystem link:
+       ```cmake
+       target_link_libraries(corefine_faults PRIVATE stdc++fs)
+       ```
+       (Required only on Frontera's gcc 9.1 / intel 19.1.1;
+       no-op on macOS clang and conda-forge gcc 13+.  Wrap in
+       `if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)` for
+       safety.)
+
+    d. **Acceptance**: the same `corefine_faults_hello` smoke
+       binary builds and prints the CGAL version on a Frontera
+       login node before any compute-node work begins.
+
+    e. **Failure mode**: if the build dies on missing
+       fftw3 / LD_LIBRARY_PATH / Boost — copy the module list
+       VERBATIM from the seas_tpv102 sbatch (per project memory
+       `feedback_sbatch_modules.md`); do NOT trim modules even
+       if they look unrelated.
```

**Test case:**
```bash
# tools/tests/frontera_build_smoke.sh (run only by user, not in CI)
set -euo pipefail
module reset
module load intel/19.1.1 impi/19.0.9 boost/1.71 gcc/9.1.0 \
            cmake/3.24.2 fftw3/3.3.10 gmp mpfr
test -d "$CGAL_DIR" || { echo "CGAL_DIR not set"; exit 1; }
cmake -S miniapps/seas/safs/tools -B build_frontera \
      -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR="$CGAL_DIR"
cmake --build build_frontera --target corefine_faults_hello
./build_frontera/corefine_faults_hello | grep "CGAL version"
```

---

### [P-011] [MODERATE] [Phase 1 §Edge cases] — Coincident-triangle case (fault A and fault B sharing 3 identical vertices) is not specified

**Category:** EDGE_CASE

**Description:**
CGAL's `corefine` documentation states that the input meshes
"do not share a triangle." If fault A's mesh contains a triangle
whose 3 vertex coordinates are bit-identical to a triangle in
fault B's mesh (possible at parallel-overlap regions of adjacent
SAF segments — see REVIEW.md R-009 about the coplanar-skip in
the Python implementation), the result is **undefined**.

The plan does not say what `corefine_pair` should do.

**Trigger:** Phase 5 all-8 input where Mojave-SAF and SBMT-SAF
have triangulation overlap (real CFM data per the implementer's
report).

**Suggested fix:**

```diff
@@ Phase 1 §Edge cases
+ - **Coincident triangles between fault A and fault B**
+   (3 bit-identical vertex coords): pre-flight detect with a
+   hash of `(min(v0, v1, v2), median, max)` over each face;
+   if any pair of faces collides across faults, exit 2 with
+   the per-fault triangle indices and a message: "fault {a}
+   tri {i} and fault {b} tri {j} share all 3 vertices.  CGAL
+   corefine() has undefined behaviour on coincident triangles.
+   Investigate the CFM source (this typically indicates a
+   triangulation that should not have been split into two
+   faults — see PLAN_multifault_intersections.md §Caveats)."
+   Implement in `corefine.cpp` as a pre-corefine pass per pair.
```

**Test case:**
```cpp
TEST(P011) {
  Mesh A = build_one_triangle({0,0,0}, {1,0,0}, {0,1,0});
  Mesh B = build_one_triangle({0,0,0}, {1,0,0}, {0,1,0}); // identical
  EdgeConstrainedMap ecmA = ..., ecmB = ...;
  // Before fix: undefined behaviour — may return 0 segments
  // silently or crash.
  // After fix: exit code 2 with the offending indices.
  ASSERT_EXIT(corefine_pair(A, B, ecmA, ecmB, "A", "B"),
              testing::ExitedWithCode(2),
              "share all 3 vertices");
}
```

---

### [P-012] [LOW] [Phase 3] — `OUTPUT_SUFFIX` modification to existing script is mentioned but not specified at line level

**Category:** QUALITY

**Description:**
Phase 3 §1 says:
> Add `OUTPUT_SUFFIX=_cgal` (defaulted) so the two runs land in
> `output/two_crossing_2000m_cgal/` vs `output/two_crossing_2000m/`.

But `run_two_crossing_2000m.sh:41` is `OUTDIR=output/two_crossing_${RES}m`
hardcoded. The plan should specify the exact line edit.

**Suggested fix:**

```diff
@@ Phase 3 §1 — at the end of the diff
+    The OUTDIR change is also needed inside the cloned script.
+    Line 41 of the existing `run_two_crossing_2000m.sh` is
+    hardcoded:
+    ```diff
+    -OUTDIR=output/two_crossing_${RES}m
+    +OUTDIR=output/two_crossing_${RES}m${OUTPUT_SUFFIX:-}
+    ```
+    The `_cgal` script sets `OUTPUT_SUFFIX=_cgal` near the
+    top.  When Phase 4 lands, the legacy script keeps the
+    empty default and the production output goes to
+    `output/two_crossing_2000m/` as today.
```

---

### [P-013] [LOW] [Phase 4] — Deprecation gating doesn't specify how the legacy tests are kept fresh

**Category:** QUALITY

**Description:**
Phase 4 says: "annotate the test module with `pytest.mark.legacy`;
add a `pytest.ini` filter so the legacy tests run only when
explicitly requested (`pytest -m legacy …`)."

This is fine, but if no CI job ever runs `-m legacy`, the
legacy tests bit-rot. The plan should specify either:
- A weekly cron that runs `-m legacy`, OR
- A retention policy: "legacy tests are removable in 90 days
  if the CGAL path has been the production default that long."

**Suggested fix:**

```diff
@@ Phase 4 §Files to modify
+ - Document in `tests/README.md` (or `pytest.ini` comment) the
+   legacy-test policy:
+   - `pytest -m "not legacy"` (the default) skips legacy tests.
+   - `pytest -m legacy` runs the Python conformalizer regression
+     suite explicitly.  Required to pass at every release.
+   - 90 days after Phase 4 lands, if no CGAL regression has
+     surfaced, the legacy tests + sources are removable in a
+     follow-up ticket.  The clock starts at Phase-4 merge.
```

---

## Summary

- Critical issues: **4** (P-001, P-002, P-003, P-004)
- Moderate issues: **6** (P-005, P-006, P-007, P-008, P-009, P-010, P-011)
- Low issues: **2** (P-012, P-013)
- Plan compliance: **N/A** (this is a plan review, not a code review)
- Verdict: **PASS WITH FIXES** — the plan structure is sound,
  but the 4 critical findings each independently can break the
  Phase-1 acceptance gate (P-004 STL recipe), the Phase-2
  output quality (P-001 split_long_edges, P-002 ECM sharing),
  or the Phase-3 acceptance bar (P-003 11/11 unachievable).
  Apply the suggested-fix diffs before launching `/code-implement`.

## Why this review pushed back hard on a plan you wrote

The previous round's REVIEW.md correctly identified that the
Python pipeline silently failed because nobody pushed back on
its own plan. The CGAL plan is a step up in robustness, but
inherits two patterns from the previous one that this review
catches:

1. **Schema-equivalence claims that were never verified against
   the consumers** (P-005). The previous Python plan claimed
   `triangle_to_fault.json` was consumed by
   `write_fault_provenance.py`; in fact it isn't.
2. **Acceptance criteria padded to look strong without proof
   that they are achievable** (P-003). The previous plan said
   "10/10 validator checks"; the actual baseline was 7/11 and
   no analysis was done to identify which of the failing 4
   were owned by Phase-1+2 vs by gmsh size-field tuning.

Both patterns are addressed in the suggested fixes above.

## Unreviewed Areas

- **Specific CGAL 5.6.1 ABI/API guarantees**: I trusted the
  user-supplied example code (REVIEW.md R-010) and CGAL's
  online docs.  The `protect_constraints(true)` semantics for
  `isotropic_remeshing` are documented as "preserved exactly"
  but if a future CGAL version changes this (and CGAL has
  backwards-incompatible ABI changes between major versions),
  the tool's correctness regresses silently.  Pin the CGAL
  version in CMake (`find_package(CGAL 5.6 EXACT REQUIRED)`)
  and fail loud if the host CGAL is a different version.
- **`PMP::repair_polygon_soup` flag set**: I specified
  `erase_all_duplicates(true)` and `require_same_orientation(false)`
  but did not verify these are the optimal flags for STL input
  produced by `ts_to_stl.py` specifically.  An alternate flag
  set (`merge_duplicate_points(true)` instead of
  `erase_all_duplicates`) may better preserve degenerate
  triangles that `ts_to_stl.py` already handled.  Cross-check
  during Phase 1 implementation.
- **HXT acceptance** of CGAL output on Mill Creek × SBMT-SAF:
  the plan asserts this will work because corefine produces a
  valid PLC; I have not actually built the binary and run HXT
  on its output.  If HXT still rejects (e.g., for adjacency
  reasons unrelated to corefine), Phase 3 must add a
  diagnostic step to surface which gmsh check fired.
- **`tools/compare_meshes.py`**: marked optional in Phase 3.
  I did not specify its API; if it's added, Phase 3 should
  pin its acceptance contract too.
- **Volume meshing with `make_mesh_3` / TetGen** (REVIEW.md
  R-010 §3): explicitly out of scope per the plan, not
  re-evaluated here.
