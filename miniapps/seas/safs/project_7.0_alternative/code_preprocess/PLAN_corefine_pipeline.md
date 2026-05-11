# Implementation Plan: Multi-fault corefine + Delaunay-remesh pipeline with quality control

## Overview

Build a CGAL-based preprocessor that ingests the SAFS multi-fault `.ts`
collection from `project_7.0_preferred/raw_data/`, computes pairwise
conformal corefinements, runs CGAL's surface Delaunay remesher with hard
quality controls (minimum element edge ≥ 100 m, target ~1.5 km, no
needles/caps), and emits a set of conformal STL files that the existing
`code_meshing/safs_fault_box*.geo` style ingests via the proven
`Merge "*.stl"; Surface{tag} In Volume{1};` pattern.

## Locked-in design decisions (from planning session 2026-05-07)

The "Open questions" section at the bottom of the previous draft has
been resolved by the user. The pipeline implements:

1. **Clipping at z = 0 is OUT OF SCOPE for this plan.** The user will
   ship a separately-prepared, already-cleaned `.ts` (or `.off` /
   `.stl`) for each fault. This pipeline takes that input verbatim
   and does not include any z>0 filter or top-edge handling.
   Phase 3 therefore drops `--keep-positive-z` from its CLI; Phase 1
   does no z-direction logic.
2. **Triple junctions: pure pairwise (option i)** as designed in
   Phase 2, with an explicit triple-junction sentinel test in
   Phase 5. The fallback options (ii) "snap pre-extracted triple
   junctions" and (iii) "single autorefine on the union" are *not*
   implemented in this plan; if the Phase 5 sentinel detects a
   T-junction, escalate via a follow-up plan.
3. **Validation: 2000 m only.** Phase 5 fixture is the
   `*_2000m.ts` (or equivalent) set. Higher-resolution validation
   is deferred to a separate plan/run; no `--runlong` gate or 1000 m
   /500 m fixtures in this iteration.

## Answer to the BooleanFragments question (informational, see Phase 4)

The current `code_meshing/` files (`safs_fault_box.geo`,
`safs_fault_box_buried.geo`, `safs_fault_box_freesurface.geo`) **do not
use `BooleanFragments`**. They use `Merge` + `Surface{tag} In Volume{1}`
+ `Mesh.Algorithm3D = 1`. `BooleanFragments` was empirically tested
earlier on a single-fault, single-box configuration and produced
"No tetrahedra in region 1" — the OCC fragmentation operation does not
work on discrete (STL-imported) surfaces; it requires parametric OCC
surfaces. This plan therefore extends the proven `Merge + Surface{} In
Volume{}` path to the multi-fault case. A short Phase 4.5 records the
BooleanFragments retest as a controlled experiment so the file at
`code_meshing/` documents *why* this is not the chosen path.

## Constraints

### Hard constraints (verified by acceptance criteria)
- **Minimum element edge size ≥ 100 m** in every output triangle (surface)
  and every tet (bulk). This is a hard floor; sub-100 m edges trigger
  build failure.
- **Surface triangle quality** `q_tri = 4·√3·A / Σ e_i²` (1 = equilateral,
  0 = degenerate): `min(q_tri) ≥ 0.3` over every output fault.
- **Bulk tet quality** `q_tet = isoperimetric ratio` (per
  `code_preprocess/msh_to_vtu.py`): `min(q_tet) ≥ 0.1`; at least 99.9 %
  of tets satisfy `q_tet ≥ 0.3`.
- **Conformality**: for every pair of corefined faults `(A, B)`, every
  vertex on the shared intersection polyline appears in *both* meshes at
  identical coordinates. Tolerance: bit-identical in OFF I/O,
  ≤ `1e-6 × bbox_diag` after STL ASCII round-trip with `%.10g` writer.

### Interface constraints
- The output of this pipeline is a set of `.stl` files consumable by
  the existing `code_meshing/` `.geo` family. **No change** to the
  `Merge + Surface{} In Volume + Algorithm3D = 1` recipe.
- `code_preprocess/ts_to_stl.py` and `code_preprocess/msh_to_vtu.py`
  are unchanged; the new tools sit alongside them.

### Dependency constraints
- C++ stage built against `cgal-61` conda env: `cgal-cpp 6.1.1`,
  `clang 19.1.7`, `cmake 4.3`, `eigen 3.4`, `boost`.
- Python orchestration uses only the deps already in `pythonenv`:
  `numpy`, `meshio`, `pymeshlab` (optional, for sanity checks),
  `gmsh` (for the `.msh` reader if needed).
- Build artefacts go in `code_preprocess/corefine_cgal/build/`
  (gitignored).

### Convention constraints
- C++ source style: match CGAL examples — `Surface_mesh<Point_3>` mesh
  type, `Epick` kernel, `namespace PMP = CGAL::Polygon_mesh_processing`
  alias.
- File naming: input `*.ts` → intermediate `*.off` (full precision,
  in `code_preprocess/work/`) → output `*_corefined.stl` (high-precision
  ASCII, in `preprocess_data/`).
- Python: argparse CLI, `main() -> int`, prints stats matching the
  existing `ts_to_stl.py` output style.
- Filename for the constraint-edge property map in C++:
  `using ECM = boost::dynamic_edge_property<bool>;`
  `auto ecmA = get(boost::dynamic_edge_property_t<bool>(), A);`

### Numerical constraints
- Kernel: `Exact_predicates_inexact_constructions_kernel` (Epick).
  `corefine` uses lazy exact predicates internally; constructions are
  done in double. This matches what `igl::copyleft::cgal::mesh_boolean`
  uses, so coordinates remain numerically clean.
- Geographic coordinates (UTM ~6-digit easting, 7-digit northing).
  Vertex-coordinate magnitude ~10⁶ m; need ≥10 significant decimal
  digits in any ASCII I/O to keep sub-millimeter precision (otherwise
  conformality breaks during Merge dedup in gmsh).

## Phase 1: Single-pair `corefine_pair` C++ binary, with quality control

### Goal
A standalone C++ binary that takes two `.off` files, applies the
canonical corefine + surface_Delaunay_remesh + degeneracy_repair
pipeline with explicit quality-control parameters, and writes two
`.off` files that are conformal at the shared intersection polyline.
After this phase, on the SBMT-Garnet_Hill × SBMT-Mission_Creek 2000 m
pair, every output triangle has edge length ≥ 100 m.

### Files to Create
- `code_preprocess/corefine_cgal/CMakeLists.txt` — cmake project that
  finds `CGAL` via the `CGAL_DIR` env var (set by `conda activate
  cgal-61`), enables C++17, links `CGAL::CGAL`, `Boost::boost`,
  `Eigen3::Eigen`.
- `code_preprocess/corefine_cgal/corefine_pair.cpp` — CLI entry point
  and corefine-pair driver.
- `code_preprocess/corefine_cgal/polyline_resample.h` — header-only
  polyline resampler (arc-length-uniform).
- `code_preprocess/corefine_cgal/quality_repair.h` — header-only
  quality-control wrapper around `PMP::remove_almost_degenerate_faces`
  + edge-length validation.
- `code_preprocess/corefine_cgal/io_helpers.h` — wrappers for
  `CGAL::IO::read_polygon_mesh` / `write_polygon_mesh` with explicit
  `OFF`/`STL` format detection by extension and a high-precision ASCII
  STL writer (`%.10g`).

### Files to Modify
- None (Phase 1 is purely additive).

### Detailed Requirements

#### CLI
```
corefine_pair INPUT_A.off INPUT_B.off OUT_A.off OUT_B.off
              [--mesh-edge-size SIZE]       (default: 1500.0 m)
              [--min-edge SIZE]             (default: 100.0 m)
              [--polyline-spacing SPACING]  (default: 0.5 * mesh-edge-size; clamped to >= min-edge)
              [--features-angle-bound DEG]  (default: 60.0)
              [--max-iterations N]          (default: 5; for repair pass)
              [--verbose]
```

#### Type aliases (top of `corefine_pair.cpp`)
```cpp
using K        = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point    = K::Point_3;
using Mesh     = CGAL::Surface_mesh<Point>;
using ed       = boost::graph_traits<Mesh>::edge_descriptor;
using vd       = boost::graph_traits<Mesh>::vertex_descriptor;
namespace PMP  = CGAL::Polygon_mesh_processing;
namespace pp   = PMP::parameters;
```

#### `polyline_resample.h` interface
```cpp
namespace polyline_resample {

// Re-sample a single polyline by arc length so that consecutive
// points are exactly `target_spacing` apart along the polyline,
// keeping the first vertex fixed and stopping at the last vertex
// (the final segment may be shorter than `target_spacing`).
//
// If the polyline has length < target_spacing, return [first, last].
//
// For a *closed* polyline (first == last with bit-identical coords),
// keep the first vertex fixed and ensure the last sample wraps back
// to it; the closing edge length will be in [target_spacing/2,
// 3·target_spacing/2].
//
// Vertices on the polyline that are < min_edge_along_polyline apart
// after the initial uniform sampling are merged (averaged); the
// resulting polyline's minimum edge is therefore guaranteed
// >= min_edge_along_polyline.
std::vector<Point_3>
resample(const std::vector<Point_3>& polyline,
         double                      target_spacing,
         double                      min_edge_along_polyline);

// Apply `resample` to every polyline in a set.
std::vector<std::vector<Point_3>>
resample_all(const std::vector<std::vector<Point_3>>& polylines,
             double target_spacing,
             double min_edge_along_polyline);

}  // namespace polyline_resample
```

The arc-length parameterization is the standard
$s_i = \sum_{k=1}^{i} \|p_k - p_{k-1}\|$, $s \in [0, L]$, then sample
at $s = j \cdot \mathrm{target\_spacing}$, $j = 0, 1, \ldots, \lfloor
L / \mathrm{target\_spacing} \rfloor$, plus $s = L$.

#### `quality_repair.h` interface
```cpp
namespace quality_repair {

struct Stats {
    std::size_t n_verts, n_faces;
    double      edge_len_min, edge_len_p1, edge_len_med, edge_len_max;
    double      tri_q_min,    tri_q_p1,    tri_q_med;
    std::size_t n_edges_below_floor;     // edges with length < min_edge_size
    std::size_t n_tri_q_below_threshold; // q_tri < 0.3
};

// Run PMP::remove_almost_degenerate_faces with parameters tuned to
// the user's quality floor; iterate until either nothing changes or
// `max_iterations` is reached.  Returns the post-repair stats.
//
// Guarantees on return:
//   * No interior face has q_tri < q_floor (default 0.3) ... unless
//     the only fix was forbidden by the constrained edge map.
//   * No interior edge has length < min_edge_size ... same caveat.
// (Caveat: constrained edges are never collapsed; if a constrained
//  edge happens to be < min_edge_size, it survives. The polyline
//  resampler in Phase 1 ensures this does not happen for our inputs.)
template <class ECM>
Stats run(Mesh& m,
          ECM   ecm,
          double min_edge_size,
          double q_floor                  = 0.3,
          double cap_threshold_cos        = std::cos(160. * CGAL_PI / 180.),
          double needle_threshold_ratio   = 4.0,
          unsigned max_iterations         = 5);

// Per-edge length statistics (no modification).
Stats compute_stats(const Mesh& m, double min_edge_size, double q_floor);

// Validate post-pipeline quality.  Returns true iff the mesh
// satisfies the user's hard floor (`min_edge_size`, `q_floor`).
template <class ECM>
bool validate(const Mesh& m, ECM ecm,
              double min_edge_size, double q_floor);

}  // namespace quality_repair
```

The triangle quality metric used:
$$q_{tri} = \frac{4\sqrt{3}\, A}{e_a^2 + e_b^2 + e_c^2}$$
This matches `code_preprocess/msh_to_vtu.py`.

#### `corefine_pair.cpp::main` algorithm

1. **Parse CLI** via `<getopt.h>` or hand-rolled; default values as
   documented above.
2. **Read inputs** with `CGAL::IO::read_polygon_mesh(args.in_a, A)` and
   same for B. `read_polygon_mesh` auto-detects format by extension;
   `.off`, `.stl`, `.ply`, `.obj` all work.
3. **Self-intersection guard** for each input:
   ```cpp
   if (PMP::does_self_intersect(A)) {
       if (args.verbose) std::cerr << "A has self-intersections; running PMP::autorefine\n";
       PMP::autorefine(A);
   }
   ```
4. **Compute intersection polylines** (read-only):
   ```cpp
   std::vector<std::vector<Point>> polylines;
   PMP::surface_intersection(A, B, std::back_inserter(polylines));
   ```
   If `polylines.empty()`, write A and B verbatim to outputs and return
   exit code 2 with message "no intersection between A and B".
5. **Resample polylines**:
   ```cpp
   double spacing = std::max(args.polyline_spacing, args.min_edge);
   auto resampled = polyline_resample::resample_all(polylines, spacing,
                                                    args.min_edge);
   ```
6. **Corefine** with constraint maps:
   ```cpp
   auto ecmA = get(boost::dynamic_edge_property_t<bool>(), A);
   auto ecmB = get(boost::dynamic_edge_property_t<bool>(), B);
   PMP::corefine(A, B,
       pp::edge_is_constrained_map(ecmA),
       pp::edge_is_constrained_map(ecmB));
   ```
7. **Surface Delaunay remesh** A and B independently, both seeing the
   *same* resampled polyline as `polyline_constraints`:
   ```cpp
   Mesh A_remeshed = PMP::surface_Delaunay_remeshing<Mesh>(A,
       pp::polyline_constraints(resampled)
         .protect_constraints(true)
         .mesh_edge_size(args.mesh_edge_size)
         .features_angle_bound(args.features_angle_bound));
   Mesh B_remeshed = PMP::surface_Delaunay_remeshing<Mesh>(B,
       pp::polyline_constraints(resampled)
         .protect_constraints(true)
         .mesh_edge_size(args.mesh_edge_size)
         .features_angle_bound(args.features_angle_bound));
   ```
8. **Quality repair** on each remeshed mesh, with the resampled
   polyline marked as constrained (re-derive the constraint map by
   matching vertex coordinates against `resampled`; tolerance 1e-9):
   ```cpp
   auto ecmA2 = mark_polyline_as_constrained(A_remeshed, resampled);
   auto stA = quality_repair::run(A_remeshed, ecmA2, args.min_edge);
   // same for B
   ```
9. **Conformality check**: every vertex of `A_remeshed` whose coords
   match a polyline vertex (tolerance 1e-9) must also exist in
   `B_remeshed` at the same coords. Count the matches; if any
   polyline vertex is missing from either mesh, return exit code 3
   with message "conformality lost during remesh: K of N polyline
   vertices missing".
10. **Hard-floor validation**: each output mesh must pass
    `quality_repair::validate(...)`. If not, return exit code 4 with
    a list of the first 10 offending edges and faces (vertex coords).
11. **Write outputs**:
    - `OUT_A.off` and `OUT_B.off`: full-precision OFF.
    - If outputs end in `.stl`, write a high-precision ASCII STL via
      our custom writer (`%.10g`); CGAL's default STL writer uses
      `<<` with default precision and is too lossy.
12. **Print stats** in the same style as `ts_to_stl.py`:
    ```
    A:  V_in=2350  F_in=4137  V_out=...  F_out=...
        edge:  min=...  median=...  max=...
        q_tri: min=...  p1=...  median=...
    B:  ... same ...
    Polyline: count=K  resampled_count=K'  min_seg=...  max_seg=...
    Conformality: K' / K' polyline vertices match in both meshes.
    Quality floor passed: yes
    ```

#### `CMakeLists.txt`

```cmake
cmake_minimum_required(VERSION 3.20)
project(safs_corefine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(corefine_pair corefine_pair.cpp)
target_link_libraries(corefine_pair PRIVATE CGAL::CGAL Eigen3::Eigen)
target_compile_options(corefine_pair PRIVATE -O2 -Wall -Wno-deprecated-copy)
```

Build invocation (used in CI and locally):
```bash
conda activate cgal-61
cd code_preprocess/corefine_cgal
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j 4
```

### Edge Cases to Handle
- **Empty polyline output** from `surface_intersection`: meshes don't
  intersect. Return exit 2; do not modify inputs.
- **Single-point intersection**: `surface_intersection` may emit a
  zero-length polyline. The resampler must filter polylines with < 2
  points.
- **Closed polyline**: at corefinement junctions where one fault loops
  through another, the polyline may be closed. The resampler must
  detect closure (`first.coords == last.coords` exactly, or distance
  < 1e-9) and produce a closed resampled polyline.
- **Coplanar overlap**: two faults that share a coplanar patch
  produce a 2-D intersection (a disk, not a polyline). `corefine`
  may or may not handle this gracefully. Detect via
  `PMP::do_intersect` returning a polyline with degree-3+ vertices;
  if found, log a warning and abort with exit 5.
- **Polyline vertex appears twice (degree-4 junction)**: where three
  or more faults meet at a single curve, an interior polyline vertex
  has degree > 2. Surface_intersection requires "all vertices but
  endpoints being of degree 2" (precondition documented at
  `intersection.h:1684`). For such inputs, run `surface_intersection`
  on N(N-1)/2 pairs separately and merge the resulting polyline sets;
  Phase 2 handles this dispatch.
- **OFF read failure** (malformed input file): exit 6 with the I/O
  error message.

### Acceptance Criteria
- [ ] Builds cleanly with `cgal-61` conda env (`conda activate cgal-61
      && cmake .. && make`); zero compiler warnings except for
      `-Wno-deprecated-copy`.
- [ ] On synthetic test (two perpendicular 10×10 rectangles
      triangulated at 4×4 each, intersecting at y=5, z=0):
      output meshes have ≥ 1 vertex on the intersection line at every
      multiple of `polyline-spacing`, and the per-mesh per-vertex
      coordinates on the polyline are bit-identical between A and B.
- [ ] On the SBMT-Garnet_Hill × SBMT-Mission_Creek 2000 m pair (largest
      XY overlap): exit 0; both outputs satisfy `min(edge) >= 100 m`,
      `min(q_tri) >= 0.3`, `n_edges_below_floor == 0`.
- [ ] STL ASCII round-trip preserves polyline conformality: writing
      output as ASCII STL with `%.10g`, reading back via `meshio`,
      every polyline vertex coordinate diff between A and B < 1e-6.

### Dependencies
- Depends on: cgal-cpp 6.1.1 (`PMP::corefine`, `PMP::surface_intersection`,
  `PMP::surface_Delaunay_remeshing`, `PMP::remove_almost_degenerate_faces`,
  `PMP::autorefine`, `PMP::does_self_intersect`).
- Required by: Phase 2 (multi-fault driver), Phase 4 (gmsh integration
  consumes the STL output).

## Phase 2: Multi-fault driver `corefine_set` C++ binary

### Goal
A C++ binary that takes a list of N input meshes, identifies which
pairs intersect, processes them in a deterministic order, and outputs
N meshes that are pairwise conformal at every shared intersection.

### Files to Create
- `code_preprocess/corefine_cgal/corefine_set.cpp` — the multi-fault
  driver.
- `code_preprocess/corefine_cgal/intersection_graph.h` — pair-discovery
  helper using `PMP::intersecting_meshes`.

### Files to Modify
- `code_preprocess/corefine_cgal/CMakeLists.txt` — add
  `add_executable(corefine_set corefine_set.cpp)` and the same link
  list.

### Detailed Requirements

#### CLI
```
corefine_set IN_DIR OUT_DIR
             [--mesh-edge-size SIZE]    (default: 1500.0)
             [--min-edge SIZE]          (default: 100.0)
             [--polyline-spacing S]     (default: 0.5*mesh-edge-size, >= min-edge)
             [--ext .off]               (file extension to scan in IN_DIR)
             [--manifest path.json]     (optional; if given, output a JSON manifest)
             [--verbose]
```

The driver scans `IN_DIR` for files matching `*<ext>` and treats each
as an input mesh.

#### Algorithm
1. **Load all N meshes** into a `std::vector<Mesh>`. Record file names
   for output naming.
2. **Bbox-prune pair candidates** with
   `PMP::intersecting_meshes(meshes, std::back_inserter(pairs))`
   (returns a list of pair indices `(i, j)` whose surfaces actually
   intersect — this does the geometric test, not just bbox).
3. **Build the intersection graph** as an undirected graph
   `G = (V, E)` with `V = {0..N-1}` and `E = pairs`. Compute the
   **edge order** via:
   - Greedy heuristic: sort pairs by the **bbox-overlap-volume of the
     pair**, descending. Process the largest intersection first; this
     fixes the most "load-bearing" polylines into the meshes early.
4. **Per-pair processing loop**:
   ```
   For (i, j) in sorted_pairs:
       polylines = surface_intersection(meshes[i], meshes[j])
       resampled = resample_all(polylines, spacing, min_edge)
       corefine(meshes[i], meshes[j], ecms[i], ecms[j])
       all_polylines[i].push_back(resampled);  // accumulate for mesh i
       all_polylines[j].push_back(resampled);
   End
   ```
   The corefine result is stored in-place; subsequent corefines for
   mesh `i` see the already-corefined version, which is exactly what
   we want for transitive conformality (vertices shared with mesh j
   stay in mesh i for the (i,k) corefine too).
5. **Per-mesh remesh**:
   ```
   For i in 0..N-1:
       polyline_set = flatten(all_polylines[i])
       output_meshes[i] = surface_Delaunay_remeshing(meshes[i],
           polyline_constraints = polyline_set,
           protect_constraints = true,
           mesh_edge_size = mesh-edge-size)
   End
   ```
   Each mesh sees ALL polylines that pass through it (one per
   intersecting partner). Because the polyline_constraints are
   resampled to the SAME points before remeshing, the output meshes
   remain conformal pairwise.
6. **Quality repair** on each output, exactly as Phase 1.
7. **Pairwise conformality re-check**: for every (i, j) pair in
   `sorted_pairs`, count the polyline vertices in `output_meshes[i]`
   that match a polyline vertex in `output_meshes[j]`. Both counts
   should equal `|resampled[(i,j)]|`. If not, log which pairs failed
   and return exit code 7.
8. **Write outputs**: for each input file `IN_DIR/<basename>.off`,
   write `OUT_DIR/<basename>_corefined.off` and (if the input was
   `.stl`) also `OUT_DIR/<basename>_corefined.stl` via the
   high-precision writer.
9. **Manifest JSON** (if `--manifest` given): list, for each output
   mesh:
   - input filename, output filename
   - vertex / face count before and after
   - bounding box
   - list of (partner index, polyline vertex count) tuples
   - quality stats from `quality_repair::compute_stats`

#### `intersection_graph.h` interface
```cpp
namespace intersection_graph {

// Returns the list of pairs (i, j), i < j, for which mesh[i] and
// mesh[j] have a non-empty surface intersection.  Internally uses
// PMP::intersecting_meshes.
std::vector<std::pair<std::size_t, std::size_t>>
find_intersecting_pairs(const std::vector<Mesh>& meshes);

// Sort pairs by descending bbox-overlap volume (a heuristic; see
// requirement 3 above).
void sort_by_overlap_volume(std::vector<std::pair<std::size_t, std::size_t>>& pairs,
                             const std::vector<Mesh>& meshes);

}  // namespace intersection_graph
```

### Edge Cases to Handle
- **Triple junction** (3 faults meeting at one polyline): handled
  naturally by the per-pair loop — pair (A,B) corefines a polyline
  into both; pair (A,C) corefines a polyline that may pass through
  the existing A∩B polyline, creating a degree-3 junction vertex.
  CGAL's `surface_intersection` precondition (degree-2 polyline
  vertices) holds *per pair*, not globally, so this is fine.
- **Disjoint subsets**: the intersection graph may have multiple
  connected components. Process each independently; report per
  component.
- **Empty `IN_DIR`**: return exit 0 with message "no inputs".
- **Single mesh with no intersections**: copy through unchanged.

### Acceptance Criteria
- [ ] Builds with `cgal-61`.
- [ ] On the 6 SAFS faults at 2000 m resolution, finds the 8
      intersecting pairs documented in
      `code_preprocess/EXPLORE_pymeshlab_intersections.md`.
- [ ] Outputs 6 `_corefined.off` files; pairwise conformality check
      passes for all 8 pairs.
- [ ] All 6 outputs satisfy `min(edge) >= 100 m`,
      `min(q_tri) >= 0.3`.
- [ ] Manifest JSON, when requested, parses with `python -m json.tool`
      and contains the expected per-mesh structure.

### Dependencies
- Depends on: Phase 1 (`quality_repair.h`, `polyline_resample.h`,
  `io_helpers.h`).
- Required by: Phase 3 (Python orchestration).

## Phase 3: Python orchestration `corefine_faults.py`

### Goal
A Python driver in `code_preprocess/` that reads `.ts` files, converts
them to `.off` (preserving per-vertex precision), invokes the
`corefine_set` C++ binary, and reports the result back as a set of
`.stl` files in `preprocess_data/`.

### Files to Create
- `code_preprocess/corefine_faults.py` — the Python driver.

### Files to Modify
- None.

### Detailed Requirements

#### CLI
```
corefine_faults.py [--in-dir DIR]                 (default: ../preprocess_data/cleaned_ts)
                   [--out-dir DIR]                (default: ../preprocess_data)
                   [--res {2000m}]                (only resolution validated by this plan)
                   [--mesh-edge-size SIZE]        (default: 1500.0)
                   [--min-edge SIZE]              (default: 100.0)
                   [--cgal-bin PATH]              (default: corefine_cgal/build/corefine_set)
                   [--workdir DIR]                (default: ../work/corefine_2000m)
                   [--keep-intermediate]          (don't delete OFF intermediates)
                   [--verbose]
```

The `--in-dir` default points at the *cleaned* TS directory the
upstream tool produces; this plan does not depend on that tool's
internal layout, only on the fact that each fault is a single
manifold-or-near-manifold `.ts`/`.off`/`.stl` in the input dir.

#### Algorithm
1. **Locate inputs**: glob `in-dir/SAFS-SAFZ-*_<res>.{ts,off,stl}`. If
   empty, exit 1 with helpful message. The pipeline accepts any of
   the three formats; selection is by extension.
2. **Convert each `.ts` to `.off`** in `workdir/` (skipped if input
   is already `.off` or `.stl`):
   - Reuse the parser from `code_preprocess/ts_to_stl.py:parse_ts()`.
   - **Pass through every vertex**. No z-filter; the input is
     assumed to already be cleaned by the upstream free-surface tool
     per the locked-in design decision in the plan header.
   - Write OFF in full double precision.
3. **Invoke `cgal-bin`**:
   ```python
   subprocess.run([str(cgal_bin), str(workdir), str(workdir),
                   "--mesh-edge-size", str(mesh_edge_size),
                   "--min-edge", str(min_edge),
                   "--manifest", str(workdir / "manifest.json")],
                  check=True)
   ```
4. **Convert each `*_corefined.off` to `*_corefined.stl`**:
   - Read OFF with `meshio`; write ASCII STL with `meshio.write(stl,
     mesh, file_format="stl-ascii")`.
   - Inspect the manifest for any non-conformal pair; if any, abort
     with exit 7 and print the offending pair.
5. **Move outputs to `out-dir`**: `*_corefined.stl` and the
   `manifest.json`.
6. **Print summary**: per-fault, per-pair, and overall quality stats
   read from `manifest.json`.
7. **Cleanup**: unless `--keep-intermediate`, delete `workdir`.

### Edge Cases to Handle
- **`cgal-bin` does not exist**: print build instructions
  (`conda activate cgal-61 && (cd code_preprocess/corefine_cgal &&
  mkdir build && cmake .. && make)`) and exit 1.
- **Subprocess returns non-zero**: capture stderr, print the C++ tool's
  error message, return the same exit code.
- **OFF / STL conversion failure**: caught by `meshio` exceptions;
  print which file failed.

### Acceptance Criteria
- [ ] On `project_7.0_preferred/raw_data/` at 2000 m: produces 6
      `_corefined.stl` files in `preprocess_data/`.
- [ ] Each STL passes `meshio.read(...)` without warnings.
- [ ] Manifest JSON conformality check shows all 8 pairs conforming.

### Dependencies
- Depends on: Phase 2 (binary).
- Required by: Phase 4 (gmsh consumes the STL output).

## Phase 4: gmsh multi-fault `.geo`

### Goal
A new `safs_multi_fault_box.geo` that mirrors `safs_fault_box_buried.geo`
but merges N corefined STL faults and embeds each as a constraint in a
single bounding-box volume, using the proven `Merge` + `Surface{} In
Volume{1}` pattern. Confirms that the conformal STLs from Phase 3 are
ingested correctly by gmsh and that bulk tet quality satisfies the hard
constraint (`min q_tet >= 0.1`, `>= 99.9% with q_tet >= 0.3`).

### Files to Create
- `code_meshing/safs_multi_fault_box.geo` — multi-fault `.geo`.

### Files to Modify
- None (existing single-fault `.geo` files unchanged).

### Detailed Requirements

#### `.geo` structure (mirror of `safs_fault_box_buried.geo`)
```gmsh
// Header comments same style as buried.geo.
SetFactory("OpenCASCADE");

// 1. Domain bounds: union of all corefined STL bboxes
//    + PAD_XY (50 km), PAD_TOP (5 km), PAD_BOTTOM (25 km).
//    Hardcoded after running Phase 3 (the manifest gives the union).
xmin = ...; xmax = ...;  // from manifest
ymin = ...; ymax = ...;
zmin = ...; zmax = ...;
dx = xmax - xmin; dy = ymax - ymin; dz = zmax - zmin;

// 2. Mesh size controls (HARD floor: LC_MIN must be >= 100.0).
LC_MIN     =  100.0;     // hard floor (matches Phase-1 --min-edge)
LC_NEAR    = 1500.0;     // target near-fault element edge
LC_FAR     = 10000.0;
DIST_INNER =  3000.0;
DIST_OUTER = 40000.0;

Mesh.MeshSizeMin               = LC_MIN;
Mesh.MeshSizeMax               = LC_FAR;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;

// 3. Box volume. OCC tags 1..6 for box faces; volume tag 1.
Box(1) = {xmin, ymin, zmin, dx, dy, dz};

// 4. Merge each corefined fault.
//    The first fault gets discrete-surface tag 7 (one past 1..6).
//    The k-th fault (1-indexed) gets tag 6+k.
//    HARD-CODED list -- driver writes this section using the manifest.
Merge "<safs_filename_1>_corefined.stl";
Merge "<safs_filename_2>_corefined.stl";
...
Merge "<safs_filename_N>_corefined.stl";

n_box_surfs = 6;
N_FAULTS    = N;     // filled in from the manifest
fault_tag_start = n_box_surfs + 1;
fault_tag_end   = n_box_surfs + N_FAULTS;
fault_surfs[] = {};
For k In {0 : N_FAULTS - 1}
    fault_surfs[] += {fault_tag_start + k};
EndFor

// 5. Embed each fault.
For k In {0 : #fault_surfs[]-1}
    Surface{fault_surfs[k]} In Volume{1};
EndFor

// 6. Distance-based size field that uses ALL fault surfaces.
Field[1] = Distance;
Field[1].SurfacesList = {fault_surfs[]};
Field[1].Sampling     = 100;

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = LC_NEAR;
Field[2].SizeMax  = LC_FAR;
Field[2].DistMin  = DIST_INNER;
Field[2].DistMax  = DIST_OUTER;
Background Field = 2;

// 7. Physical groups.
Physical Volume("rock", 1) = {1};
For k In {0 : N_FAULTS-1}
    Physical Surface(Sprintf("fault_%g", k+1), 100 + k + 1) = {fault_surfs[k]};
EndFor
Physical Surface("top", 200)    = {6};
Physical Surface("bottom", 201) = {5};
Physical Surface("sides", 202)  = {1, 2, 3, 4};

// 8. Mesh algorithms (proven on single-fault buried geometry).
Mesh.Algorithm   = 6;     // 2-D: Frontal-Delaunay
Mesh.Algorithm3D = 1;     // 3-D: Delaunay (HXT fails on embedded discrete surfaces)
```

The `<safs_filename_k>` and `N_FAULTS` are filled in by the driver
from `manifest.json`. The driver may emit the `.geo` directly so
that the file is regenerable.

#### Conformality preservation through STL Merge

Gmsh's `Merge` for STL deduplicates vertices within
`Geometry.Tolerance`, default `1e-8 * bbox_diag`. For a box ~500 km
across, that's ~5 mm. The Phase 1 high-precision STL writer (`%.10g`
or 15-digit `printf`) produces coordinates with sub-mm precision.
Therefore, after merging N corefined STLs, polyline vertices that are
shared between meshes A and B will be deduplicated into a single
vertex by gmsh, automatically.

The `.geo` should set the tolerance explicitly to the floor of the
Phase 1 writer's precision, e.g.:
```gmsh
Geometry.Tolerance = 1e-3;  // meters; well above STL ASCII precision
```
This protects against future precision-loss changes in the writer.

#### Verification

After `gmsh -3 safs_multi_fault_box.geo -o safs_multi_fault_box.msh`:
- Use `code_preprocess/msh_to_vtu.py` to split into bulk + per-fault
  VTU files.
- Use the quality analysis function from earlier work
  (`safs_fault_box_buried_bulk.vtu` analysis) to verify:
  - `min(q_tet) >= 0.1`
  - `>= 99.9 %` of tets satisfy `q_tet >= 0.3`
  - Shortest edge anywhere `>= LC_MIN = 100 m`

### Edge Cases to Handle
- **Disjoint fault subsets** (intersection graph has multiple
  connected components): each component meshes independently as a set
  of constraints inside the same volume. No special handling.
- **Fault touching box boundary**: solved as in
  `safs_fault_box_buried.geo` — set `PAD_TOP`, `PAD_BOTTOM`, `PAD_XY`
  large enough so all faults are interior to the box. The driver
  (Phase 3) computes the correct paddings from the manifest.
- **Triple junction at fault-fault-box-edge**: should not occur if
  paddings are sufficient; if the user shrinks the box, gmsh may
  reject. Fail-loud: assert via a `Geometry.Tolerance` check.

### Acceptance Criteria
- [ ] Generated by Phase 3 driver from the 6-fault 2000 m manifest.
- [ ] `gmsh -3 safs_multi_fault_box.geo` runs without errors and
      reports `0 ill-shaped tets`.
- [ ] Output `.msh` parses via `meshio`, has 1 volume + 6 fault
      surface physical groups + the 3 boundary groups.
- [ ] Bulk tet quality: `min q_tet >= 0.1`, `>= 99.9 % q_tet >= 0.3`,
      shortest edge >= 100 m.
- [ ] Each fault Physical Surface has the same triangle count as the
      corresponding corefined STL — i.e., gmsh preserved the
      triangulation as 2-D constraint without re-meshing.

### Dependencies
- Depends on: Phase 3 (corefined STLs in `preprocess_data/`).
- Required by: nothing in this plan; downstream MFEM consumes the
  resulting `.msh`.

## Phase 4.5: Optional — `BooleanFragments` retest on a corefined fault

### Goal
Document, in code, why `BooleanFragments` is not used for fault
embedding, by retrying it now that the input is a clean conformal
triangulation rather than a clipped+snapped STL.

### Files to Create
- `code_meshing/safs_multi_fault_box_boolfrag.geo` — same setup as
  Phase 4 but using `BooleanFragments` instead of `Surface{} In
  Volume{}`.
- `code_meshing/PHASE_4_5_BOOLFRAG_NOTES.md` — short note recording
  the result.

### Detailed Requirements
1. Same domain box and Merge as Phase 4.
2. After Merge, run `ClassifySurfaces{40 * Pi/180, 1, 0, 60 * Pi/180}; CreateGeometry;`
   to convert each merged STL into geometric (parametric) surfaces.
3. Track which surfaces came from each STL by snapshotting tags before
   and after each Merge.
4. Run `BooleanFragments` between the volume and the union of fault
   surfaces:
   ```gmsh
   v_out[] = BooleanFragments{ Volume{1}; Delete; }
                              { Surface{fault_surfs[]}; Delete; };
   ```
5. Mesh and capture the gmsh output.

### Acceptance Criteria
- [ ] One of: (a) `BooleanFragments` succeeds and produces a valid
      multi-region mesh with the expected per-fault Physical Surface
      groups, in which case Phase 4 should be reconsidered (use
      `BooleanFragments` instead of `Surface{} In Volume`); or
      (b) `BooleanFragments` fails (timeout, "No closed volume",
      empty regions, or wrong topology), in which case the note
      records the exact failure mode and we keep Phase 4's design.

### Dependencies
- Depends on: Phase 4.
- Required by: nothing.

## Phase 5: Validation harness

### Goal
End-to-end automated test that runs the full pipeline on the SAFS
6-fault 2000 m set, verifies all hard constraints, and produces a
reproducibility log.

### Files to Create
- `code_preprocess/test_corefine_pipeline.py` — pytest-style
  end-to-end test.
- `code_preprocess/VALIDATION_corefine_2000m.md` — generated report
  with quality histograms, conformality matrices, and timings.

### Files to Modify
- None.

### Detailed Requirements
1. Test fixture: the 6 cleaned `.ts` (or `.off`/`.stl`) files at
   2000 m resolution, supplied by the user from the upstream cleaning
   tool. **Only 2000 m is in scope** for this plan (per locked-in
   decision 3); higher-resolution validation belongs to a separate
   plan.
2. Test 1: build the C++ binaries via `cmake`. Skip with `pytest.skip`
   if `cgal-61` env is not active.
3. Test 2: run `corefine_faults.py` end-to-end into a temp directory.
4. Test 3: load each `_corefined.stl`, assert `min(edge) >= 100 m`,
   `min(q_tri) >= 0.3` per file.
5. Test 4: load each pair-wise polyline from the manifest, assert
   conformality (vertex-coord match in both meshes, max diff < 1e-6).
6. **Test 4b — triple-junction sentinel** (per locked-in decision 2):
   for every triplet `(A, B, C)` of meshes such that all three
   pairwise corefines are non-empty, locate the polyline endpoints
   shared between pairs `(A,B)` and `(A,C)` in mesh A. Assert that
   each such endpoint also exists at identical coords in B and in C
   (max diff < 1e-6 m). If any triple junction fails, log all
   offending triplets and fail the test. **This is the gate that
   would trigger escalation to options (ii) or (iii) in a follow-up
   plan.**
7. Test 5: run gmsh on the generated multi-fault `.geo`, parse the
   `.msh`, assert tet-quality criteria.
8. Test 6: emit the validation report (markdown) into the
   `code_preprocess/` dir, with per-fault stats, per-triplet
   conformality results, and the timing breakdown (corefine, remesh,
   gmsh).

### Acceptance Criteria
- [ ] All tests pass on the 2000 m fixture in under 5 min.
- [ ] The validation report contains tables and quality histograms;
      humans can read it without running the code.
- [ ] The triple-junction sentinel test (Test 4b) passes on all
      detected triplets, OR the report explicitly lists the failing
      triplets so a follow-up plan can address them.

### Dependencies
- Depends on: Phases 1, 2, 3, 4.
- Required by: nothing.

## Testing Strategy

### Unit tests
- `polyline_resample.h::resample`: synthetic open and closed
  polylines, asserting (a) first/last preserved, (b) all interior
  segments within `[0.5, 1.5] × target_spacing`, (c) total length
  preserved within target_spacing.
- `quality_repair::compute_stats`: hand-crafted mesh with known min
  edge, q_tri, n_below_floor; verify exact match.

### Integration tests (Phase 1)
- Two perpendicular rectangles (synthetic): exact geometric truth
  for the intersection line; verify polyline coords match analytically.
- Two slightly-misaligned planes (1° dihedral): tests the
  `features_angle_bound = 60°` default — should treat the
  intersection as a feature.
- Two near-coplanar planes (0.1°): expected to abort with the
  coplanar-overlap warning (Phase 1 edge case).

### End-to-end tests (Phase 5)
- 6-fault SAFS 2000 m: covered by acceptance criteria.
- 6-fault SAFS 1000 m: scale test (more triangles → more polyline
  points). May relax the < 5 min timing.

### Reference solutions
- For a single corefine pair, the analytic intersection line is
  known when the inputs are planes; use this as ground truth.
- For multi-fault, no analytic reference exists; cross-check by
  asserting:
  1. The total polyline-vertex count in mesh A across all its
     intersection partners equals the union (no double-counting).
  2. The Euler characteristic of each output mesh equals
     V − E + F (genus-preserving check).

## Risk Assessment

### High-risk areas
- **Surface_Delaunay_remeshing instability on near-coplanar inputs.**
  The Mesh_3 surface remesher uses a 3-D Delaunay refinement
  internally; near-coplanar surfaces produce flat tetrahedra in the
  Delaunay, which can cause non-termination or geometric drift.
  Mitigation: detect via the Phase 1 coplanar-overlap edge case;
  reject inputs with dihedral < ~5°.
- **STL-precision-induced conformality loss.** Default ASCII STL
  writers use ~6-digit precision, which is insufficient for 7-digit
  UTM coordinates. Mitigation: custom writer in `io_helpers.h` with
  explicit `%.10g`. Acceptance criterion 4 of Phase 1 verifies this.
- **Multi-pair corefine ordering.** Corefine for (A, B) modifies A
  in-place; subsequent corefine for (A, C) sees the modified A. If
  the (A, B) corefine introduces vertices very close to A∩C, the
  exact predicates may produce a different polyline than expected.
  Mitigation: process pairs in deterministic order by descending
  bbox-overlap volume (Phase 2); document that re-ordering may give
  bit-different but equally valid output.

### Medium-risk areas
- **`PMP::corefine` self-intersection precondition.** Inputs with
  near-coincident faces (which the SAFS data may have at fault
  perimeters) can cause `corefine` to throw `Self_intersection_exception`.
  Mitigation: pre-clean each input with `PMP::autorefine` if
  `PMP::does_self_intersect` returns true (Phase 1 step 3).
- **`surface_Delaunay_remeshing` polyline-resampling determinism.**
  Two calls to the function with the same `polyline_constraints`
  *should* produce the same resampled polyline vertex sequence, but
  this is not explicitly documented. Mitigation: pre-resample once
  at the polyline-resample step (Phase 1 step 5) and pass the
  pre-resampled polyline; the Delaunay remesher's internal resampling
  becomes a no-op on already-uniform input.
- **gmsh `Geometry.Tolerance` mismatch.** Default is bbox-diag-relative;
  if user changes domain size, the tolerance scales but the STL
  precision does not. Mitigation: Phase 4 sets
  `Geometry.Tolerance = 1e-3` explicitly.

### Low-risk areas
- **C++ build reproducibility**: `cgal-61` env is pinned;
  `cgal-cpp 6.1.1` API stable.
- **Triangulation density blowup on 500 m resolution**: the polyline
  is resampled at `mesh-edge-size / 2 = 750 m` regardless of source
  density, so polyline complexity does not explode. The total mesh
  size scales like the input triangle count, which is bounded.

## Resolved questions (locked in)

- **Q1: Pre-clipping at z=0** — RESOLVED: out of scope. A separate
  upstream tool (not part of this plan) emits cleaned `.ts` files;
  this pipeline takes them verbatim. No z-filter, no top-edge logic.
- **Q2: Triple-junction handling** — RESOLVED: pure pairwise corefine
  (option i). Phase 5 includes an explicit triple-junction sentinel
  test (Test 4b). If the sentinel fails on the SAFS 2000 m fixture,
  escalate via a separate follow-up plan that implements either
  option (ii) — pre-extracted snapped junctions — or option (iii) —
  single autorefine on the union.
- **Q3: Validation scope** — RESOLVED: 2000 m fixture only. No
  `--runlong` gate. Higher resolutions are a separate plan.
