# pymeshlab for intersecting fault-trace co-refinement — exploration report

## Question

Can [pymeshlab](https://pymeshlab.readthedocs.io/) take the multiple
intersecting fault-trace surfaces in
`miniapps/seas/safs/project_7.0_preferred/raw_data/` (six fault families ×
three resolutions, all open triangulated `.ts` sheets) and produce a single
**conformal triangulation** — i.e., a triangulation in which every
fault-fault intersection is represented by a polyline of edges that is
*shared between both faults*? If yes, gmsh can ingest the result via
`Surface{} In Volume{}` and we keep our existing pipeline. If no, what's
missing and what's the alternative?

## TL;DR

**No** — pymeshlab in its current form cannot produce a conformal
intersection of open fault-trace surfaces. The screenshots the user shared
describe MeshLab's GUI workflow ("Select Self Intersecting Faces", "Uniform
Mesh", "Boolean Union via Flatten Layers") sourced from a Google AI
overview. Tested against the actual library: the boolean filter rejects
open surfaces; the "Uniform Mesh" filter is a voxel/marching-cubes
resampler that destroys the source triangulation; and merge-and-clean
produces a non-conformal concatenation with T-junctions where the surfaces
pass through each other.

The right algorithm for this job is **CGAL's
`Polygon_mesh_processing::corefine()`**, which exists in the user's
`cgal-61` conda env as the C++ library `cgal-cpp 6.1.1` but is **not**
exposed through pymeshlab. Recommended path: a small C++ tool built
against `cgal-61` that corefines the fault pair (or set), emits a pair of
STLs sharing the intersection polyline, and feeds them into the existing
gmsh `.geo`.

## Scope

- Target: confirm whether pymeshlab can replace a hand-written CGAL
  preprocessor for SAFS multi-fault meshing.
- Entry points:
  - pymeshlab installation:
    `/Users/chunhuizhao/miniforge/envs/pythonenv/lib/python3.13/site-packages/pymeshlab/`
  - Plugin DSOs:
    `…/pymeshlab/PlugIns/libfilter_mesh_booleans.so`,
    `libfilter_mesh_alpha_wrap.so`, `libfilter_meshing.so`,
    `libfilter_clean.so`, `libfilter_plymc.so`, `libfilter_select.so`.
  - Data:
    `miniapps/seas/safs/project_7.0_preferred/raw_data/SAFS-SAFZ-*.ts`
    (six fault families: COAV-Mission_Creek, MJVS-San_Andreas,
    MULT-Southern_San_Andreas_and_Banning, SBMT-Garnet_Hill,
    SBMT-Mission_Creek, SBMT-San_Andreas).
- User's `cgal-61` env: build-time C++ environment (clang 19 + cmake 4.3
  + cgal-cpp 6.1.1 + eigen 3.4.0 + boost). **No Python interpreter**, so
  the env is for *compiling* CGAL programs, not for using CGAL from
  Python.

## Spatial overlap of the SAFS fault traces

Bounding boxes (UTM meters; printed for the 2000 m resolution; smaller
resolutions are denser but cover the same regions):

```
fault                                                   verts    x_lo    x_hi    y_lo    y_hi   z_lo  z_hi
COAV-Mission_Creek_fault_strand-CFM4                      201  541475  570749 3737744 3762850 -16813   486
MJVS-San_Andreas_fault-CFM6                               295  365054  428026 3808335 3838998 -13272  1538
MULT-Southern_San_Andreas_fault_and_Banning-CFM6          354  557182  622329 3692278 3750692 -13521   134
SBMT-Garnet_Hill_fault-CFM6                               535  504574  560211 3744545 3769519 -18073   683
SBMT-Mission_Creek_fault_strand-CFM4                      260  492981  543146 3761457 3772891 -15557  2067
SBMT-San_Andreas_fault-CFM6                               750  427968  538743 3756115 3808445 -17019  2061
```

Pairwise XY-bbox-AND-Z-bbox overlap (necessary condition for actual 3-D
intersection):

```
COAV-Mission_Creek          ×  MULT-Southern_SAF_Banning      13.6 km × 12.9 km × 13.7 km
COAV-Mission_Creek          ×  SBMT-Garnet_Hill               18.7 km × 18.3 km × 17.3 km
COAV-Mission_Creek          ×  SBMT-Mission_Creek              1.7 km ×  1.4 km × 16.0 km
MJVS-SAF                    ×  SBMT-SAF                        0.1 km ×  0.1 km × 14.8 km
MULT-Southern_SAF_Banning   ×  SBMT-Garnet_Hill                3.0 km ×  6.1 km × 13.7 km
SBMT-Garnet_Hill            ×  SBMT-Mission_Creek             38.6 km ×  8.1 km × 16.2 km
SBMT-Garnet_Hill            ×  SBMT-San_Andreas               34.2 km × 13.4 km × 17.7 km
SBMT-Mission_Creek          ×  SBMT-San_Andreas               45.8 km × 11.4 km × 17.6 km
```

Eight pairs have triple overlap; physically these are exactly the splay
junctions of the SAFS where fault planes meet (Garnet Hill / Mission
Creek / San Andreas / Banning). The intersection-curve preprocessor
must produce one or more shared polylines per overlapping pair.

## Architecture of pymeshlab's relevant filters

`MeshSet` exposes ~360 filters (each filter is also a MeshSet method).
The ones relevant to intersection handling, with backends identified
from the plugin's exported symbols:

| Filter (pymeshlab name)                                  | Plugin DSO                       | Backend (from symbols)            | What it does |
|----------------------------------------------------------|----------------------------------|-----------------------------------|--------------|
| `compute_selection_by_self_intersections_per_face`       | `libfilter_select.so`            | VCG `tri::Clean::SelfIntersections` | Marks self-intersecting faces; diagnostic only. |
| `generate_boolean_union/intersection/difference/xor`     | `libfilter_mesh_booleans.so`     | **`igl::copyleft::cgal::mesh_boolean`** (CGAL `Lazy_exact_nt` + libigl winding-number CSG) | True CSG boolean on **closed manifold solids**. |
| `generate_resampled_uniform_mesh`                        | `libfilter_plymc.so`             | VCG `PlyMC` marching-cubes (`MCMesh`, `MCVertex`, `MCFace`) | Voxel rebuild of input geometry. |
| `meshing_close_holes`                                    | `libfilter_meshing.so`           | VCG hole-filling                  | Patches free borders only. |
| `meshing_isotropic_explicit_remeshing`                   | `libfilter_meshing.so`           | VCG split/collapse/swap/smooth, `reproject` flag | Topology-preserving remesh of one mesh. |
| `meshing_merge_close_vertices`                           | `libfilter_meshing.so`           | VCG vertex weld                   | Welds verts within a distance threshold. |
| `meshing_snap_mismatched_borders`                        | `libfilter_meshing.so`           | VCG border snap                   | Snaps free-border verts that are nearly coincident. |
| `meshing_repair_non_manifold_edges/vertices`             | `libfilter_meshing.so`           | VCG topology repair               | Removes non-manifold edges/verts. |
| `generate_by_merging_visible_meshes`                     | core                             | naive concat                      | Concatenates triangle lists; vertex dedup at exact coincidence only. |
| (alpha-wrap-of-3D-points filter)                         | `libfilter_mesh_alpha_wrap.so`   | CGAL `Polygon_mesh_processing` (`Surface_mesh<Point_3<Epick>>`) | Shrink-wraps a manifold around inputs at chosen alpha. |

Confirmed concretely from the boolean DSO: symbols include
`igl::copyleft::cgal::SelfIntersectMesh<CGAL::Epeck, …>`,
`igl::copyleft::cgal::BinaryWindingNumberOperations<MeshBooleanType, K>`,
`igl::copyleft::cgal::WindingNumberFilter<KeeperType>`. The dependency
list pulls in `libgmpxx`, `libmpfr`, `libgmp` for CGAL's exact
arithmetic. So `generate_boolean_*` is the *exact* libigl-CGAL
boolean — no precision concerns, but it requires closed inputs.

## Empirical test on synthetic open intersecting surfaces

Two perpendicular triangulated rectangles (each 4×4 grid → 25 verts /
32 tris): rectangle A in the xy-plane, rectangle B in the xz-plane
through y=5, intersecting along the line y=5, z=0, x ∈ [0, 10].

Filter chain → result `(V, F, verts_on_intersection_line, edges_on_line)`:

| Chain                                              | V    | F    | verts on line | edges on line | conformal? |
|----------------------------------------------------|------|------|---------------|---------------|------------|
| `generate_boolean_union(A, B)`                     |   —  |   —  |       —       |       —       | **errors out**: `Mesh inputs must induce a piecewise constant winding number field. Make sure that both the input mesh are watertight (closed).` |
| `generate_by_merging_visible_meshes`               |  45  |  64  |       5       |      16       | no — naive concat; verts on line are pre-existing grid verts, not new intersection verts. T-junctions persist. |
| `merge` → `meshing_isotropic_explicit_remeshing(target=2%)` | 2145 | 4096 |      33       |     128       | no — local subdivide/swap/smooth; doesn't compute A∩B. New verts come from edge splits, not intersections. |
| `merge` → `generate_resampled_uniform_mesh(2%)`     | 2968 | 5558 |       0       |       0       | no — voxel marching cubes. **Original triangulation is destroyed**; the output is a manifold envelope whose verts have no relation to A or B. |

The `verts_on_intersection_line` and `edges_on_line` columns count
analytic colinearities with the geometric A∩B line, not new vertices
introduced by an intersection algorithm. The boolean filter is the only
filter that even *attempts* to handle two-mesh intersection, and it
rejects open inputs immediately.

A "thicken each fault to a slab → boolean union → extract a midsurface"
work-around is theoretically possible but introduces (a) an arbitrary
thickness parameter, (b) a slab-collapse step that's not built in, and
(c) drift away from the original fault triangulation. Not recommended.

## Why pymeshlab can't do this

The two-surface co-refinement algorithm we need is:

1. Compute `intersection_polyline = A ∩ B` exactly.
2. Insert the polyline as new edges in **both** A and B (each new
   intersection vertex appears in *both* meshes; each intersection
   segment becomes a triangle edge in *both* meshes' triangulations).
3. Optionally retriangulate the affected faces of A and B locally so the
   result remains a triangulation with reasonable quality.

This is `CGAL::Polygon_mesh_processing::corefine(A, B, params...)` —
documented as "computes the corefinement of two triangulated surface
meshes" and works on open or closed surfaces. It is implemented in
header `CGAL/Polygon_mesh_processing/corefinement.h` and lives inside
`cgal-cpp` (which the `cgal-61` env provides at version 6.1.1).

pymeshlab does not wrap PMP::corefine. Its boolean plugin uses libigl's
**volumetric** boolean (winding-number CSG on closed solids), not
PMP::corefine. The alpha-wrap plugin uses PMP, but only to wrap a
manifold around point clouds — not to corefine two meshes.

**No combination of currently-installed pymeshlab filters produces a
conformal intersection of two open meshes.**

## Consumption by the existing pipeline

Our existing pipeline is:

```
.ts source      ts_to_stl.py                      mesh/safs_fault_box.geo
    │     →  (clip @ z=0, snap top)  →  .stl  →    Merge "*.stl";          → gmsh -3 → .msh
    │                                              Surface{7} In Volume{1};
```

For multiple intersecting faults, the gmsh side already supports
embedding *several* `Surface{} In Volume{1}` constraints — but each
constraint must be a triangulation that **agrees with all the others on
their intersections**. If two faults pass through each other but their
triangulations don't share edges along the intersection line, gmsh's
3-D mesher reports "No closed volume" (same failure mode we hit earlier
with `ClassifySurfaces`). So the corefined STL pair *is* the missing
input for a multi-fault `.geo`.

If we had a co-refined STL set
`fault_A_corefined.stl, fault_B_corefined.stl, fault_C_corefined.stl, …`
the existing `safs_fault_box.geo` extends to:

```gmsh
Merge "fault_A_corefined.stl";  // → discrete surface 7
Merge "fault_B_corefined.stl";  // → discrete surface 8
Merge "fault_C_corefined.stl";  // → discrete surface 9
…
For s In {7, 8, 9, …}
    Surface{s} In Volume{1};
EndFor
```

with no other change. We've already verified `Algorithm3D = 1`
(Delaunay) without `ClassifySurfaces` is the right choice for embedding
discrete STL surfaces.

## Recommended path forward

**Build a small C++ corefine tool against `cgal-61`.** The tool is ~80
lines of CGAL code:

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/IO/STL.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>                      Mesh;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv) {
    Mesh A, B;
    CGAL::IO::read_polygon_mesh(argv[1], A);
    CGAL::IO::read_polygon_mesh(argv[2], B);
    PMP::corefine(A, B);
    CGAL::IO::write_polygon_mesh(argv[3], A);  // refined A
    CGAL::IO::write_polygon_mesh(argv[4], B);  // refined B
}
```

Build with cmake using `cgal-61`'s `cgal_create_CMakeLists` helper. For
N fault surfaces, run pairwise corefine in topological order (a single
`PMP::corefine(A, B)` modifies both A and B in-place, so subsequent
`PMP::corefine(A, C)` sees the already-refined A, and so on). CGAL also
exposes `PMP::corefine_and_compute_union/difference` if a CSG result is
ever wanted, but for our use case the plain `corefine` is what we need.

Pros:
- Uses the C++ env the user already provisioned.
- Exact predicates (no welding/snapping tolerance to tune).
- Output is two STLs that are conformal at A∩B, ready for the existing
  `Merge` + `Surface{} In Volume` pipeline.
- ~100 LoC + a CMakeLists.

Cons:
- One C++ tool to build (and rebuild after CGAL bumps).

**Alternative (if a Python-only solution is preferred):** install
`compas_cgal` or `cgal-swig-bindings` in `pythonenv`. Both expose
`PMP::corefine`. They add a runtime dep but avoid the C++ build.

**Not recommended:** any pymeshlab-only chain. The closest viable hack is
"thicken each fault to a slab → boolean → extract midsurface", which
needs ~150 LoC of bespoke geometry handling and still does not match the
exactness of CGAL's corefine.

## Conventions noted

- **OpenMP duplicate-init**: pymeshlab loads `libomp.dylib` from its
  own `Frameworks/`; if the host process already has libomp (numpy /
  scipy via OpenBLAS), `import pymeshlab` after numpy crashes with
  `OMP: Error #15`. Workaround:
  `KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1`. Worth wrapping in any
  pymeshlab-using script.
- **`PercentageValue` vs absolute**: pymeshlab parameters that have
  units are wrapped with `pml.PercentageValue(p)` for fractions of the
  bbox-diagonal, or raw floats for absolute values. Confused these
  silently if you pass a bare float to e.g. `targetlen`.
- **Filter names**: pymeshlab maintains a translation table from
  MeshLab GUI names to Python filter names
  (`pml.replace_pymeshlab_filter_names()`). The GUI names in the
  screenshots ("Uniform Mesh", "Merge Close Vertices", "Boolean Union")
  map to `generate_resampled_uniform_mesh`,
  `meshing_merge_close_vertices`, `generate_boolean_union`.

## Gotchas

- The screenshots' "Uniform Mesh" advice is **actively harmful** for
  fault meshing: it discards the source triangulation in favor of a
  voxel marching-cubes envelope. Geological surfaces should never go
  through this filter.
- `generate_by_merging_visible_meshes` does **not** call
  `meshing_remove_duplicate_vertices` afterwards. The 5-vertex
  reduction in the synthetic test (50 → 45) was from MeshLab's
  exact-coordinate dedup at merge time, not a topological union.
- `meshing_isotropic_explicit_remeshing` with `reprojectflag=True`
  reprojects new vertices onto the source mesh, but doesn't notice
  the *other* mesh's surface, so the two faults remain independent
  triangulations even after remeshing.
- `meshing_snap_mismatched_borders` only operates on **free borders**
  (open boundaries) — not on the interior of the surface. Two
  faults that pass through each other's *interior* triangles are
  invisible to this filter.

## Open questions

- Is there a downstream MFEM constraint that two co-refined fault
  surfaces, embedded as two `Physical Surface` groups in the same
  `Volume`, both work as fault interfaces in `seas_bp5_full`?
  (Probably yes, since the BP5 driver already handles separate fault
  groups, but verify before scaling beyond 2 faults.)
- For N > 2 faults, the order of pairwise PMP::corefine calls can
  matter for vertex stability. CGAL has
  `PMP::corefine_with_constrained_edges()` and a multi-mesh corefine
  variant in newer releases (≥ 6.0?) — worth checking what 6.1.1
  exposes.

## Action items if we proceed with the C++ corefine path

1. New directory: `code_preprocess/corefine_cgal/` with
   `corefine_pair.cpp`, `CMakeLists.txt`, and a small Python driver
   `code_tools/corefine_faults.py` that builds a topological-sort plan
   over the pairwise-overlap graph and invokes the C++ binary.
2. Output convention: each input
   `*.stl` produces a sibling `*_corefined.stl`; the Python driver
   reports which intersection pairs added vertices/edges.
3. Update `safs_fault_box.geo` to merge a list of `*_corefined.stl`
   files and embed each as `Surface{tag} In Volume{1}`.
4. Test on the SBMT-Garnet_Hill × SBMT-Mission_Creek pair first (the
   largest XY overlap, 38.6 × 8.1 × 16.2 km), then on the chain
   GH × MC × SAF.
