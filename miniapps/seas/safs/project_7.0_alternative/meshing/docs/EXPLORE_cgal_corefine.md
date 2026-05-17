# CGAL `Polygon_mesh_processing::corefine` for fault-trace intersections — exploration report

## Question

User context: previous attempts at corefine-based multi-fault meshing
produced *conformal* output (intersection polyline shared between
meshes) but with **very bad triangle quality at the intersection** —
sliver/needle/cap triangles that downstream gmsh 3-D meshing chokes on.

Question: does CGAL itself provide tools to recover quality after
corefine while keeping the intersection conformal? If yes, which
function(s) and named parameters?

## TL;DR

**Yes, CGAL has the right tools.** What was missing in the prior
attempt is the *post-processing* stage. `PMP::corefine()` only inserts
the intersection edges — it makes no quality guarantee. CGAL ships
*three* downstream filters that, used together with the constrained
edge map produced by corefine, recover near-equilateral triangles
while preserving the conformal intersection:

1. `PMP::remove_almost_degenerate_faces()` — explicit cap/needle
   collapse, parameterized by `cap_threshold` (cosine of max angle)
   and `needle_threshold` (longest/shortest edge ratio). Designed
   exactly for the slivers that corefine introduces.
2. `PMP::isotropic_remeshing()` with
   `edge_is_constrained_map = ecm; protect_constraints = false;
   collapse_constraints = true` — VCG-style split/collapse/swap/smooth.
   Constrained sub-edges that are too short are *collapsed*, fixing
   the slivers along the intersection without losing the polyline as
   a topological feature.
3. `PMP::surface_Delaunay_remeshing()` (CGAL 6.0+, present in
   `cgal-61`) — runs Mesh_3's Delaunay refinement on the input
   surface, with `polyline_constraints = …` (the exact corefine
   output polyline) and `protect_constraints = true`. This is the
   highest-quality option: the polyline is held exact, every other
   triangle becomes Delaunay with size bounded by `mesh_edge_size`.

The probable cause of the user's prior bad-quality output: corefine
was used either standalone, or with `isotropic_remeshing(..., protect_constraints=true)`,
which *locks the intersection slivers in place* instead of collapsing
them. Switching to `protect_constraints=false; collapse_constraints=true`
or to `surface_Delaunay_remeshing` recovers quality.

## Scope

- Target: post-corefine quality recovery in CGAL 6.1.1
  (`cgal-cpp` package in the user's `cgal-61` conda env).
- Headers explored:
  `/Users/chunhuizhao/miniforge/envs/cgal-61/include/CGAL/Polygon_mesh_processing/`
  - `corefinement.h` — `corefine`, `corefine_and_compute_*`
  - `intersection.h` — `surface_intersection` (polyline output, no
    mesh modification)
  - `remesh.h` — `isotropic_remeshing`, `split_long_edges`
  - `surface_Delaunay_remeshing.h` — Delaunay refinement remesher
  - `repair_degeneracies.h` — `remove_almost_degenerate_faces`,
    cap/needle predicates
  - `tangential_relaxation.h` — vertex smoothing
  - `autorefinement.h` — `autorefine` for self-intersections
  - `Adaptive_sizing_field.h`, `Uniform_sizing_field.h` — sizing fields

## Architecture: the relevant CGAL functions and what they do

### `PMP::corefine(tm1, tm2, np1, np2)` — the intersection inserter

Header: `corefinement.h` (line 752 in v6.1.1).

What it does (read literally from doxygen + body):
- Computes A ∩ B as a polyline using exact predicates (`Lazy_exact_nt`).
- For every triangle of A or B that is cut by the polyline, replaces
  the cut triangle with a sub-triangulation that includes the new
  intersection vertices and edges.
- Output: A and B are mutated in-place; both now contain the
  intersection polyline as triangle edges. **No vertex on the polyline
  in A is geometrically distinct from the corresponding vertex in B**
  — they share point coordinates exactly.
- Quality of the new triangles: **no guarantee**. The new sub-
  triangulation can have arbitrarily thin slivers when the polyline
  cuts triangles near a vertex or grazes an edge.

Critical named parameters (verified in `corefinement.h:719-748`):

| Parameter (np1 / np2) | Type | Default | Meaning |
|---|---|---|---|
| `vertex_point_map` | RPM<vd, Point_3> | `vertex_point` | how to read vertex coordinates |
| `edge_is_constrained_map` | RPM<ed, bool> | constant `false` | **OUT**: every edge of `tm1` (resp. `tm2`) created by the intersection is set to `true` here |
| `visitor` | model of `PMPCorefinementVisitor` | `Default_visitor` | hook to track new face creation |
| `throw_on_self_intersection` | bool | `false` | check for self-intersections in the strip near A∩B |
| `do_not_modify` | bool | `false` | only modify the *other* mesh; useful when one input is the canonical fault and the other is being refined to it |

The `edge_is_constrained_map` parameter is the key handoff to the
quality-recovery step. After `corefine` returns, `ecm[e] == true`
for **every edge `e` that was inserted by the corefinement**, and
`ecm` is a *read/write* property map that downstream filters
(`isotropic_remeshing`, `remove_almost_degenerate_faces`,
`surface_Delaunay_remeshing`) accept verbatim.

### `PMP::surface_intersection(tm1, tm2, polyline_output, np)` — polyline-only

Header: `intersection.h:1731`.

Computes A ∩ B as a sequence of `std::vector<Point_3>` polylines
(each polyline open or closed, all interior vertices degree 2). Does
**not** modify `tm1` or `tm2`. Useful for:
- Pre-computing the polyline so it can be **resampled** at a uniform
  spacing before being handed to the remesher (this is the trick
  for guaranteed conformality of multiple Delaunay-remeshed meshes;
  see "Pipeline" below).
- Detecting which pairs in a fault list actually intersect (faster
  than full corefine when most pairs don't).

### `PMP::isotropic_remeshing(faces, target_edge_length, pmesh, np)` — quality remesher

Header: `remesh.h:217`.

Algorithm (from CGAL doxygen):
1. **Split** every edge longer than `4/3 · target_edge_length`.
2. **Collapse** every edge shorter than `4/5 · target_edge_length`
   (subject to `collapse_constraints` flag, see below).
3. **Flip** edges to improve valence and shape.
4. **Tangential relaxation** of free (non-constrained) vertices.
5. **Reproject** moved vertices onto the input surface (controlled by
   `do_project`, default `true`).
6. Repeat for `number_of_iterations` (default 1; for fault data, 5–10
   is typical).

Critical named parameters (verified in `remesh.h:94-196`):

| Parameter | Type | Default | What it does, and the trap |
|---|---|---|---|
| `edge_is_constrained_map` | RWPM<ed, bool> | constant `false` | constrained edges are **never flipped** and their endpoints are **never moved by smoothing**. Sub-edges generated by splits inherit the flag. |
| `vertex_is_constrained_map` | RWPM<vd, bool> | constant `false` | constrained vertices are never moved |
| **`protect_constraints`** | bool | **`false`** | if `true`, constrained edges are *neither split nor collapsed* — held exactly. **This is the trap**: any sliver edge that corefine produced will stay tiny; the remesher then has to mesh around it with similarly tiny triangles, or fail. |
| **`collapse_constraints`** | bool | **`true`** | if `true` (default), constrained edges *can* be collapsed when too short, then the merged vertex inherits the constraint. Sub-edges that are still on the polyline keep the constraint flag. **Use this to remove corefine-induced slivers without losing the polyline.** Ignored when `protect_constraints = true`. |
| `relax_constraints` | bool | `false` | endpoints of constrained edges may slide ALONG the polyline (1-D smoothing on the constraint). Useful for redistributing nodes uniformly along the intersection. |
| `number_of_iterations` | unsigned | `1` | typical: 5–10 |
| `number_of_relaxation_steps` | unsigned | `1` | tangential smoothing iters per outer iter |
| `do_split` / `do_collapse` / `do_flip` / `do_project` | bool | `true` each | pipeline stage toggles |
| `face_patch_map` | RWPM<fd, ID> | components-of-constrained | patches are maintained; remesher won't merge patches across constrained edges |

Pre-condition (from `remesh.h:60`): if `protect_constraints = true`,
the constrained edges must already be ≤ `4/3 · target_edge_length`.
Otherwise the algorithm cascades vertex insertions and may not
terminate.

### `PMP::surface_Delaunay_remeshing(tmesh, np)` — Delaunay refinement remesher

Header: `surface_Delaunay_remeshing.h:180`. Newer than
`isotropic_remeshing`; lives in the `Mesh_3` directory of CGAL but is
exposed under PMP. Uses Mesh_3's tetrahedral Delaunay refinement
algorithm restricted to the surface — the same engine that produces
high-quality tet meshes for `Mesh_3`.

Critical named parameters (verified in `surface_Delaunay_remeshing.h:87-129`):

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `mesh_edge_size` | scalar / sizing field | required | uniform or space-varying upper bound for triangle edge length |
| `features_angle_bound` | scalar (degrees) | `60` | edges with sharper dihedral than this are auto-detected as features and protected; `180` means only border edges are protected |
| `edge_is_constrained_map` | RWPM<ed, bool> | constant `false` | mark explicit constraints; takes priority over `features_angle_bound` |
| `polyline_constraints` | range of polylines | empty | set of `MeshPolyline_3` objects that are **resampled** and protected; takes priority over `features_angle_bound` (but yields to `edge_is_constrained_map` when both given) |
| `protect_constraints` | bool | `false` | enable protection of features / constrained edges / polyline constraints |
| `face_patch_map` | RWPM<fd, ID> | connected components | preserves face-patch structure |

Output: a *new* `TriangleMesh` (the output type can differ from the
input). The original mesh is **read-only** to this function — it just
serves as the surface that gets re-meshed.

This is the highest-quality option. The Delaunay refinement engine
guarantees:
- All triangles are Delaunay (no edge flippable to improve circumradius/edge ratio).
- Triangle aspect ratios bounded by the refinement algorithm's quality criterion.
- Constrained polylines are resampled at the target edge length.

The cost is: it's slower than `isotropic_remeshing` (Delaunay
refinement is more compute-intensive) and the algorithm is sometimes
fragile on near-coplanar geometry — a fault dipping at a few degrees
should be fine; a flat-on-flat overlap is not.

### `PMP::remove_almost_degenerate_faces(face_range, tmesh, np)` — explicit needle/cap killer

Header: `repair_degeneracies.h:643`.

What it does (verified in `repair_degeneracies.h:557-563`): walks the
mesh, classifies each triangle as a needle (longest/shortest edge
ratio above `needle_threshold`) or a cap (an angle above
`cap_threshold` in degrees), then:
- **Needles**: collapse the shortest edge (subject to
  `collapse_length_threshold`).
- **Caps**: flip the edge opposite to the largest angle (subject to
  `flip_triangle_height_threshold`); if the cap is on the boundary,
  delete the face.

Critical named parameters:

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `cap_threshold` | double | `cos(160°)` | cosine of min cap angle; values in `[-1, 0]` |
| `needle_threshold` | double | `4` | edge-ratio threshold for needle |
| `collapse_length_threshold` | double | `0` (no upper limit) | refuse to collapse edges longer than this; useful to avoid collapsing the polyline itself |
| `flip_triangle_height_threshold` | double | `0` (no limit) | refuse to flip if the resulting triangle would be too "tall" |
| `edge_is_constrained_map` | RPM<ed, bool> | none | constrained edges are **never collapsed nor flipped** |
| `vertex_is_constrained_map` | RPM<vd, bool> | none | constrained vertices are guaranteed to remain in `tmesh` |
| `filter` | functor | always-`true` | per-operation veto |

Returns `true` iff every almost-degenerate face was successfully
removed; some faces may resist removal due to topological constraints,
in which case the function returns `false` and leaves them.

This filter is the **specific cure** for the slivers corefine
introduces. Run it immediately after corefine, with
`edge_is_constrained_map = ecm` and
`collapse_length_threshold = 1.5 × target_edge_length` so that the
intersection polyline is preserved as a 1-D feature but its
sub-edges that are absurdly short get collapsed.

### `PMP::tangential_relaxation(vertices, tm, np)` — pure smoothing

Header: `tangential_relaxation.h:130`.

Vertex-only smoother: each vertex moves to the centroid of its
1-ring projected back onto the local tangent plane. Constrained
edges and constrained vertices are honored. No connectivity changes.

Useful as a final pass after isotropic_remeshing or
remove_almost_degenerate_faces to clean up local quality without
touching the topology.

### `PMP::autorefine(tm, np)` — self-intersection autorefiner

Header: `autorefinement.h:1787`.

Computes self-intersections of a single mesh and inserts the
intersection polyline as new edges (analogous to `corefine` but for
one mesh against itself). This is what you use if you concatenate
multiple fault triangulations into one mesh first and want CGAL to
sort out the inter-fault intersections in one pass.

Pros: handles N faults in one call; the constraint graph captures
*all* pairwise intersections.

Cons: any pre-existing self-intersection in the input (e.g., a fault
that already self-intersects in the source `.ts`) gets resolved at
the same time, which may not be the geological behavior you want.

## Recommended pipeline

```cpp
//
//  --- Pre-flight: check inputs are clean (no self-intersections) ---
//
if (PMP::does_self_intersect(A)) PMP::autorefine(A);
if (PMP::does_self_intersect(B)) PMP::autorefine(B);

//
//  --- Step 1: compute the intersection polyline (read-only) ---
//
std::vector<std::vector<K::Point_3>> polylines;
PMP::surface_intersection(A, B, std::back_inserter(polylines));
//   `polylines` is a list of open/closed sequences of points; for two fault
//   surfaces meeting in a line, expect one polyline per connected
//   intersection curve.
//
//  --- Step 1.5: resample polylines at the target spacing ---
//   (Critical for multi-mesh conformality: both A and B must see the SAME
//    resampled polyline so they agree on the vertex set after remeshing.)
//
auto resampled = resample_polylines(polylines, target_edge_length);

//
//  --- Step 2: corefine A and B, recording new edges into ecmA, ecmB ---
//
auto ecmA = make_ed_pmap(A);   // boost dynamic edge property map<bool>, default false
auto ecmB = make_ed_pmap(B);
PMP::corefine(A, B,
    PMP::parameters::edge_is_constrained_map(ecmA),
    PMP::parameters::edge_is_constrained_map(ecmB));
//   After this call: ecmA[e] == true for every edge of A that was inserted
//   by the corefinement (and similarly ecmB for B). The two meshes are
//   now conformal at A∩B but probably contain slivers.

//
//  --- Step 3 (option A, recommended): surface_Delaunay_remeshing each ---
//   Highest quality. Pass the resampled polyline as polyline_constraints
//   so the same vertices appear in both A_out and B_out's intersection.
//
TriangleMesh A_out, B_out;
A_out = PMP::surface_Delaunay_remeshing(A,
    PMP::parameters::polyline_constraints(resampled)
                   .protect_constraints(true)
                   .mesh_edge_size(target_edge_length));
B_out = PMP::surface_Delaunay_remeshing(B,
    PMP::parameters::polyline_constraints(resampled)
                   .protect_constraints(true)
                   .mesh_edge_size(target_edge_length));

//
//  --- Step 3 (option B, faster): remove_degenerate + isotropic_remeshing ---
//   Less expensive, slightly lower quality, conformal as long as both
//   meshes collapse the same slivers consistently.
//
PMP::remove_almost_degenerate_faces(faces(A), A,
    PMP::parameters::edge_is_constrained_map(ecmA)
                   .needle_threshold(4.0)
                   .cap_threshold(std::cos(160. * CGAL_PI/180.))
                   .collapse_length_threshold(1.5 * target_edge_length));
PMP::isotropic_remeshing(faces(A), target_edge_length, A,
    PMP::parameters::edge_is_constrained_map(ecmA)
                   .protect_constraints(false)        // <-- KEY
                   .collapse_constraints(true)        // <-- KEY (default)
                   .relax_constraints(true)
                   .number_of_iterations(5));
// (and the same for B)

//
//  --- Step 4: write conformal STLs for gmsh ---
//
CGAL::IO::write_polygon_mesh("fault_A_corefined.stl", A_out);
CGAL::IO::write_polygon_mesh("fault_B_corefined.stl", B_out);
```

**For N > 2 faults**, use `PMP::intersecting_meshes(range, out)` from
`intersection.h:1664` to find which pairs intersect in O(N²)
bbox-pruned tests, then walk the intersection graph and corefine each
pair, accumulating the polyline set. The polyline set is then
resampled once and passed to a single `surface_Delaunay_remeshing`
call per fault (each fault sees only its own subset of the polyline
set).

Alternatively, concatenate all faults into one polygon soup, repair
the soup, and call `PMP::autorefine` once — but you lose per-fault
identity, which we'll need for assigning Physical Surface tags in
gmsh.

## Why the prior attempt produced bad quality

Without seeing the prior code, the most likely cause is one of:

1. **Used corefine alone**, then exported. Corefine guarantees
   conformality, **not** quality. Slivers are inherent.
2. **Used `isotropic_remeshing` with `protect_constraints = true`**,
   which is the obvious-looking option but locks the slivers. The
   pre-condition documented in `remesh.h:60` —
   *"if constraints protection is activated, the constrained edges
   must not be longer than 4/3·target_edge_length"* — is the only
   pre-condition CGAL checks; it does *not* check that constrained
   edges are not too *short*, because by-design the user is supposed
   to have pre-cleaned the constraint polyline at uniform spacing.
   With raw corefine output as the constraint, this pre-condition is
   not satisfied (constraint edges are wildly varying length).
3. **Skipped polyline resampling**: the corefine intersection has
   variable-length sub-edges (from sub-millimeter to several km).
   Without resampling at a chosen target length, both
   `protect_constraints` and `collapse_constraints` paths struggle.
4. **Skipped `remove_almost_degenerate_faces`**: in the option-B
   pipeline, this is the explicit slivers-killer. Without it, the
   subsequent isotropic remesher inherits the corefine's needles
   and caps as starting state and can take many iterations (or
   diverge) trying to flip/collapse them back to good shape.

Any combination of (1)–(4) explains "conformal but bad quality".

## Conventions noted

- **Kernel choice**: `corefine` requires that vertex coordinates be
  comparable with exact predicates. CGAL's
  `Exact_predicates_inexact_constructions_kernel` (Epick) is the
  default and is what cgal-cpp 6.1.1 ships with. The Lazy_exact_nt
  symbols already seen in libigl-CGAL (libfilter_mesh_booleans.so)
  are the same machinery — exact comparison, double constructions —
  so coordinates stay numerically clean.
- **Property maps**: most filters take `boost::dynamic_*_property_map`
  or boost-graph internal maps. For `edge_is_constrained_map` the
  recommended idiom is
  `auto ecm = boost::get(CGAL::dynamic_edge_property_t<bool>(), mesh);`
  or `using ECM = std::unordered_map<edge_descriptor, bool>;` wrapped
  in `boost::make_assoc_property_map(ecm_map)`.
- **`Surface_mesh<Point_3>`** is the recommended `TriangleMesh` model
  for these workflows. CGAL also accepts `Polyhedron_3` but
  `Surface_mesh` is faster and what newer PMP examples use.
- **CGAL 6.1.1 specifics** (relative to 5.x):
  - `surface_Delaunay_remeshing` was added in 6.0.
  - `Adaptive_sizing_field` and `Uniform_sizing_field` for
    space-varying target edge length: 6.0+.
  - `relax_constraints` flag in `isotropic_remeshing`: 5.4+.
  - `do_not_modify` flag in `corefine`: 5.5+, useful when one fault
    is the canonical reference and the other should be refined to it
    without modifying the reference.

## Gotchas

- **Polyline resampling is NOT automatic in `isotropic_remeshing`.**
  Only `surface_Delaunay_remeshing` (with `polyline_constraints`)
  resamples the polyline. `isotropic_remeshing` either preserves the
  constraint exactly (`protect_constraints=true`) or allows collapse
  (`protect_constraints=false`). Neither **resamples uniformly** —
  if you want that, do it yourself before passing the polyline.
- **Both meshes' constraints must match for conformality after
  remeshing.** The safe path is: compute the polyline once
  (`surface_intersection`), resample to uniform spacing, and pass the
  *same* resampled polyline as `polyline_constraints` to both
  `surface_Delaunay_remeshing` calls. If you instead pass the
  per-mesh `edge_is_constrained_map` from corefine and let each
  remesher collapse independently, the two meshes can drift apart
  along the polyline and gmsh will reject them as non-conformal.
- **`corefine` is sensitive to self-intersections** in the inputs.
  Geological `.ts` files with overlapping triangles (common at fault
  perimeters) will throw `Self_intersection_exception` if
  `throw_on_self_intersection = true` is set, or silently produce a
  malformed result if not. Always pre-clean each input with
  `PMP::does_self_intersect` followed by `PMP::autorefine` if
  needed.
- **Fault perimeter (border) edges** are constrained by default in
  `isotropic_remeshing` (`patch boundary edges are always considered
  as constrained` per `remesh.h:101`). This is good — the perimeter
  of the fault stays put. But if you ALSO have an
  `edge_is_constrained_map` for the intersection, the union of the
  two constraint sets is what gets protected.
- **`face_patch_map` is auto-built from constrained edges if not
  provided.** This means the corefined intersection polyline
  automatically partitions A into "left of intersection" and "right
  of intersection" patches, and the remesher won't move triangles
  across the intersection. Usually what you want.
- **Floating-point drift on read-back**: when writing/reading STL
  through `CGAL::IO::write_polygon_mesh` / `read_polygon_mesh`,
  vertex coordinates are stored as `float`. If the remeshed
  polyline vertices are at coordinates with > 6-digit decimals,
  reading back through STL will produce slightly different points
  and break conformality. Use OFF or PLY for high-precision I/O,
  or pin coordinates to some known grid before writing STL.

## Open questions

- Does `surface_Delaunay_remeshing` produce identical polyline
  vertex sets when called twice with the same `polyline_constraints`
  but different input surfaces? The doxygen says polylines are
  *resampled* — verify empirically that two calls with the same
  polyline produce the same resampled vertex sequence.
- For the SAFS use case (6 faults, 8 intersecting pairs), is the
  pairwise corefine + per-fault Delaunay-remesh approach faster
  than concatenating everything into a single polygon soup and using
  `PMP::autorefine` + a single Delaunay remesh? Worth benchmarking
  on the 2000m resolution before committing to a pipeline structure.
- Does CGAL 6.1.1 expose `corefine_with_constrained_edges()` (a
  variant noted in some 6.x changelogs) that takes pre-existing
  constrained edges and preserves them through corefine? Not
  visible in `corefinement.h` of this version — search shows only
  the standard `corefine` and `corefine_and_compute_*` family.

## Next-step action items if we adopt this pipeline

1. Add `code_preprocess/corefine_cgal/` containing:
   - `corefine_pair.cpp` — single-pair corefine + Delaunay remesh,
     reads two STLs, writes two STLs.
   - `corefine_set.cpp` — N-fault driver that calls
     `intersecting_meshes` to discover pairs, walks the dependency
     graph, and produces a corefined set.
   - `CMakeLists.txt` configured for `cgal-61` env.
2. Update `code_tools/` Python driver to invoke the C++ binaries and
   stitch the output back into the existing `safs_fault_box.geo` flow.
3. Test plan:
   1. Build `corefine_pair` and run on the SBMT-Garnet_Hill ×
      SBMT-Mission_Creek 2000m pair (largest XY overlap, simplest
      geometry).
   2. Inspect the output STLs in ParaView; verify:
      (a) intersection polyline is shared (same coordinate values
      in both files at every polyline vertex);
      (b) triangle quality near the polyline is good (min angle
      > 20° as a rough target);
      (c) `ts_to_stl.py` no longer needs to run — the corefined STL
      is already at the right z-range and orientation.
   3. Run gmsh on a `.geo` that merges both STLs, embeds each as a
      `Surface{} In Volume{1}` with `Algorithm3D = 1`. Confirm "No
      ill-shaped tets" in gmsh output.
   4. Repeat for the SBMT-Garnet_Hill × SBMT-Mission_Creek ×
      SBMT-San_Andreas triple.
4. If the per-pair approach scales, extend to all 8 SAFS pairs at the
   2000m resolution before moving to 1000m / 500m.
