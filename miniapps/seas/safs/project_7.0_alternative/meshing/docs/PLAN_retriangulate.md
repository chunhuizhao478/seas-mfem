# Implementation Plan: Drop z > 0 vertices and re-triangulate `.ts` fault surfaces

## Overview

The CFM/Fuis raw `.ts` fault surfaces include vertices with positive elevation
(z > 0), which is incompatible with our simulation convention of a flat free
surface at z = 0. The existing tool `tools/ts_to_stl.py` handles this by
*per-triangle linear interpolation* — clipping each straddling triangle at the
z = 0 plane and emitting a chain of new edge-aligned vertices. The user wants
a different, cleaner approach for the input data itself: **drop every vertex
with z > 0 entirely, then re-triangulate the surviving point cloud** so the
top boundary is reconstructed by Delaunay rather than by linear edge cuts.

This plan adds a new tool `tools/retriangulate_ts.py` that emits a cleaned
`.ts` file as a drop-in replacement for the raw input. The downstream pipeline
(`ts_to_stl.py` → gmsh → MFEM) continues to work unchanged on the cleaned
file. The first acceptance target is

  `raw_data/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m.ts`

(2646 VRTX, 4729 TRGL, 296 vertices with z > 0). The implementation must
generalize to the other 1000m / 500m sibling files and to the
`project_7.0_preferred/raw_data/` family without code changes.

## Constraints

- **Input format**: GoCAD TSurf (`.ts`) ASCII. Only `VRTX <id> x y z` and
  `TRGL <a> <b> <c>` are present in the target file (no `ATOM`, no `PVRTX`),
  but the parser must still tolerate `ATOM <new> <ref>` and
  `PVRTX <id> x y z ...` because sibling files in `project_7.0_preferred/` may
  contain them. Reuse the parsing logic from
  `tools/ts_to_stl.py:parse_ts()`; do not re-derive it.
- **Output format**: Same GoCAD TSurf ASCII. Header (lines 1–13 of the
  source: `GOCAD TSurf 1`, `HEADER { … }`,
  `GOCAD_ORIGINAL_COORDINATE_SYSTEM … END_ORIGINAL_COORDINATE_SYSTEM`,
  `TFACE`) is copied verbatim from the input. Vertices are renumbered
  sequentially starting at 1. The file ends with `END`.
- **No external mesh libraries**. `numpy` and `scipy.spatial.Delaunay` only,
  matching the dependency footprint already used by other safs python tools.
  These are present in `conda activate pythonenv` (verified: scipy 1.16.2,
  numpy 2.3.3).
- **No hardcoded magic numbers**. Sample spacing, edge-length pruning
  thresholds, and PCA-quality warnings must derive from the input mesh's
  edge-length statistics, not from filename or constants. A single
  user-tunable multiplier (CLI flag) is allowed.
- **Z-cut convention**: keep vertex iff `z <= +EPS` (same `EPS = 1e-9` used by
  `ts_to_stl.py`). Vertices exactly at z = 0 are kept.
- **Triangle orientation**: each output triangle must carry an outward normal
  consistent with the *original* surface orientation. Concretely, if the
  Delaunay triangle's vertices in (u, v) order produce a normal n_uv in PCA
  coords, the lifted 3D triangle's vertex order must be reversed when
  necessary so that `dot(n_3d, n_pca_basis_w) > 0`, where w is the
  smallest-singular-vector direction (the fault normal) chosen so that the
  *majority* of original input triangle normals dot positively with it.
- **Style**: match the existing `ts_to_stl.py` — module docstring at top,
  argparse CLI, `main() -> int` returning shell exit code, helpful prints
  reporting counts and bounding boxes.

## Phase 1: Parsing and z > 0 filtering

### Goal
After this phase the tool reads a `.ts` file, drops vertices with `z > EPS`
along with every triangle that references one of them, prints diagnostics,
and writes a `.ts` file containing only the surviving vertices and triangles
(same connectivity, no re-triangulation yet). This is a verifiable
intermediate output: a "naive deletion" of the cap, with the jagged top
boundary that motivates Phase 2.

### Files to Create
- `miniapps/seas/safs/project_7.0_alternative/tools/retriangulate_ts.py`
  — new standalone CLI script. Phase 1 implements the parse → filter → emit
  pipeline; Phase 2 inserts the re-triangulation step before emission.

### Files to Modify
- None.

### Detailed Requirements

1. **CLI interface** (argparse):
   ```
   retriangulate_ts.py INPUT.ts OUTPUT.ts
       [--edge-length-multiplier FLOAT]   # Phase 2; ignored in Phase 1
       [--no-retriangulate]               # Phase 1 behavior even after Phase 2
       [--verbose]
   ```
   `INPUT` and `OUTPUT` are `pathlib.Path`. `--no-retriangulate` is a flag
   that short-circuits to the Phase 1 emit path (drop-only, original
   connectivity). `--edge-length-multiplier` defaults to `1.8` and is unused
   until Phase 2.

2. **Parser** — function:
   ```python
   def parse_ts(path: Path) -> tuple[dict[int, tuple[float, float, float]],
                                     list[tuple[int, int, int]],
                                     list[str]]:
   ```
   Returns `(verts, tris, header_lines)` where `header_lines` is the list of
   raw lines (with trailing newlines preserved) from the start of the file
   up to and including the first `TFACE` line. Vertex/triangle parsing is
   identical in semantics to `ts_to_stl.py:parse_ts()` (handle `VRTX`,
   `PVRTX`, `ATOM`, `TRGL`).

3. **Filter** — function:
   ```python
   def filter_above_z0(verts: dict, tris: list, eps: float = 1.0e-9
                      ) -> tuple[dict, list, dict]:
   ```
   Returns `(kept_verts, kept_tris, stats)`. A vertex is dropped iff
   `z > eps`. A triangle is dropped iff *any* of its three vertex IDs were
   dropped. `stats` is a dict with keys: `n_verts_in`, `n_verts_dropped`,
   `n_verts_kept`, `n_tris_in`, `n_tris_dropped_missing_vertex`,
   `n_tris_kept`, `z_min_kept`, `z_max_kept`, `bbox_kept` (a 6-tuple
   `(xmin, xmax, ymin, ymax, zmin, zmax)`).

4. **Renumber + emit** — function:
   ```python
   def write_ts(path: Path,
                header_lines: list[str],
                verts: dict[int, tuple[float, float, float]],
                tris: list[tuple[int, int, int]]) -> None:
   ```
   Writes a valid GoCAD TSurf file:
   - Copy `header_lines` verbatim.
   - Renumber vertices: build `id_old → id_new` where `id_new` runs 1..N in
     iteration order over `sorted(verts.keys())`.
   - Emit `VRTX <id_new>  <x> <y> <z>` lines using the same float format as
     the input (use `repr()` is acceptable; preferred is `f"{x:.7f}"` which
     preserves typical CFM precision).
   - Emit `TRGL <a_new> <b_new> <c_new>` for each triangle. Skip any
     triangle whose vertices are no longer in the renumber map (defensive;
     should be empty after `filter_above_z0`).
   - Final line `END\n`.

5. **`main()`**:
   - Parse CLI; verify input exists.
   - `verts, tris, header_lines = parse_ts(args.input)`.
   - `kept_verts, kept_tris, stats = filter_above_z0(verts, tris)`.
   - If `args.no_retriangulate` (or in Phase 1, always):
     `write_ts(args.output, header_lines, kept_verts, kept_tris)`.
   - Print: input file name, counts (in / dropped / kept) for vertices and
     triangles, kept bounding box, and the output path with file size.

6. **Header line copy** must include the `TFACE` line so the output is
   immediately valid as a single-face TSurf. The parser must stop appending
   to `header_lines` *after* it sees a line whose first token is `TFACE` —
   subsequent `VRTX`/`TRGL` lines are content, not header.

### Interfaces
- Public functions exposed at module level for unit tests:
  `parse_ts`, `filter_above_z0`, `write_ts`, `main`.
- Module-level constant `EPS = 1.0e-9`, matching `ts_to_stl.py`.

### Edge Cases to Handle
- Vertex exactly at z = 0 → kept (`z > eps` is strict).
- Triangle with one or more vertices not present in `verts` (e.g., dangling
  reference) → dropped, counted in `n_tris_dropped_missing_vertex`.
- Empty surviving set (every vertex above z = 0) → print error to stderr,
  return exit code 1, do not write output.
- Duplicate vertex IDs in input → last one wins (matches Python dict
  semantics; same behavior as `ts_to_stl.py`).
- `ATOM` aliases — already resolved at parse time into `verts`; treat
  resolved entries identically.

### Acceptance Criteria
- [ ] `python retriangulate_ts.py raw_data/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m.ts /tmp/out_phase1.ts --no-retriangulate`
      runs without error and prints `n_verts_dropped: 296`,
      `n_verts_kept: 2350`, `n_tris_kept: 4137` (4729 minus the 592
      triangles that reference at least one of the 296 dropped vertices —
      exact count to be verified at runtime; the test asserts
      `n_tris_dropped_missing_vertex == n_tris_in - n_tris_kept` and
      `n_tris_kept > 0.8 * n_tris_in`).
- [ ] `awk '/^VRTX/ {print $5}' /tmp/out_phase1.ts | sort -g | tail -1`
      reports a value `<= 1.0e-9`.
- [ ] Vertex IDs in the output are dense and sequential 1..N (verifiable
      via `awk '/^VRTX/{print $2}' /tmp/out_phase1.ts | uniq -c | awk '{if($1!=1)bad++}END{print bad+0}'` returning `0`,
      and `awk '/^VRTX/{c++}/^VRTX/{if($2!=c){print "gap at",c;exit}}' /tmp/out_phase1.ts` printing nothing).
- [ ] Every TRGL index in the output is in range [1, N].
- [ ] Header section (everything up to and including `TFACE`) is byte-equal
      to the corresponding lines in the input.
- [ ] The existing `ts_to_stl.py` runs successfully on the cleaned file:
      `python tools/ts_to_stl.py /tmp/out_phase1.ts /tmp/out_phase1.stl`
      reports `Vertices with z > 0: 0` and writes a non-empty STL.

### Dependencies
- Depends on: nothing.
- Required by: Phase 2 (re-triangulation) extends the same script.

## Phase 2: Re-triangulate the surviving point cloud

### Goal
After this phase the default (no `--no-retriangulate`) behavior emits a `.ts`
file in which the top boundary is no longer the jagged remnant of dropping
cap-touching triangles, but a fresh Delaunay triangulation of the surviving
point cloud projected onto the fault's best-fit plane, with sliver and
long-edge artifacts pruned. The new triangulation preserves the original
surface orientation (outward normal direction).

### Files to Create
- None (extends the script from Phase 1).

### Files to Modify
- `miniapps/seas/safs/project_7.0_alternative/tools/retriangulate_ts.py`
  — add the re-triangulation pipeline and wire it into `main()` as the
  default path.

### Detailed Requirements

1. **PCA parametrization** — function:
   ```python
   def fit_fault_plane(verts: dict[int, tuple[float, float, float]]
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
   ```
   Returns `(centroid, u_axis, v_axis, w_axis, sigma)` where:
   - `centroid` is the mean of the input points (shape `(3,)`).
   - `u_axis`, `v_axis`, `w_axis` are unit 3-vectors forming a right-handed
     orthonormal basis.
   - The basis is the right singular vectors V of `X - centroid` from
     `numpy.linalg.svd(..., full_matrices=False)`. `u_axis = V[0]` (largest
     singular value), `v_axis = V[1]`, `w_axis = V[2]` (smallest, the fault
     normal). After computation, flip `v_axis` so `dot(v_axis, [0,0,-1]) >= 0`
     (in-plane "down" direction), and recompute `u_axis = cross(v_axis,
     w_axis)` with the appropriate sign to keep the basis right-handed.
   - `sigma` is the length-3 array of singular values.

2. **Quality warning** — if `sigma[2] / sigma[0] > 0.20`, print a warning
   to stderr: the surface deviates significantly from a single best-fit
   plane and the (u, v) parametrization may produce overlapping triangles
   when lifted back to 3D. Continue anyway. Threshold `0.20` is a module
   constant `PLANARITY_WARN_RATIO`.

3. **Project to (u, v)** — function:
   ```python
   def project_to_plane(verts: dict, centroid, u_axis, v_axis
                       ) -> tuple[list[int], np.ndarray]:
   ```
   Returns `(ids, uv)` where `ids` is the list of vertex IDs in iteration
   order over `sorted(verts.keys())` and `uv` is an `(N, 2)` array with
   `uv[k] = ((p_k - centroid) · u_axis, (p_k - centroid) · v_axis)`.

4. **Compute reference edge length** — function:
   ```python
   def reference_edge_length(verts, tris) -> float:
   ```
   Returns the *median* edge length of the input mesh (post-filter) using
   3D Euclidean distance on the surviving triangles. Build the edge set
   `frozenset({a, b})` for each pair to avoid double-counting. This number
   is the "natural sample spacing" of the data and is used as the pruning
   threshold base in step 6.

5. **Delaunay** — call `scipy.spatial.Delaunay(uv)`. The result `.simplices`
   is an `(M, 3)` int array indexing rows of `uv`; map back through `ids`
   to original vertex IDs to obtain candidate triangles.

6. **Prune long-edge triangles** — drop any triangle whose maximum 3D edge
   length exceeds `args.edge_length_multiplier * reference_edge_length`.
   Default multiplier `1.8` chosen so that a regularly-spaced 2000 m grid
   admits diagonal edges of `~2828 m` (multiplier 1.414+) while rejecting
   long shortcuts across the convex hull's concavities. Print the pruning
   ratio (kept / candidate) and the threshold used.

7. **Sliver pruning** — additionally drop any triangle with 2D area in
   `(u, v)` below `1e-6 * (median_2d_edge_length ** 2)`. These arise when
   four nearly-collinear vertices produce a near-degenerate Delaunay cell.

8. **Orientation fix** — for each surviving triangle:
   - Compute the candidate 3D normal `n_cand = cross(p1 - p0, p2 - p0)`
     (un-normalized).
   - Compute the *reference* fault outward normal once, before the
     pruning loop, as the sign-corrected `w_axis`. Sign correction:
     for each input triangle `(a, b, c)` whose vertices all survive,
     accumulate `sum(sign(dot(cross(p_b - p_a, p_c - p_a), w_axis)))`.
     If the sum is negative, replace `w_axis ← -w_axis`. This makes
     `w_axis` agree with the *majority* of input-triangle outward
     normals.
   - If `dot(n_cand, w_axis) < 0`, swap two vertex indices to flip the
     orientation of the output triangle.
   This guarantees the new mesh's per-triangle normal direction agrees with
   the input's, which is required for downstream STL/gmsh tooling.

9. **Wire into `main()`**: when `--no-retriangulate` is not set, after
   `filter_above_z0`:
   - `centroid, u_axis, v_axis, w_axis, sigma = fit_fault_plane(kept_verts)`.
   - Sign-correct `w_axis` per step 8.
   - Print sigma values and the planarity ratio.
   - `ids, uv = project_to_plane(...)`.
   - `ref_edge = reference_edge_length(kept_verts, kept_tris)`.
   - Run Delaunay, prune, fix orientation.
   - Build a new triangle list using *original* vertex IDs from
     `kept_verts` (not new IDs — `write_ts` will renumber).
   - `write_ts(args.output, header_lines, kept_verts, new_tris)`.
   - Print: `n_tri_delaunay`, `n_tri_pruned_long_edge`, `n_tri_pruned_sliver`,
     `n_tri_orientation_flipped`, `n_tri_out`.

10. **Diagnostics flag** — if `--verbose`, additionally print:
    - Histogram (10 bins) of input edge lengths and output edge lengths.
    - The PCA basis vectors (u, v, w) as 3-vectors.
    - The min / max / median 2D edge length.

### Interfaces
- New module-level functions: `fit_fault_plane`, `project_to_plane`,
  `reference_edge_length`, `delaunay_retriangulate`, `prune_triangles`,
  `fix_orientation`.
- Module-level constants: `EPS = 1.0e-9`,
  `PLANARITY_WARN_RATIO = 0.20`,
  `DEFAULT_EDGE_MULTIPLIER = 1.8`,
  `SLIVER_AREA_FRACTION = 1.0e-6`.

### Edge Cases to Handle
- Highly non-planar surface (`sigma[2] / sigma[0] > 0.5`): still produce
  output, but emit a stronger warning that lifted triangles may
  self-intersect; recommend using `--no-retriangulate` instead.
- Fewer than 3 surviving vertices: error out with message
  `not enough vertices to triangulate (N=<n>)`, exit 1.
- Delaunay failure (degenerate point cloud): catch
  `scipy.spatial.qhull.QhullError`, print the message, exit 2.
- Duplicate (u, v) points (two different 3D points project to the same
  plane location): scipy handles this by raising; pre-deduplicate by
  rounding (u, v) to 6 decimals and keeping the first; report dedup count.
- All input triangles end up dropped by long-edge pruning (multiplier too
  tight): print error suggesting a larger multiplier, exit 1.

### Acceptance Criteria
- [ ] `python tools/retriangulate_ts.py raw_data/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m.ts /tmp/out_phase2.ts`
      runs without error.
- [ ] Output has `n_verts == 2350` (same vertex set as Phase 1, just
      re-renumbered).
- [ ] Every output VRTX has `z <= 1e-9`.
- [ ] Every output TRGL references vertex indices in `[1, 2350]`.
- [ ] Output triangle count is within `[0.7, 1.3] × 4137` (i.e.,
      re-triangulation does not radically change connectivity density).
- [ ] No output triangle has any edge longer than
      `args.edge_length_multiplier * reference_edge_length`.
- [ ] At most 5 % of input triangles change orientation (those whose
      original normal was inconsistent with the majority); the per-triangle
      normal direction agrees with `w_axis` for ≥ 95 % of output triangles.
- [ ] `python tools/ts_to_stl.py /tmp/out_phase2.ts /tmp/out_phase2.stl`
      reports `Vertices with z > 0: 0`, `Triangles dropped: 0`,
      `Triangles clipped: 0`. (All work done by re-triangulation; nothing
      left for the clipper.)
- [ ] Visual check: open `/tmp/out_phase2.stl` in ParaView/Gmsh; the top
      edge at z = 0 is closed and continuous; no triangles cross or
      self-intersect noticeably; the global surface looks like the
      original San Andreas fault Fuis-ALT6 model with the cap removed.

### Dependencies
- Depends on: Phase 1 (parser, filter, writer).
- Required by: Phase 3.

## Phase 3: Apply to remaining sibling files and document

### Goal
The cleaned `.ts` files for all three resolutions of the
San_Andreas_fault_Fuis-ALT6 model exist in
`raw_data_clipped/`, and a one-line README in
`raw_data_clipped/README.md` records the source command and date. The
existing downstream pipeline (`safs_fault_box.geo`, gmsh meshing) is shown
to consume the cleaned file without modification.

### Files to Create
- `miniapps/seas/safs/project_7.0_alternative/raw_data_clipped/` — new
  directory holding cleaned `.ts` outputs:
  - `SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m.ts`
  - `SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_1000m.ts`
  - `SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_500m.ts`
- `miniapps/seas/safs/project_7.0_alternative/raw_data_clipped/README.md` —
  short note: how the files were generated (the exact command line),
  the multiplier used, and the date.

### Files to Modify
- None.

### Detailed Requirements

1. Run the tool on each of the three resolutions, defaulting the
   multiplier to 1.8. If the 500 m file produces unacceptable artifacts
   at default settings (visible in ParaView), tune the multiplier per
   file and record the chosen value in the README.

2. Verify the gmsh pipeline still works: regenerate the STL with
   `tools/ts_to_stl.py` from the cleaned 2000 m `.ts`, then run gmsh on
   `mesh/safs_fault_box.geo` (or the appropriate `.geo` file) referencing
   the new STL. Confirm the volume mesh produces the same topology as the
   pre-existing run.

3. Write the README. Include:
   ```
   Generated: <YYYY-MM-DD>
   Source: ../raw_data/<filename>
   Command: python ../tools/retriangulate_ts.py <in> <out>
   Multiplier: <value>
   Output stats: N_verts, N_tris, bbox
   ```

### Acceptance Criteria
- [ ] All three cleaned `.ts` files are present and pass Phase 1 acceptance
      (no z > 0 vertices, dense numbering, etc.).
- [ ] `tools/ts_to_stl.py` on each cleaned file reports zero clipped /
      dropped triangles.
- [ ] gmsh successfully meshes the resulting volume on at least the 2000 m
      input.
- [ ] README is present and accurate.

### Dependencies
- Depends on: Phase 2.
- Required by: nothing.

## Testing Strategy

Phase-1 testable behaviors:
- Parser round-trip: parse → write_ts (no filter, all vertices kept) →
  re-parse should give an isomorphic vertex/triangle set.
- Filter unit test on a hand-crafted 6-vertex / 4-triangle synthetic
  surface where vertices 5 and 6 have z > 0; expected: kept_verts has IDs
  {1, 2, 3, 4}, kept_tris has only the triangles with all indices in that
  set.
- Renumber test: input vertex IDs {3, 7, 11} should renumber to {1, 2, 3}
  in sorted order, and the TRGL referring to ID 7 should now refer to 2.

Phase-2 testable behaviors:
- PCA on a synthetic planar point cloud (z = 0 plane plus tiny noise)
  returns `w_axis ≈ ±[0, 0, 1]` and `sigma[2] / sigma[0]` near zero.
- Orientation fix: a synthetic point cloud whose canonical Delaunay
  triangle ordering produces normals opposite to the reference is flipped;
  after the fix, ≥ 95 % of triangle normals align with w_axis.
- Long-edge pruning: a synthetic point cloud with one outlier far from
  the rest produces Delaunay shards that the pruner removes; without the
  outlier the pruner removes 0 triangles.

Integration / end-to-end:
- The full pipeline applied to
  `SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m.ts` produces an output
  that, when fed into `ts_to_stl.py`, reports zero clipping work — the cap
  is already gone.
- ParaView visual diff: original (with cap) vs cleaned (no cap) should
  show the fault body identical except for the missing top sliver.

Suggested test files (created together with the script):
- `miniapps/seas/safs/project_7.0_alternative/tools/test_retriangulate_ts.py`
  using `pytest` with the synthetic fixtures listed above. Run via
  `conda activate pythonenv && pytest tools/test_retriangulate_ts.py`.

## Risk Assessment

- **Non-planar fault surfaces** (e.g., curved-along-strike SAF segments)
  are the dominant failure mode. The PCA fit smears curvature into
  apparent thickness in the third singular direction; lifted Delaunay
  triangles may self-intersect when the fault folds back on itself in
  the (u, v) projection. Mitigation: the planarity warning in step 2 of
  Phase 2; if it fires for a real input, fall back to `--no-retriangulate`
  and accept the jagged top, or switch to an arclength-along-strike
  parametrization (`s, z`) — out of scope for this plan but a natural
  Phase 4 if needed.

- **Edge-length pruning with non-uniform sampling**: `*.500m.ts` and
  `*.2000m.ts` have very different sample densities. The reference edge
  length is computed *per file*, so the absolute threshold differs
  appropriately, but if a single file mixes resolutions (e.g., refined
  near the surface trace and coarse at depth), the median may be too tight
  for the coarse region or too loose for the refined region. Mitigation:
  expose `--edge-length-multiplier`; document tuning in the Phase 3
  README.

- **Coordinate precision loss in `.ts` writer**: GoCAD `.ts` files use
  ~7 decimal digits of precision in the original CFM data. Use
  `f"{x:.7f}"` in `write_ts`; do not use `repr()` (varies by Python
  version) or `f"{x}"` (which can drop trailing zeros). Verify by
  diffing a parse → write_ts → parse round-trip on a 100-vertex slice.

- **Header semantics**: the GoCAD coordinate-system block carries axis
  directions and units. We copy it verbatim, so as long as the input is
  valid, the output remains valid. We do not parse or transform
  coordinates; the cleaned file inherits the same coordinate system.

- **Concavity at the new top edge**: Delaunay's convex hull will fill any
  bay or notch in the original surface trace at z = 0 with a long-edge
  triangle, which is then dropped by the long-edge pruner. The user should
  verify visually that the resulting top edge is what they expect; in
  particular, if two segments of the SAF are separated by a gap in the
  raw data, the pruner should leave that gap open.

- **Tooling: `tools/ts_to_stl.py` already does its own clipping**. After
  Phase 2 the cleaned `.ts` has no z > 0 content, so `ts_to_stl.py` will
  pass everything through unchanged. This is *not* a regression — the
  acceptance criterion checks that the clipping counts are zero. We
  do *not* need to modify `ts_to_stl.py`.
