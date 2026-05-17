# Implementation Plan: NW Hard-Cut of the Alternative ALT6 Long-Strip Fault STLs

## Overview

The alternative SAF representation
`SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_*_clean_clip.stl` is a long strip
that extends ~236 km further to the north-west than the preferred SAF
representation
`SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6_2000m_clean_clip.stl`. We will write a
single Python driver that (a) reads the preferred mesh and computes the
position of its NW-most surface point, (b) defines a vertical cutting plane
whose horizontal normal is the NW direction and which passes through that
point, and (c) hard-cuts every alternative ALT6 STL in
`project_7.0_alternative/data_cleanfreesurf/`, removing all geometry on the
NW side of the plane and writing a new `_nwcut.stl` next to each input.
Triangles that straddle the plane are clipped via linear interpolation along
each crossed edge so the new free boundary lies exactly on the plane (this
mirrors the existing `clip_triangle_at_z0` idiom in
`code_preprocess/ts_to_stl.py`).

## Constraints

### Interface / file constraints

- **Do not modify the preferred file** at any time; it is read-only.
- **Do not modify** the existing alternative `_clean_clip.stl` files; outputs
  go to *new* files at the same directory.
- **Do not modify shared BP5 / fault solver / friction code** (memory rule
  C2). All work lives under
  `project_7.0_alternative/code_preprocess/` and `…/document/`.
- The `data_cleanfreesurf/` filename pattern `<base>_clean_clip.stl` is
  consumed downstream by `code_meshing/safs_fault_box*.geo`. Outputs of this
  task must be a *new* basename (`<base>_clean_clip_nwcut.stl`) so existing
  meshing pipelines that look up the un-cut filename keep working unchanged.

### Convention constraints (mirror existing preprocess scripts)

- ASCII STL only (matches what `clean_freesurface_mesh.py` writes via
  `ms.save_current_mesh(..., binary=False)` and what
  `code_meshing/safs_fault_box*.geo` consumes).
- Use the project's existing conda env: `conda activate pythonenv` for
  any execution. Required libs: `numpy` (mandatory), `pymeshlab` (optional
  cleanup pass; only used if available — see Phase 3, step 5).
- `argparse`-driven CLI with a `--batch` mode, mirroring
  `clean_freesurface_mesh.py`.
- File header docstring in the same explanatory style as
  `clean_freesurface_mesh.py`.
- Use the same triangle-clip idiom as `ts_to_stl.clip_triangle_at_z0`
  (1-below / 2-below cases, vertex order preserved). We will *generalise* it
  to a generic plane, not duplicate it.

### Numerical constraints (measured on the actual data, see Phase 1)

- The preferred mesh
  `SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6_2000m_clean_clip.stl` has bbox
  `x ∈ [365072.3, 428025.6]`, `y ∈ [3808374.0, 3838985.0]`, `z ∈ [-13272.1, 0.0]`.
- Its NW-most surface point (max of `−x + y` over all vertices) is at
  `(x, y, z) = (365072.3, 3838982.0, 0.0)`. This is the **anchor**.
- The alternative ALT6 strip (all three resolutions: 500 m, 1000 m, 2000 m)
  extends to `(x, y, z) = (179309.0, 3987261.0, -9582.6)` and
  `y_max ≈ 3988004` at the surface — i.e., the NW projection
  `s = (−x + y) / √2` reaches ≈ 2,692,628.7 vs. the preferred's
  ≈ 2,456,425.1, a ~236 km overshoot along the NW axis.
- `EPS = 1.0e-9` (consistent with `ts_to_stl.EPS`) for float comparisons of
  the signed plane distance.

## Definitions used throughout

- **NW direction (horizontal unit vector):** `n_h = (−1, +1, 0) / √2` in
  UTM zone 11 N (the file is in metres of UTM Easting/Northing). "Left,
  upper" in plan view = West + North = `+n_h` direction.
- **NW projection of a point `p = (x, y, z)`:**
  `s(p) = (−x + y) / √2`. (We deliberately ignore Z so the cutting plane is
  vertical — see "Cutting plane" below.)
- **Anchor point `p_anchor`:** vertex of the *preferred* mesh that maximises
  `s(·)`. Numerically this is `(365072.3, 3838982.0, 0.0)` (Phase 1
  computes it from the file at run time, not as a hard-coded constant).
- **Cutting plane Π:** vertical plane (parallel to the Z axis) passing
  through `p_anchor` with horizontal normal `n_h`. Equivalently:
  `Π = { p : n_h · (p − p_anchor) = 0 }`, i.e.
  `(−x + y) − (−x_anchor + y_anchor) = 0` (the `√2` divisor cancels).
- **Signed distance from Π:** `d(p) = (−x + y) − (−x_anchor + y_anchor)`.
  Sign convention:
  - `d(p) > 0`: point lies NW of the plane → **DROP** (further left + upper
    than the preferred mesh's NW edge).
  - `d(p) ≤ 0`: point lies SE of the plane → **KEEP**.
- **Half-space to keep:** `H_keep = { p : d(p) ≤ 0 }`. The clipped output
  is the intersection of the input mesh with `H_keep`.

The plane is purely vertical (its third row in normal form has `n_z = 0`),
which is what the user asked for ("the most north-west location … cut at
that location"). A dipping cut would couple Z to the cut location and is
not requested.

## Phase 1 — NW-anchor extraction utility

### Goal

After this phase a single function exists that, given any STL/PLY/OBJ path,
returns the NW-most vertex (the "anchor") and the plane equation that hard-
cuts everything NW of it. The function is import-safe (no side effects) so
Phase 3 can reuse it programmatically.

### Files to Create

- `project_7.0_alternative/code_preprocess/nw_cut_strip.py` — single new
  module that owns Phases 1–3. (We do *not* split into one file per phase;
  the existing preprocess scripts are also single-file.)

### Files to Modify

- *(none)*

### Detailed Requirements

1. Add a module-level constant
   `NW_DIRECTION_XY: np.ndarray = np.array([-1.0, 1.0]) / np.sqrt(2.0)`
   used everywhere so the convention is defined once.
2. Add a function
   ```python
   def load_vertices(path: pathlib.Path) -> np.ndarray:
       """Return an (N, 3) float64 array of vertex coordinates from any
       mesh file readable by `meshio.read`. ASCII STL, binary STL, PLY, OBJ,
       and `.ts` (via the existing `ts_to_stl.parse_ts`) must all work; raise
       FileNotFoundError if the file is missing and ValueError if no
       triangles are present."""
   ```
   - **Orphan-vertex filter (applies to all branches):** the docstring
     promises an `(N, 3)` array of *triangle-referenced* vertices, so
     orphan vertices that no triangle references must be stripped. This
     is needed because pymeshlab-cleaned STLs and raw `.ts` files can
     both carry orphaned points (or ATOM definitions) that drift outside
     the mesh and would otherwise pollute the NW anchor.
   - `.ts` branch: import `parse_ts` from `ts_to_stl` (same `sys.path.insert`
     trick already used at the top of `clean_freesurface_mesh.py`).
     `parse_ts` returns `(verts: dict, tris: list)` with 1-indexed IDs.
     Build the set of triangle-referenced IDs
     `referenced = {i for ia, ib, ic in tris for i in (ia, ib, ic)}`,
     then stack `verts[i]` for each `i in referenced` (in sorted order
     for determinism), as a `(N, 3)` float64 array.
   - All other suffixes: `meshio.read(str(path))`. Build the set of
     triangle-referenced indices by gathering the unique vertex indices
     used by every `m.cells` block of type `"triangle"`; slice `m.points`
     with that index array (cast to `np.float64`).
   - The output is always a contiguous `(N, 3)` `np.float64` array. If
     no triangle cells are present, raise `ValueError(f"no triangles in
     {path}")`. If the path does not exist, let the underlying loader
     raise `FileNotFoundError`.
3. Add a function
   ```python
   def nw_anchor(verts: np.ndarray) -> tuple[np.ndarray, float]:
       """Return (p_anchor, s_max) where p_anchor is the (3,) vertex of the
       input that maximises s(p) = (-x + y) / sqrt(2), and s_max is that
       maximum projection. If two vertices tie within EPS, the lower-Z one
       wins (deterministic + matches 'top trace' intuition)."""
   ```
   Implementation: compute `s = verts[:, :2] @ NW_DIRECTION_XY` (shape
   `(N, 2) @ (2,) → (N,)`); let `s_max = float(s.max())`;
   `mask = s >= s_max - EPS`; among the masked rows, pick `argmin(z)`.
   Return that row (as a `(3,)` float64 array) and `s_max`.

   **Note on operand order:** the converse expression
   `NW_DIRECTION_XY @ verts[:, :2]` raises `ValueError` because numpy
   promotes the 1-D LHS to `(1, 2)` and then requires the inner dim of
   `(N, 2)` to equal 2, which fails for any real mesh (N = 1479 / 13722
   / 54168 / 214938). Whenever this plan re-computes the projection (e.g.
   the post-cut `s_max_after` in Phase 3 step 3.9 and the smoke-test
   assertion in Phase 4 §5), use the same `verts[:, :2] @
   NW_DIRECTION_XY` form.
4. Add a function
   ```python
   def cutting_plane(p_anchor: np.ndarray) -> tuple[np.ndarray, float]:
       """Return (n, c) such that the cutting plane is { p : n · p = c },
       with n = (-1, +1, 0)/sqrt(2) (horizontal NW normal) and
       c = n · p_anchor. The 'keep' half-space is { p : n · p <= c }."""
   ```
   This is the canonical form Phase 2 will consume.
5. Make the module runnable as `python nw_cut_strip.py --print-anchor PATH`
   so the user can sanity-check the anchor without invoking the cut.
   Output a single line like
   `anchor: x=3.650723e+05  y=3.838982e+06  z=0.000000e+00  s=2.456425e+06`.

### Interfaces

- `load_vertices(path) -> np.ndarray (N, 3)`
- `nw_anchor(verts) -> (p_anchor: np.ndarray (3,), s_max: float)`
- `cutting_plane(p_anchor) -> (n: np.ndarray (3,), c: float)`
- Module-level constants `NW_DIRECTION_XY`, `EPS = 1.0e-9`.

### Edge Cases to Handle

- File missing → `FileNotFoundError` with the path.
- File present but contains no triangles (e.g. point cloud) → `ValueError`.
- File contains orphan vertices (in `m.points` but no triangle uses them) →
  ignore them in the anchor computation (see step 2).
- Multiple vertices tie for the max NW projection (within `EPS`) → break the
  tie by smallest Z (step 3). This is rare but happens when a mesh edge is
  exactly NW-aligned.
- All vertices coplanar in NW direction (e.g. user passes a line mesh) →
  `s_max` is well-defined; `nw_anchor` still returns one vertex; downstream
  Phase 2 will then drop the entire alternative mesh, which is the correct
  behaviour and *not* an error.

### Acceptance Criteria

- [ ] `load_vertices` of `…/SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6_2000m_clean_clip.stl`
      returns an array of shape `(1479, 3)` with bbox
      `x ∈ [365072.3, 428025.6]`, `y ∈ [3808374.0, 3838985.0]`,
      `z ∈ [-13272.1, 0.0]` (all to 1 m).
- [ ] `nw_anchor(verts)` of the same file returns
      `p_anchor ≈ (365072.3, 3838982.0, 0.0)` (each component within 1 m)
      and `s_max ≈ 2456425.1` (within 1 m).
- [ ] `cutting_plane(p_anchor)` returns `n = (-0.7071…, 0.7071…, 0.0)`
      (within 1e-12 of `±1/√2`) and `c = s_max` (within 1e-9).
- [ ] `python nw_cut_strip.py --print-anchor <preferred-stl>` exits 0 and
      prints the line shown in step 5 above.

### Dependencies

- Depends on: nothing.
- Required by: Phase 2, Phase 3.

## Phase 2 — Generic vertical-plane triangle clipper

### Goal

After this phase a clipper exists that, given any triangle and the plane
returned by Phase 1, emits the (zero, one, or two) sub-triangles of the
triangle that lie in `H_keep`. Behaviour mirrors
`ts_to_stl.clip_triangle_at_z0` but the plane is generic, not the z = 0
plane.

### Files to Modify

- `project_7.0_alternative/code_preprocess/nw_cut_strip.py` — append the new
  helpers (do not import them from `ts_to_stl`; we want a clean dependency
  for the generic plane).

### Detailed Requirements

1. Add helper
   ```python
   def signed_dist(p: tuple, n: np.ndarray, c: float) -> float:
       """Return n · p − c. > 0 means 'NW of plane' (drop); ≤ 0 means
       'SE of plane' (keep)."""
   ```
2. Add helper
   ```python
   def lerp_to_plane(p_keep: tuple, p_drop: tuple,
                     n: np.ndarray, c: float) -> tuple[float, float, float]:
       """Linearly interpolate between p_keep (d ≤ 0) and p_drop (d > 0) to
       find the point where the segment crosses the plane n · p = c."""
   ```
   Implementation: `d_keep = signed_dist(p_keep, n, c)`,
   `d_drop = signed_dist(p_drop, n, c)`. If
   `abs(d_drop − d_keep) < EPS` (segment is parallel to or in the plane
   to within float64 noise — both endpoints have signed distances within
   `EPS` of each other regardless of sign; the segment cannot meaningfully
   cross the plane), return `p_keep` snapped onto the plane (subtract
   `d_keep * n` from `p_keep`; correctness check:
   `signed_dist(p_keep − d_keep·n, n, c) = (n·p_keep − d_keep·||n||²) − c
   = (c + d_keep) − d_keep·1 − c = 0` since `||n|| = 1`).
   Else `t = -d_keep / (d_drop - d_keep)` (matches `ts_to_stl.lerp_to_z0`'s
   `t = -zb / (za - zb)`), return
   `p_keep + t * (p_drop - p_keep)` as a 3-tuple. The interpolated point
   *should* satisfy `signed_dist(., n, c) ≈ 0` to within `EPS * (|c| + 1)`;
   assert this in a debug build (see Phase 4 tests).
3. Add the per-triangle clipper
   ```python
   def clip_triangle_at_plane(p0, p1, p2, n: np.ndarray,
                              c: float) -> list[tuple]:
       """Clip a triangle to the half-space n · p <= c. Returns 0, 1, or 2
       sub-triangles preserving original CCW vertex order (so outward normals
       are unchanged)."""
   ```
   Algorithm copied verbatim from `ts_to_stl.clip_triangle_at_z0`, with the
   *only* differences being:
   - replace `z[i] = pts[i][2]` with `d[i] = signed_dist(pts[i], n, c)`,
   - replace `below = [zi <= EPS for zi in z]` with
     `keep = [di <= EPS for di in d]`,
   - replace `lerp_to_z0(...)` calls with `lerp_to_plane(..., n, c)`,
   - the `n_above == 1` (= `n_keep == 2`) and `n_above == 2` (= `n_keep == 1`)
     branches stay byte-for-byte identical otherwise (vertex labelling and
     output triangle ordering).
4. Provide a sibling array form for batch use:
   ```python
   def clip_mesh_at_plane(verts: np.ndarray, faces: np.ndarray,
                          n: np.ndarray, c: float
                          ) -> tuple[np.ndarray, np.ndarray, dict]:
       """verts: (N, 3) float64. faces: (M, 3) int32 indices into verts.
       Returns (V_out, F_out, stats) where V_out is a possibly-duplicated
       vertex array (caller welds with pymeshlab if available), F_out is
       (M', 3), and stats contains:
         - n_in (= M)
         - n_out (= M')
         - n_kept_whole, n_dropped, n_clipped_2to1, n_clipped_1to2
       (same five counts that `clean_freesurface_mesh.clip_to_arrays`
       reports). The case names track *number-above-the-plane*
       (= number on the drop side, where d > 0), exactly as in
       `clean_freesurface_mesh.clip_to_arrays`:
         - n_above == 0 (all on keep side): kept_whole
         - n_above == 1 (2 keep + 1 drop): emit 2 triangles -> '1to2'
         - n_above == 2 (1 keep + 2 drop): emit 1 triangle  -> '2to1'
         - n_above == 3 (all on drop side): dropped."""
   ```
   This matches the tally produced by
   `clean_freesurface_mesh.clip_to_arrays` so the user gets familiar
   diagnostic output.

### Interfaces

- `signed_dist(p, n, c) -> float`
- `lerp_to_plane(p_keep, p_drop, n, c) -> (x, y, z)`
- `clip_triangle_at_plane(p0, p1, p2, n, c) -> list[(p, p, p)]`
- `clip_mesh_at_plane(verts, faces, n, c) -> (V_out, F_out, stats)`

### Edge Cases to Handle

- All three triangle vertices have `d ≤ EPS` (entire triangle SE of or on
  the plane): return `[(p0, p1, p2)]` unchanged.
- All three `d > EPS` (entire triangle NW): return `[]`.
- One vertex exactly on plane (`|d| ≤ EPS`) and two on the keep side: treat
  the on-plane vertex as keep (this is what the `≤ EPS` test does); no
  interpolation needed.
- One vertex exactly on plane and two on drop side: treat on-plane as keep
  → produces a degenerate triangle with an interpolation point coincident
  with the on-plane vertex. The downstream `meshing_remove_null_faces`
  pymeshlab call (Phase 3 step 5) cleans these up. Document this.
- `lerp_to_plane` denominator near zero (segment lies in the plane): snap
  `p_keep` onto the plane and return it (step 2). This degenerate case is
  the `clip_triangle_at_z0` analogue's behaviour.
- Numerical cancellation when `s_max` is large (~2.5e6): `(-x + y) - c` is
  evaluated in float64, which loses ≤ 1 mm of precision for our inputs
  (max coord ≈ 4e6 m → rel. eps ≈ 4e6 × 2.22e-16 ≈ 1e-9 m). Acceptable for
  a 500 m–2000 m mesh.

### Acceptance Criteria

- [ ] Unit test: triangle entirely SE of plane → `clip_triangle_at_plane`
      returns the single input triangle unchanged.
- [ ] Unit test: triangle entirely NW of plane → returns `[]`.
- [ ] Unit test: triangle with 2 verts SE, 1 NW → returns 2 triangles whose
      union covers the SE quadrilateral; the two new vertices satisfy
      `|signed_dist(., n, c)| < EPS_TEST = 1e-6`; CCW orientation is
      preserved (cross product has the same sign as the input).
- [ ] Unit test: triangle with 1 vert SE, 2 NW → returns 1 triangle whose
      two new vertices both lie on the plane (within `EPS_TEST`) and whose
      orientation matches the input.
- [ ] Numerical-equivalence test: when `n = (0, 0, 1)` and `c = 0`,
      `clip_triangle_at_plane(p0, p1, p2, n, c)` returns the same triangles
      as `ts_to_stl.clip_triangle_at_z0(p0, p1, p2)` (up to vertex order
      within each output triangle).

### Dependencies

- Depends on: Phase 1 (`EPS`, `cutting_plane`).
- Required by: Phase 3.

## Phase 3 — Driver: hard-cut every ALT6 long-strip STL

### Goal

After this phase, running
`python nw_cut_strip.py --batch` (with no other flags) reads
`SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6_2000m_clean_clip.stl` from the
*preferred* directory, computes the cutting plane once, applies it to every
`SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_*_clean_clip.stl` in the
*alternative* `data_cleanfreesurf/`, and writes a sibling
`<base>_clean_clip_nwcut.stl` for each. Originals are not touched; the
existing meshing pipeline (which references `_clean_clip.stl` by name) is
unaffected unless someone re-points it at the `_nwcut` files.

### Files to Modify

- `project_7.0_alternative/code_preprocess/nw_cut_strip.py` — add CLI and
  the batch driver.

### Detailed Requirements

1. Hard-coded path defaults (overridable via flags):
   ```python
   DEFAULT_REFERENCE = (
       Path(__file__).resolve().parents[2]
       / "project_7.0_preferred/data_cleanfreesurf"
       / "SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6_2000m_clean_clip.stl"
   )
   DEFAULT_ALT_DIR = (
       Path(__file__).resolve().parent.parent / "data_cleanfreesurf"
   )
   DEFAULT_ALT_GLOB = "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_*_clean_clip.stl"
   ```
   `parents[2]` lands at `safs/`; verify in the implementation that the
   reference path resolves to a real file and abort with a clear message
   if not.
2. CLI:
   ```
   nw_cut_strip.py [--reference PATH] [--alt-dir DIR] [--glob PAT]
                   [--suffix SFX] [--batch | INPUT OUTPUT]
                   [--print-anchor PATH] [--keep-orphan-verts]
                   [--verbose | -q]
   ```
   - `--reference`: path to the file used to compute the NW anchor;
     default `DEFAULT_REFERENCE`.
   - `--alt-dir`: directory glob'd for inputs in `--batch` mode; default
     `DEFAULT_ALT_DIR`.
   - `--glob`: glob pattern; default `DEFAULT_ALT_GLOB`.
   - `--suffix`: appended to each input stem before `.stl`; default
     `_nwcut`. Output path = `<input.stem><suffix>.stl` in the same
     directory as the input.
   - `--batch`: run over `--alt-dir / --glob`; otherwise positional
     `INPUT OUTPUT` is required.
   - `--print-anchor`: from Phase 1, exits 0 immediately after printing.

   **Mode mutual exclusion.** The three top-level modes —
   `--print-anchor PATH`, `--batch`, and the positional `INPUT OUTPUT`
   form — are mutually exclusive. The implementer must enforce this with
   `argparse.add_mutually_exclusive_group()` (or, because positional
   pairs cannot live inside a mutex group, with explicit post-parse
   validation). Exit 2 with a message such as
   `"--print-anchor / --batch / INPUT OUTPUT are mutually exclusive"` if
   the user supplies more than one. Treat `INPUT` and `OUTPUT` as a
   single mode (both required if either is given). If none of the three
   modes is supplied, print the CLI help and exit 2.
3. Top-level pipeline (single file, no `--batch`):
   1. Load reference vertices (`load_vertices`).
   2. Compute `p_anchor, s_max = nw_anchor(verts_ref)`. Print it
      (verbose only).
   3. Compute `n, c = cutting_plane(p_anchor)`.
   4. Load the alternative file as `(verts_alt, faces_alt)` using `meshio`
      (or pymeshlab MeshSet → `current_mesh().vertex_matrix()` /
      `face_matrix()`; either is fine). Use the *vertex-matrix-and-face-
      matrix* form, not the dict form, because Phase 2 expects ndarrays.
      Immediately after loading, coerce dtypes:
      `verts_alt = np.ascontiguousarray(verts_alt, dtype=np.float64)`,
      `faces_alt = np.ascontiguousarray(faces_alt, dtype=np.int32)`.
      This is required because `meshio.read` of an ASCII or binary STL
      can return `float32` on some versions, and `clip_mesh_at_plane`'s
      contract specifies `float64` (precision loss on a ~4e6-magnitude
      UTM coordinate base would otherwise cause non-deterministic
      mis-classification of vertices within ~1 m of the cutting plane).
   5. Call `clip_mesh_at_plane(verts_alt, faces_alt, n, c)`; capture
      `(V_out, F_out, stats)`.
   6. If `len(F_out) == 0`: in the **CLI path**, print the message
      `"refusing to write empty STL: every triangle was clipped away"`
      to stderr and exit with a non-zero status. In the **programmatic
      path** (`nw_cut_file()`), raise
      `ValueError("refusing to write empty STL: every triangle was
      clipped away")`. In **both** paths, do not write the empty STL —
      the downstream gmsh `.geo` files would silently load it and emit
      an unmeshable input. (This shouldn't happen on the real data —
      the SE portion of the alternative is enormous — but a typo in
      `--reference` or a swapped reference path could trigger it.)
   7. Build a `pymeshlab.MeshSet`, `add_mesh(ml.Mesh(V_out, F_out))`, then
      run the same cleanup chain as
      `clean_freesurface_mesh.clean_freesurface_clip` does after clipping:
      `meshing_remove_duplicate_vertices`,
      `meshing_remove_null_faces`,
      `meshing_remove_unreferenced_vertices`,
      `meshing_repair_non_manifold_edges(method="Remove Faces")`. We do
      *not* call `meshing_isotropic_explicit_remeshing` — the input is
      already remeshed and we don't want to alter the bulk of the strip.
   8. Save with `ms.save_current_mesh(str(out_path), binary=False)`.
   9. Print the post-cut bbox, the number of triangles in / out, and the
      max-NW-projection of the surviving mesh, **measured after the
      pymeshlab cleanup chain in step 7 has run** (i.e. read it from
      `ms.current_mesh().vertex_matrix()`, not from the pre-cleanup
      `V_out`). Compute it with the same operand order as Phase 1 step 3:
      `s_after = V_clean[:, :2] @ NW_DIRECTION_XY`, then
      `s_max_after = float(s_after.max())`. By construction this must
      satisfy `s_max_after ≤ s_max + EPS`. The acceptance criteria in
      this Phase, and the `s_max_out` field of `nw_cut_file`'s return
      dict, refer to this same post-cleanup quantity.
4. Batch mode iterates over all matching files. The cutting plane is
   computed *once* (constant across the batch).
5. **pymeshlab is optional.** Wrap step 7 in `try: import pymeshlab as ml`
   inside the function. If the import fails, fall back to a pure-NumPy
   weld using `np.unique` along the vertex rows (round to nearest
   `1e-6` m before unique-ing to avoid float-equality misses). Do *not*
   require pymeshlab just for the cleanup; the user runs this in
   `pythonenv` which already has pymeshlab, but the unit tests in Phase 4
   should be runnable without it (CI).

### Interfaces

CLI shown in step 2. Importable Python entry point:

```python
def nw_cut_file(input_path: Path, output_path: Path,
                plane: tuple[np.ndarray, float],
                cleanup: bool = True, verbose: bool = True
                ) -> dict:
    """Apply (n, c) = plane to input_path and write output_path. Returns
    the stats dict from clip_mesh_at_plane plus 'bbox_in', 'bbox_out',
    's_max_in', 's_max_out'."""
```

### Edge Cases to Handle

- Reference file missing → exit 2 with
  `"reference file not found: <path>"`.
- Reference file is itself an ALT6 long strip (user accidentally points
  `--reference` at the alternative): the anchor is the alternative's own
  NW corner, the cut keeps everything → output is a copy of the input
  (no triangles dropped). Detect "no triangles dropped" and print a
  warning so the user notices.
- Output file already exists → overwrite without asking (matches the rest
  of the preprocess pipeline; user keeps git-committed originals).
- `--alt-dir` glob matches zero files → exit 1 with a clear message,
  matching `_batch` in `clean_freesurface_mesh.py`.
- `--alt-dir` matches a file whose stem already ends in `_nwcut` → skip
  with a warning (we don't want to re-cut an already-cut output if the
  user accidentally globs it back in).
- pymeshlab cleanup fails on a non-manifold edge (rare) → print the
  exception and continue without the repair step (matches existing
  `try/except` around `meshing_repair_non_manifold_edges` in
  `clean_freesurface_mesh.py`).

### Acceptance Criteria

- [ ] On the real data, `python nw_cut_strip.py --batch -v` produces three
      new files in `…/project_7.0_alternative/data_cleanfreesurf/`:
      - `SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m_clean_clip_nwcut.stl`
      - `SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_1000m_clean_clip_nwcut.stl`
      - `SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_500m_clean_clip_nwcut.stl`
- [ ] For each output, the post-cut max NW projection satisfies
      `s_max_after ≤ 2456425.1 + 1e-6` (i.e. nothing further NW than the
      preferred mesh's NW corner survives).
- [ ] For each output, the post-cut bbox satisfies `x_min ≥ 365072.3 −
      ε`, where ε ≤ 5 m (small `ε` allowed because the cutting plane is
      not aligned with the X axis — a vertex with `y` slightly less than
      `y_anchor` can have `x` slightly less than `x_anchor` and still
      satisfy `−x + y ≤ −x_anchor + y_anchor`).
- [ ] No vertex in the output has `(−x + y) > (−x_anchor + y_anchor) +
      1e-6`.
- [ ] Output STL is ASCII (first line `solid ...`).
- [ ] The total triangle count in the output is strictly less than the
      input (we expect ~30–40% removal on the 2000 m mesh; a mesh with no
      removal flags a configuration error).
- [ ] All three originals (`*_clean_clip.stl`) are unchanged
      (`md5sum` before vs after).
- [ ] The existing test suite in
      `project_7.0_alternative/code_preprocess/` (if `pytest` is run there
      pre-task) still passes.

### Dependencies

- Depends on: Phase 1, Phase 2.
- Required by: Phase 4 (tests use the same driver).

## Phase 4 — Tests

### Goal

After this phase, regressions in either the geometric clipping or the
pipeline glue are caught by `pytest` before any real run.

### Files to Create

- `project_7.0_alternative/code_preprocess/test_nw_cut_strip.py`

### Detailed Requirements

1. Synthetic test fixtures (no I/O against the real data):
   - `make_triangle_with_signed_dists(d0, d1, d2, n, c, *, rng=None)`:
     build a triangle whose three vertices have signed distances
     `d0, d1, d2` from the plane `n·p = c`. Implementation: sample three
     random in-plane points `q0, q1, q2` (e.g.
     `qi = rng.uniform(-10, 10, 3)`, then project onto the plane via
     `qi -= signed_dist(qi, n, c) * n`), then translate each `qi` along
     `n` by `di` to land at the requested signed distance:
     `pi = qi + di * n`. Returns `(p0, p1, p2)` as 3-tuples of floats.
     `rng` defaults to `np.random.default_rng(0)` for deterministic
     tests.
2. Unit tests for `clip_triangle_at_plane` (run all four with one
   randomly-rotated cutting plane to break any axis-aligned coincidence;
   use only "keep" / "drop" terminology — the test plane is *not*
   geographically NW vs SE, it is just a random orientation):
   - `test_all_keep`: all 3 vertices on the keep side (`d ≤ 0`) → 1
     output triangle = input.
   - `test_all_drop`: all 3 vertices on the drop side (`d > 0`) → 0
     output triangles.
   - `test_one_keep_two_drop`: 1 keep + 2 drop → 1 output triangle whose
     two new vertices satisfy
     `|signed_dist(v, n, c)| < 1e-6 * (|c| + 1)`; the surviving keep-
     side vertex must be one of the input triangle's vertices.
   - `test_two_keep_one_drop`: 2 keep + 1 drop → 2 output triangles
     whose union (as a polygon in 3-D) equals the **keep-side
     quadrilateral** (the original triangle ∩ `H_keep`). Verify by
     polygon area sum: the sum of the two output triangle areas must
     equal the area of the keep-side quadrilateral computed
     independently from the parametric clip, to within
     `1e-9 * area_input`. Do not refer to compass directions here — the
     plane is randomly oriented.
3. Equivalence test with the existing `clip_triangle_at_z0`:
   - For 100 random triangles spanning `z = 0`, assert
     `clip_triangle_at_plane(p0, p1, p2, n=(0,0,1), c=0)` returns the same
     triangles (up to vertex permutation within each tri) as
     `ts_to_stl.clip_triangle_at_z0`.
4. Plane API test: `cutting_plane(np.array([0, 0, 0]))` returns
   `n = (-1, 1, 0)/sqrt(2)`, `c = 0`. For
   `p_anchor = np.array([1, 1, 7])`, `c = 0/sqrt(2) = 0`. (Z is
   ignored by construction.)
5. Real-data smoke (skipped under
   `@pytest.mark.skipif(not REFERENCE.exists(), ...)`):
   - Compute anchor from the preferred 2000 m STL; assert
     `p_anchor` matches `(365072.3, 3838982.0, 0.0)` to 1 m and
     `s_max ≈ 2456425.1`.
   - Apply the cut to the *2000 m* alternative only (the 500 m mesh has
     ~215 k vertices and would slow the test; 2000 m has ~14 k). Assert
     post-cut `s_max ≤ s_max_pref + 1e-6` and post-cut triangle count is
     in the range `(0.5 * n_in, 0.95 * n_in)` (we expect substantial but
     not total clipping).
6. Run via `cd project_7.0_alternative/code_preprocess && pytest -q`.
   Document this in the file's docstring.

### Acceptance Criteria

- [ ] All non-real-data tests pass under `pytest -q` in either
      `mfem-dev` or `pythonenv` envs (no pymeshlab required for these
      tests).
- [ ] The real-data smoke test, when not skipped, passes within 30 s
      wall-clock on the 2000 m mesh.
- [ ] No test writes anything outside `tmp_path`.

### Dependencies

- Depends on: Phase 1, Phase 2, Phase 3.
- Required by: nothing.

## Testing Strategy

- **Phase 1 / Phase 2** are pure functions over numpy arrays → unit-test
  exhaustively with synthetic triangles (Phase 4 §1–4). No real data
  needed.
- **Phase 3** is mostly orchestration; the only unique correctness claim
  is "the post-cut max NW projection ≤ the preferred's max NW projection".
  This is checked as a `pytest` smoke (Phase 4 §5) and again as the last
  printed line of the verbose driver output (Phase 3 §3.9).
- **Visual sanity** (recommended but not required): after a successful
  batch run, open both the preferred 2000 m STL and the
  `_nwcut` 2000 m STL in MeshLab or paraview; the NW edge of the
  alternative should align with the NW corner of the preferred. The
  user already has this workflow set up for the related
  `clean_freesurface_mesh.py` outputs.
- **Mesh validity** post-cut is enforced by reusing the pymeshlab cleanup
  steps from `clean_freesurface_mesh.clean_freesurface_clip` (weld dups,
  drop nulls, drop orphans, repair non-manifold edges). We do *not*
  remesh: the input was already remeshed at the target edge length and
  re-running the remesh would shift interior vertices.

## Risk Assessment

- **Wrong NW direction sign.** The codebase already uses `(−1, +1)/√2` for
  NW because positive X is east and positive Y is north in UTM. Plan
  references this in two places and the unit test in Phase 4 §4 catches
  any flip.
- **Float precision near the cutting plane.** Coords are ~1e6 m with
  float64 ≈ 1e-9 m relative precision; on a 500 m–2000 m mesh this is
  irrelevant. The `EPS = 1e-9` threshold in Phase 2 is on the *signed
  distance* (also ~1e6 metres scale), which is loose; the unit test in
  Phase 4 §2 verifies new vertices land on the plane to 1e-6.
- **Triangles with one vertex *exactly* on the plane.** The `≤ EPS` rule
  in `clip_triangle_at_plane` treats them as keep. In rare cases this
  produces a degenerate triangle (two output vertices coincident); the
  pymeshlab `meshing_remove_null_faces` cleanup removes them. Documented
  in Phase 2 step 3 and Phase 3 step 7.
- **Anchor mis-identification on a "noisy" preferred mesh.** If the
  preferred mesh has stray orphan vertices outside the SAF trace, the
  anchor would be wrong. Phase 1 step 2 mitigates this by ignoring
  orphan vertices not used by any triangle. The current preferred
  2000 m STL is clean (its NW-most vertex is the same in raw and
  triangle-referenced sets), but this guard makes the driver robust to
  future re-runs of `clean_freesurface_mesh.py`.
- **Downstream gmsh pipeline assumes the un-cut filename.** Outputs are
  written under a *new* name (`_nwcut.stl`) so existing `.geo` files keep
  working. Switching the meshing pipeline to consume the cut version is
  out of scope for this plan and should be a separate PR.
- **Did the user actually want a vertical cut in (X, Y), or a fault-strike-
  perpendicular cut in 3-D?** The user said "left, upper" (plan view) and
  "the most NW location," both of which are horizontal-plane concepts.
  This plan therefore commits to a vertical cutting plane with horizontal
  normal `(-1, +1, 0)/√2`. If the user later wants a 3-D plane (e.g.
  perpendicular to the local fault strike at the anchor), Phases 1–3
  generalise trivially: only `cutting_plane` changes.
- **What if the preferred mesh's "NW corner" is not the natural cut you
  want?** The script accepts `--reference` so the user can swap in any
  mesh whose NW corner defines the desired cut. We do not hard-code the
  filename anywhere except in `DEFAULT_REFERENCE`.
