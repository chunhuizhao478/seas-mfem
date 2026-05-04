# Implementation Plan: San Andreas Fault System (SAFS) Realistic-Geometry Gmsh Mesh

## Overview
Build a 3D conformal tetrahedral finite-element mesh of eight SCEC CFM fault surfaces in the
greater San Andreas Fault System using Gmsh, starting from GOCAD TSurf (`.ts`) point/triangle
data shipped under `~/Documents/Earthquake Cycle Modeling of San Andreas Fault System/CFM_data/`.
The output is a tagged `.msh` file (one physical attribute per fault, plus boundary tags for
free surface and far-field) suitable for downstream consumption by the existing SEAS-MFEM
domain/fault infrastructure used by BP5.

This plan is scoped strictly to **mesh generation**: parsing CFM data → cleaned per-fault
STL → assembled Gmsh BREP/GEO model → conformal tetrahedral `.msh`. Solver integration
(boundary-condition assignment in C++, friction-law setup, fault-attribute mapping in
`fault/`, driver wiring) is explicitly **out of scope** here and will be planned separately
once the mesh exists and is validated.

## Input Data Inventory

CFM data directory contains 24 `.ts` files (8 faults × {500m, 1000m, 2000m} target edge
length). All files are GOCAD TSurf 1 format, X/Y/Z in meters, `ZPOSITIVE Elevation` (Z up,
depth negative). Triangle counts at the 2000m level:

| Fault file (CFM ID)                                             | VRTX | TRGL | Strike | Dip |
| --------------------------------------------------------------- | ---- | ---- | ------ | --- |
| SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6 (Mojave)                  |  295 |  504 |   116  |  90 |
| ETRA-PMFZ-MULT-Pinto_Mountain_fault-CFM5                        |  371 |  639 |   265  |  90 |
| SAFS-SAFZ-SBMT-Mill_Creek_fault_strand-CFM4                     |  444 |  794 |   286  |  83 |
| SAFS-SAFZ-SBMT-Mission_Creek_fault_strand-CFM4 (San Bern)       |  260 |  449 |   282  |  82 |
| SAFS-SAFZ-COAV-Mission_Creek_fault_strand-CFM4 (Coachella)      |  201 |  346 |   311  |  85 |
| SAFS-SAFZ-MULT-Banning_fault-CFM6                               |  176 |  234 |   109  |  90 |
| SAFS-SAFZ-MULT-Southern_San_Andreas_fault_and_Banning-CFM6      |  354 |  610 |   310  |  62 |
| SAFS-SAFZ-SBMT-San_Andreas_fault-CFM6 (San Bern)                |  750 | 1346 |   115  |  51 |

Combined spatial extent (UTM-like meters, NAD83 / UTM Zone 11N inferred from coordinate
ranges and Southern California map):

```
X: [365 054, 622 329]   →   ΔX ≈ 257.3 km   (along-strike NW→SE corridor)
Y: [3 692 278, 3 838 998] →  ΔY ≈ 146.7 km
Z: [   −18 014,    +2 593] →  ΔZ ≈  20.6 km   (≈18 km depth + ≈2.6 km topography)
```

The 500m files are dense enough to over-constrain the bulk mesh; the 2000m files are the
preferred default (they are still finer than the 1–2 km bulk resolution typical for regional
SEAS runs). The 500m and 1000m sets are kept for refinement studies.

## Constraints

- **Coordinate convention and origin.** Specified separately in
  [`PLAN_origin.md`](PLAN_origin.md). Summary: model frame is local-Cartesian metres
  aligned with UTM Zone 11N E/N (no rotation), with origin at the fixed UTM grid point
  `(E0, N0, Z0) = (500 000, 3 765 000, 0)`. So `v_local = v_utm − (E0, N0, 0)`. Z is not
  translated — z = 0 is sea level (CFM `ZPOSITIVE Elevation` matches the project's
  Z<0=depth convention from `miniapps/seas/CLAUDE.md`). Free-surface clamp (z>0 → 0) is
  applied in the local frame during Phase 2. A schema-validated `transform.json` records
  the choice; downstream tooling reads it to recover UTM. **Read `PLAN_origin.md` before
  implementing Phase 2.**

- **Numerical conditioning.** Raw UTM coordinates (10^5–10^6 m) cause loss of precision in
  the Gmsh CAD kernel (OpenCASCADE) tolerance handling. The fixed-grid-origin translation
  (above) keeps all coordinates < 1.4×10^5 m in magnitude — see `PLAN_origin.md` for the
  detailed conditioning argument.

- **Gmsh kernel.** Use `SetFactory("OpenCASCADE")`. The built-in geo kernel cannot handle
  triangulated freeform surfaces. OpenCASCADE supports STL import via `Merge` and surface
  classification via `ClassifySurfaces`/`CreateGeometry`.

- **Conformal fault interfaces.** The downstream SEAS-MFEM fault infrastructure assumes the
  fault is an internal interface where adjacent tetrahedra on the two sides share faces
  (the InterfaceTransformation pattern used in `fault/`). The mesher must therefore embed
  each fault surface in the bulk volume so that bulk tets on either side conform to the
  fault triangulation, and so that fault–fault intersections (e.g., Mission Creek meets
  San Andreas at the San Gorgonio Pass knot) are represented by a shared 1-D edge curve in
  the mesh, not by stairstepped/independent triangulations.

- **No fault-line clipping inside the cleanup step.** We do not modify the geometry of CFM
  surfaces (no decimation, no smoothing) beyond (a) removing duplicate vertices below
  tolerance, (b) flattening z>0 vertices to z=0, (c) deleting degenerate (zero-area)
  triangles. Anything more changes scientific content and must be approved separately.

- **Reproducibility.** Every step is driven by a script (Python or `.geo`) checked into
  `miniapps/seas/safs/`. No interactive Gmsh GUI editing in the pipeline. Random elements
  (e.g., Gmsh meshing seed) are fixed via `Mesh.RandomSeed = 1`.

- **Project conventions.**
  - All scripts live under `miniapps/seas/safs/mesh/` mirroring the BP5 layout.
  - Python tooling runs in the `pythonenv` conda env (per project `CLAUDE.md`).
  - Gmsh invocation lives in a `generate_safs_mesh.py` driver mirroring
    `bp5/mesh/generate_bp5_mesh.py`.
  - All numerical constants (resolution, buffer size, depth) come from CLI flags or
    `DefineConstant[…]` blocks — no hardcoded magic numbers (per user feedback memory
    `feedback_no_hardcoded_numbers.md`).

- **Module scope.** New top-level miniapp folder `miniapps/seas/safs/`. The `bp5/`,
  `bp2/`, `bp1/` siblings define the layout convention.

- **Single-fault smoke test (Mill Creek strand).** Specified in
  [`PLAN_smoke_millcreek.md`](PLAN_smoke_millcreek.md). Run this end-to-end on the
  Mill Creek fault strand before exercising the eight-fault pipeline; it adds an
  `--include-fault` flag to every driver and a 10-check validation suite tuned to
  catch infrastructure-level bugs that would otherwise hide behind multi-fault noise.

- **Domain box, mesh sizing, physical-tag scheme.** Specified separately in
  [`PLAN_domain.md`](PLAN_domain.md). Summary: bulk box = `floor/ceil` of fault-union
  bbox + 50 km lateral / 50 km depth buffers, snapped to a 10 km grid (default
  370 × 260 × 50 km box). Size field = single `Distance` field over all fault patches +
  linear `Threshold` ramping `res_f = 1000 m` on-fault to `res_ff = 20 km` at far-field
  over 30 km. **Physical tags exactly match `bp5_v2.geo`**: 1–4 = lateral Dirichlet,
  5 = ztop Natural (free surface), 6 = zbot Natural, 10 = bulk volume, **100 = fault
  (all 8 CFM surfaces unified under a single tag)**. Per-fault distinction is preserved
  in `output/fault_provenance.json`. **Read `PLAN_domain.md` before implementing
  Phases 3 and 5.**

## Layout

```
miniapps/seas/safs/
├── PLAN.md                              ← this document
├── README.md                            ← (created in Phase 6) usage + provenance
└── mesh/
    ├── ts_to_stl.py                     ← Phase 1+2: parse .ts, clean, write STL per fault
    ├── audit_ts_quality.py              ← Phase 1: per-file QC report (CSV)
    ├── safs.geo                         ← Phase 3+4+5: Gmsh assembly + bulk meshing
    ├── generate_safs_mesh.py            ← Phase 5: thin Python driver around gmsh CLI
    ├── validate_msh.py                  ← Phase 6: post-mesh validation
    ├── stl/                             ← (generated) per-fault cleaned STL
    ├── transform.json                   ← (generated) UTM-↔-local translation metadata
    └── output/                          ← (generated) safs_*.msh + audit CSVs
```

---

## Phase 1: Parse and Audit CFM TSurf Data

### Goal
After this phase, every `.ts` file has a row in a QC CSV (vertex count, triangle count,
bounding box, manifold-edge violations, duplicate-vertex count, zero-area triangle count,
top-edge z-extrema) so we know exactly what defects to expect before any meshing.

### Files to Create
- `mesh/audit_ts_quality.py` — Python script. Walks the CFM data directory, parses each
  `.ts`, and writes `output/cfm_audit_2000m.csv` (and similarly for 500m/1000m if `--res`
  is passed). No external dependencies beyond the Python standard library and `numpy`.

### Files to Modify
None.

### Detailed Requirements

1. **TSurf parser** — `parse_tsurf(path) → (V: ndarray[N,3], T: ndarray[M,3], header: dict)`.
   - Reads line by line. Recognized line prefixes: `VRTX`, `PVRTX`, `ATOM`, `TRGL`, plus the
     header block. `PVRTX` and `ATOM` are vertex aliases used by GOCAD; treat all three as
     producing one vertex entry. `ATOM id ref` shares coordinates of vertex `ref`.
   - Vertex IDs are 1-based in the file. Convert to 0-based internally.
   - The triangle indices in `T` reference the **file-order** vertex array; ATOMs collapse
     to the referenced vertex *position* but keep their own slot in the array (so the
     index map stays consistent with `TRGL` references).
   - Multiple `TFACE` blocks per file: concatenate; record block boundaries in `header`.

2. **QC computations** for each surface:
   - `n_vertices`, `n_triangles` (raw counts).
   - Bounding box `(xmin, ymin, zmin, xmax, ymax, zmax)` in meters.
   - **Duplicate-vertex count** — vertices within `1e−3 m` (1 mm) of another vertex.
     Implementation: snap to a `1 mm` grid (`np.round(V * 1000).astype(int64)`),
     `np.unique` along axis 0, count duplicates.
   - **Zero-area triangle count** — triangles whose two edge cross-product magnitude is
     `< 1e−6 m²`.
   - **Non-manifold edge count** — edges shared by `> 2` triangles. Build the edge
     multiset as `frozenset({a, b})` for each triangle edge; count edges with multiplicity
     `> 2`.
   - **Boundary-edge count** — edges shared by exactly 1 triangle (these become the fault
     trace at the free surface and the bottom/lateral fault edges).
   - **Top-edge clipping** — count vertices with `z > 0`. (Some CFM files have a few
     vertices a few meters above sea level due to topography sampling; these must be
     flattened in Phase 2.)
   - **Estimated mean edge length** (sanity check against the 500m/1000m/2000m label).

3. **CLI**:
   ```
   python audit_ts_quality.py \
       --cfm-dir "/Users/.../CFM_data" \
       --res 2000 \
       --out output/cfm_audit_2000m.csv
   ```

4. **Output CSV columns**: `file, n_vrtx, n_trgl, xmin, xmax, ymin, ymax, zmin, zmax,
   n_dup_vrtx, n_zero_area, n_nonmanifold_edges, n_boundary_edges, n_z_above_zero,
   mean_edge_len_m`.

### Interfaces
Public Python API (used by Phase 2):
```python
parse_tsurf(path: str | Path) -> dict  # keys: V, T, name, color
audit(V: np.ndarray, T: np.ndarray) -> dict  # all metrics above
```

### Edge Cases to Handle
- File contains `ATOM id ref` aliasing a previously declared vertex → duplicate vertex
  position is expected; do NOT count as a defect.
- File contains multiple `TFACE` blocks → concatenate triangles; note in header that the
  surface is multi-patch (relevant if any fault file is delivered as branches).
- Vertex `id` field is not contiguous (CFM does not guarantee 1..N) → build an explicit
  `id → row` lookup.
- File ends without `END` token → emit warning, do not fail.
- Coordinates with `>10` significant digits (CFM uses `%.13g`-like formatting) → use
  Python `float`, do not assume `float32`.

### Acceptance Criteria
- [ ] `python audit_ts_quality.py --res 2000` produces a CSV with one row per `.ts` file
      and exits 0.
- [ ] CSV `n_vrtx`/`n_trgl` columns match the values in the inventory table above
      (295/504 for MJVS-SAF, 750/1346 for SBMT-SAF, etc., to ±0).
- [ ] CSV `n_nonmanifold_edges` is 0 for every fault at every resolution (we expect CFM
      surfaces to be 2-manifold). If not, the offending file is flagged for manual review
      before proceeding to Phase 2.
- [ ] Combined bounding box matches `X[365054, 622329] Y[3692278, 3838998] Z[−18014, 2593]`
      (already measured during inventory) to ±1 m.

### Dependencies
- Depends on: nothing.
- Required by: Phase 2.

---

## Phase 2: Per-Fault Cleanup and STL Export

### Goal
After this phase we have one cleaned STL file per fault in `mesh/stl/`, plus a
`transform.json` recording the UTM-to-local translation, suitable for direct `Merge` into
Gmsh.

### Files to Create
- `mesh/ts_to_stl.py` — converts each `.ts` file at the chosen resolution into a binary STL.
  Produces `mesh/stl/<short_name>.stl` and writes `mesh/transform.json`.

### Files to Modify
None.

### Detailed Requirements

1. **Translation to local frame.** Specified in `PLAN_origin.md`. Use
   `safs_origin.utm_to_local(...)`; do NOT compute a centroid. Origin is the fixed UTM
   grid point `(E0, N0, Z0) = (500 000, 3 765 000, 0)`, identical for every fault and for
   every resolution. Z is not translated. `safs_origin.write_transform_json(...)` writes
   the metadata file with the schema specified in `PLAN_origin.md`.

2. **Free-surface flattening.** Any vertex with `z > 0` after translation is clamped to
   `z = 0`. Rationale: the SEAS half-space domain is bounded above by `z = 0`, and CFM
   fault traces extending into topography (a few meters at most per audit) would create
   spurious mesh elements above the free surface. This is a *modelling choice*, not a
   geometric correction; the count of clamped vertices is logged.

3. **Duplicate-vertex collapse.** Snap vertices to a `1 mm` grid and rebuild the triangle
   index list against the unique set. This removes the file-level VRTX duplicates noted
   in Phase 1 audit.

4. **Zero-area triangle removal.** Drop triangles whose area is `< 1e−6 m²` (these survive
   when CFM authors leave collapsed micro-triangles after editing). Re-emit a renumbered
   triangle list.

5. **STL output.** Write binary STL via a manual writer (Python struct), with one solid
   per file, normals computed from CCW vertex order:
   ```
   <80-byte header containing "SAFS:<fault_short_name>:CFM6">
   <uint32 n_triangles>
   <50-byte triangle records>
   ```
   Filenames use a stable short name table (declared at the top of the script):
   ```
   SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6              → safs_mjvs_saf.stl
   ETRA-PMFZ-MULT-Pinto_Mountain_fault-CFM5           → safs_pmfz_pinto.stl
   SAFS-SAFZ-SBMT-Mill_Creek_fault_strand-CFM4        → safs_sbmt_millcreek.stl
   SAFS-SAFZ-SBMT-Mission_Creek_fault_strand-CFM4     → safs_sbmt_missioncreek.stl
   SAFS-SAFZ-COAV-Mission_Creek_fault_strand-CFM4     → safs_coav_missioncreek.stl
   SAFS-SAFZ-MULT-Banning_fault-CFM6                  → safs_mult_banning.stl
   SAFS-SAFZ-MULT-Southern_San_Andreas_fault_and_Banning-CFM6 → safs_mult_ssaf_banning.stl
   SAFS-SAFZ-SBMT-San_Andreas_fault-CFM6              → safs_sbmt_saf.stl
   ```
   These short names also become the basis for Gmsh physical-group names in Phase 5.

6. **CLI**:
   ```
   python ts_to_stl.py \
       --cfm-dir "/Users/.../CFM_data" \
       --res 2000 \
       --out-dir mesh/stl
   ```

### Interfaces
- Reads from Phase 1's `parse_tsurf`.
- Writes:
  - `mesh/stl/safs_*.stl` (8 files at the chosen `--res`).
  - `mesh/transform.json`.
  - `mesh/output/cleanup_log_<res>.csv` — per-fault row: `file, n_clamped_z, n_dup_collapsed,
    n_zero_area_dropped, final_n_vrtx, final_n_trgl, final_bbox`.

### Edge Cases to Handle
- A fault triangle becomes degenerate after duplicate-vertex collapse (two of its three
  indices map to the same unique vertex) → drop the triangle.
- A fault becomes disconnected after zero-area drop (very unlikely at 2000m, but possible
  at native resolution) → emit warning and continue; Phase 5 will treat the surface as a
  single multi-patch entity anyway.
- The cleaned STL has zero triangles → hard error.

### Acceptance Criteria
- [ ] All 8 STL files exist under `mesh/stl/` after one `python ts_to_stl.py --res 2000`
      invocation.
- [ ] Each STL re-imports into Gmsh (`gmsh -open mesh/stl/safs_mjvs_saf.stl -0`) without
      "non-manifold edge" warnings.
- [ ] `transform.json` is written and round-trips: re-applying `+translation_m` to a
      known vertex recovers its CFM-UTM coordinate to ±1 m.
- [ ] `cleanup_log_2000m.csv` shows `n_clamped_z` ≤ 50 per fault (small topographic
      bleed) and `n_zero_area_dropped` = 0 for all faults at 2000m. Otherwise stop and
      investigate before Phase 3.

### Dependencies
- Depends on: Phase 1.
- Required by: Phase 3.

---

## Phase 3: Regional Domain and Free-Surface Definition

### Goal
After this phase we have a Gmsh model that contains (a) a rectangular half-space bulk
volume sized to enclose all fault surfaces with a configurable buffer, and (b) the eight
fault STLs imported and *located in space* (but not yet stitched, intersected, or
embedded). The model can be inspected in the Gmsh GUI to confirm fault positions look
geographically correct.

### Files to Create
- `mesh/safs.geo` — Gmsh script (geo + OCC). This phase establishes only the bounding box
  + STL imports; Phase 4 and 5 extend the same file.

### Files to Modify
None.

### Detailed Requirements

1. **Domain box.** Read fault bounding box from `transform.json` + the per-fault STL extents
   computed in Phase 2 (write a side-car `bbox.json` listing per-fault extents). In
   `safs.geo`, accept these as `DefineConstant` parameters with sensible defaults:
   ```
   DefineConstant[ buf_x  = {30000, Name "Buffer along ±X (m)"} ];
   DefineConstant[ buf_y  = {30000, Name "Buffer along ±Y (m)"} ];
   DefineConstant[ z_top  = {0,     Name "Free surface elevation (m)"} ];
   DefineConstant[ z_bot  = {-40000, Name "Domain bottom (m, negative = depth)"} ];
   DefineConstant[ res_f  = {2000,  Name "Target fault edge length (m)"} ];
   DefineConstant[ res_ff = {10000, Name "Far-field edge length (m)"} ];
   ```
   The combined fault bbox is `X∈[365054, 622329], Y∈[3692278, 3838998]` in UTM, which
   becomes (after centroid translation Δx=−493691, Δy=−3765638) approximately
   `X∈[−128637, +128638], Y∈[−73360, +73360]`. The box is then
   `X∈[xmin−buf_x, xmax+buf_x], Y∈[ymin−buf_y, ymax+buf_y], Z∈[z_bot, z_top]`.

   The exact `xmin/xmax/ymin/ymax` are emitted into `safs.geo` by the Phase 5 driver
   (`generate_safs_mesh.py`) as `-setnumber x_min …` etc., so the `.geo` itself stays
   independent of the specific data.

2. **Bulk volume.** `Box(1) = {x_min−buf_x, y_min−buf_y, z_bot,  Lx, Ly, z_top−z_bot};`
   where `Lx = (x_max−x_min)+2*buf_x`, `Ly = (y_max−y_min)+2*buf_y`. This produces a
   single OCC volume whose top face lies in `z = z_top = 0` (free surface).

3. **STL imports.** For each of the 8 cleaned STLs:
   ```
   Merge "stl/safs_mjvs_saf.stl";
   ```
   In OCC mode this loads them as discrete surfaces. We then call
   `ClassifySurfaces[ angle, force_param, include_boundary, curve_angle ]` followed by
   `CreateGeometry;` on the merged discrete entities so OCC can treat them as
   parametrizable surfaces for booleans (Gmsh ≥ 4.10 supports this; otherwise fall back to
   keeping them as discrete and using `MeshSize Field`-driven embedding without booleans
   — see Phase 4 risk note).

4. **No booleans yet.** This phase only verifies that all 8 surfaces import, sit in the
   correct location relative to the bulk box, and survive `ClassifySurfaces`. Visual
   inspection: open the resulting `.geo` in Gmsh GUI (`gmsh safs.geo`) and confirm the
   fault traces project sensibly onto a Southern California map.

5. **Phase-3 output target.** The `.geo` produced after Phase 3 should mesh to a
   surface-only `.msh` (`gmsh -2 safs.geo -o output/safs_phase3_check.msh`) without errors,
   even if the bulk is not yet decomposed. This is a sanity gate before continuing.

### Interfaces
None new — `safs.geo` is consumed only by `gmsh` CLI.

### Edge Cases to Handle
- An STL file is missing → Gmsh `Merge` errors; the script aborts. Add a `Printf` listing
  the expected files at the top of the `.geo`.
- `ClassifySurfaces` rejects a fault (returns 0 patches) → fall back to discrete-only
  treatment for that surface (Phase 4 will mesh it as a non-CAD embedded surface). Log the
  fault name as a warning.

### Acceptance Criteria
- [ ] `gmsh -2 safs.geo -o output/safs_phase3_check.msh` exits 0.
- [ ] Visual check (manual): in the Gmsh GUI, all 8 fault patches appear inside the bulk
      box, with the Mojave SAF segment to the NW, the Coachella Mission Creek segment to
      the SE, and Pinto Mountain crossing roughly E–W.
- [ ] Each fault patch reports `> 0` triangles in the Gmsh statistics window.

### Dependencies
- Depends on: Phase 2 (STL files), implicit on Gmsh ≥ 4.10.
- Required by: Phase 4.

---

## Phase 4: Fault–Fault Intersection Resolution

### Goal
After this phase, all eight fault surfaces are stitched into a single set of OCC faces
that share **conformal 1-D edge curves** wherever two fault surfaces touch (most notably:
Mission Creek meets San Andreas (San Bernardino) at the San Gorgonio Pass; Banning meets
Southern San Andreas; Pinto Mountain meets Mission Creek). This is the precondition for
producing a tetrahedral bulk mesh in Phase 5 where each fault remains a manifold internal
interface.

### Files to Create
None — extends `mesh/safs.geo`.

### Files to Modify
- `mesh/safs.geo` — append intersection logic after the imports of Phase 3.

### Detailed Requirements

1. **Strategy choice — primary path: OCC `BooleanFragments`.** Once `ClassifySurfaces` +
   `CreateGeometry` from Phase 3 has promoted each STL to an OCC surface, run
   ```
   BooleanFragments{ Volume{1}; Delete; }
                   { Surface{f1, f2, …, f8}; Delete; };
   ```
   This computes all pairwise intersections (volume↔surface and surface↔surface) and
   produces a non-overlapping cell complex: one or more bulk volumes plus a set of internal
   faces along each fault, with shared edges where faults cross. This is the same
   technique the BP5 `bp5_v2.geo` uses to create internal fault faces, scaled up to
   freeform inputs.

2. **Intersection tolerance.** Set `Geometry.OCCFixDegenerated = 1;
   Geometry.OCCFixSmallEdges = 1; Geometry.OCCFixSmallFaces = 1;
   Geometry.Tolerance = 1.0;` (1 m tolerance, well below the 500m–2000m mesh size). Larger
   values cause CFM surfaces that are physically distinct (parallel strands ≈100 m apart)
   to be glued; smaller values miss legitimate touches.

3. **Strategy fallback — discrete embedding.** If `ClassifySurfaces` fails for one or more
   faults (Phase 3 logged it), fall back to using `Mesh.Algorithm3D = 10` (HXT) with
   discrete-surface embedding via `Mesh.MeshSizeFromCurvature = 0` and constraint surfaces
   declared with `Surface{…} In Volume{1};`. This skips OCC booleans for the affected
   surface and lets HXT honor the discrete triangulation as a *Steiner constraint*. The
   trade-off: HXT will not enforce fault-fault edge conformity if both intersecting faults
   are discrete-only — this case must be flagged as known limitation.

4. **Identification of post-fragment fault tags.** After `BooleanFragments`, the surface
   tags change. We need to recover, for each fault, the new surface tags. Use coordinate
   queries:
   ```
   safs_mjvs_saf_surfaces[] = Surface In BoundingBox{ <bbox + 100 m pad of original STL> };
   ```
   Cross-reference with each fault's pre-fragment bounding box (loaded from `bbox.json`).
   This is fragile when faults overlap in 2-D projection but distinct in Z (San Andreas SBMT
   dips at 51° while the Banning fault is vertical, so vertical separation is significant).
   Where bounding boxes overlap, also filter by surface-normal direction (compare the average
   normal of the post-fragment patch against the pre-fragment STL's average normal; require
   `dot > 0.95`).

5. **Diagnostic output.** Emit a `safs_intersection_report.txt` with:
   - Number of OCC volumes after fragments (expected: 1 or 2, depending on whether faults
     fully partition the bulk — they don't, they are slivers).
   - For each fault: number of surface patches found, their summed area, and the
     pre-fragment STL area (must agree to ±5%).
   - List of intersection curves and the two faults each curve belongs to.

### Interfaces
None new — extension of `safs.geo`.

### Edge Cases to Handle
- **Two faults nearly parallel and within OCC tolerance** (e.g., Mill Creek and Mission
  Creek strands ≈500 m apart in the San Bernardino mountains). With `Geometry.Tolerance = 1`
  this is fine, but if a user passes a coarser `--res` we must not silently raise the
  tolerance. Hard-code the OCC tolerance independent of mesh resolution.
- **Fault dips below domain bottom** (`z_bot = −40 km`). All CFM fault z-extents are
  shallower than 18 km, so `z_bot = −40 km` is safely below; if the user changes
  `z_bot > -20 km`, a fault may extend below it and `BooleanFragments` will leave a
  dangling fault face below the box. Detect by checking whether the post-fragment fault
  has any vertex with `z < z_bot + 100 m` and erroring out.
- **Fault trace exits the buffered domain laterally.** Same risk; `buf_x` and `buf_y`
  default to 30 km which keeps all CFM faults strictly interior, but enforce by the same
  bbox check.

### Acceptance Criteria
- [ ] `gmsh -2 safs.geo -o output/safs_phase4_check.msh` exits 0.
- [ ] `safs_intersection_report.txt` shows summed-area-ratio in `[0.95, 1.05]` for every
      fault (no surface lost or duplicated to fragmentation).
- [ ] Manual GUI inspection: the San Gorgonio knot region shows a single shared 1-D edge
      where Mission Creek strand and San Andreas SBMT cross, NOT two separated edges.

### Dependencies
- Depends on: Phase 3.
- Required by: Phase 5.

---

## Phase 5: Volumetric Tetrahedral Meshing

### Goal
After this phase we have a 3-D conformal tetrahedral mesh (`output/safs_<res>m.msh`) where
each fault is an internal face set with its own Gmsh physical group, the free surface and
far-field box faces have boundary physical groups, and mesh size is locally refined around
each fault.

### Files to Create
- `mesh/generate_safs_mesh.py` — driver that emits per-fault bbox numbers into `safs.geo`
  and invokes `gmsh -3`. Mirrors `bp5/mesh/generate_bp5_mesh.py`.

### Files to Modify
- `mesh/safs.geo` — add Mesh size field, physical groups, and 3-D meshing block.

### Detailed Requirements

1. **Mesh size field — distance-driven refinement around faults.**
   ```
   Field[1] = Distance;
   Field[1].SurfacesList = { <all post-fragment fault surface tags> };
   Field[1].Sampling = 100;

   Field[2] = Threshold;
   Field[2].InField = 1;
   Field[2].SizeMin = res_f;        // edge length on the fault
   Field[2].SizeMax = res_ff;       // edge length far away
   Field[2].DistMin = 0;
   Field[2].DistMax = 4 * res_ff;   // ramp transition over 4 far-field cells

   Background Field = 2;
   Mesh.MeshSizeExtendFromBoundary = 0;
   Mesh.MeshSizeFromPoints = 0;
   Mesh.MeshSizeFromCurvature = 0;
   ```

2. **3-D meshing algorithm.** `Mesh.Algorithm3D = 10;` (HXT, parallel-friendly,
   robust on freeform constraint surfaces). Set `Mesh.OptimizeNetgen = 1;` to clean up
   slivers post-mesh.

3. **Physical groups.**
   ```
   // Bulk
   Physical Volume("bulk", 1) = { <all post-fragment volumes> };

   // Faults — one tag per CFM fault. Tag values 100-series for clarity.
   Physical Surface("safs_mjvs_saf",         101) = { safs_mjvs_saf_surfaces[] };
   Physical Surface("safs_pmfz_pinto",       102) = { safs_pmfz_pinto_surfaces[] };
   Physical Surface("safs_sbmt_millcreek",   103) = { safs_sbmt_millcreek_surfaces[] };
   Physical Surface("safs_sbmt_missioncreek",104) = { safs_sbmt_missioncreek_surfaces[] };
   Physical Surface("safs_coav_missioncreek",105) = { safs_coav_missioncreek_surfaces[] };
   Physical Surface("safs_mult_banning",     106) = { safs_mult_banning_surfaces[] };
   Physical Surface("safs_mult_ssaf_banning",107) = { safs_mult_ssaf_banning_surfaces[] };
   Physical Surface("safs_sbmt_saf",         108) = { safs_sbmt_saf_surfaces[] };

   // External boundaries
   Physical Surface("free_surface", 200) = { <top face of bulk: z = 0> };
   Physical Surface("bottom",       201) = { <bottom face: z = z_bot> };
   Physical Surface("xmin",         202) = { <face: x = x_min - buf_x> };
   Physical Surface("xmax",         203) = { <face: x = x_max + buf_x> };
   Physical Surface("ymin",         204) = { <face: y = y_min - buf_y> };
   Physical Surface("ymax",         205) = { <face: y = y_max + buf_y> };
   ```
   Tag identification for the box faces uses `Surface In BoundingBox{ … }` queries with a
   1 m pad.

4. **Driver script `generate_safs_mesh.py`.** CLI mirrors BP5:
   ```
   python generate_safs_mesh.py \
       --res 2000 \
       --buf-xy 30000 \
       --z-bot -40000 \
       -o output/safs_2000m.msh
   ```
   Internally:
   - Reads `mesh/transform.json` and `mesh/bbox.json` (per-fault, written in Phase 2).
   - Computes domain bbox.
   - Invokes `gmsh -3 mesh/safs.geo -setnumber res_f <res> -setnumber buf_x <buf-xy> …`.

5. **Mesh format.** `.msh` v2.2 ASCII (the format consumed by MFEM). Add
   `Mesh.MshFileVersion = 2.2;` to `safs.geo`. Also emit a `.vtk` for ParaView inspection:
   the driver runs Gmsh a second time with `-format vtk` (or uses `Save "...";` from inside
   `safs.geo`).

6. **Resolution tiers.**
   - `res_f = 2000 m, res_ff = 10 km`: development tier, ≈O(few hundred k tets), runs in
     minutes on a laptop.
   - `res_f = 1000 m, res_ff = 5 km`: science tier.
   - `res_f =  500 m, res_ff = 2.5 km`: production tier (use 500m STLs as input).

   The `--res` flag selects the tier; the driver picks the matching STL set automatically:
   `2000` → uses 2000m STLs, etc. (Phase 2 must be re-run for the chosen STL set.)

### Interfaces
- `generate_safs_mesh.py` invocation contract documented in `README.md` (Phase 6).

### Edge Cases to Handle
- HXT runs out of memory at 500m tier → log peak memory; suggest re-running on larger node.
- A fault becomes very thin sliver after fragmentation (edge < res_f / 10) → HXT may
  refuse. Mitigation: bump `Mesh.AlgorithmSwitchOnFailure` and re-run with `Algorithm3D = 1`
  (Delaunay) as fallback, behind a `--fallback-mesh-algo` flag. Do not silently switch.
- Internal-face orientation inconsistency → MFEM's mesh reader is tolerant, but the SEAS
  fault interface code expects oriented normals. The fault attribute scheme above gives
  one tag per fault; in C++ post-processing, a per-fault canonical normal is computed
  from the average of triangle normals (out of scope for this plan).

### Acceptance Criteria
- [ ] `python generate_safs_mesh.py --res 2000` produces `output/safs_2000m.msh` and exits 0.
- [ ] `gmsh -info output/safs_2000m.msh` reports `nb tetrahedra > 0`, `nb triangles > 0`,
      and lists all 14 physical groups (1 volume + 8 fault + 5 box-face + 1 top = 14;
      adjust if `bottom`/`xmin`/`xmax`/`ymin`/`ymax` collapse to fewer faces).
- [ ] Per-fault triangle count in the `.msh` matches the post-fragment count from Phase 4
      to ±2% (small tolerance for HXT's local re-triangulation of constraint surfaces).
- [ ] Free surface (`Physical Surface 200`) is the entire top face of the bulk and contains
      the surface trace of every fault that reaches `z = 0` (visual check via
      `validate_msh.py` plotting traces).
- [ ] Mesh quality: `gmsh -info` reports `min Gamma > 0.1` (no degenerate tets).

### Dependencies
- Depends on: Phase 4.
- Required by: Phase 6.

---

## Phase 6: Mesh Validation and Provenance

### Goal
After this phase we have a `validate_msh.py` script that confirms the mesh is suitable for
SEAS-MFEM consumption, plus a `README.md` documenting how to regenerate the mesh and what
each physical tag means.

### Files to Create
- `mesh/validate_msh.py` — Python script using `meshio` (already a project dependency in
  the `pythonenv` env, per BP5 scripts) to inspect the `.msh`.
- `safs/README.md` — usage + provenance documentation.

### Files to Modify
None.

### Detailed Requirements

1. **`validate_msh.py` checks.**
   - All expected physical groups present (by name and tag).
   - For each fault tag: triangulation is connected (no orphan patches), no duplicate
     triangles, no triangles with all 3 vertices on the free surface (those would be
     degenerate fault traces, not fault interior).
   - For each fault triangle: both adjacent tetrahedra exist (i.e., it really is internal).
     This is the critical SEAS assumption: if a fault triangle has only one adjacent tet,
     the fault has been turned into a boundary by mistake.
   - Free-surface triangle count matches `nb triangles in Physical Surface 200`.
   - Print summary: per-fault triangle count, per-fault area, per-fault depth extent,
     per-fault dip-angle histogram (computed from triangle normals; should match the CFM
     average dip values from the inventory table to ±5°).

2. **`README.md` content.**
   - One-paragraph overview.
   - The fault inventory table (copy from this PLAN.md).
   - Coordinate convention note + pointer to `transform.json`.
   - Regeneration instructions (3 commands: `audit_ts_quality.py`, `ts_to_stl.py`,
     `generate_safs_mesh.py`).
   - Physical-tag table mapping name ↔ ID ↔ CFM source file ↔ description.
   - Known limitations (any discrete-only fallbacks from Phase 4; any
     fault near the buffer edge; any HXT slivers post-optimization).

3. **Reproducibility seal.** `validate_msh.py` writes
   `output/safs_<res>m.provenance.json` containing:
   - Git commit hash of `seas-mfem` repo.
   - SHA-256 of every input `.ts` file used.
   - Gmsh version (`gmsh -version`).
   - Driver-script CLI flags as parsed.
   - Output mesh SHA-256.

### Interfaces
None new.

### Edge Cases to Handle
- `meshio` cannot be imported → emit a clear error pointing the user to
  `conda activate pythonenv && pip install meshio`.
- `.msh` file lacks a fault physical group (Phase 5 mis-tagged) → `validate_msh.py` reports
  it and exits non-zero, blocking commit.

### Acceptance Criteria
- [ ] `python validate_msh.py output/safs_2000m.msh` exits 0 for a freshly generated mesh.
- [ ] `output/safs_2000m.provenance.json` exists, contains all 8 input-file SHAs, and
      reproduces (i.e., the same inputs + same git commit + same Gmsh version reproduce
      the same output mesh SHA bit-for-bit modulo HXT non-determinism — accept ±0.5% tet
      count drift).
- [ ] `safs/README.md` references this PLAN.md and includes a one-paragraph "what is in
      scope vs out of scope" section so the next agent knows the solver wiring is still
      pending.

### Dependencies
- Depends on: Phase 5.
- Required by: nothing (terminal phase).

---

## Testing Strategy

- **Phase 1**: golden-file test — run the audit on the 2000m set, check that the produced
  CSV equals a checked-in `tests/golden_audit_2000m.csv` (regenerate when `audit_ts_quality.py`
  changes intentionally).
- **Phase 2**: round-trip test — for each cleaned STL, parse it back with `meshio`,
  re-translate to UTM via `transform.json`, and check that the recovered bounding box
  matches the original `.ts` bbox to ±2 m.
- **Phase 3**: smoke test — `gmsh -2 safs.geo` succeeds; the resulting surface mesh has
  at least N triangles where N = sum of input STL triangles.
- **Phase 4**: area-conservation test as in Phase 4 acceptance.
- **Phase 5**: connectivity test — every fault triangle has 2 adjacent tets. This is the
  same property that BP5 relies on; the project's existing `face_audit_*.csv` machinery
  (visible in the repo root) already exercises this assumption against BP5 meshes and can
  be pointed at the SAFS mesh as an integration check.
- **Phase 6**: validation script return code as the canonical "is the mesh good" signal.

There are no analytical/benchmark comparisons in this plan because the SAFS mesh is purely
a geometric artifact — scientific validation comes later when a SEAS run on it is compared
against, e.g., GPS-inferred slip rates.

---

## Risk Assessment

1. **OCC `BooleanFragments` on freeform discrete surfaces is the central technical risk.**
   `ClassifySurfaces` succeeds easily on simple shapes but can fail or produce many tiny
   patches on noisy CFM triangulations. Mitigation: Phase 3 already runs `ClassifySurfaces`
   in isolation as a gate; Phase 4 has an explicit fallback. If the fallback is needed for
   more than two faults, the discrete-only path may not produce a mesh with conformal
   fault-fault intersections — at that point we stop and consider an alternative pipeline
   (e.g., TetGen with PLC input, or pre-conforming the surfaces in MeshLab/PyMesh before
   Gmsh ever sees them). That alternative is **out of scope** for this plan and will be
   planned separately if triggered.

2. **CFM coordinate-system assumption.** This plan assumes UTM Zone 11N / NAD83 based on
   the X/Y ranges. If CFM's published metadata says otherwise (e.g., NAD27, or a different
   projection), the `transform.json` `utm_zone`/`datum` fields will be wrong — but this
   does not affect mesh geometry (it is local-frame) and only affects future tooling that
   wants to re-georeference the mesh. We will verify against the SCEC CFM documentation
   page (linked from the original download portal) before merging Phase 2.

3. **Manifold violations from CFM source.** Phase 1 audits non-manifold edges; if any
   appear, Phase 2 cleanup (currently only duplicate-vertex collapse and zero-area drop)
   will not fix them. We do not blindly add manifold repair — that would alter scientific
   content. Instead, the affected fault is flagged and the user is asked to choose between
   (a) accepting the higher-resolution CFM file (which sometimes lacks the defect),
   (b) excluding the fault from the model, or (c) accepting a manifold repair (e.g.,
   `pymeshlab` `meshing_repair_non_manifold_edges`). No automatic answer.

4. **Mesh size at the SAF/Mission Creek/Banning knot.** The San Gorgonio Pass region has
   three faults converging within ~5 km. With `res_f = 2000 m` this is fine; at `res_f =
   500 m` the local tet count near the knot may explode. Phase 5's `Threshold` field uses
   a smooth ramp; if memory becomes a problem, switch the field to a per-fault `Threshold`
   so different faults get different `res_f` values (this is a Phase 5+ extension, not
   in this plan).

5. **Pre-existing project conventions for fault tagging.** BP5 uses a single fault tag
   (= 3 in `bp5.geo`). The SAFS mesh introduces **eight** fault tags. The existing
   `fault/` infrastructure may assume a single fault attribute in places. Reviewing and
   extending that infrastructure is **out of scope** for this plan but is the most
   important downstream dependency. The plan's Phase 6 README explicitly flags this so
   the solver-side plan is written with eyes open.

6. **The user's note "Always follow SCEC benchmark specifications exactly" (memory).**
   This is a SAFS realistic-geometry build, not a benchmark, so there is no SCEC SEAS
   spec to follow exactly. The CFM data themselves are the spec. Cite CFM version (CFM4,
   CFM5, CFM6 per fault — see inventory) in the README so any future "did we use the
   right data" question has a clear answer.
