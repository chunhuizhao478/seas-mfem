# Implementation Plan: On-Fault Stress Projection from Regional Bulk Stress

## Overview

Build an end-to-end pipeline that takes the regional bulk Cauchy stress
tensor σ⁰ from
`data_stress_frictionlaw/regional_stress_projection/hickman_and_zoback_regional_stress_projection.py`,
delivers it as a per-DOF pre-stress field into the SEAS quasi-dynamic and
dynamic-rupture solvers (`tau0`, `sigma_n`, `tau1_0`/`tau2_0`/`sigma_n0`),
and produces ParaView artefacts for human verification along the way.
The pipeline mirrors the **velocity-model pipeline** end-to-end:

```
                preprocess (Python)                     runtime (C++)
   ┌──────────────────────────────────────┐   ┌────────────────────────────┐
   σ⁰ (H&Z) ───► project_to_fault_stress  ───► stress_safs.h5 ──► StressField3D
                ├── VTU (debug/QA)            (schema-v1 sidecar)   │
                └── JSON summary                                    ▼
                                                          FieldProjector
                                                          ├── ::ProjectStress (bulk)
                                                          └── ::ProjectFaultPreStress
                                                                    │
                                                                    ▼
                                              FaultGeometry::ComputeSAFSParams
                                                  tau_pre_  (per-DOF tau1, tau2)
                                                  sigma_n_  (per-DOF |σ_n|)
                                                                    │
                                              ┌─────────────────────┴───────┐
                                              ▼                             ▼
                                      QD: RateStateFaultOperator    Dyn: DOFData
                                          .tau0 = tau_pre_              .tau1_0, tau2_0
                                                                        .sigma_n0
```

The plan now has **nine phases**:

| Phase | What is delivered |
|-------|---|
| 0 | H&Z dump entry point (`--dump <path>`); the only edit to the H&Z file |
| 1 | Mesh I/O + per-cell geometry (centroids, normals, areas, volumes) |
| 2 | Per-triangle fault basis with global orientation harmonisation; **cell→node averaging helpers** for the continuous fault field |
| 3 | Stress projection: constant σ⁰ rotation onto fault frame + depth-dependent option; **single sign-flip site** (geomechanics → SEAS, compression positive) |
| 4 | ParaView VTU + JSON summary writers — **continuous fault field as point-data** (cell-data kept as `_cell` channel for traceability) |
| 5 | **schema-v1 HDF5 sidecar writer** (`stress_safs.h5`) — canonical machine-readable artefact, mirrors `velocity_safs.h5`; sign convention compression-positive Pa |
| 6 | **C++ consumer**: `StressField3D` loader + `FieldProjector::ProjectStress` / `::ProjectFaultPreStress` + `FaultGeometry::ComputeSAFSParams` seam — pure pass-through of the sign convention |
| 7 | **Verification driver** `seas_project_stress_to_mesh` — analog of `seas_project_velocity_to_mesh` |
| 8 | **VTU postprocess verifier** `verify_onfault_stress.py` — reads every written VTU back and checks against the analytic prediction within an error tolerance |

Phases 0–5 and 8 are **Python**, run in the `pythonenv` conda
environment. Phase 6 is **C++**, built in `mfem-dev`. Phase 7 is a
thin C++ driver that exercises Phase 6 end-to-end.

---

## Constraints

### Interface constraints (must not change)
- `hickman_and_zoback_regional_stress_projection.py` **may be modified
  in two narrow ways only**:
  1. Add an entry point that dumps the resolved σ⁰ tensor and the
     scalar input parameters (`SHmax`, `Shmin`, `Sv`, `P_p`,
     `SHmax_az_deg`, and the regime label) to a sidecar file (default
     `hickman_zoback_sigma0.json`) so the new module can read those
     values without re-running the demo or importing private state.
  2. Add a CLI flag `--dump <path>` that invokes the new entry point.

  No other changes are permitted; the function signatures
  `build_bulk_stress_tensor`, `fault_basis_vectors`, `resolve_traction`,
  and `compute_fault_stress` keep their current names, argument lists,
  and numerical semantics. The new module imports them as a library
  in addition to reading the dumped sidecar.
- `msh_to_vtu.py` is **read-only**; we consume its `*_fault.vtu` and
  `*_bulk.vtu` outputs as-is. They share a common point cloud, both
  contain `gmsh:physical` cell-data, and the fault VTU is pure-triangle
  while the bulk VTU is pure-tet (verified on
  `code_meshing/safs_fault_box_nwcut_2000m_*.vtu`).
- Output directory **must** be
  `/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_alternative/code_preprocess/data_projection_onfaultstress/`
  (user-explicit override of the original plan, which had the
  output directory at the project-root level: outputs are
  co-located with the plan document under `code_preprocess/` per
  R-701).
  Subfolder per mesh resolution (`500m/`, `1000m/`, `2000m/`).

### Dependency constraints
- **Phases 1–5 (Python)** run in the `pythonenv` conda environment
  (same as `nw_cut_strip.py`, `clean_freesurface_mesh.py`,
  `msh_to_vtu.py`, `data_projection/sidecar.py`).
  Allowed packages: `numpy`, `meshio`, `h5py`, `pytest`. **Do not**
  pull in `pymeshlab`, `pyvista`, or any C++ dependency — meshio writes
  legacy unstructured VTU which ParaView reads natively.
- **Phase 5 sidecar writer** must reuse the existing schema-v1 writer
  `code_preprocess/data_projection/sidecar.py:write_sidecar`. Do not
  fork a new HDF5 writer — schema-v1 already declares the
  `stress` field set (`sigma_xx`, `sigma_yy`, `sigma_zz`, `sigma_xy`,
  `sigma_yz`, `sigma_xz`; units Pa each — see
  `miniapps/seas/document/features_dev/data_projection_schema_v1.md`
  §5 table).
- **Phases 6–7 (C++)** run in `mfem-dev`. Allowed dependencies are
  the dependencies already used by `io/data_field_3d.hpp` and
  `io/field_coefficient.hpp` (MFEM, HDF5, MPI). Do not introduce new
  third-party libraries.
- No new top-level data products outside
  `data_projection_onfaultstress/`. Per-mesh subfolders only.

### Convention constraints
- File naming follows existing pipeline:
  `safs_fault_box_nwcut_<lc>m_fault_stress.vtu`, `..._bulk_stress.vtu`,
  `..._summary.json`, `..._stress.h5`. Stem is the input-mesh stem
  with the suffix `_fault_stress` / `_bulk_stress` / `_summary` /
  `_stress`. The HDF5 sidecar shares the naming pattern of
  `data_projected/velocity_safs.h5` — i.e. one canonical SAFS-scale
  file (`stress_safs.h5`) that drops in next to the velocity sidecar.
- CLI layout follows `nw_cut_strip.py`: `argparse`, three mutually
  exclusive modes (single, `--batch`, `--print-info`), `python -m`
  invocable, module-level constants for defaults.
- C++ class naming follows `DataField3D` / `FieldCoefficient` /
  `FieldProjector`: capitalised camel-case, header-only where
  possible, namespace `mfem::seas`.
- **Unified sign convention on emitted artefacts: compression
  POSITIVE everywhere.** Every output written to disk — fault VTU
  point-data and cell-data, bulk VTU cell-data, JSON summary,
  schema-v1 HDF5 sidecar — and every C++ consumer
  (`StressField3D`, `FieldProjector::ProjectStress`,
  `FieldProjector::ProjectFaultPreStress`,
  `FaultGeometry::tau_pre_`, `sigma_n_per_dof_`,
  `DOFData.sigma_n0` / `tau1_0` / `tau2_0`) carries the SEAS
  internal "compression positive, σ_n > 0 = compression" convention,
  exactly as documented in `miniapps/seas/CLAUDE.md` ("Sign
  Conventions" → "Normal stress: sigma_n > 0 = compression
  (geology convention)").
- **The sign flip is a SINGLE full-tensor flip at the source.**
  `bulk_stress_tensor_field` (Phase 3 §2) multiplies H&Z's output
  by −1 once at the source, returning σ_seas = −σ_HZ in
  compression-positive SEAS convention. Every downstream consumer
  (resolver, fault VTU, bulk VTU, JSON summary, schema-v1 sidecar,
  C++ `StressField3D`, `ProjectStress`, `ProjectFaultPreStress`,
  `FaultGeometry::tau_pre_`, `sigma_n_per_dof_`, DOFData slots,
  Phase 8 verifier) consumes compression-positive Pa/MPa as a
  pass-through. **No further sign manipulation anywhere.**

  The R-501/R-502 derivation: under Phase 2's locked Tandem basis
  (`s = up × n`, `d = s × n`, `n` harmonised to the H&Z half-space),
  the basis vectors satisfy `s_T = −s_HZ` and `d_T = −d_HZ`. The
  source σ-flip combined with the basis-flip composes as:

  - σ_n = nᵀ σ_seas n = nᵀ (−σ_HZ) n = −σ_n_HZ
    → compression positive (one flip, propagates).
  - τ_strike = s_Tᵀ σ_seas n = (−s_HZ)ᵀ (−σ_HZ) n = +τ_strike_HZ
    → unchanged from H&Z (two flips cancel symmetrically).
  - τ_dip = d_Tᵀ σ_seas n = (−d_HZ)ᵀ (−σ_HZ) n = +τ_dip_HZ
    → unchanged from H&Z (two flips cancel symmetrically).

  Result: SEAS-internal sign convention holds simultaneously for
  σ_n (compression positive), τ_strike (right-lateral positive),
  and τ_dip (reverse positive). No special-case logic anywhere.

  (Historical note: round-1 R-001 originally argued for an
  asymmetric "split-site flip" because the plan then used H&Z's
  basis where the basis-flip does NOT happen. Round-2 R-102
  switched to Tandem basis but R-001's split-site reasoning was
  not re-evaluated. Round-3 R-501/R-502 restored the simpler
  single-site flip now that the basis-flip provides the
  cancellation.)
- The continuum-mechanics convention is preserved **only inside** the
  H&Z module (its docstring states the convention explicitly) and
  inside the Phase 0 dump file (whose `convention` string is
  `"compression negative; continuum-mechanics; MPa; frame (east,
  north, up)"`). The JSON summary records
  `"convention": "compression POSITIVE (SEAS internal); ..."`
  so a reader cannot mistake the active convention on the emitted
  artefact.
- A separate parity test against H&Z
  (`test_resolve_traction_against_HZ_safod`) confirms the
  sign-flip-and-basis-flip cancellation for SAFOD:
    σ_n_seas      = −σ_n_HZ      = +104.78 MPa
    τ_strike_seas = +τ_strike_HZ = +21.41  MPa  (right-lateral +)
    τ_dip_seas    = +τ_dip_HZ    ≈ 0       MPa  (pure strike-slip)
  Anchors the "positive τ_strike ↔ right-lateral on the fault"
  relationship to the CLAUDE.md "Slip rate direction" rule.
- Coordinate frame **must** be the H&Z convention `(x = east, y = north,
  z = up)`. The SAFS mesh is in UTM Zone 11 metres
  (verified bbox: X ∈ [303 000, 697 500] m, Y ∈ [3.612e6, 3.903e6] m,
  Z ∈ [−66 607, 100] m). The local UTM frame is already
  (east, north, up) — no axis flip needed.
- Pore pressure `P_p` is reported positive in every convention;
  under SEAS (compression-positive), `σ_n_eff = σ_n_total − P_p`
  (P_p reduces the effective normal load). This sign sense
  matches Terzaghi effective stress and the friction
  `μ · σ_n_eff` convention.
- Units: input MPa (matching H&Z); written-out fields stay in MPa
  (annotated in the JSON summary and as a `units` cell-data attribute
  if meshio supports it; otherwise documented in the summary). The
  downstream C++ consumer converts to Pa as needed; see CSM
  equilibration plan §1.
- **No hardcoded numerics.** Every physical constant (SHmax, Shmin, Sv,
  P_p, SH_max azimuth, default dip, rake sense, lithostatic ρ̄ if used,
  g if used) must come from CLI flags with defaults pulled from the
  H&Z `demo_safod` block — and named so a reader of the code does not
  need to open the H&Z file to understand the meaning.

### Numerical constraints
- Per-triangle rotation must be numerically stable for **near-vertical**
  faults (the SAF is approximately vertical). The Tandem
  `facetBasis(up, n)` algorithm (`tandem/src/geometry/Curvilinear.cpp:260`)
  throws when `up ∥ n`; we cannot afford that. For a triangle whose
  normal has `|n_z| > 1 − tol` (effectively horizontal triangle, dip ≈ 0),
  the strike direction is degenerate — emit a NaN-marked row and warn
  in the summary; do not crash.
- Triangle normal orientation is ambiguous (CCW vs CW). The H&Z basis
  function chooses the normal to **point to the left of strike in map
  view** for a right-lateral fault. We must enforce this globally:
  after building the raw normal from the triangle vertex order, flip
  it if `n_raw · n_consensus < 0`, where `n_consensus` is the mean
  raw normal over the trace-projected mid-strike triangles (see
  Phase 2 §Normal-orientation harmonisation).
- All depth coordinates Z are negative in this mesh frame. When
  depth-dependent σ_v is requested, use `depth = max(0, -z_centroid)`
  to clamp the (rare) above-free-surface triangles.

---

## Phase 0 — H&Z dump entry point

### Goal
After Phase 0 the H&Z script can be invoked as
`python hickman_and_zoback_regional_stress_projection.py --dump <path>`
to write a single JSON file containing σ⁰ (3×3, MPa, geomechanics
convention) and the scalar input parameters used to build it.
Downstream phases read this JSON and never need to re-import private
state from the H&Z module beyond `build_bulk_stress_tensor`,
`fault_basis_vectors`, and `resolve_traction`.

### Files to Create
- None.

### Files to Modify
- `data_stress_frictionlaw/regional_stress_projection/hickman_and_zoback_regional_stress_projection.py`
  — add two functions and a CLI flag (described below).

### Detailed Requirements

1. Add the function

   ```python
   def dump_safod_sigma0(
       out_path: Path,
       *,
       SHmax: float = 113.0,
       Shmin: float = 49.0,
       Sv: float = 45.0,
       P_p: float = 16.0,
       SHmax_azimuth_deg: float = 23.0,
       fault_strike_azimuth_deg: float = 314.0,
       fault_dip_deg: float = 90.0,
       rake_sense: str = "right-lateral",
   ) -> None
   ```

   that calls `compute_fault_stress(...)` with the supplied arguments
   and writes a JSON file with this exact schema:

   ```json
   {
     "schema": "hickman_zoback_sigma0_v1",
     "params": {
       "SHmax_MPa": 113.0, "Shmin_MPa": 49.0, "Sv_MPa": 45.0,
       "P_p_MPa": 16.0,
       "SHmax_azimuth_deg": 23.0,
       "fault_strike_azimuth_deg": 314.0,
       "fault_dip_deg": 90.0,
       "rake_sense": "right-lateral"
     },
     "sigma0_MPa": [[...3x3 row-major...]],
     "convention": "compression negative; continuum-mechanics; MPa; frame (east, north, up)",
     "regime": "<value of compute_fault_stress.last_regime>"
   }
   ```

   The defaults match `demo_safod()` exactly (lines 378–387 of the H&Z
   file) so that the call with no kwargs reproduces the SAFOD reference
   case used by all downstream tests.

2. Add a `--dump <path>` CLI argument that calls
   `dump_safod_sigma0(path)`. If the user combines `--dump` with other
   flags from a future Phase, those flags pass through to
   `dump_safod_sigma0` as kwargs; with this Phase 0 only `--dump
   <path>` is supported.
   The recommended canonical path is
   `data_stress_frictionlaw/regional_stress_projection/hickman_zoback_sigma0.json`,
   declared as `DEFAULT_DUMP_PATH = Path(__file__).resolve().parent / "hickman_zoback_sigma0.json"`
   at module top so Phase 3 (`--hz-dump-file`) and Phase 5 can
   import it for their parameter-echo cross-checks without
   per-test configuration.

3. **No other behaviour changes.** Running the script with no
   arguments must still execute `demo_safod()` and `sweep_strike()` as
   it does today.

### Acceptance Criteria
- [ ] `python hickman_and_zoback_regional_stress_projection.py --dump
      /tmp/safod.json` writes a file whose
      `json.load(open(...))["sigma0_MPa"]` matches the printed σ⁰
      from `demo_safod()` to 1e-9.
- [ ] `json.load(open(...))["params"]["SHmax_azimuth_deg"]` equals 23.0.
- [ ] Running the script with no arguments still prints the demo
      output (no regression on the existing behaviour).

### Dependencies
- Depends on: nothing.
- Required by: Phase 3 (reads the dumped JSON to populate CLI
  defaults when `--hz-dump-file <path>` is passed), Phase 5 (sidecar
  writer uses the same JSON as its parameter echo), and the
  Phase 0 / 8 cross-validation tests.

---

## Phase 1 — Mesh I/O and per-cell geometry

### Goal
After Phase 1 the new module can load a `<base>_fault.vtu` and
`<base>_bulk.vtu` pair, expose triangle centroids/normals/areas and
tetrahedron centroids/volumes as plain NumPy arrays, and round-trip
the meshes back to VTU with arbitrary new cell-data attached.

### Files to Create
- `code_preprocess/project_to_fault_stress.py` — new module skeleton
  with mesh-I/O functions only (Phase 2/3/4 fill in the rest).

### Files to Modify
- None.

### Detailed Requirements

1. Module header, mirroring `nw_cut_strip.py` (lines 1–46):
   - Triple-quoted module docstring covering purpose, conventions,
     CLI synopsis, and pytest invocation.
   - Standard imports: `argparse`, `json`, `sys`, `pathlib.Path`,
     `dataclasses.dataclass`, `numpy as np`, `meshio`.
   - Sibling import of the H&Z module via `sys.path.insert(0, ...)`,
     exactly as `clean_freesurface_mesh.py` does for `ts_to_stl`
     (clean_freesurface_mesh.py:56–59). Imports needed:
     `build_bulk_stress_tensor`, `fault_basis_vectors`,
     `resolve_traction`, `ResolvedTraction`.

2. Module-level constants (no magic numbers; names match `nw_cut_strip.py`
   style):
   - `EPS: float = 1.0e-9`
   - `UP_VECTOR: np.ndarray = np.array([0.0, 0.0, 1.0])`
   - `NEAR_HORIZONTAL_NZ_TOL: float = 1.0 - 1.0e-6`   # |n_z| above this
     means dip ≈ 0 and strike is degenerate
   - `DEFAULT_BULK_NAME: str = "rock"`
   - `DEFAULT_FAULT_NAME: str = "fault"`
   - `DEFAULT_OUT_DIR: Path = Path(__file__).resolve().parent /
     "data_projection_onfaultstress"`  # user-explicit override
     (R-701): outputs land in `code_preprocess/data_projection_onfaultstress/`
     so they are co-located with the plan document, not at the
     project-root level.
   - `DEFAULT_MESH_GLOB: str = "safs_fault_box_nwcut_*m.msh"`  # used
     by `--batch`
   - `DEFAULT_INPUT_BASES: tuple[str, ...] = (
        "safs_fault_box_nwcut_500m",
        "safs_fault_box_nwcut_1000m",
        "safs_fault_box_nwcut_2000m",
     )`  # what `--batch` iterates over; derived from
     `ls code_meshing/`, never hardcoded MB sizes or counts.

3. Add a dataclass `FaultCellGeometry`:
   ```python
   @dataclass
   class FaultCellGeometry:
       centroids: np.ndarray   # (N_tri, 3)
       normals: np.ndarray     # (N_tri, 3), unit, orientation-harmonised
       strikes: np.ndarray     # (N_tri, 3), unit; NaN for degenerate rows
       dips: np.ndarray        # (N_tri, 3), unit; NaN for degenerate rows
       areas: np.ndarray       # (N_tri,)
       n_degenerate: int       # count of rows where strike was NaN'd
   ```

4. Add a dataclass `BulkCellGeometry`:
   ```python
   @dataclass
   class BulkCellGeometry:
       centroids: np.ndarray   # (N_tet, 3)
       volumes: np.ndarray     # (N_tet,)
   ```

5. Functions to add (all pure, no I/O side-effects):

   ```python
   def load_fault_mesh(path: Path) -> meshio.Mesh
   def load_bulk_mesh(path: Path) -> meshio.Mesh

   def triangle_geometry(points: np.ndarray,
                         tri_conn: np.ndarray) -> tuple[
                             np.ndarray,  # centroids (N, 3)
                             np.ndarray,  # raw normals (N, 3), unit
                             np.ndarray,  # areas (N,)
                         ]
   def tet_geometry(points: np.ndarray,
                    tet_conn: np.ndarray) -> tuple[
                        np.ndarray,  # centroids (N, 3)
                        np.ndarray,  # signed volumes (N,)
                    ]
   ```

   Implementation specifics:
   - `triangle_geometry` computes `n = (p1-p0) × (p2-p0)`; area =
     `0.5 * |n|`; centroid = `(p0+p1+p2)/3`; unit normal = `n / |n|`.
     Mirror the algebra in `msh_to_vtu.py:43-58` and `ts_to_stl.py:154-160`.
   - `tet_geometry` computes `V_signed = (1/6) * det([p1-p0, p2-p0, p3-p0])`
     and uses `|V_signed|` as the cell weight; centroid = mean of the four
     vertices. Mirror `msh_to_vtu.py:tet_quality:40-45`.

6. Helper to extract a single physical-group from a meshio.Mesh, factored
   out of `msh_to_vtu.py:120-145`:

   ```python
   def extract_cells_by_physical(
       mesh: meshio.Mesh,
       cell_type: str,        # "triangle" or "tetra"
       phys_name: str,
   ) -> np.ndarray   # connectivity (N, 3) or (N, 4)
   ```

   This must read `mesh.field_data` and `mesh.cell_data["gmsh:physical"]`
   exactly as `msh_to_vtu.py:90-145` does. Re-implementing here keeps
   the dependency one-way (we read VTU directly so the helper actually
   keys on the `gmsh:physical` cell-data array that meshio attaches to
   VTU outputs of `msh_to_vtu.py`).

   **Implementation note**: because we will consume `*_fault.vtu`
   directly, the field_data block may be missing (meshio's VTU writer
   doesn't always carry `field_data`). In that case the physical tag
   is the unique value in `cell_data["gmsh:physical"]` for that cell
   block, and the function returns the entire block. Guard with:
   `if len(np.unique(phys_arr)) == 1: return all rows`.

### Interfaces

| Function | Inputs | Outputs |
|---|---|---|
| `load_fault_mesh(path)` | path to `*_fault.vtu` | meshio.Mesh |
| `load_bulk_mesh(path)` | path to `*_bulk.vtu` | meshio.Mesh |
| `triangle_geometry(pts, conn)` | (N_pt,3), (N_tri,3) | centroids, normals, areas |
| `tet_geometry(pts, conn)` | (N_pt,3), (N_tet,4) | centroids, volumes |

### Edge Cases to Handle
- VTU file missing: `FileNotFoundError` with the absolute path in the
  message, matching `nw_cut_strip.py:load_vertices` (lines 91–106).
- VTU file has no triangle/tetra cells: `ValueError` with a clear
  message (same idiom as `msh_to_vtu.py:138-145`).
- Degenerate triangle (zero area): set its normal/centroid to NaN
  vectors and exclude from the degenerate-row count (Phase 2 picks
  up NaN propagation).
- Mesh in non-UTM coordinates (e.g. someone pipes in a different mesh):
  detect via `points[:, 0].max() < 10_000`, emit a warning suggesting
  CRS check. Do **not** fail; the math is purely linear-algebraic so
  it still works.

### Acceptance Criteria
- [ ] `load_fault_mesh("safs_fault_box_nwcut_2000m_fault.vtu")` returns
      a `meshio.Mesh` with `len(mesh.cells[0].data) == 2685` and a single
      triangle block. (Cross-checks the earlier bash probe.)
- [ ] `triangle_geometry` on the same mesh returns centroids whose
      bounding box matches `[X=303 000..697 500, Y=3.612e6..3.903e6,
      Z=−66 607..100]`.
- [ ] `tet_geometry` on the same bulk VTU returns 145 300 cells with
      positive volumes; sum of volumes equals the gmsh box volume within
      0.1 %.
- [ ] `pytest -q` passes the new geometry unit tests
      (`test_triangle_geometry_against_ts_to_stl`,
      `test_tet_geometry_matches_paraview_quality` — both compare
      against precomputed analytic values for a synthetic mesh).

### Dependencies
- Depends on: nothing (Phase 1 is the foundation).
- Required by: Phases 2, 3, 4.

---

## Phase 2 — Per-triangle fault basis with global orientation harmonisation

### Goal
After Phase 2 the module exposes a `build_fault_basis(...)` function
that returns a `FaultCellGeometry` whose `(strikes, dips, normals)`
form a right-handed orthonormal frame at every triangle, with
**consistent global orientation**: for a right-lateral SAF with NW
strike (azimuth 314°), the normal points to the SW half-space
everywhere (H&Z convention: n is "to the left of strike in map
view"; see `hickman_and_zoback_regional_stress_projection.py`
lines 159-166).

### Files to Create
- (none — extends `project_to_fault_stress.py`)
- `code_preprocess/PLAN_onfaultstress_notes.md` — *optional* short
  appendix capturing the orientation-harmonisation derivation; only
  create if Phase 2 requires more discussion than fits in code
  comments.

### Files to Modify
- `code_preprocess/project_to_fault_stress.py` — add basis-builder.

### Detailed Requirements

1. Add the function

   ```python
   def per_triangle_basis_raw(
       normals: np.ndarray,                 # (N, 3), unit
       up: np.ndarray = UP_VECTOR,
   ) -> tuple[
       np.ndarray,  # strikes (N, 3), unit; rows may be NaN if degenerate
       np.ndarray,  # dips (N, 3), unit; rows may be NaN if degenerate
       np.ndarray,  # bool mask (N,) True = degenerate
   ]
   ```

   following the Tandem `facetBasis` algorithm
   (`tandem/src/geometry/Curvilinear.cpp:281-292`):

   ```
   for each i:
       s_i = up × n_i
       if |s_i| < tol:   # n_i ~ ±up: horizontal triangle
           strikes[i] = NaN, dips[i] = NaN, degen[i] = True
       else:
           strikes[i] = s_i / |s_i|
           dips[i]   = strikes[i] × n_i   (Tandem / SEAS convention:
                                            d = s × n = DOWN-dip;
                                            matches `facetBasis` at
                                            tandem/src/geometry/Curvilinear.cpp:291
                                            and CLAUDE.md "Fault-local
                                            tangent frame" rule
                                            `can_t1 = (0, 0, -1)`.
                                            The H&Z `fault_basis_vectors`
                                            uses the OPPOSITE convention
                                            `d_HZ = n × s = UP-dip`; we
                                            deliberately diverge from
                                            H&Z to match the production
                                            C++ consumer.)
   ```

   Tolerance: `np.cross(up, n)` norm < `np.sqrt(1.0 - NEAR_HORIZONTAL_NZ_TOL**2)`.

   Math justification: with `up = ẑ` and `n` unit, `s = ẑ × n =
   (−n_y, n_x, 0)` lies in the horizontal plane and is perpendicular
   to both `ẑ` and `n` — i.e., it is along strike of the fault plane
   at this triangle. `d = s × n` then completes the right-handed
   triple and is the **down-dip** direction (Tandem / SEAS canonical
   frame; CLAUDE.md "Fault-local tangent frame" rule). The H&Z
   module uses the opposite (up-dip) sign convention by virtue of
   `d_hat = np.cross(n_hat, s_hat)` at lines 178-180; the parity
   test `test_resolve_traction_against_HZ_safod` (Phase 3 §1
   cross-check) accounts for this by comparing tau_dip magnitudes
   while inverting the sign relative to H&Z.

2. Add the function

   ```python
   def harmonise_normal_orientation(
       normals: np.ndarray,                 # (N, 3) unit
       rake_sense: str,                     # "right-lateral" or "left-lateral"
       trace_centroids: np.ndarray,         # (N, 3) — for picking the global ref
       fault_strike_azimuth_hint_deg: float | None,
   ) -> np.ndarray   # (N, 3) flipped to consistent orientation
   ```

   Algorithm (right-lateral case; left-lateral is symmetric):

   a. Compute the median triangle position `p_med` (median over each
      coordinate). This is robust to outlier surface anchors that
      `nw_cut_strip.py` already accommodates.

   b. Determine the **expected global outward direction** `n_global`
      from the fault-strike hint:
      - If `fault_strike_azimuth_hint_deg` is provided, compute
        `s_hat_global, _, n_hat_global = fault_basis_vectors(
              fault_strike_azimuth_hint_deg, dip_deg=90.0,
              rake_sense=rake_sense)`
        and use that `n_hat_global` as the reference.
      - If `None`, fit a best-plane PCA to all triangle centroids,
        take the smallest-eigenvalue eigenvector as the global fault
        normal, and pick its sign by majority vote against the raw
        triangle normals.

   c. Flip per-triangle normals so that `dot(n_i, n_global) > 0`. This
      sign-flip then propagates through the basis: after the flip we
      re-run `per_triangle_basis_raw` so that `s` and `d` use the
      corrected normal.

   d. The function must be **idempotent**: calling it twice yields the
      same output.

   **Why this is the right approach** (not just lifting the H&Z
   global-fault construction): the ALT6 Fuis–Seih SAF representation
   the user is targeting is *non-planar*. A single global strike will
   misorient the basis at any triangle whose local strike deviates by
   more than ~10°. The PCA-anchored sign flip preserves the global
   half-space convention without forcing a global strike on every
   triangle.

3. Add a top-level helper:

   ```python
   def build_fault_basis(
       fault_mesh: meshio.Mesh,
       fault_phys_name: str = DEFAULT_FAULT_NAME,
       rake_sense: str = "right-lateral",
       fault_strike_azimuth_hint_deg: float | None = 314.0,  # SAF default
       up: np.ndarray = UP_VECTOR,
   ) -> FaultCellGeometry
   ```

   that wires Phase 1 geometry + the two new functions. Default
   strike hint `314.0` matches the H&Z demo (SAF strike N46W);
   passing `None` triggers PCA fallback.

4. Add the cell→node averager described in Phase 4 §"Cell→Node
   averager":

   ```python
   def cell_to_node_average(
       points: np.ndarray,            # (N_pt, 3)
       tri_conn: np.ndarray,          # (N_tri, 3)
       cell_values: np.ndarray,       # (N_tri, ...) — scalar or vector
       areas: np.ndarray,             # (N_tri,)
   ) -> np.ndarray                    # (N_pt, ...)  — node-averaged
   ```

   - Build the per-vertex area sum `w[v] = Σ_{t ∈ T_v} A_t` via
     `np.add.at(w, tri_conn.flatten(), np.repeat(areas, 3))`.
     Branch-free.
   - Accumulate weighted values
     `acc[v] = Σ_{t ∈ T_v} A_t · cell_values[t]` as follows. Let
     `K = cell_values.shape[1:]` (the trailing shape: `()` for
     scalar, `(3,)` for vector). The single correct broadcast is:

     ```python
     weighted = (
         areas.reshape((-1,) + (1,) * (1 + len(K)))      # (N_tri, 1, *K)
         * cell_values.reshape((N_tri, 1) + K)           # (N_tri, 1, *K)
     ).repeat(3, axis=1).reshape((3 * N_tri,) + K)       # (3*N_tri, *K)
     np.add.at(acc, tri_conn.flatten(), weighted)
     ```

     This works uniformly for `K = ()` and `K = (3,)` (or higher
     trailing dimensions). Document this precondition in the
     docstring: `cell_values.shape[0]` must equal `N_tri`, and
     `cell_values.shape[1:]` is the "trailing" broadcast partner.
   - Final: `node = acc / w.reshape((-1,) + (1,) * len(K))` with
     the divisor reshaped to broadcast against `acc.shape[1:]`.
     Vertices with `w == 0` (orphans) get `np.nan`.

   The function must work for both scalar `cell_values` (shape
   `(N_tri,)`) and vector `cell_values` (shape `(N_tri, K)` for any
   trailing K). Test both.

5. Add the orthonormal-basis-at-vertex helper (used by the fault
   VTU writer to project the per-triangle `(s, d, n)` basis onto
   per-vertex slots):

   ```python
   def basis_to_node(
       points: np.ndarray,
       tri_conn: np.ndarray,
       strikes_cell: np.ndarray,      # (N_tri, 3)
       dips_cell:    np.ndarray,
       normals_cell: np.ndarray,
       areas:        np.ndarray,
   ) -> tuple[np.ndarray, np.ndarray, np.ndarray]
       # strikes_node, dips_node, normals_node, each (N_pt, 3) and unit
   ```

   Implementation:
   1. Compute `normals_node = cell_to_node_average(..., normals_cell)`
      then `normals_node /= |normals_node|`.
   2. Compute `s_avg = cell_to_node_average(..., strikes_cell)`.
   3. Gram–Schmidt onto the plane perpendicular to `normals_node`:
      `s_proj = s_avg - (s_avg · n_node) · n_node`.
   4. Normalise: `strikes_node = s_proj / |s_proj|`.
   5. `dips_node = cross(strikes_node, normals_node)` (down-dip,
      Tandem; matches `per_triangle_basis_raw`'s `d = s × n` and
      CLAUDE.md "Fault-local tangent frame" rule `can_t1 = (0, 0,
      -1)`). Already unit. **NOT** `n × s` (H&Z up-dip), which
      would invert the dip sign relative to the cell-data dip.
   6. Mark degenerate vertices (`|s_proj| < EPS`) with NaN rows.

### Interfaces

The public new API surface from this phase is `build_fault_basis`,
`cell_to_node_average`, and `basis_to_node`. Internal helpers
`per_triangle_basis_raw` and `harmonise_normal_orientation` are
exported (no leading underscore) so that tests can drive them
directly.

### Edge Cases to Handle
- Horizontal triangle (dip ≈ 0): NaN strike/dip, counted in
  `FaultCellGeometry.n_degenerate`.
- All triangles co-planar (typical for the box-cut SAF): PCA gives a
  clean normal; the third eigenvalue is much smaller than the other
  two. Assert this ratio < 1e-2 and warn otherwise.
- Mesh with multiple fault physical groups (future multi-fault case):
  the function only sees one group at a time; the **caller** loops
  over groups. Document this in the docstring; do not silently merge.

### Acceptance Criteria
- [ ] On a synthetic mesh of one vertical NW-striking triangle (strike
      azimuth 314°), `build_fault_basis` returns `strikes[0]` aligned
      with H&Z's `s_hat` for strike 314°, vertical dip, right-lateral,
      to within 1e-10.
- [ ] Idempotency: `harmonise_normal_orientation(harmonise_normal_orientation(...))`
      equals `harmonise_normal_orientation(...)` exactly.
- [ ] On the real 2000 m SAFS fault VTU, `n_degenerate / N_tri < 0.001`
      (we expect zero degenerate rows for a sub-vertical fault).
- [ ] At every non-degenerate row: `np.dot(s_i, d_i) < 1e-12`,
      `np.dot(s_i, n_i) < 1e-12`, `np.dot(d_i, n_i) < 1e-12`,
      `|s_i| = |d_i| = |n_i| = 1` to within 1e-12.

### Dependencies
- Depends on: Phase 1.
- Required by: Phase 3.

---

## Phase 3 — Stress projection (constant background + depth-dependent option)

### Goal
After Phase 3 the module can resolve the H&Z bulk stress σ⁰ onto every
non-degenerate fault triangle and produce a `ResolvedTraction` (or a
vectorised equivalent) per triangle, with the option to evaluate σ⁰ at
the triangle centroid's depth instead of using a single global tensor.

### Files to Create
- (none — extends `project_to_fault_stress.py`)

### Files to Modify
- `code_preprocess/project_to_fault_stress.py` — add the rotator.

### Detailed Requirements

1. Define a vectorised rotator:

   ```python
   def resolve_traction_per_cell(
       sigma_global: np.ndarray,            # (3, 3) OR (N_tri, 3, 3)
       strikes: np.ndarray,                 # (N_tri, 3)
       dips:    np.ndarray,                 # (N_tri, 3)
       normals: np.ndarray,                 # (N_tri, 3)
       P_p_per_cell: np.ndarray | float = 0.0,
   ) -> dict[str, np.ndarray]
   ```

   **`sigma_global` is in compression-POSITIVE SEAS convention here**
   (σ_seas = −σ_HZ). The single-site sign flip from H&Z
   continuum-mechanics to SEAS happens at the source in
   `bulk_stress_tensor_field` (Phase 3 §2 below); the resolver itself
   is a pure rotation onto Phase 2's Tandem basis with **no further
   sign flips**. Under the unified σ_seas convention with Tandem
   basis (s = up × n, d = s × n; locked in by Phase 2 R-301 / R-102),
   all three emitted scalar components carry the SEAS-internal sign
   convention without special-case logic. See R-501 in `REVIEW.md`
   for the algebraic justification (the round-1 R-001 split-site
   flip was correct under H&Z basis but became wrong when R-102
   switched the basis to Tandem; under Tandem the σ-flip and the
   basis-flip cancel symmetrically on shear, leaving correct signs
   for σ_n, τ_strike, and τ_dip after a single bulk flip).

   Returns a dict with keys (each value is `(N_tri,)` except `traction_vec`):
   - `"sigma_n_total"` — `n · σ_seas · n`, MPa, **compression positive**
   - `"sigma_n_eff"`   — `sigma_n_total - P_p`, MPa, compression positive,
     effective normal stress reduced by pore pressure (consistent
     with `friction = μ · σ_n_eff` where σ_n_eff ≥ 0).
   - `"tau_strike"`    — `s · σ_seas · n`, MPa (positive = right-lateral
     under the Tandem basis; matches BP5 / TPV102 production
     convention and the C++ Phase 6 `ProjectFaultPreStress` output)
   - `"tau_dip"`       — `d · σ_seas · n`, MPa (positive = reverse,
     consistent with BP5 t1 = down-dip convention)
   - `"tau_magnitude"` — `sqrt(tau_s² + tau_d²)`, MPa
   - `"rake_deg"`      — `atan2(tau_d, tau_s)` in degrees
   - `"mu_apparent"`   — `tau_magnitude / sigma_n_eff`, NaN where
     `sigma_n_eff <= 0` (effective tension — frictional Mohr-Coulomb
     model not applicable). Tracks effective-tension DOFs and logs
     the count via the JSON summary's `stats.mu_apparent.nan_count`.
   - `"traction_vec"`  — `(N_tri, 3)`, the full global-frame traction
     in **SEAS (compression-positive) convention**: `t = σ_seas · n`.
     Consistent with `sigma_n_total`, `tau_strike`, `tau_dip` (all
     computed on σ_seas). A downstream consumer can verify the
     rotation by checking `n · traction_vec ≈ sigma_n_total`,
     `s · traction_vec ≈ tau_strike`, `d · traction_vec ≈ tau_dip`
     — all consistent without sign-asymmetry caveats.

   Math (the contract for the implementer — implementer must use these
   exact formulas):

   For a single triangle with Tandem basis `(s, d, n)` and σ_seas
   (compression positive; the bulk-path flip happens once in
   `bulk_stress_tensor_field` Phase 3 §2):
   ```
       t      = σ_seas · n            # 3-vector, SEAS convention
       σ_n    = nᵀ σ_seas n           # scalar, compression positive
       τ_s    = sᵀ σ_seas n           # scalar, right-lateral positive
       τ_d    = dᵀ σ_seas n           # scalar, reverse positive
       σ_n,eff = σ_n - P_p            # P_p REDUCES the effective normal load
       |τ|     = √(τ_s² + τ_d²)
       rake   = atan2(τ_d, τ_s)
       μ_app  = |τ| / σ_n,eff         (NaN if σ_n,eff ≤ 0)
       traction_vec = t               # 3-vector emitted to VTU directly
   ```

   **Cross-check against the H&Z scalar path**: H&Z's
   `resolve_traction` returns `sigma_n_eff = sigma_n_total + P_p`
   in continuum-mechanics convention (compression negative).
   Phase 3 output, computed on σ_seas = −σ_HZ with Tandem basis
   (s_T = up × n = −s_HZ, d_T = s × n = −d_HZ, n_T = n_HZ), is
   related to H&Z's output by **two simultaneous sign flips** on
   each rotated component (σ → −σ_HZ AND basis vector →
   −its-H&Z-counterpart for s and d), which **cancel** for the
   shear components but compose for the normal:

   - σ_n_seas = n_T · σ_seas · n_T = n_HZ · (−σ_HZ) · n_HZ
              = −(n_HZ · σ_HZ · n_HZ) = −σ_n_HZ
     → flips sign (compression-negative → compression-positive)
   - τ_strike_seas = s_T · σ_seas · n_T = (−s_HZ) · (−σ_HZ) · n_HZ
                   = s_HZ · σ_HZ · n_HZ = τ_strike_HZ
     → SAME sign and magnitude as H&Z (both flips cancel)
   - τ_dip_seas = d_T · σ_seas · n_T = (−d_HZ) · (−σ_HZ) · n_HZ
                = d_HZ · σ_HZ · n_HZ = τ_dip_HZ
     → SAME sign and magnitude as H&Z (both flips cancel)

   The test `test_resolve_traction_against_HZ_safod` (Phase 3 tests)
   asserts:
   - `σ_n_seas       = -σ_n_HZ`       to 1e-12
   - `τ_strike_seas  =  τ_strike_HZ`  to 1e-12
   - `τ_dip_seas     =  τ_dip_HZ`     to 1e-12

   For SAFOD parameters (verified live in Phase 0):
     σ_n_seas ≈ +104.78 MPa,
     τ_strike_seas ≈ +21.41 MPa (right-lateral positive ✓),
     τ_dip_seas ≈ 0 (pure strike-slip).

   Vectorised: when `sigma_global` is `(3, 3)`,
   ```
       t = np.einsum('ij,kj->ki', sigma_global, normals)   # (N, 3)
   ```
   When `sigma_global` is `(N_tri, 3, 3)` (depth-dependent),
   ```
       t = np.einsum('kij,kj->ki', sigma_global, normals)
   ```
   Then component projections:
   ```
       # sigma_global is σ_seas (compression positive, SEAS
       # convention) — all three components rotated uniformly.
       sigma_n_total =  np.einsum('ki,ki->k', normals, t)
       tau_strike    =  np.einsum('ki,ki->k', strikes, t)
       tau_dip       =  np.einsum('ki,ki->k', dips,    t)
   ```
   The H&Z `resolve_traction` function (lines 216–269) is the
   single-triangle reference; the new vectorised path **must** agree
   with the scalar H&Z path to within 1e-12 on a random sample of
   triangles (unit test enforces this).

2. Define a depth-dependent σ⁰ builder. **This is the SINGLE
   sign-flip site for the entire pipeline**: the function calls
   H&Z's `build_bulk_stress_tensor` (which returns compression-
   negative tensors per continuum-mechanics convention) and
   multiplies by −1 once at the source, returning σ_seas in
   compression-positive SEAS convention. Every downstream consumer
   (resolver, fault VTU, bulk VTU, sidecar, C++ projector, verifier)
   sees compression-positive Pa/MPa without any further sign flips.

   ```python
   def bulk_stress_tensor_field(
       centroids_z: np.ndarray,            # (N,) z-coordinate (≤ 0 in subsurface)
       *,
       SHmax_az_deg: float,
       depth_model: str = "constant",      # "constant" | "lithostatic_sv"
       SHmax_top: float, Shmin_top: float, Sv_top: float,
       SHmax_grad: float = 0.0,            # MPa/m
       Shmin_grad: float = 0.0,
       Sv_grad: float = 0.0,
   ) -> np.ndarray  # (N, 3, 3) in compression-POSITIVE SEAS convention
   ```

   Implementation:
   ```python
   def bulk_stress_tensor_field(centroids_z, *, SHmax_az_deg, ...):
       sigma_HZ = _build_HZ_field(centroids_z, ...)            # (N, 3, 3), compression NEG
       sigma_seas = -sigma_HZ                                  # one-shot source-site flip
       _assert_symmetric(sigma_seas)                           # within 1e-9 MPa
       return sigma_seas                                       # SEAS convention henceforth
   ```

   - `"constant"` (default): every row equals
     `-build_bulk_stress_tensor(SHmax_top, Shmin_top, Sv_top, SHmax_az_deg)`
     (H&Z output flipped to compression positive).
   - `"lithostatic_sv"`: principal magnitudes at depth `d = max(0, -z)`
     are `SHmax(d) = SHmax_top + SHmax_grad * d` (and analogously for
     Shmin and Sv). Each row rotates the (depth-evaluated)
     principal-frame tensor about the vertical exactly as in
     `build_bulk_stress_tensor`, then the whole stack is flipped
     uniformly to compression-positive. The user's later "static
     solve for heterogeneous stress field" plugs in here.
   - All gradients default to 0.0, so depth-dependence is opt-in by
     CLI flag. **Do not** infer a default ρ̄·g lithostatic gradient
     without an explicit CLI value (no hardcoded numerics rule).
   - **Unified sign convention**: σ_seas is compression-positive
     throughout the pipeline downstream of this function. The
     resolver (§1), the fault VTU and bulk VTU writers (Phase 4),
     the schema-v1 sidecar writer (Phase 5), the C++ consumer
     (Phase 6), and the verifier (Phase 8) all consume
     compression-positive Pa/MPa with no further sign manipulation.
     Under Phase 2's Tandem basis (s = up × n, d = s × n), this
     unified σ_seas convention gives correct signs for σ_n
     (compression positive), τ_strike (right-lateral positive), and
     τ_dip (reverse positive) simultaneously — the basis-flip and
     the σ-flip cancel symmetrically on shear (see R-501 algebra in
     `REVIEW.md`).

3. Pore-pressure field builder:

   ```python
   def pore_pressure_field(
       centroids_z: np.ndarray,            # (N,)
       *,
       P_p_top: float,                     # MPa at z=0 (default = 0)
       P_p_grad: float = 0.0,              # MPa/m
   ) -> np.ndarray  # (N,)
   ```

   Returns `P_p_top + P_p_grad * max(0, -z)` per cell. Default keeps
   the H&Z demo behaviour (single P_p value applied uniformly).

4. Wire-up function:

   ```python
   def project_stress_onto_fault(
       fault_geom: FaultCellGeometry,
       *,
       SHmax: float, Shmin: float, Sv: float, P_p: float,
       SHmax_az_deg: float,
       depth_model: str = "constant",
       SHmax_grad: float = 0.0, Shmin_grad: float = 0.0,
       Sv_grad:    float = 0.0, P_p_grad: float = 0.0,
   ) -> dict[str, np.ndarray]
   ```

   Returns the same dict as `resolve_traction_per_cell` plus a key
   `"sigma_field"` `(N_tri, 3, 3)` so the bulk path (Phase 4) can reuse
   the rotation results without re-evaluating the depth model.

### Interfaces
- New public functions:
  `resolve_traction_per_cell`, `bulk_stress_tensor_field`,
  `pore_pressure_field`, `project_stress_onto_fault`.

### Edge Cases to Handle
- Degenerate row (NaN basis): every output for that row is NaN. The
  rotation math naturally propagates NaN; no special-casing needed
  beyond making sure `mu_apparent` is NaN-safe.
- `sigma_n_eff ≤ 0` (effective tension; compression-positive
  convention is violated either because the input σ⁰ is small or
  because P_p exceeds σ_n_total): `mu_apparent[i] = np.nan` and the
  per-cell event is counted in the JSON summary's
  `stats.mu_apparent.nan_count`. Document in the docstring; under
  the unified compression-positive convention the H&Z `|·|` guard
  is not appropriate because a negative σ_n_eff is unphysical for
  Mohr-Coulomb friction and must NOT be silently masked into a
  positive μ.
- A user passing a non-symmetric σ⁰ accidentally: assert
  `np.allclose(sigma_global, sigma_global.swapaxes(-1, -2), atol=1e-9)`
  at the entry point of `resolve_traction_per_cell`. If not symmetric,
  raise `ValueError` with the maximum off-diagonal asymmetry.

### Acceptance Criteria
- [ ] On the SAFOD parameters
      (`SHmax=113, Shmin=49, Sv=45, P_p=16, SHmax_az=23`) and a single
      synthetic triangle with strike azimuth 314° / dip 90° / right-lateral,
      `project_stress_onto_fault` returns `sigma_n_eff ≈ +89` MPa
      (compression-positive SEAS convention; **note the sign flip
      relative to H&Z's `-89` MPa**) and `mu_apparent ≈ 0.24` (sign
      flip cancels because μ uses magnitudes). The test must accept
      values within 1e-3 MPa / 1e-4.
- [ ] Trace invariance: for any orthonormal frame `(s, d, n)`,
      `σ_n + sᵀσs + dᵀσd = trace(σ)` to within 1e-9 MPa
      (unit test sweeps 100 random frames and random σ).
- [ ] Vectorised result equals 100 independent scalar
      `resolve_traction(...)` calls on the same per-triangle bases,
      element-wise within 1e-12.
- [ ] With `depth_model="lithostatic_sv"` and non-zero gradients, the
      σ⁰ field at z = −1671 m equals (depth-extrapolated) H&Z values
      to within 1e-9.

### Dependencies
- Depends on: Phases 1, 2.
- Required by: Phase 4.

---

## Phase 4 — VTU + JSON writers and the CLI

### Goal
After Phase 4, running the script writes the three artefacts under
`data_projection_onfaultstress/<lc>m/`, both via a single-file invocation
and via `--batch`. Outputs are loadable in ParaView and the JSON
summary captures every input parameter plus min/median/max for each
output field.

### Files to Create
- (none — extends `project_to_fault_stress.py`)

### Files to Modify
- `code_preprocess/project_to_fault_stress.py` — add writers and CLI.

### Detailed Requirements

1. Fault VTU writer. **The fault VTU emits a continuous (point-data,
   piecewise-linear) field as its primary output**, plus the raw
   per-triangle (cell-data) field as a secondary `_cell` channel so
   the user can spot-check the source values that fed the averaging.

   ```python
   def write_fault_vtu(
       fault_mesh: meshio.Mesh,           # original input mesh (point cloud reused)
       tri_conn: np.ndarray,              # (N_tri, 3) — fault-phys triangles
       fault_geom: FaultCellGeometry,
       resolved: dict[str, np.ndarray],   # output of project_stress_onto_fault (cell-centred)
       out_path: Path,
       convention_str: str,               # for metadata
   ) -> None
   ```

   Inside the writer, **the cell-centred arrays are first projected
   onto the fault vertices** using the area-weighted node-averaging
   helper `cell_to_node_average` (Phase 2 extension below). The
   resulting point-data is the continuous field. Cell-data is kept
   alongside for traceability.

   ### Continuous fault field — math

   Let `T_v = { t : v ∈ vertices(t) }` be the triangles incident on
   vertex `v`, and `A_t` be the area of triangle `t`. The continuous
   value at `v` for any scalar field `f` known per-cell is

   ```
       f_node(v) = ( Σ_{t ∈ T_v}  A_t · f_cell(t) ) / ( Σ_{t ∈ T_v} A_t )
   ```

   This is the classical area-weighted vertex average. The output VTU
   carries the per-vertex `f_node`; ParaView's default rendering
   then interpolates linearly across each triangle (piecewise-linear
   on the fault surface). The resulting field is **C⁰ continuous**
   across triangle edges by construction — no isolated triangle
   values.

   Vector fields (`strike_vec`, `dip_vec`, `normal_vec`,
   `traction_vec`) follow the same recipe component-wise. For the
   unit basis vectors `strike_vec_node` and `dip_vec_node`, the
   area-weighted average is **re-orthonormalised at each vertex**:

   ```
       n_node = normalize( Σ A_t · n_cell )
       s_node = normalize( Σ A_t · s_cell - (Σ A_t · s_cell · n_node) · n_node )   # Gram–Schmidt onto plane ⊥ n_node
       d_node = s_node × n_node                                                     # Tandem down-dip (matches Phase 2 §5; R-704)
   ```

   This guarantees `|n_node| = |s_node| = |d_node| = 1` and
   `s_node ⊥ d_node ⊥ n_node` at every fault vertex. The
   re-orthonormalisation is what distinguishes a "geometrically
   meaningful" continuous basis from a naive componentwise average
   (which generally yields a non-unit, non-orthogonal triple).

   ### Edge vertices (fault boundary)

   Vertices that lie on the fault outer boundary (free-surface trace,
   bottom edge, NW/SE end-caps) have fewer incident triangles than
   interior vertices. The area-weighted recipe still applies; no
   special case is needed mathematically. The `n_degenerate` count
   from Phase 2 must propagate: if any incident cell carries a NaN
   basis, that NaN propagates into the node value. The writer logs
   the count of NaN-tainted vertices to stdout and the JSON summary.

   ### Cell-data and point-data fields written

   | Field name                  | Where         | Shape          | Notes |
   |-----------------------------|---------------|----------------|---|
   | `sigma_n_total_MPa`         | **point-data**| `(N_pt,)`      | compression POSITIVE |
   | `sigma_n_eff_MPa`           | **point-data**| `(N_pt,)`      | after P_p |
   | `tau_strike_MPa`            | **point-data**| `(N_pt,)`      | + for right-lateral |
   | `tau_dip_MPa`               | **point-data**| `(N_pt,)`      | + for reverse |
   | `tau_magnitude_MPa`         | **point-data**| `(N_pt,)`      | ≥ 0 |
   | `rake_deg`                  | **point-data**| `(N_pt,)`      | atan2(τ_d, τ_s) on the node-averaged components |
   | `mu_apparent`               | **point-data**| `(N_pt,)`      | NaN allowed |
   | `strike_vec`                | **point-data**| `(N_pt, 3)`    | unit, re-orthonormalised |
   | `dip_vec`                   | **point-data**| `(N_pt, 3)`    | unit, re-orthonormalised |
   | `normal_vec`                | **point-data**| `(N_pt, 3)`    | unit, area-weighted then renormalised |
   | `traction_vec_MPa`          | **point-data**| `(N_pt, 3)`    | global-frame traction `σ_seas · n_cell` (area-weighted from the cell-data via `cell_to_node_average`; **NOT** recomputed from `σ·n_node` at the vertex). Convention: compression POSITIVE (consistent with `sigma_n_total_MPa`). |
   | `sigma_n_total_MPa_cell`    | cell-data     | `(N_tri,)`     | raw, per-triangle (traceability) |
   | `sigma_n_eff_MPa_cell`      | cell-data     | `(N_tri,)`     | raw |
   | `tau_strike_MPa_cell`       | cell-data     | `(N_tri,)`     | raw |
   | `tau_dip_MPa_cell`          | cell-data     | `(N_tri,)`     | raw |
   | `gmsh:physical`             | cell-data     | `(N_tri,)`     | carried through |
   | `gmsh:geometrical`          | cell-data     | `(N_tri,)`     | if present |

   **`rake_deg` is recomputed from the node-averaged `tau_strike` and
   `tau_dip`, not averaged directly**: averaging an angle wraps the
   periodic discontinuity and produces wrong values near
   ±180°. Same for `tau_magnitude` (recomputed
   `sqrt(τ_s_node² + τ_d_node²)`) and `mu_apparent`. Document this
   in the writer's docstring.

   ### Cell→Node averager — public helper

   Add the function (in Phase 2 to make it available to the verifier
   in Phase 8 as well):

   ```python
   def cell_to_node_average(
       points: np.ndarray,            # (N_pt, 3)
       tri_conn: np.ndarray,          # (N_tri, 3)
       cell_values: np.ndarray,       # (N_tri, ...) — scalar or vector
       areas: np.ndarray,             # (N_tri,)
   ) -> np.ndarray                    # (N_pt, ...)  — node-averaged
   ```

   that builds the incident-cell list once, then accumulates
   `A_t * cell_values[t]` into each vertex slot, dividing by the
   per-vertex area sum at the end. Vertices that are not referenced
   by any triangle (orphans — `nw_cut_strip.py` already strips
   these) get `NaN`.

   The implementation must be O((N_tri · 3) + N_pt) using a single
   numpy `np.add.at(...)` accumulation; do NOT loop in Python over
   triangles.

   ### meshio specifics

   - Point-data is passed as the `point_data=` kwarg to
     `meshio.Mesh(...)`. Vector fields go in as `(N_pt, 3)` numpy
     arrays.
   - Cell-data follows the dict-of-list-per-block convention
     (matches `msh_to_vtu.py:170-198`).
   - The convention string is attached as
     `mesh.info = {"stress_convention": convention_str, "units": "MPa",
     "frame": "(east, north, up) / UTM Zone 11N"}`. If meshio drops
     `info` on VTU write (round-trip test in Phase 4 acceptance
     covers this), fall back to encoding the convention as a 0-D
     `field_data` entry per meshio's existing handling.

2. Bulk VTU writer (gated by `--write-bulk`):

   ```python
   def write_bulk_vtu(
       bulk_mesh: meshio.Mesh,
       tet_conn: np.ndarray,              # (N_tet, 4)
       bulk_geom: BulkCellGeometry,
       sigma_field: np.ndarray,           # (N_tet, 3, 3), compression-POSITIVE SEAS
       out_path: Path,
       convention_str: str,
   ) -> None
   ```

   Cell-data:
   | Field | Shape | Notes |
   |---|---|---|
   | `sigma_xx_MPa` | `(N_tet,)` | compression POSITIVE (SEAS) |
   | `sigma_yy_MPa` | `(N_tet,)` | compression POSITIVE (SEAS) |
   | `sigma_zz_MPa` | `(N_tet,)` | compression POSITIVE (SEAS) |
   | `sigma_xy_MPa` | `(N_tet,)` | compression POSITIVE (SEAS) |
   | `sigma_xz_MPa` | `(N_tet,)` | compression POSITIVE (SEAS) |
   | `sigma_yz_MPa` | `(N_tet,)` | compression POSITIVE (SEAS) |
   | `sigma_tensor_MPa` | `(N_tet, 9)` | row-major 3x3 (`σ_xx, σ_xy, σ_xz, σ_yx, σ_yy, σ_yz, σ_zx, σ_zy, σ_zz`), so ParaView's TensorGlyph filter can consume it. Compression POSITIVE. |
   | `gmsh:physical`/`gmsh:geometrical` | carried through |

   **Pass-through writer (no sign flip).** Under the unified
   single-site flip contract (R-501/R-502), `sigma_field` arrives
   already in compression-positive SEAS convention from Phase 3's
   `bulk_stress_tensor_field`. This writer just splits the
   (N_tet, 3, 3) tensor into six scalar cell-data fields plus the
   9-component flat tensor; no sign manipulation. The docstring
   should explicitly state "expects σ in compression-positive SEAS
   convention; no flip applied" so a future caller doesn't
   accidentally double-flip.

   **Important note re: dynamic-rupture use** — the docstring of this
   function MUST cite the user-stated rule: *"the bulk stress field is
   not consumed by the dynamic-rupture run because that run assumes
   ∇·σ⁰ ≈ 0 (equilibrium); the artefact is for the static-solve
   step downstream"*. This anchors the rationale in the file so future
   readers don't try to read it into the production DG operator.

3. JSON summary writer:

   ```python
   def write_summary_json(
       out_path: Path,
       *,
       input_mesh_path: Path,
       fault_geom: FaultCellGeometry,
       bulk_geom: BulkCellGeometry | None,
       sigma_global_at_z0: np.ndarray,    # (3, 3) for reference
       params: dict[str, float | str],    # full CLI input echo
       resolved: dict[str, np.ndarray],
   ) -> None
   ```

   JSON schema (no nesting except `params` and `stats`):
   ```json
   {
     "input_mesh": "...",
     "fault_n_cells": 2685,
     "fault_n_degenerate_basis": 0,
     "bulk_n_cells": 145300,
     "params": {
        "SHmax_MPa": 113.0, "Shmin_MPa": 49.0, "Sv_MPa": 45.0,
        "P_p_MPa": 16.0, "SHmax_azimuth_deg": 23.0,
        "rake_sense": "right-lateral", "depth_model": "constant",
        "SHmax_grad_MPa_per_m": 0.0, "Shmin_grad_MPa_per_m": 0.0,
        "Sv_grad_MPa_per_m": 0.0, "P_p_grad_MPa_per_m": 0.0,
        "fault_strike_azimuth_hint_deg": 314.0
     },
     "sigma_global_at_z0_MPa": [[...3x3...]],
     "stats": {
        "sigma_n_total_MPa":  {"min": 10.2, "median": 68.3, "max": 120.5, "nan_count": 0},
        "sigma_n_eff_MPa":    {...},
        "tau_strike_MPa":     {...},
        "tau_dip_MPa":        {...},
        "tau_magnitude_MPa":  {...},
        "rake_deg":           {...},
        "mu_apparent":        {...}
     },
     "convention": "compression POSITIVE (SEAS internal); P_p positive; sigma_n_eff = sigma_n_total - P_p; tau_strike + = right-lateral (H&Z natural sign, preserved through the rotation); frame (east, north, up) UTM Zone 11N; units MPa"
   }
   ```

   Use `json.dump(..., indent=2, allow_nan=True)`. NaN serialisation
   matters when `mu_apparent` is NaN; `allow_nan=True` is the default
   in CPython but document it.

4. CLI (mirroring `nw_cut_strip.py:main`, with three mutually exclusive
   modes):

   ```
   python project_to_fault_stress.py <INPUT_MESH> <MESH_BASE>
       [--bulk-name rock] [--fault-name fault]
       [--hz-dump-file PATH]
       [--SHmax 113.0] [--Shmin 49.0] [--Sv 45.0] [--P_p 16.0]
       [--SHmax-az 23.0] [--rake-sense right-lateral]
       [--strike-hint-az 314.0] [--no-strike-hint]
       [--depth-model {constant,lithostatic_sv}]
       [--SHmax-grad 0.0] [--Shmin-grad 0.0] [--Sv-grad 0.0]
       [--P_p-grad 0.0]
       [--write-bulk]
       [--out-dir DATA_PROJECTED_ONFAULTSTRESS]

   python project_to_fault_stress.py --batch [same flags above; iterates
       over DEFAULT_INPUT_BASES, reading code_meshing/<base>_fault.vtu
       and (if --write-bulk) code_meshing/<base>_bulk.vtu]

   python project_to_fault_stress.py --print-info INPUT_MESH
       (prints the bulk σ⁰, fault cell count, strike-azimuth statistics,
        then exits without writing)
   ```

   When `--hz-dump-file <path>` is passed, the file is loaded via
   `json.load`, validated to match `schema == "hickman_zoback_sigma0_v1"`,
   and its `params` block is used as the source of `SHmax`, `Shmin`,
   `Sv`, `P_p`, and `SHmax_az_deg`. Per-flag CLI overrides still
   win (an explicit `--SHmax 200.0` on the command line beats the
   dump-file value). This wires Phase 0's dump into Phase 3's
   parameter pipeline so the H&Z module remains the single source
   of truth.

   Defaults are **identical to H&Z `demo_safod` SAFOD parameters**.
   The defaults are loaded by *importing the constants the H&Z module
   uses for demo_safod*, **not** by retyping the numbers. Implement
   by reading the docstring constants once at module load and exposing
   them as `DEFAULTS = dict(SHmax=113.0, Shmin=49.0, Sv=45.0, P_p=16.0,
   SHmax_az=23.0)`. If the H&Z file is later parameter-swept, the user
   will edit both — that's fine, but document the coupling.

5. Mode logic, replicating `nw_cut_strip.py:main`:
   - Reject combining `--batch` with positional args (`argparse`
     mutually-exclusive groups).
   - In `--batch`, glob `DEFAULT_MESH_GLOB` under `code_meshing/`; for
     each match, derive `out_base = <data_projection_onfaultstress>/<lc>m/`
     and recurse into the single-file path.
   - In single-file mode, the user passes the mesh's `*_fault.vtu`
     path. The corresponding bulk VTU is the same stem with `_fault`
     swapped for `_bulk`. If `--write-bulk` and the bulk VTU does not
     exist, error out with a clear pointer to `msh_to_vtu.py`.

6. Logging / stdout style: mirror `msh_to_vtu.py:84-204`. Print one
   line per phase (`Reading ...`, `Loaded ... triangles ...`, `Computed
   per-triangle basis ...`, `Rotated stress ...`, `Wrote ...`).

### Interfaces

The CLI is the public surface. The module is also importable from
other Python files (`from project_to_fault_stress import
project_stress_onto_fault, write_fault_vtu` etc.).

### Edge Cases to Handle
- Output directory does not exist: `out_path.parent.mkdir(parents=True,
  exist_ok=True)` (matches `msh_to_vtu.py:202`).
- `--write-bulk` requested but bulk VTU missing: `FileNotFoundError`
  with a one-line hint suggesting to run `msh_to_vtu.py` first.
- `--depth-model lithostatic_sv` requested with all gradients zero:
  warn (`f"depth_model=lithostatic_sv with all gradients = 0 is
  equivalent to 'constant'"`).
- Mesh has no degenerate triangles vs has some: both should pass
  through; `n_degenerate` is logged to stdout and written to the
  summary.

### Acceptance Criteria
- [ ] Running
      ```
      python project_to_fault_stress.py
          ../code_meshing/safs_fault_box_nwcut_2000m_fault.vtu
          safs_fault_box_nwcut_2000m
          --write-bulk
      ```
      produces
      `data_projection_onfaultstress/2000m/safs_fault_box_nwcut_2000m_fault_stress.vtu`,
      `..._bulk_stress.vtu`, and `..._summary.json`, all readable by
      meshio.read(...) with the expected cell-data field names.
- [ ] `--batch` produces all three output bundles for the 500 m, 1000 m,
      and 2000 m meshes.
- [ ] `--print-info` prints σ⁰ and exits without creating any output
      files.
- [ ] ParaView opens the fault VTU and renders `tau_magnitude_MPa`
      as a colour-mapped cell scalar (manual check; document the
      verification recipe in the file docstring).
- [ ] Summary JSON for the SAFOD-parameter SAF 2000 m run shows
      `mu_apparent.median` in the range `[0.20, 0.35]` (loose sanity
      bound based on H&Z 2004 Fig. 4 / the file's `sweep_strike`
      output).
- [ ] Total runtime on the 500 m mesh (the densest at ~10⁵ triangles
      if upscaled — current is 2685 tris × 4 for 500m ~ 10⁴) is under
      30 seconds wall clock.

### Dependencies
- Depends on: Phases 1, 2, 3.
- Required by: Phase 5 (the sidecar writer reuses the bulk and fault
  centroid+geometry products, plus the resolved-stress dict).

---

## Phase 5 — schema-v1 HDF5 sidecar writer (`stress_safs.h5`)

### Goal
After Phase 5 the preprocessor produces a single canonical
`stress_safs.h5` schema-v1 sidecar (analog of
`data_projected/velocity_safs.h5`) that the C++ runtime can load with
the existing `DataField3D` class **without** any C++-side changes.
This is the artefact the dynamic-rupture / quasi-dynamic solvers will
consume.

### Files to Create
- `code_preprocess/build_stress_safs.py` — new module that calls the
  Phase 3 σ⁰ field builder on a *regular* (Nx, Ny, Nz) UTM grid that
  encloses the SAFS mesh, then hands the six tensor-component arrays
  to `data_projection.sidecar.write_sidecar`.
- `code_preprocess/test_build_stress_safs.py` — pytest companion.

### Files to Modify
- None.

### Detailed Requirements

1. Module header mirrors `data_projection/build_velocity_cvmh.py`:
   docstring; CLI synopsis; sibling-import of
   `data_projection.sidecar.write_sidecar` via the same
   `sys.path.insert(0, str(Path(__file__).resolve().parent /
   "data_projection"))` idiom.

2. Module-level constants (no magic numbers — derive from the H&Z
   `demo_safod` set and the mesh bbox; the user feedback memory
   `feedback_no_hardcoded_numbers.md` applies):
   - `DEFAULT_GRID_DX_M: float = 1000.0`
   - `DEFAULT_GRID_PAD_M: float = 2000.0`    # one mesh cell of padding
   - `DEFAULT_BOUNDS_FACTOR: float = 1.2`    # min/max declared in sidecar
   - `DEFAULT_OUT: Path = Path(__file__).resolve().parent /
     "data_projection_onfaultstress" / "stress_safs.h5"`
     # R-802: outputs land under `code_preprocess/data_projection_onfaultstress/`,
     # co-located with Phase 4 outputs (the post-R-701 user-override
     # location).

3. Add the function

   ```python
   def build_uniform_utm_grid(
       mesh_bbox: tuple[float, float, float, float, float, float],
       *,
       dx_m: float,
       pad_m: float,
   ) -> tuple[np.ndarray, np.ndarray, np.ndarray]
   ```

   that returns `(x, y, z)` 1-D arrays satisfying schema-v1's
   strict-monotone-increasing requirement, padded by `pad_m` on each
   face of the mesh bbox so that any DOF of any candidate SAFS mesh
   resolution (500 m / 1000 m / 2000 m) is strictly inside the grid.
   The mesh bbox must be supplied by the caller (Phase 1 already
   computes it from the input VTU); **do not** hardcode UTM extents.

4. Add the function

   ```python
   def evaluate_stress_field_on_grid(
       x: np.ndarray, y: np.ndarray, z: np.ndarray,
       *,
       SHmax: float, Shmin: float, Sv: float,
       SHmax_az_deg: float,
       depth_model: str = "constant",
       SHmax_grad: float = 0.0, Shmin_grad: float = 0.0,
       Sv_grad: float = 0.0,
   ) -> dict[str, np.ndarray]   # six (Nx, Ny, Nz) arrays in Pa
   ```

   that:
   - Builds a `(Nx, Ny, Nz, 3, 3)` σ_seas field
     (**compression-positive SEAS convention**, already sign-flipped
     once at the source by Phase 3 §2's `bulk_stress_tensor_field`)
     on a vectorised z-only stack (constant in x and y for the
     homogeneous-σ⁰ case).
   - **Unit conversion only** (no sign flip — R-501/R-502 contract):
     `sigma_seas_Pa = 1e6 * sigma_seas_MPa`. The single bulk-path
     flip from H&Z to SEAS happens once in
     `bulk_stress_tensor_field`; this function is a pure
     MPa → Pa pass-through.
   - Returns six independent `(Nx, Ny, Nz)` arrays
     `{"sigma_xx", "sigma_yy", "sigma_zz", "sigma_xy", "sigma_yz",
     "sigma_xz"}` in **Pa, compression POSITIVE SEAS convention**.
   - Asserts the symmetry constraint `σ_xy == σ_yx`, `σ_xz == σ_zx`,
     `σ_yz == σ_zy` holds within 1e-9 MPa (= 1e-3 Pa, matching
     Phase 3's `_SIGMA_SYMMETRY_TOL_MPA` after the MPa → Pa unit
     conversion; R-801) across the entire grid before return
     (single-source-of-truth check; the H&Z rotation should produce
     symmetric tensors).
   - Asserts every array contains no NaN before return (defence in
     depth — `sidecar.write_sidecar` also re-checks).

5. Add the function

   ```python
   def derive_field_bounds(
       fields: dict[str, np.ndarray],
       *,
       safety_factor: float = DEFAULT_BOUNDS_FACTOR,
   ) -> dict[str, tuple[float, float, str]]
   ```

   that returns `{name: (min_value, max_value, "Pa")}` for each input
   field with `min_value = min(field) * safety_factor` if the min is
   negative else `min(field) / safety_factor`, and analogously for
   max. The bounds must enclose the actual field values strictly so
   that schema-v1's guard `G-2` (every cell in
   `[min_value, max_value]`) passes. **Do not** hardcode example
   bounds.

6. Top-level wrapper:

   ```python
   def build_stress_safs(
       mesh_paths: list[Path],
       out_path: Path,
       *,
       SHmax: float, Shmin: float, Sv: float, SHmax_az_deg: float,
       depth_model: str = "constant",
       SHmax_grad: float = 0.0, Shmin_grad: float = 0.0,
       Sv_grad: float = 0.0,
       dx_m: float = DEFAULT_GRID_DX_M,
       pad_m: float = DEFAULT_GRID_PAD_M,
       mesh_tag: str = "safs_fault_box_nwcut",
   ) -> None
   ```

   that:
   1. Reads every mesh in `mesh_paths` (the 500 m / 1000 m / 2000 m
      bulk VTUs from `code_meshing/`) and computes the union bbox.
   2. Calls `build_uniform_utm_grid`.
   3. Calls `evaluate_stress_field_on_grid`.
   4. Calls `derive_field_bounds`.
   5. Calls
      `sidecar.write_sidecar(out_path, x, y, z, fields, attrs={"schema_version": SCHEMA_VERSION, "crs": CANONICAL_CRS, "units": CANONICAL_UNITS, "z_positive": CANONICAL_Z_POSITIVE, "source": ",".join(m.name for m in mesh_paths), "source_crs": "EPSG:32611", "mesh_tag": mesh_tag}, field_bounds={...})`.

7. CLI mirroring `build_velocity_cvmh.py`:
   ```
   python build_stress_safs.py
       --meshes <vtu1> <vtu2> ...
       --out    DATA_PROJECTED_ONFAULTSTRESS/stress_safs.h5
       --SHmax 113.0 --Shmin 49.0 --Sv 45.0 --SHmax-az 23.0
       [--depth-model {constant,lithostatic_sv}]
       [--SHmax-grad 0.0 --Shmin-grad 0.0 --Sv-grad 0.0]
       [--dx 1000.0 --pad 2000.0]
       [--mesh-tag safs_fault_box_nwcut]
   ```

### Interfaces
The public API is `build_stress_safs(...)` and the CLI. The sidecar is
loadable by any schema-v1 reader, including the existing
`DataField3D` C++ class — Phase 6 imposes no additional constraints
on the file.

### Edge Cases to Handle
- Mesh bbox unioned across resolutions: pad once on the union, not on
  each mesh — keeps the sidecar generation idempotent across
  reruns.
- All gradients zero with `depth_model="lithostatic_sv"`: warn (same
  as Phase 3 §Edge Cases).
- `out_path.parent` does not exist: created via
  `sidecar.write_sidecar` (which already does this — see
  `sidecar.py:78`). Do not re-mkdir here; doing so duplicates a
  defensive sibling in the writer.

### Acceptance Criteria
- [ ] Running
      ```
      python build_stress_safs.py
          --meshes ../code_meshing/safs_fault_box_nwcut_500m_bulk.vtu
                   ../code_meshing/safs_fault_box_nwcut_1000m_bulk.vtu
                   ../code_meshing/safs_fault_box_nwcut_2000m_bulk.vtu
          --out    data_projection_onfaultstress/stress_safs.h5
      ```
      produces a sidecar that `h5dump -A` reports with the required
      schema-v1 attributes (`schema_version`, `crs`, `units`,
      `z_positive`).
- [ ] `python -c "import h5py; f=h5py.File('stress_safs.h5'); print(list(f['fields']))"`
      lists exactly the six `sigma_*` names.
- [ ] Loading the same file via the C++ `DataField3D("...", "sigma_xx")`
      ctor in a unit test (Phase 6 stub) succeeds — does not abort —
      with no Phase-6 code changes.
- [ ] The sidecar's grid bbox strictly encloses the 500 m mesh bbox by
      at least `pad_m` on every face.

### Dependencies
- Depends on: Phase 3 (the σ⁰ field builder), Phase 1 (bbox).
- Required by: Phase 6 (the C++ consumer reads this file).

---

## Phase 6 — C++ consumer: `StressField3D`, `ProjectStress`, `ComputeSAFSParams`

### Goal
After Phase 6 the SEAS runtime loads `stress_safs.h5` and uses it to
populate per-DOF pre-stress for the quasi-dynamic and dynamic-rupture
solvers. The seam matches the velocity pipeline: a sidecar reader
backs a per-coefficient Coefficient, which a projector pushes onto
either a bulk H1 space (for visualisation / future static-solve) or
into the existing `FaultGeometry::tau_pre_` / `sigma_n` slots (for
production runs).

### Files to Create
- `io/stress_field_3d.hpp` — bulk 6-component sidecar reader.
- `io/stress_field_coefficient.hpp` — `mfem::VectorCoefficient` (6
  components) backed by six `DataField3D` instances; a thin wrapper.

### Files to Modify
- `io/field_coefficient.hpp` — add
  `FieldProjector::ProjectStress` (bulk, six GFs) and
  `FieldProjector::ProjectFaultPreStress` (fault, per-DOF) static
  methods. Mirror the existing `ProjectVelocity` signature
  (`field_coefficient.hpp:143-150`).
- `io/field_coefficient.cpp` — implementations.
- `fault/fault_geometry.hpp` — add a parallel
  `ComputeSAFSParams(const StressField3D& field, ...)` next to
  `ComputeBP5Params()` (`fault_geometry.hpp:703-757`) that populates
  `tau_pre_` and a new `sigma_n_per_dof_` member array from the
  sidecar instead of the analytic `bp5_params_.tau0_vec / .sigma_n`.
- `fault/rate_state_fault.hpp` — extend the per-DOF prestress
  consumption path (currently `tau0_ = params_.tau0()` at lines 255 /
  691) to optionally read from `FaultGeometry::tau_pre_per_dof()` and
  `sigma_n_per_dof()` when SAFS mode is on (new bool flag on the
  operator).
- `drivers/seas_driver.cpp` — wire a new TOML option
  `[stress]` with `sidecar = "..."` and a `use_sidecar` boolean,
  defaulting to `false`. When `use_sidecar = true`, load the sidecar,
  call `FaultGeometry::ComputeSAFSParams`, and pass the resulting
  per-DOF tau/sigma_n into the operator.
- `config/seas_config.hpp`, `seas_config_parser.hpp`,
  `seas_config_bridge.hpp` — add the `[stress]` section to the
  TOML schema.

### Detailed Requirements

1. `StressField3D` — header-only class in `mfem::seas`:

   ```cpp
   class StressField3D
   {
   public:
       /// Load all six tensor components from a schema-v1 sidecar.
       /// Each component is a separate DataField3D; the ctor instantiates
       /// them in order (sigma_xx, sigma_yy, sigma_zz, sigma_xy,
       /// sigma_yz, sigma_xz) and shares the schema-version /
       /// CRS / bbox checks with DataField3D.
       StressField3D(const std::string& sidecar_path,
                     OOBPolicy oob = OOBPolicy::Abort);

       /// Evaluate the symmetric stress tensor at (x, y, z) in canonical
       /// CRS. Returns a (3, 3) dense matrix in Pa, **compression
       /// POSITIVE (SEAS internal convention)**.  The full-tensor sign
       /// flip from geomechanics happens in the Python preprocessor at
       /// the SOURCE site — Phase 3
       /// `code_preprocess/project_to_fault_stress.py:bulk_stress_tensor_field`
       /// (R-906 reconciliation; the earlier plan revision misattributed
       /// the flip to Phase 5).  Phase 5's
       /// `build_stress_safs.py:evaluate_stress_field_on_grid` is a
       /// pure MPa → Pa unit conversion (no further sign manipulation
       /// — R-501 / R-502).  This C++ reader returns the on-disk
       /// values unchanged.
       mfem::DenseMatrix Evaluate(real_t x, real_t y, real_t z) const;

       /// Access the underlying six DataField3D instances by index 0..5
       /// in the order {xx, yy, zz, xy, yz, xz}.
       const DataField3D& Field(int component_index) const;

       const std::array<real_t, 6>& BBox() const;
       bool ContainsBBox(real_t xmin, real_t xmax,
                         real_t ymin, real_t ymax,
                         real_t zmin, real_t zmax,
                         real_t eps = 0.0) const;

       /// Match the InterpMode propagation pattern from DataField3D
       /// (data_field_3d.hpp:88-90).  Sets all six underlying readers'
       /// InterpMode in lock-step.
       void SetInterpMode(InterpMode mode);
   };
   ```

   Implementation: construct six `DataField3D` members
   `{"sigma_xx", ..., "sigma_xz"}`. `Evaluate(x, y, z)` returns:
   ```cpp
       DenseMatrix S(3);
       const real_t xx = sigma_xx_.Evaluate(x, y, z);
       const real_t yy = sigma_yy_.Evaluate(x, y, z);
       const real_t zz = sigma_zz_.Evaluate(x, y, z);
       const real_t xy = sigma_xy_.Evaluate(x, y, z);
       const real_t yz = sigma_yz_.Evaluate(x, y, z);
       const real_t xz = sigma_xz_.Evaluate(x, y, z);
       S(0,0) = xx; S(1,1) = yy; S(2,2) = zz;
       S(0,1) = xy; S(1,0) = xy;
       S(0,2) = xz; S(2,0) = xz;
       S(1,2) = yz; S(2,1) = yz;
       return S;
   ```
   `BBox()` returns the intersection of all six component bboxes — they
   are identical by construction (same grid), but the intersection is
   the schema-conformant value.

2. `StressFieldCoefficient` (header-only):

   ```cpp
   class StressFieldCoefficient : public mfem::VectorCoefficient
   {
   public:
       /// 6-component vector coefficient: returns (σ_xx, σ_yy, σ_zz,
       /// σ_xy, σ_yz, σ_xz) at the evaluation point.  This is the
       /// schema-v1 canonical order
       /// (`data_projection_schema_v1.md` §5) — **NOT** standard Voigt
       /// notation, which uses (σ_xx, σ_yy, σ_zz, σ_yz, σ_xz, σ_xy).
       /// The order is the one expected by the static-equilibration
       /// plan §1 (sigma_csm VectorCoefficient).
       StressFieldCoefficient(const StressField3D& field,
                              real_t scale = 1.0,
                              real_t offset = 0.0)
           : VectorCoefficient(6), field_(field), scale_(scale),
             offset_(offset) {}

       void Eval(mfem::Vector& v, mfem::ElementTransformation& T,
                 const mfem::IntegrationPoint& ip) override;
   };
   ```
   `Eval` calls `T.Transform(ip, x)` and then six
   `DataField3D::Evaluate` calls in the schema-v1 canonical order
   `(σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz)`, applying the shared
   `scale` / `offset`.

3. `FieldProjector::ProjectStress` — analog of `ProjectVelocity`:

   ```cpp
   struct StressFields
   {
       std::shared_ptr<mfem::ParGridFunction> sigma_xx;
       std::shared_ptr<mfem::ParGridFunction> sigma_yy;
       std::shared_ptr<mfem::ParGridFunction> sigma_zz;
       std::shared_ptr<mfem::ParGridFunction> sigma_xy;
       std::shared_ptr<mfem::ParGridFunction> sigma_yz;
       std::shared_ptr<mfem::ParGridFunction> sigma_xz;
       real_t min_sigma_xx, max_sigma_xx, ..., max_sigma_xz;   // six pairs
   };

   static StressFields ProjectStress(
       const std::string&           sidecar_path,
       mfem::ParFiniteElementSpace& target_fes,
       InterpMode interp = InterpMode::Trilinear);
   ```

   Behaviour: six independent `Project(...)` calls (one per
   component), each on the same `target_fes`. Each call increments
   `call_count_` exactly once (so the one-time-load contract is six,
   not one — but it is a single sidecar load thanks to `StressField3D`
   pre-loading every component at construction; we project the
   already-in-memory fields rather than re-reading HDF5 six times).
   The min/max bounds are forwarded from the sidecar's per-field
   declarations.

### Phase 6.A — `FaultGeometry` per-DOF 3-D coordinate and basis accessors

`ProjectFaultPreStress` (Phase 6 §4) needs three things per fault
DOF that the existing `FaultGeometry` class does not currently
expose: (a) the 3-D global UTM coordinate `(x_i, y_i, z_i)` to
query `StressField3D::Evaluate`, (b) the unit normal `n_i`, and
(c) the unit tangents `t1_i` (dip) / `t2_i` (strike). The
existing `coords_x2_` / `coords_x3_` accessors hold 2-D BP5-frame
along-strike / depth scalars (`fault_geometry.hpp:715-716`), which
do not generalise to a curved SAF, and `fault_basis_q_` exists per
quadrature point, not per DOF.

Phase 6.A adds the missing per-DOF accessors **before** Phase 6 §4
can be implemented.

#### Files to Modify
- `fault/fault_geometry.hpp` — add private members + initialiser +
  public const accessors.

#### Detailed Requirements

1. Add private members to `FaultGeometry`:

   ```cpp
   Vector       dof_coords_3d_;   // 3 * num_fault_dofs_; (x_i, y_i, z_i) interleaved
   DenseMatrix  dof_basis_;       // (9, num_fault_dofs_); col i = [n_i; t1_i; t2_i]
                                  // unit, orthonormal (Gram–Schmidt re-orthonormalised)
   ```

2. Initialiser (called from the same constructor site that fills
   `coords_x2_` / `coords_x3_`):
   - For each fault DOF `i`, fetch the parent face's element
     transformation `T_F` and the corresponding reference-point
     coordinate `ip` from the fault DG space. Evaluate
     `T_F.Transform(ip, x3d)` to get `(x_i, y_i, z_i)`; write into
     `dof_coords_3d_` at offsets `3*i .. 3*i+2`.
   - Evaluate the fault-face basis at each DOF's reference
     coordinate directly: invoke the `FaultBasis` evaluator (the
     same one that fills `fault_basis_q_` at quadrature points) at
     the per-DOF reference coordinate `ip_dof` on the parent face.
     The basis is per-face and discontinuous across face
     boundaries; this is correct and matches the per-face storage
     of `tau_pre_`. **Do NOT** project the quadrature-point values
     via H1 shape functions — the fault basis is not in H1.
   - Re-orthonormalise the projected basis per DOF (Gram–Schmidt,
     mirroring Phase 2 `basis_to_node`):
     `n_i ← n_i / |n_i|`,
     `t1_i ← (t1_i − (t1_i·n_i) n_i)`, then `t1_i ← t1_i / |t1_i|`,
     `t2_i ← n_i × t1_i` (already unit).
   - Store columns `[n_i; t1_i; t2_i]` of `dof_basis_` at column `i`.

3. Public const accessors:

   ```cpp
   const Vector&      fault_dof_coords_3d() const { return dof_coords_3d_; }
   const DenseMatrix& fault_dof_basis()     const { return dof_basis_;     }
   ```

#### Edge Cases to Handle
- Curved fault with sub-vertical triangles producing a near-degenerate
  projected basis: when `|t1_i_projected − (t1_i·n_i)·n_i| < 1e-12`,
  emit an MFEM_WARNING and fall back to deriving `t1_i` from the
  reference up-vector `(0, 0, 1)` via `t1 = n × (up × n)`, then
  re-normalise. The fallback site count is reported to stdout at
  init.

#### Acceptance Criteria
- [ ] `geom.fault_dof_coords_3d().Size() == 3 * geom.num_fault_dofs()`.
- [ ] `geom.fault_dof_basis().Height() == 9`,
      `geom.fault_dof_basis().Width() == geom.num_fault_dofs()`.
- [ ] For every DOF `i` in the SAFS 2 km test mesh:
      `|n_i| = |t1_i| = |t2_i| = 1` to within 1e-12;
      `n_i · t1_i = n_i · t2_i = t1_i · t2_i = 0` to within 1e-12.
      A new unit test `test_fault_dof_basis_orthonormal` asserts
      this on the real 2 km mesh.
- [ ] Existing BP5 verification (`make test-bp5-integration`) passes
      bit-exact: the new members are write-once at init and never
      read by BP5/BP2/TPV102 code paths.

#### Dependencies
- Depends on: nothing in this plan (extends existing `FaultGeometry`).
- Required by: Phase 6 §4 (`ProjectFaultPreStress`) and Phase 6 §5
  (`ComputeSAFSParams`).

---

4. `FieldProjector::ProjectFaultPreStress` — pure rotation onto the
   per-DOF fault basis. **No sign flip lives here.** The sidecar
   already carries compression-positive Pa (Phase 3 / Phase 5 did
   the flip once at the source); this function only rotates and
   subtracts pore pressure:

   ```cpp
   /// Project a sidecar stress tensor onto the per-DOF fault basis
   /// to produce SEAS-internal pre-stress slots.
   ///
   /// Inputs from the sidecar are already in SEAS internal sign
   /// convention (compression positive, Pa).  No sign flip is
   /// performed in this function; the rotation is purely linear
   /// algebra:
   ///
   ///   sigma_n_per_dof(i)       = n_i^T σ(x_i) n_i
   ///                              - (P_p_pa + P_p_grad * max(0, -z_i))
   ///   tau_pre_per_dof(2*i)     = t1_i^T σ(x_i) n_i   (dip,    BP5 t1)
   ///   tau_pre_per_dof(2*i+1)   = t2_i^T σ(x_i) n_i   (strike, BP5 t2)
   ///
   /// where (t1, t2, n) is the canonical SEAS fault-local frame
   /// (CLAUDE.md "fault-local tangent frame" rule: t1 = dip, t2 = strike).
   ///
   /// The sign convention contract is enforced upstream via a
   /// single-site flip (R-501/R-502):
   ///   - Phase 0 dump file carries compression-negative σ⁰ (H&Z).
   ///   - Phase 3 `bulk_stress_tensor_field` performs the SINGLE
   ///     full-tensor flip (σ_seas = −σ_HZ) at the source and
   ///     returns σ_seas in compression-positive SEAS convention.
   ///   - Phase 5 sidecar writer is a pass-through (MPa → Pa unit
   ///     conversion only; no further sign manipulation).
   ///   - Phase 6 `StressField3D` is a pass-through reader.
   ///   - This function receives σ_seas (compression-positive) and
   ///     rotates it onto the per-DOF Tandem fault basis (Phase 6.A).
   ///     Under Tandem basis with σ_seas, all three rotated
   ///     components (σ_n, τ_strike, τ_dip) carry the SEAS-internal
   ///     sign: σ_n > 0 = compression; τ_strike > 0 = right-lateral;
   ///     τ_dip > 0 = reverse. No additional sign manipulation here.
   ///   - QD/dyn consumers expect σ_n > 0 = compression and
   ///     positive τ_strike for right-lateral slip; read this output
   ///     directly into tau_pre_ / sigma_n_per_dof / DOFData slots.
   ///
   /// If a future caller introduces an inverted σ, the
   /// `test_field_projector_project_fault_prestress_basic` parity test
   /// (matches the Phase 4 JSON medians) detects it.
   static void ProjectFaultPreStress(
       const StressField3D&             field,
       const FaultGeometry&             fault_geom,
       mfem::Vector&                    sigma_n_per_dof,   // (num_fault_dofs,)
       mfem::Vector&                    tau_pre_per_dof,   // (2 * num_fault_dofs,)
       real_t                           P_p_pa = 0.0,
       real_t                           P_p_grad_pa_per_m = 0.0);
   ```

   The implementation loops over `i = 0 .. num_fault_dofs - 1`, fetches
   `(x_i, y_i, z_i)` from `FaultGeometry::fault_dof_coords_3d()` and
   `(n_i, t1_i, t2_i)` from `FaultGeometry::fault_dof_basis()` —
   both surfaced by Phase 6.A. It then evaluates the bulk tensor
   `S = field.Evaluate(x_i, y_i, z_i)` (compression-positive,
   sidecar convention) and computes the three components above.
   Pore pressure is **subtracted** from the normal stress
   (compression-positive convention: P_p reduces the effective
   normal load): `sigma_n_per_dof(i) = nᵀSn - (P_p_pa + P_p_grad *
   max(0, -z_i))`.

5. `FaultGeometry::ComputeSAFSParams` — parallel slot to
   `ComputeBP5Params`:

   ```cpp
   /// SAFS-mode pre-stress initialisation: same outputs as
   /// ComputeBP5Params (a_values_, eta_values_, dc_values_, tau_pre_,
   /// V_init_vec_) but with `tau_pre_` and the new `sigma_n_per_dof_`
   /// sourced from a sidecar.  `a`, `Dc`, `eta`, and `V_init_vec_` keep
   /// their existing analytic spatial dependence (BP5 a(x2, x3),
   /// Dc(x2, x3), etc.).  Note: SAFS-specific a/Dc spatial functions
   /// are out of scope here; the user can plug in the production
   /// SAFS-params analytic forms when ready.
   void ComputeSAFSParams(const StressField3D& field,
                          real_t P_p_pa = 0.0,
                          real_t P_p_grad_pa_per_m = 0.0);
   ```

   New private member: `Vector sigma_n_per_dof_;` of size
   `num_fault_dofs_`. New public const accessor:
   `const Vector& sigma_n_per_dof() const { return sigma_n_per_dof_; }`.

6. `RateStateFaultOperator` consumer side (vector path only — SAFS
   is a 3-D fault, `SlipComponents == 2`, so the BP2 scalar branch
   at `rate_state_fault.hpp` lines 267 / 390 / 806 / 872 is **NEVER
   taken** under SAFS instantiation and must **NOT** be modified):
   - Add private members
     `const Vector* tau_pre_per_dof_ = nullptr;`,
     `const Vector* sigma_n_per_dof_ = nullptr;`, and
     `bool safs_mode_ = false;`.
   - In the BP5 vector `Init` path (`rate_state_fault.hpp:303-304`,
     where `tau_pre_(2*i)` / `tau_pre_(2*i+1)` are read), replace
     the read with `(*tau_pre_per_dof_)(2*i)` / `(*tau_pre_per_dof_)(2*i+1)`
     under `safs_mode_`. NB `tau_pre_` is *already* the per-DOF
     vector for BP5; the only change is whether it is sourced from
     `bp5_params_.tau0_vec` or from `StressField3D` via
     `ComputeSAFSParams`. The same applies to the BP5 vector
     `ComputeRHS` path (`rate_state_fault.hpp:825-826`).
   - Replace the scalar `sigma_n_bp5_` at the
     `rate_state_fault.hpp:831, 838` `SolveSlipRateVectorPsi` call
     sites with
     `safs_mode_ ? (*sigma_n_per_dof_)(i) : sigma_n_bp5_` so the
     heterogeneous normal stress reaches the friction solver. The
     same substitution applies to the equilibrium-verification path
     around `rate_state_fault.hpp:856` and to any other call site
     whose audit (per CLAUDE.md "Files Requiring Extreme Care")
     shows it reads `params_.sigma_n` directly. Enumerate every
     site in the Phase 6 commit message; do not leave any
     unaudited (per the user-feedback memory
     `feedback_complete_sign_sites.md`).
   - The BP5 scalar `Init` (line 255 `tau0_ = params_.tau0()`) and
     BP2 init (line 691) remain **untouched** — both are scalar
     paths not exercised by SAFS (`SlipComponents == 1`).
   - Audit `friction/dieterich_ruina.hpp` for any direct reads of
     `params_.sigma_n`; the SAFS-mode dispatch must thread the
     per-DOF value through the same call chain to the friction
     solver. Enumerate the touched sites in the commit message.

7. TOML schema extension. Add to `config/seas_config.hpp`:

   ```cpp
   struct StressConfig
   {
       bool use_sidecar = false;
       std::string sidecar_path;
       real_t P_p_pa = 0.0;
       real_t P_p_grad_pa_per_m = 0.0;
   };
   ```
   Add `StressConfig stress;` to `SeasConfig`. Parser updates in
   `seas_config_parser.hpp`. Bridge updates in `seas_config_bridge.hpp`.
   Driver: when `stress.use_sidecar`, instantiate `StressField3D`,
   call `FaultGeometry::ComputeSAFSParams`, mark the operator
   `safs_mode_ = true`, attach the per-DOF vectors.

### Edge Cases to Handle
- Sidecar mesh bbox does not contain the SAFS fault DOFs: the
  per-DOF projector aborts via `DataField3D::Evaluate` (interpolation-
  only contract). This is the desired behaviour; the Phase 5 grid
  padding is what prevents it from triggering in practice.
- `use_sidecar = true` but no sidecar path: parser aborts with a
  clear message.
- `use_sidecar = false`: every existing BP5 / BP2 / TPV102 / TPV205
  code path remains untouched. The new behaviour is **strictly
  opt-in**.
- Per-DOF normal stress is negative or zero (rare for SAF geometry
  but possible in the upper few hundred metres where σ_v is small):
  emit a warning the first time and clamp via `sigma_n_per_dof(i) =
  max(sigma_n_per_dof(i), 1.0e6)` (1 MPa floor) only when the value
  is < 1 MPa. **This clamp is the only numeric default**; document
  it in the docstring of `ProjectFaultPreStress`. Make the clamp
  configurable via a `min_sigma_n_pa` argument that defaults to 0 (no
  clamp) so the user must opt in.

### Acceptance Criteria
- [ ] `StressField3D("stress_safs.h5")` constructs successfully against
      the Phase-5 output and exposes a `(3, 3)` symmetric tensor at
      a sample interior point that matches the H&Z `demo_safod` σ⁰
      values to 1e-3 Pa.
- [ ] `FieldProjector::ProjectStress` returns six `ParGridFunction`s
      whose pointwise values at a sample of cell centroids match the
      sidecar values to within the trilinear interpolation noise
      (≤ 1 Pa for a constant field).
- [ ] `FieldProjector::ProjectFaultPreStress` produces
      `sigma_n_per_dof` whose median equals
      `-median(sigma_n_total)` from the Phase 4 JSON summary
      (sign-flipped), to within 1e-3 MPa.
- [ ] `FieldProjector::ProjectFaultPreStress` produces
      `tau_pre_per_dof(2*i+1)` (strike) whose median matches the
      `tau_strike` median from the Phase 4 summary to within 1e-3 MPa.
- [ ] Running the existing BP5 verification driver with
      `--config bp5_example.toml` (which has `stress.use_sidecar =
      false` by construction) produces bit-exact output relative to
      the pre-Phase-6 baseline. (The `safs_mode_` guard must compile
      out cleanly when not in use.)
- [ ] Running a new minimal SAFS TOML
      (`drivers/safs_smoke.toml`, added in Phase 6) with
      `stress.use_sidecar = true` and a homogeneous H&Z σ⁰ sidecar
      reaches t = 1 yr without nucleating an event (sanity smoke;
      should be stable under the equilibrium assumption).

### Dependencies
- Depends on: Phase 5 (sidecar exists).
- Required by: Phase 7.

---

## Phase 7 — `seas_project_stress_to_mesh` verification driver

### Goal
After Phase 7, a user can visually verify the stress projection in
ParaView before launching a long SAFS simulation, exactly as
`project_velocity_to_mesh` already does for velocity. The driver
loads a mesh + a stress sidecar, projects all six tensor components
onto an H1(p) space, and writes a `.pvd` + per-rank `.vtu` set.

### Files to Create
- `drivers/project_stress_to_mesh.cpp` — copy `project_velocity_to_mesh.cpp`
  and swap `ProjectVelocity` for `ProjectStress`; six PointData
  fields instead of five.
- `tests/unit/test_project_stress_to_mesh.cpp` — drives the driver
  via subprocess on a tiny synthetic sidecar + mesh; verifies the
  pvd / vtu files exist and round-trip-load.

### Files to Modify
- `Makefile` — add a `seas_project_stress_to_mesh` target mirroring
  the existing `seas_project_velocity_to_mesh` target.
- `tests/Makefile` (if present) — add the new unit test to the test
  suite.

### Detailed Requirements

1. CLI matches `project_velocity_to_mesh.cpp:46-95`:
   ```
   seas_project_stress_to_mesh
       --mesh    safs_fault_box_nwcut_2000m.msh
       --sidecar stress_safs.h5
       --out     projected_stress_2000m
       [--order P]
       [--interp {trilinear,catmull-rom}]
       [--ascii | --binary]
   ```
   Defaults match the velocity driver: `--order 2`, `--interp
   trilinear`, `--binary`.

2. Output: `<out>/<out>.pvd` referencing per-rank `.vtu` files with
   six PointData scalar fields `{sigma_xx, sigma_yy, sigma_zz,
   sigma_xy, sigma_yz, sigma_xz}` in Pa.

3. Pre-flight: mesh bbox vs sidecar bbox via
   `StressField3D::ContainsBBox`. Abort cleanly on failure
   with the same message format as `FieldProjector::Project`.

### Edge Cases to Handle
- Sidecar missing fields: surfaced by the `StressField3D` ctor
  (each `DataField3D` re-checks).
- `--out` path collides with an existing file: overwrite (matches
  `project_velocity_to_mesh.cpp`).

### Acceptance Criteria
- [ ] `seas_project_stress_to_mesh --mesh ... --sidecar
      stress_safs.h5 --out projected_stress_2000m` produces a `.pvd`
      readable in ParaView with six scalar PointData fields, each
      visualisable as a colour map.
- [ ] Round-trip test: load the resulting VTU back via `meshio` (the
      Phase 7 unit test does this when `--ascii`), and assert the
      per-field mean equals the sidecar's per-field mean to 1 part
      in 1e-9.
- [ ] The driver runs in under 10 s on the 2000 m mesh, single rank.

### Dependencies
- Depends on: Phase 6.
- Required by: Phase 8.

---

## Phase 8 — VTU postprocess verifier (analytical-prediction check)

### Goal
After Phase 8 the user can run a single script that reads every
written VTU artefact (fault point-data, fault cell-data, bulk
cell-data, and the Phase-7 H1-projected VTUs from
`seas_project_stress_to_mesh`) and verifies the recorded values
against the **analytic prediction** evaluated independently at the
same coordinates. The script returns a non-zero exit code if any
field exceeds the per-field error tolerance; this catches silent
corruption in the writers (cell→node averaging, meshio round-trip,
ParaView dialect bugs) without re-running the full simulation.

### Files to Create
- `code_preprocess/verify_onfault_stress.py` — verifier script + CLI.
- `code_preprocess/test_verify_onfault_stress.py` — pytest suite.

### Files to Modify
- None.

### Detailed Requirements

1. Module header following the existing pipeline style. Imports:
   `numpy`, `meshio`, `json`, `argparse`, plus a sibling import of
   `project_to_fault_stress` for the analytic evaluator
   (`project_stress_onto_fault`, `cell_to_node_average`, `basis_to_node`).

2. Module-level constants (per `feedback_no_hardcoded_numbers.md`,
   tolerances are configurable via CLI flags whose **defaults** are
   defined here once and traceable to the floating-point round-off
   limits of the affected math):
   - `DEFAULT_TOL_PA_CELL: float = 1.0`         # raw cell-data: rotation error only
   - `DEFAULT_TOL_PA_NODE: float = 1.0e3`       # node-averaged: piecewise-linear interpolation error of a smooth field
   - `DEFAULT_TOL_PA_BULK: float = 1.0`         # constant σ⁰ on bulk: exact
   - `DEFAULT_TOL_PA_H1_PROJ: float = 1.0e3`    # H1(p) projection: tied to mesh edge length
   - `DEFAULT_REL_TOL: float = 1.0e-6`          # for fields whose magnitude is well below the absolute tolerance

3. Add the dataclass:

   ```python
   @dataclass
   class VerificationResult:
       field_name: str
       l_inf_err: float          # max |observed - predicted|
       l_inf_loc: np.ndarray     # (3,) the (x, y, z) where it was reached
       l2_err: float             # sqrt(mean(diff²))
       observed_min: float
       observed_max: float
       passed: bool
       tol_abs: float
       tol_rel: float
   ```

4. Core verifier functions:

   ```python
   def verify_fault_point_data(
       fault_vtu_path: Path,
       params: dict,                       # echoed Phase 4 JSON "params" block
       *,
       tol_abs_pa: float = DEFAULT_TOL_PA_NODE,
       tol_rel: float = DEFAULT_REL_TOL,
   ) -> list[VerificationResult]

   def verify_fault_cell_data(
       fault_vtu_path: Path,
       params: dict,
       *,
       tol_abs_pa: float = DEFAULT_TOL_PA_CELL,
       tol_rel: float = DEFAULT_REL_TOL,
   ) -> list[VerificationResult]

   def verify_bulk_cell_data(
       bulk_vtu_path: Path,
       params: dict,
       *,
       tol_abs_pa: float = DEFAULT_TOL_PA_BULK,
       tol_rel: float = DEFAULT_REL_TOL,
   ) -> list[VerificationResult]

   def verify_h1_projected_vtu(
       projected_pvd_path: Path,           # output of seas_project_stress_to_mesh
       params: dict,
       *,
       tol_abs_pa: float = DEFAULT_TOL_PA_H1_PROJ,
       tol_rel: float = DEFAULT_REL_TOL,
   ) -> list[VerificationResult]
   ```

   Each function:
   1. Loads the VTU via `meshio.read`.
   2. Reads the recorded field arrays (point-data or cell-data per the
      schema written in Phase 4 / Phase 7).
   3. Evaluates the analytic prediction at the same coordinates by
      calling `bulk_stress_tensor_field` (Phase 3) for the bulk and
      `project_stress_onto_fault` for the fault, using `params`.
   4. Computes L∞ and L2 errors; flags `passed = (l_inf_err ≤
      max(tol_abs_pa, tol_rel * max(|observed|)))`.

5. Verifier-specific math — the analytic prediction at point `p`,
   under the unified single-site flip contract (R-501/R-502):

   For the **bulk** point `p = (x, y, z)`:
   ```
       σ_seas(p) = -build_bulk_stress_tensor(SHmax(z), Shmin(z),
                                             Sv(z), SHmax_az)
   ```
   (i.e., compression-positive SEAS convention, matching Phase 3 §2's
   `bulk_stress_tensor_field` output.) Direct component lookup
   against `sigma_xx_Pa` etc.

   For the **fault point** at vertex `v`:
   - Re-evaluate the node-averaged basis `(s_v, d_v, n_v)` from the
     fault mesh via Phase 2's `basis_to_node` (Tandem convention).
   - Predict directly on σ_seas (NO sign-flip asymmetry — the same
     σ_seas drives all three rotated components):
     ```
       σ_n_v        =  n_v · σ_seas(p_v) · n_v
       τ_strike_v   =  s_v · σ_seas(p_v) · n_v
       τ_dip_v      =  d_v · σ_seas(p_v) · n_v
     ```
     Matches Phase 3 §1's emission formulas exactly; the verifier
     and the writer use the same σ_seas + Tandem-basis rotation.
   - Compare to the recorded `sigma_n_total_MPa`, `tau_strike_MPa`,
     `tau_dip_MPa` at vertex `v`.
   - Tolerance: cell→node averaging is exact for *linear* fields on
     a triangulated surface — for a constant σ⁰ background the
     verifier should match to machine precision (`tol_abs_pa = 1
     Pa`). For a depth-varying σ_v the averaging is O(h · grad),
     captured by the `DEFAULT_TOL_PA_NODE` default.

   For the **fault cell** at triangle `t`:
   - Same unified-σ_seas predictor evaluated at the triangle
     centroid using the per-cell basis (no averaging). Tolerance
     should be at the rotation round-off limit (`tol_abs_pa = 1 Pa`).

6. CLI:
   ```
   python verify_onfault_stress.py
       --fault-vtu <path>
       [--bulk-vtu <path>]
       [--h1-projected-pvd <path>]
       [--summary-json <path>]    # canonical: read params from here
       [--tol-pa <float>]         # override all absolute tolerances
       [--tol-rel <float>]
       [--fail-on-warn]           # exit 1 if any warning is emitted
   ```

   The script prints a one-table-per-field summary
   (L∞, L2, observed range, tol, pass/fail) and exits 0 if all checks
   pass, 1 otherwise. The pass/fail tally is also written to a
   `verify_report.json` next to the inputs.

7. The verifier must be **resilient** to:
   - VTU files written with either point-data or cell-data alone
     (skip the missing layer with an info-level log).
   - Symmetric tensor stored as either 9-component flat or six
     `sigma_*` separate fields (Phase 4 uses the second; Phase 6
     C++ output may use either).
   - Stem-mismatched paths (mesh resolutions can be inferred from the
     stem; the verifier never re-runs the projection).

### Interfaces

The public API is `verify_fault_point_data`, `verify_fault_cell_data`,
`verify_bulk_cell_data`, `verify_h1_projected_vtu`, `VerificationResult`,
and the CLI.

### Edge Cases to Handle
- Tolerance budget for the H1-projected VTU. The error model is
  `||σ - P_h σ||_L²(K) ≤ C · h_K^{p+1} · |σ|_{p+1}` for an order-`p`
  H1 projection on element `K`. For constant σ⁰ this is zero
  analytically; the verifier should pass at `tol_abs_pa = 1`. For a
  depth-varying σ⁰ the user must pass `--tol-pa` large enough to
  cover the per-element gradient times `h`. Document this in the
  CLI docstring.
- Verifier on Phase 6 / Phase 7 output before the full pipeline is
  green: the `--fault-vtu` flag can be omitted for a bulk-only
  check; the script must not require all four artefacts to be
  present.
- Float-precision drift between Python (float64) and MFEM C++
  (real_t = float64 by default but configurable). Document that the
  `tol_abs_pa = 1 Pa` default is for the real_t = float64 build;
  bump to 1e3 Pa for the real_t = float build.

### Acceptance Criteria
- [ ] Running
      ```
      python verify_onfault_stress.py
          --fault-vtu data_projection_onfaultstress/2000m/safs_fault_box_nwcut_2000m_fault_stress.vtu
          --bulk-vtu  data_projection_onfaultstress/2000m/safs_fault_box_nwcut_2000m_bulk_stress.vtu
          --summary-json data_projection_onfaultstress/2000m/safs_fault_box_nwcut_2000m_summary.json
      ```
      exits 0 on a freshly-written Phase 4 artefact set with the
      default SAFOD parameters and reports L∞ < 1 Pa for the bulk
      cell-data and L∞ < tol for the fault point-data.
- [ ] A deliberately-corrupted VTU (one cell's `sigma_xx_Pa`
      perturbed by 1 GPa) is detected: the script exits 1 and prints
      a row identifying the offending field.
- [ ] The Phase-7 H1-projected VTU verifies with L∞ ≤ machine
      precision for the constant-σ⁰ case.
- [ ] On the lithostatic-Sv case the verifier passes with
      `tol_abs_pa = 1e3` (and fails with `tol_abs_pa = 1e-3`), which
      gives a sanity-check on the documented error model.

### Dependencies
- Depends on: Phases 1, 2, 3, 4 (Python pipeline); Phase 7 (when
  the H1-projected VTU mode is exercised).
- Required by: nothing.

---

## Testing Strategy

### Test file
`code_preprocess/test_project_to_fault_stress.py`, mirroring
`test_nw_cut_strip.py` (lines 1–40): pytest layout, `sys.path` insert
to enable sibling imports, `TestPhase1` / `TestPhase2` / ... class
groupings, R-numbered regression tests for any bugfix that lands
during implementation.

### Per-phase tests

**Phase 0 (H&Z dump entry point)**
- `test_hz_dump_demo_safod`: invoke `dump_safod_sigma0(tmp_path / "safod.json")`;
  reload via `json.load`; assert
  `data["params"]["SHmax_MPa"] == 113.0`, `data["sigma0_MPa"]` is a
  3×3 array matching `compute_fault_stress(*demo_safod_args).sigma0`
  to 1e-9 MPa, `data["convention"]` contains `"compression negative"`,
  and `data["schema"] == "hickman_zoback_sigma0_v1"`.
- `test_hz_dump_cli`: subprocess-invoke
  `python hickman_and_zoback_regional_stress_projection.py --dump
  /tmp/safod.json`; assert the file is written and JSON-valid.
- `test_hz_no_arg_still_runs_demo`: subprocess-invoke with no args
  and verify `demo_safod()` output is on stdout (regex match
  `mu_apparent\s*=\s*0\.24`). Defends the "no other behaviour
  changes" rule.

**Phase 1 (mesh I/O)**
- `test_load_fault_vtu_real_data`: loads the real 2000 m fault VTU;
  asserts triangle count = 2685, bbox matches probed values.
- `test_triangle_geometry_unit`: synthetic 3-4-5 right triangle in
  the XY plane; asserts area = 6.0, normal = ±ẑ, centroid = (mean
  of vertices). Run for 100 random orientations: |normal| = 1 to 1e-15.
- `test_tet_geometry_unit`: synthetic reference tet with vertices at
  the standard simplex; asserts volume = 1/6, centroid = (0.25, 0.25, 0.25).

**Phase 2 (fault basis)**
- `test_basis_matches_HZ_for_vertical_NW_strike`: synthetic vertical
  triangle with normal along H&Z's `n_hat(strike=314, dip=90)`; assert
  `per_triangle_basis_raw` returns `s_hat`, `d_hat` matching H&Z to
  1e-12.
- `test_basis_orthonormal`: 1000 random triangles, no degenerate cases;
  `dot(s, d), dot(s, n), dot(d, n)` all < 1e-12 and norms = 1 ± 1e-12.
- `test_horizontal_triangle_marks_degenerate`: triangle with normal
  along ẑ; assert `degen[0] == True` and `strikes[0]` is NaN.
- `test_harmonise_idempotent`: random normals; `harmonise(harmonise(x))
  == harmonise(x)`.
- `test_harmonise_flips_to_majority`: half normals flipped a priori;
  after harmonisation all dot products with reference > 0.
- `test_harmonise_pca_fallback`: pass `fault_strike_azimuth_hint=None`;
  on planar synthetic mesh, the resulting normal direction agrees with
  the analytic plane normal to within 1e-12.
- `test_cell_to_node_average_constant_scalar`: cell field is
  uniform `f = 3.14`; assert every vertex gets `f = 3.14` exactly
  (constant case: averaging is identity).
- `test_cell_to_node_average_linear_scalar`: synthetic triangulated
  square with `f_cell = x_centroid` (a linear function of position);
  assert each interior vertex receives the vertex-x to within
  1e-12 (linear-on-triangulation case: area-weighted average of
  centroid values equals the vertex value within float64 noise).
  This is the precise test of "continuous field on a planar mesh".
- `test_cell_to_node_average_vector`: stack three different scalar
  fields as a `(N_tri, 3)` vector cell-data input; assert
  componentwise that the result matches three independent scalar
  averages.
- `test_cell_to_node_average_orphan_vertex_is_nan`: synthetic mesh
  with one orphan vertex (no incident triangles); assert that
  vertex's node value is NaN.
- `test_cell_to_node_no_python_loop`: monkeypatch a slow numpy
  `np.add.at` to fail the test if the helper does more than O(1)
  numpy calls per cell — the implementation must vectorise
  (per Phase 2 §4 requirement). Or, more simply, run on a 10⁵-cell
  mesh and assert wall time < 0.5 s.
- `test_basis_to_node_unit_and_orthonormal`: random mesh of 100
  vertical triangles; assert `|n_node|, |s_node|, |d_node| = 1`
  to 1e-12 and pairwise dot products < 1e-12.
- `test_basis_to_node_matches_cell_on_uniform_fault`: all triangles
  share the same basis (perfectly planar fault); assert every
  vertex basis matches the shared cell basis exactly.

**Phase 3 (stress rotation, single sign-flip site)**
- `test_resolve_traction_against_HZ_safod`: SAFOD inputs, single
  vertical triangle aligned with H&Z's `s_hat`/`d_hat`. Assert
  `sigma_n_eff ≈ +89 MPa` (**compression positive**; flipped from
  H&Z's `-89 MPa`), `mu_apparent ≈ 0.24` (unchanged because μ takes
  magnitudes), both to 1e-3 absolute. This is the primary
  cross-check against the H&Z scalar path and serves as the
  *numerical regression anchor* for the whole module — and the
  guardrail for the one-shot sign flip in
  `bulk_stress_tensor_field`.
- `test_bulk_stress_tensor_field_sign_flip`: call
  `bulk_stress_tensor_field` and `build_bulk_stress_tensor`
  separately with identical args; assert `bulk_stress_tensor_field`
  output equals `-build_bulk_stress_tensor` to 1e-12. **This is the
  guard against an accidental missed flip.**
- `test_bulk_stress_tensor_field_symmetric`: 100 random parameter
  sets; assert each output `(3, 3)` tensor is symmetric to 1e-9.
- `test_vectorised_equals_scalar`: 100 random triangles, random
  symmetric σ⁰; vectorised result vs. 100 calls to H&Z `resolve_traction`,
  element-wise diff < 1e-12.
- `test_trace_invariance`: random `(s, d, n)` frames and random
  symmetric σ; assert `σ_n + sᵀσs + dᵀσd ≈ tr(σ)` to within 1e-9.
- `test_no_shear_on_principal_face`: σ⁰ diagonal with all three
  distinct eigenvalues, normal aligned with one eigenvector; assert
  `tau_strike ≈ tau_dip ≈ 0` (< 1e-12 MPa).
- `test_depth_dependent_lithostatic_sv`: 5 sample depths, non-zero
  gradient; assert σ_zz changes linearly with depth and σ_xx / σ_yy
  follow rotation about ẑ exactly.
- `test_nonsymmetric_sigma_raises`: deliberately non-symmetric input;
  expect `ValueError`.

**Phase 4 (writers + CLI, continuous fault field)**
- `test_write_fault_vtu_roundtrip`: write a synthetic 5-triangle
  fault VTU; reload with `meshio.read`; assert every expected
  point-data **and** cell-data field is present with the right
  shape and dtype.
- `test_fault_vtu_point_data_continuous`: synthetic 4-triangle fan
  (one interior vertex shared by all 4 triangles); the four cells
  carry distinct `sigma_n` values; assert the **interior vertex**
  point-data value equals the area-weighted average of the four
  cells exactly (1e-12). The four neighbouring boundary vertices'
  values are predictable from the same formula. **This is the
  primary continuous-field test.**
- `test_fault_vtu_rake_recomputed_at_node`: cell `rake_deg` values
  are {±175°, ∓179°} (deliberately straddling the periodic
  discontinuity); assert the point-data `rake_deg` at the shared
  vertex is recomputed from `atan2(τ_d_node, τ_s_node)`, **not**
  averaged from the cell-data angles. Detects the periodic-average
  bug.
- `test_fault_vtu_basis_orthonormal_at_node`: load the written
  fault VTU; assert at every (non-degenerate) vertex:
  - `|np.dot(strike_vec_node, normal_vec_node)| < 1e-9`
  - `|np.dot(dip_vec_node,    normal_vec_node)| < 1e-9`
  - `|np.dot(strike_vec_node, dip_vec_node)|    < 1e-9`
  - `|strike_vec_node| = 1 ± 1e-9`
  - `|dip_vec_node|    = 1 ± 1e-9`
  - `|normal_vec_node| = 1 ± 1e-9`
  All three orthogonality relations and all three unit-magnitude
  conditions must hold; partial coverage masks Gram-Schmidt
  regressions (e.g. the R-102 dip-sign bug, where a wrong-cross
  recipe would still produce a vector orthogonal to the normal
  but with the wrong sign).
- `test_fault_vtu_compression_positive`: for SAFOD parameters,
  assert every `sigma_n_total_MPa` value (point-data AND cell-data)
  is **strictly positive** — defends the sign convention contract.
- `test_summary_json_schema`: write a summary; load with
  `json.load`; assert all top-level keys present and `stats` has
  entries for every output field.
- `test_cli_print_info_no_write`: invoke `--print-info` via subprocess
  on the real 2000 m fault VTU; assert no new files appear in
  `data_projection_onfaultstress/`.
- `test_cli_single_file_smoke`: invoke the script on the real
  2000 m fault VTU into a tempdir; assert all three expected files
  appear and pass `meshio.read` / `json.load`. **This test must run
  in the pythonenv conda env**; gate with
  `pytest.importorskip("meshio")`.
- `test_vtu_units_metadata_carried`: write a fault VTU, reload, assert
  the `units` cell-data attribute is `"MPa"` for every stress field
  (or, if meshio drops it, assert the JSON summary's `units` key
  exists and equals `"MPa"`).
- `test_per_field_minmax_match_summary`: stats reported in the JSON
  summary must match `np.min` / `np.max` / `np.median` of the
  corresponding VTU cell-data array exactly (no rounding drift).
- `R-001`-style regression slots reserved for issues caught during
  implementation review.

**Phase 5 (sidecar writer)**
- `test_build_uniform_utm_grid_monotone`: random bbox, random pad,
  random dx; assert returned x/y/z are strictly monotone increasing.
- `test_build_uniform_utm_grid_contains_mesh`: pass the real 500 m
  mesh bbox; assert sidecar grid bbox encloses it by ≥ pad on every
  face.
- `test_evaluate_stress_field_constant_case`: constant depth_model;
  assert every `(Nx, Ny, Nz)` cell holds the same six values and
  those values match H&Z `build_bulk_stress_tensor` outputs to
  1e-9 Pa (after MPa → Pa conversion).
- `test_evaluate_stress_field_lithostatic_sv`: assert σ_zz changes
  monotonically with depth and σ_xx / σ_yy follow the same horizontal
  rotation at every depth slice.
- `test_derive_field_bounds_encloses_data`: every cell value lies
  strictly inside the derived `[min_value, max_value]`.
- `test_no_nan_in_emitted_fields`: every field array has no NaN
  before passing to `write_sidecar`.
- `test_sidecar_schema_v1_conformance`: write a sidecar; reopen with
  h5py; assert
  `attrs["schema_version"] == "data_projection_v1"`,
  `attrs["crs"] == "EPSG:32611"`,
  `attrs["units"] == "m"`,
  `attrs["z_positive"] == "elevation"`,
  `created_at` is ISO-8601, every field has `units`, `min_value`,
  `max_value` attributes, every dataset shape is `(Nx, Ny, Nz)`,
  every grid axis is strictly monotone increasing. (This is the
  defence-in-depth pair of the C++ reader's startup checks.)
- `test_sidecar_field_set_is_stress`: assert
  `list(f["fields"].keys()) == sorted(["sigma_xx", "sigma_yy",
  "sigma_zz", "sigma_xy", "sigma_yz", "sigma_xz"])` exactly.
- `test_sidecar_units_are_Pa`: every `/fields/<name>/units` attribute
  equals `"Pa"`.
- `test_sidecar_reloads_in_DataField3D_compatible_layout`: open the
  sidecar via the `sidecar.py:read_sidecar_attrs` helper and pass
  every required attribute check the C++ reader runs (cross-check on
  the Python side avoids a slow C++ rebuild for an attribute drift
  bug; the C++ side rechecks at load anyway).
- `test_sidecar_roundtrip_evaluator`: pick 100 random in-bbox query
  points, evaluate σ⁰ analytically from `build_bulk_stress_tensor`
  vs trilinear-interpolated from the sidecar grid; for a constant
  σ⁰ case the difference must be < 1e-9 Pa (trilinear of constant
  is exact). For a `lithostatic_sv` depth-varying case with a grid
  step of `dx_m`, the difference must be < `0.5 * dx_m * grad_pa_per_m`
  per component (linear interpolation error of a linear function is
  zero at grid points and at-most-grid-step elsewhere).
- `test_build_stress_safs_idempotent`: run `build_stress_safs` twice
  with the same args; the two output files are byte-identical except
  for the `created_at` attribute.
- `test_cli_invokes_h5py`: gate with `pytest.importorskip("h5py")`;
  invoke `build_stress_safs.py` via subprocess; verify the produced
  file passes the schema check.

**Phase 6 (C++ consumer)**

The Phase 6 unit tests live as new files under
`tests/unit/`, matching the style of `test_data_field_3d.cpp`
(custom `TEST_ASSERT` / `TEST_NEAR` macros, no gtest dependency).

- `test_stress_field_3d_load_schema_v1.cpp`: load a synthetic
  schema-v1 stress sidecar (created by the test in `setUp`);
  verify the six DataField3D members are reachable, the bbox is the
  intersection of the six per-component bboxes, and `Evaluate(x, y,
  z)` returns a symmetric `(3, 3)` tensor whose components match
  the synthetic sidecar exactly at grid corners.
- `test_stress_field_3d_evaluate_symmetric.cpp`: 1024 random
  in-bbox points; assert `S(0, 1) == S(1, 0)`, `S(0, 2) == S(2,
  0)`, `S(1, 2) == S(2, 1)` to machine precision.
- `test_stress_field_3d_ctor_aborts_on_wrong_schema.cpp`: synthetic
  sidecar with `schema_version = "data_projection_v2"`; ctor must
  abort (use the existing MFEM abort-trap pattern from
  `test_data_field_3d.cpp:257-272`).
- `test_stress_field_3d_ctor_aborts_on_missing_component.cpp`:
  synthetic sidecar with only five of the six expected fields; ctor
  must abort with a precise message naming the missing field.
- `test_stress_field_3d_ctor_aborts_on_nan.cpp`: like
  `test_data_field_3d.cpp:T-3-3` but for the stress reader.
- `test_stress_field_3d_interp_mode_lockstep.cpp`: call
  `SetInterpMode(CatmullRom)`; assert all six underlying components
  report `CatmullRom`.
- `test_field_projector_project_stress_bulk.cpp`: build a small
  serial mesh, write a synthetic constant-σ⁰ sidecar, project,
  assert every `ParGridFunction`'s pointwise value at a sample of
  DOFs equals the sidecar constant to within 1e-12 Pa.
- `test_field_projector_project_stress_call_count.cpp`: verify the
  `call_count_` advance: `ResetCallCount()` then `ProjectStress`,
  expect `CallCount() == 6` (one per component).
- `test_field_projector_project_fault_prestress_basic.cpp`: synthetic
  vertical strike-slip fault triangle; constant σ⁰ matching SAFOD;
  call `ProjectFaultPreStress`; assert the per-DOF outputs match
  the Phase 4 Python `tau_strike` / `sigma_n_eff` for the same
  triangle to 1e-3 Pa. **This is the Python↔C++ parity anchor.**
- `test_field_projector_project_fault_prestress_no_sign_flip.cpp`:
  feed in a compression-positive σ_seas (the canonical sidecar
  output); assert the emitted `sigma_n_per_dof` values are **also
  positive** (compression-positive in = compression-positive out).
  Then deliberately feed in a hand-flipped compression-negative
  tensor and assert that `sigma_n_per_dof` comes out *negative* —
  this guards that no hidden sign flip is performed inside
  `ProjectFaultPreStress`. The sign-flip contract lives upstream in
  Phase 3, not here.
- `test_field_projector_project_fault_prestress_pore_pressure.cpp`:
  call with `P_p_pa = 16e6`; assert `sigma_n_per_dof(i)` is reduced
  by 16 MPa relative to the `P_p_pa = 0` case, at every DOF.
- `test_field_projector_project_fault_prestress_pore_pressure_grad.cpp`:
  call with `P_p_grad_pa_per_m = 9.8e3 * 1e3` (hydrostatic ρ_w = 1000
  kg/m³); assert linear scaling with depth at every fault DOF.
- `test_field_projector_project_fault_prestress_min_sigma_n_clamp.cpp`:
  σ⁰ field producing negative σ_n at one DOF; call with
  `min_sigma_n_pa = 1e6`; assert that DOF receives `sigma_n_per_dof =
  1e6` and all other DOFs are unchanged.
- `test_compute_safs_params_matches_bp5_baseline.cpp`: build a synthetic
  σ⁰ sidecar whose value at every fault DOF equals the BP5 analytic
  `tau0_vec` / `sigma_n`; call `ComputeSAFSParams` and compare
  `tau_pre_` / `sigma_n_per_dof_` against the values
  `ComputeBP5Params` would produce. Match within 1e-6 Pa per
  component. **This is the BP5↔SAFS parity anchor.**
- `test_compute_safs_params_size_matches_num_dofs.cpp`:
  `tau_pre_.Size() == 2 * num_fault_dofs_` and
  `sigma_n_per_dof_.Size() == num_fault_dofs_`.
- `test_compute_safs_params_no_nan.cpp`: run on the real 2000 m mesh
  + real `stress_safs.h5`; assert no NaN in either output vector.
- `test_rate_state_fault_safs_mode_dispatch.cpp`: build a minimal
  fault operator, set `safs_mode_ = true`, plug in per-DOF vectors;
  assert that subsequent ODE residual calls use the per-DOF values
  rather than the scalar `params_.tau0()` / `params_.sigma_n`.
  (Inspect the resulting `traction` and back-solve which σ_n
  value the friction solver received.)
- `test_safs_mode_off_is_bit_exact_to_pre_phase6.cpp`: build the
  same operator with `safs_mode_ = false`; run one ODE residual;
  assert byte-exact equality to a snapshot recorded against the
  pre-Phase-6 baseline (golden vector serialised once during
  Phase 6 test authoring).
- `test_seas_config_parses_stress_section.cpp`: TOML round-trip; the
  `[stress]` section with `use_sidecar = true`, `sidecar = "..."`,
  `P_p_pa = 16e6` survives parse → bridge → operator-level state.
- `test_seas_config_default_stress_off.cpp`: TOML without a
  `[stress]` block produces `stress.use_sidecar = false`. Defends
  the "strictly opt-in" rule.
- `test_safs_smoke_runs_one_step.cpp`: TOML `safs_smoke.toml` with
  `use_sidecar = true` and a constant-σ⁰ sidecar; advance the
  operator one step; assert it does not abort and `V_max <
  V_ini * 1e6` (no immediate runaway).

**Phase 7 (verification driver)**
- `test_project_stress_to_mesh_cli_smoke.cpp`: subprocess-invoke
  `seas_project_stress_to_mesh --mesh <tiny.msh> --sidecar <synth.h5>
  --out <tmp/out> --ascii`; assert `<tmp/out>/out.pvd` exists and
  references at least one `*.vtu`. Reload the VTU via meshio in a
  helper Python script (run via `subprocess.run([sys.executable,
  ...])`), assert the six PointData fields are present.
- `test_project_stress_to_mesh_mean_matches_sidecar.cpp`: same setup;
  assert the per-field mean of each `PointData` array equals the
  sidecar's per-field mean (constant field) to 1 part in 1e-9.
- `test_project_stress_to_mesh_bbox_failure.cpp`: drive with a
  sidecar whose bbox does NOT enclose the mesh; the driver must
  exit with a non-zero code and emit the expected
  `ProjectStress: mesh bbox is NOT contained` message.

**Phase 8 (VTU postprocess verifier)**
- `test_verify_bulk_constant_field_is_exact`: feed a constant-σ⁰
  Phase-4 bulk VTU into `verify_bulk_cell_data`; expect L∞ < 1 Pa
  for every component (trilinear interpolation of a constant is
  exact) and `passed = True`.
- `test_verify_fault_point_data_constant_field_is_exact`:
  constant-σ⁰ background; assert L∞ < 1 Pa on every fault point-data
  field — the area-weighted node averaging of a constant is exact.
- `test_verify_fault_cell_data_matches_analytic_rotation`: random
  curved fault mesh + constant σ⁰; for each triangle, predict
  `(σ_n, τ_s, τ_d)` analytically and assert the recorded cell-data
  matches the prediction to 1e-9 Pa (rotation round-off only).
- `test_verify_lithostatic_sv_fault_within_tol`: depth-varying σ_v
  with `Sv_grad = 1e4 Pa/m`; verify the fault VTU with
  `tol_abs_pa = 1e3` — must pass; with `tol_abs_pa = 1e-3` — must
  fail. This sanity-checks the documented error model
  (`||σ - σ_node||_∞ ≤ O(h · grad)`).
- `test_verify_detects_single_corrupted_cell`: load a Phase-4 VTU,
  perturb one cell's `sigma_xx_Pa` by 1 GPa in-memory, re-emit,
  and re-run the verifier. Expect `passed = False` and a
  per-field row whose L∞ matches the perturbation to 1e-3.
- `test_verify_h1_projected_constant_is_machine_precision`: Phase-7
  H1-projected VTU on a constant-σ⁰ sidecar; assert L∞ ≤ 1e-9 Pa.
- `test_verify_cli_exit_code_on_pass`: subprocess-invoke
  `verify_onfault_stress.py --fault-vtu <good>`; expect exit code 0.
- `test_verify_cli_exit_code_on_fail`: subprocess-invoke on the
  corrupted VTU from above; expect exit code 1.
- `test_verify_writes_json_report`: assert a `verify_report.json`
  is produced next to the inputs with a `pass: bool` and a
  per-field section. JSON-load it and assert schema-stability.
- `test_verify_handles_missing_layer_gracefully`: write a VTU with
  only point-data (no cell-data); the verifier must skip the
  cell-data check with an info log and pass on the point-data
  check.
- `test_verify_uses_summary_json_params`: pass `--summary-json
  <phase4_summary>`; assert the verifier reads SHmax/Shmin/Sv/
  P_p/SHmax_az from the summary and produces matching analytic
  predictions (no CLI overrides needed).
- `test_verify_fault_basis_continuity`: load a real (Phase 4)
  fault VTU; for every vertex, assert
  `|strike_vec_node · normal_vec_node| < 1e-9` — i.e., the
  re-orthonormalisation actually produced an orthonormal triple
  at every node. This is the postprocess QA equivalent of the
  Phase 2 unit test on synthetic data.

### Reference data
- The H&Z `demo_safod()` numerical output is the **gold standard** for
  the SAFOD parameter case. Capture its current console-printed values
  once during Phase 3 test authoring, embed them as constants in the
  test file with a comment citing the H&Z file line numbers, and use
  them as the cross-check anchor.
- The Phase 4 JSON summary for the 2000 m mesh + SAFOD parameters
  is the **Python↔C++ parity anchor**. Generate it once at the end
  of Phase 4, freeze its per-field median values into a constant
  `tests/data/safs_2km_safod_medians.json` file checked into the
  repo, and reference it from the Phase 6 parity tests.

---

## Risk Assessment

### Numerical risks
- **Normal-orientation drift across a curved fault.** If
  `harmonise_normal_orientation` flips a sub-cluster of triangles
  incorrectly, the fault basis will have a sign inversion in that
  region, causing `tau_strike` to invert sign locally. **Detection:**
  the JSON summary's `tau_strike.min` and `tau_strike.max` should be
  same-signed for a regional-stress projection onto a single fault;
  a sign mix indicates the harmonisation failed. **Mitigation:**
  the PCA fallback + explicit strike hint + idempotency test cover
  this. Manual ParaView spot-check of `normal_vec` on the 2000 m mesh
  is the final QA.
- **Near-vertical sliver triangles** from the NW-cut pipeline can have
  `|n_z|` very close to 1 (effectively horizontal). The
  `NEAR_HORIZONTAL_NZ_TOL` threshold needs to be loose enough not to
  flag legitimate sub-vertical triangles. We log `n_degenerate`; if
  it is > 0.1 % of the mesh, investigate before treating output as
  trustworthy.
- **Symmetry in σ⁰**: small floating-point asymmetry from the H&Z
  rotation can trip the `ValueError`. Use `atol=1e-9`.

### Codebase risks
- **VTU cell-data lifetime in meshio.** meshio's VTU writer used to
  drop `field_data` silently. Confirm with a round-trip test
  (`write -> read -> compare`) at Phase 4 start; if it drops, fall
  back to writing the convention string into the JSON summary only.
- **Coordinate frame assumption.** The mesh is currently UTM Zone 11N
  metres. If a future mesh re-projects to a different CRS or to a
  rotated local frame, the H&Z assumption `(x = east, y = north,
  z = up)` breaks. **Mitigation:** the bbox check in Phase 1 (with
  warning) flags the obvious case; document the assumption
  prominently in the module docstring and the JSON summary.

### Scope risks
- **The CSM static-equilibration plan defines a fuller, MFEM-side
  pipeline** that needs σ⁰ as a `VectorCoefficient` for downstream
  L2 projection (`data_stress_frictionlaw/community_stress_model/PLAN_csm_static_equilibration.md` §1).
  Phase 5 of this plan produces the sidecar that plan's Phase 1
  consumes; we **must not** start implementing the L2-equilibration
  solve here.
- **Multi-fault future case.** The current pipeline produces a single
  fault group (`fault`, tag 101). When multi-fault becomes real, the
  module's `extract_cells_by_physical` will need a list-input form.
  Mark this as a follow-up in the module docstring; do not pre-build
  the multi-group form now (it would invite premature abstraction and
  the user's feedback memory warns against this).

### Phase-6 / Phase-7 specific risks
- **Sign-convention drift across phases.** The single sign flip
  lives in `bulk_stress_tensor_field` (Phase 3); from there
  downstream every artefact is compression positive. The risk is
  that a future change adds a second flip somewhere — e.g. a
  refactor that "fixes" what looks like a wrong sign without
  realising it has already been corrected. **Mitigations:**
  (1) `test_bulk_stress_tensor_field_sign_flip` asserts the
  Python flip is present and idempotent; (2)
  `test_fault_vtu_compression_positive` asserts every fault VTU
  value is positive in the SAFOD case; (3)
  `test_field_projector_project_fault_prestress_no_sign_flip`
  asserts the C++ rotator does not flip; (4) the parity test
  against the Phase 4 JSON medians detects any composite drift.
  The docstring of `bulk_stress_tensor_field` cites this list of
  guards.
- **Continuous-field vs. raw-cell divergence.** A future developer
  might disable the cell→node averaging "for performance" and read
  the cell-data channel directly, then be confused by the
  discontinuous field. **Mitigation:** the writer's docstring
  states that point-data is the primary channel and cell-data is
  for traceability only; Phase 8's verifier checks both layers
  independently so any quiet regression on the point-data path
  surfaces immediately.
- **One-time-load contract.** `ProjectVelocity` advances
  `call_count_` by three (one per field). `ProjectStress` advances by
  six. The existing T-4-7 unit test in
  `test_data_field_3d.cpp` will break if it asserts an exact count;
  audit that test before running Phase 6 (see
  `feedback_complete_sign_sites.md` — audit ALL occurrences).
- **`bp5_params_.sigma_n` is currently a single scalar.** Several
  consumer code paths read it directly (`rate_state_fault.hpp:806`,
  `872`, `dieterich_ruina.hpp`). The Phase 6 audit must list every
  such site and decide for each whether the SAFS-mode dispatch is
  needed. The user-feedback memory `feedback_complete_sign_sites.md`
  governs: list ALL sites in the docstring of `ComputeSAFSParams`
  before changing any one of them.
- **TOML schema additions can break existing configs.** The
  `[stress]` block must default to "off" (`use_sidecar = false`) and
  the parser must accept TOMLs that omit the block entirely. The
  `test_seas_config_default_stress_off.cpp` test enforces this.
- **C++ unit test build cost.** The Phase 6 / 7 tests require
  rebuilding the SEAS library. Schedule the test runs so each
  topic-area's unit tests can be invoked in isolation (via the
  `make test-stress-projection` target, added in Phase 6).

### Implementation-process risks
- **Easy to silently bake hardcoded numerics** (SHmax = 113 etc.).
  The user's `feedback_no_hardcoded_numbers.md` memory governs:
  *every* numerical default must come from a named CLI flag whose
  default originates from the H&Z `demo_safod` block, captured in
  one named constant near the top of the module.

---

## Files Summary

### Python (Phases 0–5 and 8)

| File | Status | Purpose |
|---|---|---|
| `data_stress_frictionlaw/regional_stress_projection/hickman_and_zoback_regional_stress_projection.py` | **modified (Phase 0 only)** | Add `dump_safod_sigma0` + `--dump` CLI flag |
| `code_preprocess/project_to_fault_stress.py`      | **new** | Library + CLI (Phases 1–4) |
| `code_preprocess/test_project_to_fault_stress.py` | **new** | pytest suite (Phases 0, 1–4) |
| `code_preprocess/build_stress_safs.py`            | **new** | schema-v1 sidecar writer (Phase 5) |
| `code_preprocess/test_build_stress_safs.py`       | **new** | pytest suite (Phase 5) |
| `code_preprocess/verify_onfault_stress.py`        | **new** | VTU postprocess verifier (Phase 8) |
| `code_preprocess/test_verify_onfault_stress.py`   | **new** | pytest suite (Phase 8) |
| `code_preprocess/PLAN_onfaultstress.md`           | **new** | This document |
| `code_preprocess/clean_freesurface_mesh.py`       | unchanged | Reference style |
| `code_preprocess/nw_cut_strip.py`                 | unchanged | Reference style |
| `code_preprocess/msh_to_vtu.py`                   | unchanged | Reference style; we read its outputs |
| `code_preprocess/ts_to_stl.py`                    | unchanged | Reference style |
| `code_preprocess/data_projection/sidecar.py`      | unchanged | Imported as library (Phase 5) |

### Python output artefacts

| File | Phase | Purpose |
|---|---|---|
| `data_projection_onfaultstress/<lc>m/<base>_fault_stress.vtu` | 4 | Fault ParaView (debug/QA) |
| `data_projection_onfaultstress/<lc>m/<base>_bulk_stress.vtu`  | 4 (opt) | Bulk ParaView (debug/QA) |
| `data_projection_onfaultstress/<lc>m/<base>_summary.json`     | 4 | Parameter + stats summary |
| `data_projection_onfaultstress/stress_safs.h5`                | 5 | **Canonical sidecar (production input)** |

### C++ (Phases 6–7)

| File | Status | Purpose |
|---|---|---|
| `io/stress_field_3d.hpp`              | **new** | Bulk 6-component schema-v1 reader |
| `io/stress_field_coefficient.hpp`     | **new** | `mfem::VectorCoefficient` adapter |
| `io/field_coefficient.hpp`            | modified | Add `ProjectStress`, `ProjectFaultPreStress` |
| `io/field_coefficient.cpp`            | modified | Implementations |
| `io/data_field_3d.hpp`                | unchanged | Reused as the per-component reader |
| `fault/fault_geometry.hpp`            | modified | Add `ComputeSAFSParams`, `sigma_n_per_dof_` |
| `fault/rate_state_fault.hpp`          | modified | `safs_mode_` dispatch on per-DOF prestress |
| `friction/dieterich_ruina.hpp`        | audited | every `sigma_n` site listed; modified only where audit requires |
| `config/seas_config.hpp`              | modified | Add `StressConfig` |
| `config/seas_config_parser.hpp`       | modified | TOML `[stress]` block parse |
| `config/seas_config_bridge.hpp`       | modified | Pass `StressConfig` to driver |
| `drivers/seas_driver.cpp`             | modified | Wire `[stress]` opt-in path |
| `drivers/project_stress_to_mesh.cpp`  | **new** | Visualisation driver (Phase 7) |
| `drivers/safs_smoke.toml`             | **new** | Minimal SAFS-mode smoke config |
| `Makefile`                            | modified | `seas_project_stress_to_mesh` target |
| `tests/unit/test_stress_field_3d_*.cpp`           | **new** | Phase 6 reader tests |
| `tests/unit/test_field_projector_project_stress_*.cpp` | **new** | Phase 6 projector tests |
| `tests/unit/test_compute_safs_params_*.cpp`       | **new** | Phase 6 fault-geometry tests |
| `tests/unit/test_rate_state_fault_safs_mode_*.cpp` | **new** | Phase 6 operator dispatch tests |
| `tests/unit/test_seas_config_*_stress*.cpp`       | **new** | Phase 6 TOML tests |
| `tests/unit/test_project_stress_to_mesh_*.cpp`    | **new** | Phase 7 driver tests |
| `tests/data/safs_2km_safod_medians.json`          | **new** | Frozen parity-anchor medians |
