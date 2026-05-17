# stress/

Community Stress Model (Johnson & Hearn) + Hickman–Zoback SAFOD regional
parameters → bulk Cauchy stress tensor σ⁰(z) → per-mesh **on-fault**
tractions and **on-bulk** stress tensor field (ParaView VTU artefacts).

The pipeline is pure Python (`meshio`, `numpy`) and produces only
visualization / postprocessing artefacts. The HDF5 sidecar
(`stress_safs.h5`, used by the C++ runtime to interpolate σ_seas via
`DataField3D`) is a separate downstream step and is **not** required to
regenerate the per-mesh VTU outputs below.

## Directory layout

```
raw/                CSM_Johnson_Hearn_raw_*km.csv   CSM principal stresses per depth
                    friction_law_params/            placeholder for upcoming
                                                    friction parameters (empty)

code/               hickman_and_zoback_regional_stress_projection.py
                                                    H&Z helpers (resolve_traction,
                                                    build_bulk_stress_tensor,
                                                    fault_basis_vectors, demo_safod)
                    project_to_fault_stress.py      fault VTU + bulk VTU →
                                                    per-mesh {bulk,fault}_stress.vtu
                                                    + summary.json
                    build_stress_safs.py            bulk meshes → stress_safs.h5
                                                    (C++ runtime sidecar; optional)
                    verify_onfault_stress.py        sanity guards (symmetry, units,
                                                    analytic re-evaluation match)
                    test_*.py                       pytest suites

docs/               PLAN_csm_static_equilibration.md
                    PLAN_onfaultstress.md (+ .pdf), PLAN_onfaultstress_fix.md
                    PHASE6_DEVIATIONS.md, PHASES_6_7_8_COMPLETION.md
                    REVIEW.md, REVIEW_fix.md

results/            stress_safs.h5                  bulk-tensor sidecar (optional)
                    {500m,1000m,2000m}/             basic NW-cut variant
                    {500m,1000m,2000m}_lcfar3000/   far-field lc = 3 km
                    {500m,1000m,2000m}_lcfar5000/   far-field lc = 5 km
                    {500m,1000m,2000m}_zgraded/     z-graded refinement
                       each holds: <base>_bulk_stress.vtu
                                   <base>_fault_stress.vtu
                                   <base>_summary.json
```

## Per-file responsibilities

### `code/hickman_and_zoback_regional_stress_projection.py`
Reference H&Z (2004) regional stress helpers — the analytic core that
everything downstream depends on.
- `azimuth_to_math_angle` / `math_angle_to_azimuth` — convert between
  geological azimuth (cw from N) and math angle (ccw from E).
- `build_bulk_stress_tensor(SHmax, Shmin, Sv, SHmax_az_deg)` — assemble
  the 3×3 bulk Cauchy σ⁰ in the (east, north, up) frame, **compression
  POSITIVE (SEAS convention)** after the single source-site sign flip.
- `fault_basis_vectors(normal, up)` — Tandem-style fault basis
  `(s = up × n̂, d = s × n̂, n̂)` (CLAUDE.md "Fault-local tangent frame").
- `resolve_traction(sigma, n̂, s, d, P_p)` → `ResolvedTraction`
  (σ_n_total, σ_n_eff, τ_strike, τ_dip, |τ|, rake, μ_apparent).
- `compute_fault_stress` / `sweep_strike` — convenience drivers for the
  analytic 2-D sweeps shown in the docs.
- `demo_safod()` — canonical SAFOD parameter set
  (SHmax=113, Shmin=49, Sv=45, P_p=16 MPa, SHmax_az=23°).
- `dump_safod_sigma0(path)` — dump the H&Z params as a JSON blob the
  Phase-3 pipeline can re-read.

### `code/project_to_fault_stress.py`
Phase 1–4 driver. Reads a `<base>_fault.vtu` (+ optional `<base>_bulk.vtu`
sibling) and writes the per-mesh artefacts. Importable as a library and
runnable as a CLI in three modes (`--print-info` | `--batch` |
positional). Key entry points:
- `load_fault_mesh` / `load_bulk_mesh` — meshio readers that enforce
  triangle / tetra cell types and warn on non-UTM-magnitude coordinates.
- `triangle_geometry`, `tet_geometry` — per-cell centroids, normals,
  areas/volumes (matches `msh_to_vtu.py` algebra).
- `build_fault_basis` → `FaultCellGeometry` (centroids, unit normals,
  strikes, dips, areas, n_degenerate). Strike orientation is harmonised
  against the configurable strike-azimuth hint (default 314° for SAF).
- `bulk_stress_tensor_field(z, …)` — vectorised σ_seas(z) in MPa
  (constant or `lithostatic_sv` depth model).
- `project_stress_onto_fault` — per-triangle resolution of σ⁰ onto the
  Tandem fault frame, returning `ResolvedTraction`s.
- `write_fault_vtu` / `write_bulk_vtu` — VTU writers with cell- and
  node-averaged sigma_n / tau_strike / tau_dip / |τ| / rake / μ_apparent.
- `write_summary_json` — per-mesh summary (parameters echo,
  σ⁰ at z=0, fault/bulk cell counts, stats min/median/max).
- `_extract_lc_tag(mesh_base)` — `..._<N>m` → `<N>m`,
  `..._<N>m_<variant>` → `<N>m_<variant>` (latter added to keep
  variant outputs in sibling subdirs alongside the basic ones).
- CLI: `--print-info`, `--batch` (iterates `DEFAULT_INPUT_BASES`),
  positional `(input_mesh, mesh_base)` for one-off runs.

### `code/build_stress_safs.py`
Phase-5 schema-v1 HDF5 sidecar writer. Reads one or more bulk meshes,
computes σ_seas (Pa) on a uniform UTM Zone 11 N rectilinear grid that
strictly encloses the union bbox, and emits `stress_safs.h5` consumable
by the C++ `DataField3D` class (no C++-side change). **Independent of
the per-mesh VTU pipeline** — only run this if the SEAS C++ runtime needs
to query σ_seas during initial-condition setup. Note: as of the
`code_preprocess → velocity/` reorg the script's
`_DATA_PROJ_DIR = parent/"data_projection"` sidecar import is stale
(sidecar.py now lives at `velocity/code/sidecar.py`); fix before the
next h5 rebuild.

### `code/verify_onfault_stress.py`
Phase-8 postprocess verifier. Re-evaluates σ_seas (and its on-fault
resolution) analytically at every recorded cell centroid and vertex
and compares against the values that `project_to_fault_stress.py`
wrote. The verifier covers the full population of cells/nodes (no
sampling), checks three layers — bulk cell-data, fault cell-data,
fault point-data — and returns non-zero exit on tolerance violation.
Required CLI args: `--fault-vtu`, `--bulk-vtu`, `--summary-json` (or
`--params`); optional `--report-json` to persist results, `--tol-pa`
to override absolute tolerance, `--tol-rel` for relative. See
workflow §D below for the canonical invocation and the all-variants
result (`12 / 12` fields pass at machine precision).

### `code/test_*.py`
Pytest suites. Cover analytic identities for the H&Z math, geometry
helpers, the Phase-3/4 writers, CLI smoke tests on the 2000m mesh, and
the sidecar contract.

## Inputs

- **Mesh** — fault + bulk VTUs from `meshing/results/vtu/`. Three NW-cut
  resolutions (500m, 1000m, 2000m) each in four variants: basic,
  `lcfar3000` (far-field lc = 3 km), `lcfar5000` (far-field lc = 5 km),
  `zgraded` (z-graded refinement). The fault triangulation is identical
  across variants of the same resolution; only the bulk tet count
  varies.
- **Regional stress (H&Z SAFOD)** — `SHmax = 113`, `Shmin = 49`,
  `Sv = 45`, `P_p = 16` MPa, `SHmax_az = 23°` cw-from-N. These are the
  H&Z `demo_safod` defaults, used as the constants in the Python file —
  the on-fault projection does **not** need them to be precomputed into
  a sidecar.
- **CSM raw csvs** (`raw/CSM_Johnson_Hearn_raw_*km.csv`) — depth slices
  of Johnson & Hearn principal stresses, reserved for future
  depth-varying models (currently the pipeline uses the constant H&Z
  σ⁰).

## Workflow

### A. Per-mesh fault + bulk projection (recommended starting point)

```bash
conda activate pythonenv
cd stress/code

# Basic --batch: the hardcoded 500m / 1000m / 2000m bases.
python project_to_fault_stress.py --batch --write-bulk

# Variant mesh (single-file mode). Repeat for each <N>m × {lcfar3000,
# lcfar5000, zgraded}. `mesh_base` is the filename stem; the output
# subdir is auto-derived as `<N>m_<variant>/`.
python project_to_fault_stress.py \
    ../../meshing/results/vtu/safs_fault_box_nwcut_500m_lcfar3000_fault.vtu \
    safs_fault_box_nwcut_500m_lcfar3000 \
    --write-bulk
```

Outputs per run land in `results/<lc_tag>/`:
- `<base>_fault_stress.vtu` — surface mesh with cell- and node-averaged
  `sigma_n_total_MPa`, `sigma_n_eff_MPa`, `tau_strike_MPa`,
  `tau_dip_MPa`, `tau_magnitude_MPa`, `rake_deg`, `mu_apparent`, plus
  the per-cell basis (`strike`, `dip`, `normal`).
- `<base>_bulk_stress.vtu` — tetra mesh with cell-data
  `sigma_tensor_MPa` (the full 6-component symmetric σ_seas evaluated
  at tet centroids).
- `<base>_summary.json` — parameter echo, σ⁰ at z=0, fault/bulk cell
  counts, and min/median/max per output field.

### B. Quick mesh diagnostic (no writes)

```bash
python project_to_fault_stress.py \
    --print-info ../../meshing/results/vtu/safs_fault_box_nwcut_500m_fault.vtu
```

Prints the H&Z params, σ⁰ at z=0, triangle count, and strike-azimuth
statistics; useful for verifying a new mesh before paying the full
projection cost.

### C. Rebuild C++ runtime sidecar (only when needed)

```bash
python build_stress_safs.py \
    --meshes ../../meshing/results/vtu/safs_fault_box_nwcut_500m_bulk.vtu \
             ../../meshing/results/vtu/safs_fault_box_nwcut_1000m_bulk.vtu \
             ../../meshing/results/vtu/safs_fault_box_nwcut_2000m_bulk.vtu \
    --out    ../results/stress_safs.h5
```

(The sidecar lookup path is stale post-reorg — see `build_stress_safs.py`
above.)

### D. Post-process validation against the analytic predictor

`verify_onfault_stress.py` re-evaluates σ_seas analytically at every
recorded centroid / vertex and compares against the values the writer
stored in the VTU. It is strictly stronger than a random-sample
spot-check — it covers the full population of cells and nodes — and
returns non-zero exit if any field exceeds its tolerance.

The verifier runs three checks per artefact set:

1. **Bulk VTU cell-data** (six `sigma_<ij>_MPa` components, converted
   to Pa). Compares to `bulk_stress_tensor_field(z)` evaluated at each
   tetra centroid. Tolerance `1.0 Pa` (machine round-off only — σ⁰ is
   constant for the `depth_model = constant` case).
2. **Fault VTU cell-data** (`sigma_n_total_MPa_cell`,
   `tau_strike_MPa_cell`, `tau_dip_MPa_cell`). Re-evaluates the
   per-triangle Tandem basis from the VTU connectivity and projects
   σ⁰ onto each triangle. Tolerance `1.0 Pa`.
3. **Fault VTU point-data** (`sigma_n_total_MPa`, `tau_strike_MPa`,
   `tau_dip_MPa`). Replays the writer's exact two-step pipeline —
   per-cell projection followed by area-weighted node averaging via
   `cell_to_node_average` — so the comparison is against the writer's
   own averaging convention, not against a node-averaged-basis
   prediction (the two differ on a faceted fault). Tolerance `1.0e3 Pa`.

Per-variant CLI:

```bash
conda activate pythonenv
cd stress

VARIANT=500m_lcfar3000
BASE=safs_fault_box_nwcut_${VARIANT}
python code/verify_onfault_stress.py \
    --fault-vtu    results/${VARIANT}/${BASE}_fault_stress.vtu \
    --bulk-vtu     results/${VARIANT}/${BASE}_bulk_stress.vtu \
    --summary-json results/${VARIANT}/${BASE}_summary.json \
    --report-json  results/${VARIANT}/verify_report.json
```

Outputs (printed and persisted to `verify_report.json`):

- Per-field `L_inf` and `L_2` error in Pa.
- The 3-D location of the L∞ violation.
- A pass/fail flag against the per-field tolerance.

**Result for all 12 mesh variants:** every variant passes 12 / 12
fields. Worst per-resolution L∞ across all fields:
`500m` 6.0e-08 Pa · `1000m` 7.5e-08 Pa · `2000m` 6.0e-08 Pa — i.e.,
machine round-off of a ~1e8 Pa stress magnitude. The `lcfar3000`,
`lcfar5000`, and `zgraded` variants are identical to the basic ones
at each resolution because they share the same fault triangulation.

Run the verifier's own unit tests alongside the writer's:

```bash
cd stress/code
pytest -q test_verify_onfault_stress.py test_project_to_fault_stress.py
```

## Conventions (carried end-to-end)

- Frame: UTM Zone 11 N metres, `(x = east, y = north, z = up)`,
  `z ≤ 0` underground.
- σ⁰ on the bulk side: **compression POSITIVE** (SEAS internal) after
  the single source-site flip in `build_bulk_stress_tensor`.
- `σ_n_eff = σ_n_total − P_p` with `P_p > 0`.
- Fault-local tangent frame: Tandem `(s, d, n̂)` with `s = up × n̂`,
  `d = s × n̂` (down-dip). `tau_strike > 0` means right-lateral.
- Units: MPa in all VTU/JSON outputs; Pa only in the C++ sidecar
  `stress_safs.h5`.

## Current results inventory (12 mesh variants)

| variant | fault tri | bulk tet | σ_n range (MPa) | \|τ\| range (MPa) |
| --- | ---:| ---:| --- | --- |
| `500m`           | 42358 | 1,004,280 | [48.9, 113.0] | [0.19, 34.0] |
| `500m_lcfar3000` | 42358 | 1,114,485 | [48.9, 113.0] | [0.19, 34.0] |
| `500m_lcfar5000` | 42358 |   528,279 | [48.9, 113.0] | [0.19, 34.0] |
| `500m_zgraded`   | 42358 | 1,457,366 | [48.9, 113.0] | [0.19, 34.0] |
| `1000m`           | 10603 |   190,831 | [50.4, 113.0] | [0.16, 34.0] |
| `1000m_lcfar3000` | 10603 |   965,127 | [50.4, 113.0] | [0.16, 34.0] |
| `1000m_lcfar5000` | 10603 |   381,757 | [50.4, 113.0] | [0.16, 34.0] |
| `1000m_zgraded`   | 10603 | 1,313,059 | [50.4, 113.0] | [0.16, 34.0] |
| `2000m`           |  2685 |   145,300 | [50.4, 113.0] | [0.21, 34.0] |
| `2000m_lcfar3000` |  2685 |   918,731 | [50.4, 113.0] | [0.21, 34.0] |
| `2000m_lcfar5000` |  2685 |   335,284 | [50.4, 113.0] | [0.21, 34.0] |
| `2000m_zgraded`   |  2685 | 1,267,727 | [50.4, 113.0] | [0.21, 34.0] |

Fault stats are identical across variants of the same resolution (same
fault triangulation → same on-fault projection); only the bulk tet
count differs, driven by the far-field sizing strategy.
