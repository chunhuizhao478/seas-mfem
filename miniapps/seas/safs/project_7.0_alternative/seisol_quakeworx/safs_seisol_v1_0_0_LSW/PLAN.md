# Implementation Plan: Wire CVM (multiscale_statewise) velocity model into the SAFS SeisSol case

## Overview

Replace the constant elastic medium (`safs_material.yaml`, ρ=2670, μ=λ=32 GPa) in the
verified SAFS SeisSol setup with the heterogeneous CVM dataset
`velocity/results/multiscale_statewise_cvm/velocity_safs.h5`. SeisSol ingests 3-D gridded
material through **easi `!ASAGI`**, which reads a **COARDS NetCDF file on a uniform
(equidistant) grid** — so the work is (1) a converter script that resamples the rectilinear
HDF5 sidecar onto a uniform grid, converts Vp/Vs/ρ → ρ/μ/λ, and writes ASAGI NetCDF;
(2) a new easi material YAML pointing at that file; (3) validation locally and via a
**QuakeWorx** (https://quakeworx.org) job. **No SeisSol source code is modified.**

Execution platform: jobs are submitted through the QuakeWorx SeisSol app, which runs a
pre-built SeisSol on gateway compute resources and consumes the case folder as-is
(parameters.par + yaml + mesh + data files) — so the entire deliverable is *files in this
folder*; there is nothing to build or install on the platform ourselves. ASAGI support
in the gateway's SeisSol is confirmed by precedent: the SeisSol Training repo's Kaikoura
QuakeWorx case runs `!ASAGI` NetCDF material on the app (see "QuakeWorx precedent"
below).

## Established facts (verified during planning — do not re-derive)

### Source data (`velocity_safs.h5`, schema `data_projection_v1`)
- Layout: `grid/x (264)`, `grid/y (195)`, `grid/z (39)` float64 axes;
  `fields/{Vp,Vs,density} (264,195,39)` float64, shape = `(nx,ny,nz)`, no NaN.
- CRS `EPSG:32611` (UTM 11N, meters, x=East y=North z=Up/elevation) — **identical to the
  SeisSol mesh CRS; no coordinate transform needed.**
- x: 303000→697500 @ uniform 1500 m. y: 3612000→3903000 @ uniform 1500 m.
- z: **non-uniform**, 39 levels from −70000 to +100:
  `[-70k,-60k,-50k,-40k,-30k]` (10 km), `[-30k..-20k]` (2 km), `[-20k..-3k]` (1 km, with
  one 500 m step at −18k→−17k... see grid/z), `[-3k..-1k]` (250–500 m), `[-1k..+100]` (100–250 m).
  Topmost level z=+100 is a synthetic clone of the surface slice (`--extend-z-top 100`).
- Units: SI (m/s, kg/m³). Field ranges: Vp 979–9171, Vs 162–5165, ρ 1244–3342.
- Derived moduli are valid everywhere: λ ∈ [1.13e9, 1.30e11] Pa (no λ≤0), μ ∈ [3.3e7, 8.1e10] Pa,
  Vp/Vs ≥ 1.469, Poisson ratio ∈ [0.068, 0.487].

### Target mesh
- `safs_mesh.puml.h5`: bbox x [314840.1, 672329.5], y [3642278, 3888985], z [−41607.6, 0].
- **Fully contained** in the CVM grid with ≥9.5 km lateral margin and the grid extending
  to +100 m above the mesh top. Truncating the grid at z=−45000 still leaves >3 km margin
  below the mesh bottom.

### SeisSol-side mechanics (file:line refs in `/Users/chunhuizhao/projects/SeisSol`)
- `MaterialFileName` (`&equations`) → easi YAML, parsed in
  `src/Initializer/Parameters/ModelParameters.cpp:59` and evaluated in
  `src/Initializer/ParameterDB.cpp:410-445`.
- Evaluation points: by default `UseCellHomogenizedMaterial = 1`
  (`ModelParameters.h:36`) → material sampled at `ConvergenceOrder³` quadrature points
  per tet and **homogenized** (`ParameterDB.cpp:467-516`): ρ arithmetic mean, μ harmonic
  mean, λ derived from averaged μ and Poisson ratio. `UseCellHomogenizedMaterial = 0` →
  single sample at the element barycenter (`ParameterDB.cpp:147-160`).
- `!ASAGI` easi map: requires **equidistant spacing per dimension** (`docs/asagi.rst:130-132`);
  NetCDF must follow COARDS: 1-D coordinate variables named like their dimensions, and a
  variable (default name `data`) of a **compound type** `material{rho;mu;lambda}` with
  dimension order **`data(z, y, x)`** — z slowest, x fastest (`docs/asagi.rst:26-37,128`).
  Reference writer: `preprocessing/science/generating_ASAGI_file.py` (float32 compound,
  variable `data`, coordinate vars float32).
- Interpolation: `interpolation: nearest` or `linear`. MFEM's projector used **trilinear**
  (`drivers/project_velocity_to_mesh.cpp:44`), so use `linear` for parity.
- Out-of-range queries: ASAGI clamps to the nearest grid point and logs
  "ASAGI: Coordinate in dimension N is out of range. Fixing." (`docs/easi.rst:255-260`).
  With our containment margins this must NOT appear; treat any occurrence as a failure.
- ASAGI is an **optional build**: CMake `-DASAGI=ON` (`CMakeLists.txt:185-194`). A YAML
  using `!ASAGI` on a non-ASAGI build fails at easi parse time with a generic exception
  (`ParameterDB.cpp:708-720`).
- All material units are SI (`docs/easi.rst:82-137`) — matches the sidecar.

### QuakeWorx precedent (github.com/SeisSol/Training, checked 2026-06-09)
- The **Kaikoura** training case ships a QuakeWorx notebook (`kaikoura/kaikoura_qwx.ipynb`)
  AND a 3-D ASAGI material model (`NZ_rhomulambda_large.yaml` referencing
  `NZ_asagi_7_4.nc`, `NZ_asagi_40_4.nc`) — i.e. **ASAGI NetCDF material runs on the
  QuakeWorx SeisSol app today**. The build-support question is settled by precedent.
- Reference `.nc` layout (inspected `sulawesi/3dvel_Sulawesi.nc`): float32 compound
  `material{rho;mu;lambda}`, variable `data(z,y,x)`, 1-D coordinate variables — exactly
  the layout this plan's converter writes. (The reference files use float32 coordinates,
  which show last-bit spacing jitter that ASAGI tolerates; we write float64 coordinates
  to avoid even that.)
- Community YAML idiom (Kaikoura & Sulawesi): an `!Any`/`!IdentityMap` component chain
  `fine !ASAGI → [coarser !ASAGI →] !ConstantMap` so any point that misses the grid falls
  through to a constant background instead of relying on ASAGI's clamp-with-warning.
  Phase 2 adopts this with our verified constant medium as the fallback.
- QuakeWorx job I/O convention (`cdb_tpv23/parameters_qwx.par`, `*_qwx.ipynb`): the app
  consumes the case folder by filename exactly as a local run would (relative
  `MaterialFileName`/`MeshFile`, `OutputFile = 'outputs/...'`); results appear under
  `/jobs/<JOB_NAME>/` on the gateway JupyterHub for post-processing.

## Constraints

- **Do not modify SeisSol source** (`/Users/chunhuizhao/projects/SeisSol` is read-only reference).
- **Do not modify the MFEM velocity pipeline** (`velocity/code/*`); the sidecar HDF5 is a
  read-only input.
- Keep the verified constant-material configuration runnable: the constant
  `safs_material.yaml` stays untouched; the CVM model goes in a NEW file and is selected
  by a one-line `MaterialFileName` change in `parameters.par` (easy A/B).
- The physics math: `μ = ρ·Vs²`, `λ = ρ·(Vp² − 2·Vs²)`. Conversion to moduli happens in
  the converter (NetCDF stores ρ/μ/λ directly) so the easi YAML stays a plain `!ASAGI`
  map matching the existing `safs_material_asagi.yaml.example`.
- All deliverable scripts/files live in `safs_seissol/` and run with
  `/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python3` (has h5py; netCDF4 must be
  verified/installed there — see Phase 1 acceptance).
- `Plasticity = 0` is unchanged → material YAML needs only `[rho, mu, lambda]`
  (no plastCo/bulkFriction/initial-stress Include).

## Phase 0: Prerequisite checks (no files changed)

### Goal
Confirm the toolchain can actually consume what we will produce, before writing anything.

### Detailed Requirements
1. Verify `netCDF4` python module availability in `pythonenv`
   (`python3 -c "import netCDF4"`); if missing, `pip install netCDF4` in that env.
2. ASAGI support in the QuakeWorx SeisSol build: **confirmed by precedent** — the
   Kaikoura training case runs `!ASAGI` NetCDF material on the QuakeWorx app (see
   "QuakeWorx precedent" above). No dedicated verification job is required; the Phase 3
   smoke run's clean init doubles as the direct confirmation on our case.
3. Record the QuakeWorx SeisSol app's `ConvergenceOrder` if discoverable from its app
   description or run logs (affects homogenization quadrature count only; informational).

### Acceptance Criteria
- [ ] `import netCDF4` succeeds in the env used for conversion (already verified
      2026-06-09 in `pythonenv`: netCDF4 imports and reads the Sulawesi reference file).

### Dependencies
- Depends on: nothing. Required by: Phase 1 (item 1), Phase 3 (item 2).

## Phase 1: Converter script `convert_cvm_to_asagi.py`

### Goal
A self-validating script exists that turns `velocity_safs.h5` into an ASAGI-ready
`safs_material_cvm.nc`, with parity proven against the source data.

### Files to Create
- `safs_seissol/convert_cvm_to_asagi.py` — sidecar HDF5 → ASAGI COARDS NetCDF converter.
- `safs_seissol/test_convert_cvm_to_asagi.py` — pytest suite (pure numpy + h5py + netCDF4).

### Detailed Requirements
1. CLI (argparse), defaults chosen so a bare invocation produces the deliverable:
   ```
   python3 convert_cvm_to_asagi.py \
     --sidecar  <path to velocity_safs.h5>          (required)
     --out      safs_material_cvm.nc                 (default)
     --z-min    -45000.0    # truncate below mesh bottom (−41608) with 3.4 km margin
     --z-max    100.0       # keep the synthetic above-surface slice
     --dz       100.0       # uniform z spacing of the output grid
     --vs-floor None        # optional: clamp Vs up to this value (m/s) BEFORE moduli
     --dtype    float32     # 'float32' (default, matches generating_ASAGI_file.py) or 'float64'
     --plot-check / --no-plot-check  # optional PNG slice comparison (matplotlib)
   ```
2. Read `grid/{x,y,z}` and `fields/{Vp,Vs,density}`; assert root attrs
   `crs == 'EPSG:32611'`, `units == 'm'`, `z_positive == 'elevation'`,
   `schema_version == 'data_projection_v1'` — hard error otherwise.
3. Assert x and y are uniform (`np.allclose(np.diff(x), dx)`); x/y are used **as-is**
   (no horizontal resampling). Assert z strictly increasing.
4. Build the uniform output z axis:
   `z_out = np.arange(z_min, z_max + 0.5*dz, dz)` (float64; 452 levels for the defaults).
   Assert `z_out[0] >= z[0]` and `z_out[-1] <= z[-1]` (interpolation only, never
   extrapolation).
5. Resample each field along z only, per (i,j) column, with
   `np.interp(z_out, z_src, field[i, j, :])` — vectorized implementation required:
   `def resample_z(field: np.ndarray, z_src: np.ndarray, z_out: np.ndarray) -> np.ndarray`
   (input `(nx,ny,nz_src)` → output `(nx,ny,nz_out)`; use `scipy.interpolate.interp1d`
   on axis=2 or a searchsorted-based weights formulation; a python loop over 264×195
   columns is acceptable as fallback if it runs < 60 s).
   Linear-in-z matches what trilinear evaluation of the source would give at the source
   levels and is the same family of interpolation MFEM used.
6. Optional Vs floor: if `--vs-floor V` is given, `Vs = np.maximum(Vs, V)` **after**
   resampling, then enforce `Vp = np.maximum(Vp, sqrt(2)*1.001*Vs)` so λ stays > 0.
   Default off (parity with MFEM first; the floor is a later scientific knob).
7. Convert to moduli (float64 math): `rho = density`, `mu = rho*Vs**2`,
   `lambda = rho*(Vp**2 - 2*Vs**2)`. Assert `lambda.min() > 0` and `mu.min() > 0`.
8. Write NetCDF (netCDF4, format `NETCDF4`):
   - Dimensions `x (264)`, `y (195)`, `z (len(z_out))`.
   - Coordinate variables `x(x)`, `y(y)`, `z(z)` — float64 (or float32 with `--dtype
     float32`; coordinate precision must represent 1500-m UTM steps exactly, float32 is
     sufficient up to 7e5 but float64 is safer: **always write coordinates as float64**).
   - Compound type `material` with members `rho`, `mu`, `lambda` of `--dtype`.
   - Variable `data(z, y, x)` of type `material` — **transpose the (nx,ny,nz) arrays with
     `np.transpose(arr, (2,1,0))`** before filling, mirroring
     `generating_ASAGI_file.py:writeNetcdf4SeisSol`.
   - Global attributes for provenance: `source_sidecar` (abspath), `source_created_at`,
     `source_mesh_tag`, `crs='EPSG:32611'`, `z_positive='elevation'`, `dz`, `vs_floor`,
     `converter='convert_cvm_to_asagi.py'`.
9. Built-in round-trip self-check (always runs, prints PASS/FAIL, nonzero exit on FAIL):
   re-open the written NetCDF, draw 1000 fixed-seed (`np.random.default_rng(20260609)`)
   random points uniform in the MESH bbox (hardcode the bbox from this plan), evaluate
   (a) trilinear interpolation of the NetCDF grid and (b) trilinear-in-x,y ×
   linear-in-rectilinear-z interpolation of the SOURCE sidecar, both for ρ/μ/λ
   (source side converts velocities→moduli per-grid-node first, i.e. the same order as
   the converter). Require max relative error ≤ 2e-3 for float32 output (≤ 1e-10 for
   float64) at points more than one source-cell away from z-levels coarser than dz...
   simplified: require **median rel. err ≤ 1e-6 and max rel. err ≤ 5e-3**; the residual
   comes only from re-gridding coarse-z regions and float32 rounding.
10. Print a summary table: output grid shape, file size, ρ/μ/λ min/max, Vs/Vp min/max
    implied, and equivalent min/max wavespeeds (for timestep expectations).

### Interfaces
- `read_sidecar(path: Path) -> tuple[np.ndarray x, y, z, dict fields, dict attrs]`
- `resample_z(field, z_src, z_out) -> np.ndarray` (see req. 5)
- `velocities_to_moduli(vp, vs, rho) -> tuple[rho, mu, lam]`
- `write_asagi_netcdf(path, x, y, z, rho, mu, lam, dtype, attrs: dict) -> None`
- `trilinear_sample(xg, yg, zg, field, pts) -> np.ndarray` (shared by self-check & tests;
  must handle non-uniform zg)
- `main(argv: list[str] | None = None) -> int`

### Edge Cases to Handle
- Output z grid not aligned with source levels (e.g. source −1250, −750 are not multiples
  of 100 offset from −45000): handled by interpolation; no special casing.
- `--z-min` below source z[0] or `--z-max` above source z[-1] → argparse-level hard error
  (req. 4 assert) with a message citing the source z range.
- Source file missing fields/attrs → hard error naming the missing item.
- `--dz` that does not divide `z_max - z_min` → last point capped at `z_max` via the
  arange epsilon in req. 4; document in --help.
- NaN anywhere after resampling → hard error (defense-in-depth; schema forbids NaN).

### Acceptance Criteria
- [ ] `python3 convert_cvm_to_asagi.py --sidecar .../multiscale_statewise_cvm/velocity_safs.h5`
      produces `safs_material_cvm.nc` ≈ 280 MB (float32, 264×195×452) with PASS self-check.
- [ ] `ncdump -h safs_material_cvm.nc` (or netCDF4 introspection in the test) shows
      compound `material{rho;mu;lambda}`, variable `data(z, y, x)`, coordinate vars x/y/z.
- [ ] pytest suite passes:
      - axes uniformity assertions fire on a fabricated non-uniform-x sidecar;
      - `velocities_to_moduli` reproduces μ=λ=32 GPa for Vp=5389.6, Vs=3461.7, ρ=2670
        (the constant-model parity point) to 1e-6 relative;
      - `resample_z` is exact at source z-levels contained in z_out;
      - round-trip check passes on a tiny synthetic sidecar (e.g. 5×4×7 grid with an
        affine field, where trilinear is exact and max rel. err must be ≤ float32 eps·10);
      - λ>0 assertion fires on a fabricated Vp/Vs violating Vp≥√2·Vs.
- [ ] Self-check median/max relative errors printed and within req. 9 bounds.

### Dependencies
- Depends on: Phase 0 item 2. Required by: Phases 2–3.

## Phase 2: easi YAML + parameter wiring

### Goal
SeisSol case files select the CVM material via one explicit switch, with the constant
model preserved for A/B.

### Files to Create
- `safs_seissol/safs_material_cvm.yaml` — community fall-through idiom (Kaikoura/Sulawesi
  training cases) with the verified constant medium as last-resort background:
  ```yaml
  !Switch
  [rho, mu, lambda]: !Any
    components:
      - !ASAGI
          file: safs_material_cvm.nc
          parameters: [rho, mu, lambda]
          var: data
          interpolation: linear
      - !ConstantMap
          map:
            rho:    2670.0
            mu:     3.2e10
            lambda: 3.2e10
  ```
  (Mesh containment makes the fallback unreachable in practice; it exists so a future
  mesh tweak degrades gracefully instead of clamp-warning. Note: implementer must
  verify the exact `!Any` nesting syntax against the bundled easi examples —
  `kaikoura/NZ_rhomulambda_large.yaml` uses `!Any` with a `- !ASAGI` list, while easi
  docs sometimes show `components:` explicitly; copy the working training-case form.)
  Header comments must state: source dataset (multiscale_statewise_cvm sidecar + its
  `created_at`), grid (264×195×452 @ 1500/1500/100 m, z ∈ [−45000, 100]), CRS, and the
  converter command line used.

### Files to Modify
- `parameters.par` — change `MaterialFileName = 'safs_material.yaml'` →
  `'safs_material_cvm.yaml'`, with a comment naming the constant file as the fallback.
  **No other key changes** (notably keep `UseCellHomogenizedMaterial` unset = default 1;
  homogenization is desirable for a heterogeneous medium).
- `safs_material_asagi.yaml.example` — either delete or update its comments to point at
  the now-real `safs_material_cvm.yaml`; it claimed "asagiconv" and a `safs_velocity.nc`
  filename that the implemented flow supersedes. Prefer: replace its body with a
  pointer comment to avoid two divergent ASAGI examples.
- `README_port.md` — add a "Heterogeneous material (CVM)" section: converter invocation,
  the QuakeWorx workflow (upload the case folder including `safs_material_cvm.nc` to the
  SeisSol app, submit), the ASAGI-support assumption and its failure signature, the A/B
  switch, the validation gates of Phase 3, and the scientific caveats (Vs_min 162 m/s vs
  mesh resolution; expected timestep reduction since Vp_max ≈ 9160 > 5390; nucleation
  overstress was tuned for the constant medium — rupture behavior will differ from the
  MFEM constant-material verification run).
- `MFEM_to_SeisSol_mapping.md` §7 — one-paragraph update: constant `!ConstantMap` is the
  verified baseline; CVM via `!ASAGI` is the heterogeneous variant (link to README
  section). Do not rewrite history elsewhere in the doc.

### Edge Cases to Handle
- Running the CVM YAML on a non-ASAGI build: README must state the expected failure mode
  (easi parse exception at init) and the resolution (QuakeWorx support request).
- Relative path resolution: `file:` in easi is resolved relative to the working directory
  of the run; keep `safs_material_cvm.nc` next to the YAML in the case folder so the
  QuakeWorx app picks both up together (same convention the existing fault YAML uses).

### Acceptance Criteria
- [ ] `parameters.par` diff is exactly one functional line (+comments).
- [ ] Constant-material run still selectable by reverting that single line.
- [ ] YAML is valid easi (validated in Phase 3 smoke run; locally at minimum
      `python3 -c "import yaml,sys; yaml.safe_load(...)"`-level syntax check is NOT
      sufficient for easi tags — note this; the real check is the smoke run).

### Dependencies
- Depends on: Phase 1. Required by: Phase 3.

## Phase 3: Validation (local + QuakeWorx)

### Goal
The CVM-material run initializes cleanly on QuakeWorx and the loaded material provably
matches the source CVM.

### Detailed Requirements
1. **Upload**: add `safs_material_cvm.nc` + `safs_material_cvm.yaml` + updated
   `parameters.par` to the case folder used by the QuakeWorx SeisSol app (same bundle
   that ran the verified constant-material job). The ~280 MB NetCDF is the only new
   large file — confirm the app's upload path handles it (use the gateway's data
   manager rather than a browser form if there is a size limit).
2. **Smoke run** submitted through the QuakeWorx SeisSol app (short `EndTime`, e.g.
   1.0 s, or reduced walltime): gates, in order:
   - init reaches time loop: no easi/ASAGI exception (this is also the ASAGI-support
     verification of Phase 0, if not front-loaded);
   - **zero** occurrences of "out of range. Fixing." in the log;
   - logged global timestep DECREASES vs the constant run (Vp_max ≈ 9160 vs 5390 m/s;
     expect roughly ×0.6 or smaller depending on where fast material meets small cells);
   - memory per node acceptable (~280 MB grid + SeisSol's usual footprint).
3. **Material parity gate** (quantitative): enable volume output for one wavefield snapshot
   OR use the `&Output` material/wavefield fields available in the build; if the build
   exposes no material output, fall back to: run with `UseCellHomogenizedMaterial = 0`
   ONCE, and compare SeisSol's reported maximum wavespeed / timestep against the
   converter's printed Vp_max at barycenters; plus visually compare the free-surface or
   wavefield snapshot pattern against the sidecar's `velocity_safs.vtr` shallow slices
   (basins must appear in the same places). Document which variant was used.
4. **A/B sanity**: rerun the constant-material configuration (one-line revert) and confirm
   it still reproduces the previously verified behavior (binary-identical parameter set
   except MaterialFileName).
5. Record outcomes (timestep ratio, parity numbers, log excerpts) in README_port.md's
   new section.

### Edge Cases to Handle
- LTS clustering (`ClusteredLTS = 2`) redistributes with heterogeneous wavespeeds —
  slower overall stepping is expected, not a bug.
- If the smoke run shows clamp warnings despite containment margins, suspect z-axis
  orientation (must be ascending in the NetCDF) — the Phase 1 round-trip check makes
  this nearly impossible, but the log gate catches it.

### Acceptance Criteria
- [ ] Smoke run completes with all gates in req. 2 green.
- [ ] Parity evidence captured (req. 3) and committed to README_port.md.
- [ ] Constant-material A/B revert verified.

### Dependencies
- Depends on: Phases 0–2. Required by: nothing (terminal).

## Testing Strategy
- Phase 1 carries the bulk of automated testing (pytest, synthetic sidecars, analytic
  parity point μ=λ=32 GPa, round-trip trilinear comparison with fixed seed).
- Phase 2 is configuration-only; its test is the Phase 3 smoke run.
- Phase 3 is empirical gating on QuakeWorx job logs plus an A/B revert; quantitative
  parity via timestep/wavespeed cross-check and slice comparison against
  `velocity_safs.vtr`.

## Design decisions and alternatives considered
- **Single uniform-z file (chosen)** at Δz=100 m, z ∈ [−45 km, +100 m] → 264×195×452,
  ~280 MB float32. Preserves the sidecar's full 100-m near-surface fidelity; deep
  oversampling is wasted bytes but harmless. **Alternative (fallback if memory/IO ever
  hurts):** 3 stacked ASAGI files at Δz = 100 m / 1 km / 2.5 km composed with easi
  `!AxisAlignedCuboidalDomainFilter` (~75 MB total) — rejected for v1 because linear
  interpolation cannot cross file boundaries (clamp seams at layer interfaces) and the
  YAML triples in complexity.
- **Store ρ/μ/λ, not Vp/Vs/ρ**: keeps the YAML a plain `!ASAGI` map (no FunctionMap
  post-processing) and matches the existing example file. Note: SeisSol then linearly
  interpolates moduli whereas MFEM interpolated velocities and derived moduli pointwise —
  differences are second-order in cell-size and far below the re-gridding error already
  accepted; additionally SeisSol's quadrature homogenization (harmonic μ) intentionally
  differs from any pointwise scheme.
- **No Vs floor by default**: first land parity with the dataset MFEM uses; clamping is a
  scientific decision (affects basins) to make explicitly later via `--vs-floor`.

## Risk Assessment
- **QuakeWorx SeisSol build lacks ASAGI** — effectively retired: the Kaikoura QuakeWorx
  training case already runs `!ASAGI` NetCDF material on the gateway. Residual exposure
  is only a regression in the deployed app build; the Phase 3 smoke run's clean init is
  the direct check, and the fallback `!ConstantMap` in the YAML does NOT mask such a
  failure (a missing-ASAGI build fails at YAML parse, not per-point lookup).
- **Upload size limits on the gateway** for the ~280 MB NetCDF — mitigation: use the
  gateway data manager; if a hard limit exists, rerun the converter with `--dz 250`
  (~112 MB) or fall back to the 3-layer multi-file variant (~75 MB; see Design
  decisions).
- **Timestep collapse / cost increase**: Vp_max rises 5390→~9160 m/s and the smallest
  cells may now sit in fast material; LTS mitigates. Detect via the smoke-run timestep
  gate; mitigation if severe: revisit `--z-min` truncation (excludes the fastest deep
  material below −45 km already) — the deep mantle Vp 9171 occurs at −70 km, outside the
  grid after truncation; recheck the converter's printed Vp_max (expect ≈ 9160 within
  [−45 km, 0], still high at −40 km).
- **Under-resolved low-Vs basins (scientific, not a crash)**: Vs ≈ 162 m/s in ~1.5 km
  far-field cells supports only very low frequencies; fine for rupture-dynamics-focused
  runs but surface-wave content in basins is unresolved. Documented caveat + `--vs-floor`
  knob.
- **Rupture behavior changes**: hypocentral Vs/μ differ from the constant medium, so the
  20 MPa Gaussian overstress tuned for μ=32 GPa may over/under-drive nucleation. This is
  the intended physics change, but flag it in README so a different rupture is not
  mistaken for a wiring bug; the §2 stress-validation gate of the mapping doc still
  applies unchanged (initial stress is material-independent).
- **Compound-type/netCDF4 pitfalls**: member order or dimension order wrong → garbage
  material that still loads. The Phase 1 round-trip check (independent trilinear read-back
  of the .nc) is specifically designed to catch transposition/member-order bugs before
  anything is uploaded to QuakeWorx.
