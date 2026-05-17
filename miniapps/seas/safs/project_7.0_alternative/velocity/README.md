# velocity/

CVM-H ASCII slices → projected velocity sidecar (and the size-field driver
that pushes velocity contrast into the mesher).

```
raw/                <version>/velocity_raw_*.bp
                                         CVM ASCII extracts (one .bp per depth,
                                         grouped by CVM version subdir:
                                         cvmh/, cvm_s4.26.m01/,
                                         multiscale_statewise_cvm/)
                    CVM_*.png            extraction screenshots from the CVM portal

code/               raw_readers.py       *.bp -> RawSliceStack
                    crs.py               geographic <-> UTM 11N
                    bbox_check.py        G-3 mesh-bbox containment guard
                    sidecar.py           schema-v1 HDF5 writer/reader
                    build_velocity_cvmh.py
                                         driver: raw/<version>/ -> results/<version>/velocity_safs.h5
                    fault_zone_metric.py per-tet L^2 by depth band
                    build_size_field.py  velocity -> .size_field.pos for Gmsh
                    plot_*.py            comparison + diagnostic plots
                    test_*.py            pytest suites

docs/               fault_zone_projection_plan_v{1,2,3}.md (+ v1 .pdf)

results/            <version>/velocity_safs.h5
                                         rectilinear UTM 11N sidecar
                                         (one per CVM version, mirroring
                                         the raw/<version>/ layout)
                    <version>/velocity_safs.vtr
                                         same data as VTK rectilinear grid
                    pngs/                compare_3x3_{Vp,Vs,density}_z*km.png
                    preview/             per-mesh proc000000.vtu previews
                                         (projected_velocity_<RES>m/Cycle000000/)
```

## Workflow

The pipeline has two preprocessing stages plus a diagnostics stage.
Stages 1 and 2 are owned by this pillar; the diagnostics stage depends
on MFEM C++ outputs that come from the SEAS driver
`drivers/project_velocity_to_mesh.cpp`, not from this directory.

```
                                                  ┌───────────────────────┐
   raw/<version>/velocity_raw_*.bp                │  meshing/code/        │
            │                                     │  run_nwcut_meshing.py │
            │  build_velocity_cvmh.py             │   (consumes the .pos) │
            ▼                                     └───────────▲───────────┘
   results/<version>/velocity_safs.h5                         │
   results/<version>/velocity_safs.vtr                        │
            │                                                 │
            │  build_size_field.py                            │
            ▼                                                 │
   ../meshing/results/msh/safs_fault_box_nwcut_<RES>m.size_field.pos
            │                                                 │
            └─────────────────────────────────────────────────┘
                       │
                       │  (MFEM driver, not this pillar)
                       │  project_velocity_to_mesh.cpp
                       ▼
            results/preview/projected_velocity_<RES>m/Cycle000000/proc000000.vtu
                       │
                       │  plot_*.py
                       ▼
            results/pngs/compare_3x3_{Vp,Vs,density}_z*km.png
```

### Stage 1 — sidecar (`build_velocity_cvmh.py`)

```bash
conda activate pythonenv
cd velocity/code

# One invocation per raw-data version.  cvmh/, cvm_s4.26.m01/, and
# multiscale_statewise_cvm/ all share the same 38 slice depths and the
# same lon/lat grid, so the same flags work for each.
for ver in cvmh cvm_s4.26.m01 multiscale_statewise_cvm; do
    python build_velocity_cvmh.py \
        --raw-dir       ../raw/${ver} \
        --out-path      ../results/${ver}/velocity_safs.h5 \
        --grid-dx       1500 \
        --mesh-msh      ../../meshing/results/msh/safs_fault_box_nwcut_500m.msh \
        --extend-z-top  100 \
        --paraview-export \
        --verbose
done
```

What this produces (per `<version>`):

- `../results/<version>/velocity_safs.h5` — schema-v1 HDF5 sidecar
  (264 × 195 × 39 UTM 11N rectilinear grid at `--grid-dx`, CRS
  EPSG:32611, z positive = elevation, z range [-70 km, +100 m])
- `../results/<version>/velocity_safs.vtr` — same data as a VTK
  rectilinear grid (when `--paraview-export` is set)

Notable flags:

- `--mesh-msh PATH` — turns on the G-3 bbox-containment guard so a
  mesh that extends past the velocity grid is caught at build time
  rather than producing silent NaNs in the projector.
- `--extend-z-top METRES` — synthetic top slice that clones the
  shallowest source slice's values. Must be `≥ meshing/pad_top`
  (currently 100 m) or the G-3 guard rejects the mesh.
- Per-field sanity bounds default to permissive values
  (`--vp-min-mps`, `--vp-max-mps`, `--vs-min-mps`, `--vs-max-mps`,
  `--rho-min-kgm3`, `--rho-max-kgm3`).

### Stage 2 — Gmsh size field (`build_size_field.py`)

```bash
# One invocation per mesh resolution (500/1000/2000 m).  The size
# field is built from the CVM-H sidecar by default; swap the
# ``cvmh`` segment for another CVM version to drive the mesher off
# that dataset instead.
for res in 500 1000 2000; do
    python build_size_field.py \
        --sidecar ../results/cvmh/velocity_safs.h5 \
        --field   Vs \
        --out-pos ../../meshing/results/msh/safs_fault_box_nwcut_${res}m.size_field.pos \
        --lc-near 1500 --lc-far 10000 --lc-min 500 --alpha 8 \
        --smooth-sigma 1.0 --voxel-stride 2
done
```

What this produces (under `meshing/results/msh/`):

- `safs_fault_box_nwcut_<RES>m.size_field.pos` — Gmsh
  `Field[PostView]` consumed by `safs_fault_box_nwcut.geo` when
  `USE_SIZE_FIELD = 1` (i.e., when `run_nwcut_meshing.py` is invoked
  with `--gen-size-field`).

The size field is `LC = clamp(LC_FAR / (1 + α ‖∇Vs‖ / Vs_med), LC_MIN,
LC_FAR)` per voxel; α controls how aggressively the mesher refines
toward high-gradient regions (defaults reproduce the v3 plan G-1
acceptance band).

### Stage 3 — Diagnostic plots (`plot_*.py`, optional)

These consume the **MFEM preview .vtu** files in `results/preview/`,
which are produced by `drivers/project_velocity_to_mesh.cpp` (built
under the SEAS C++ tree, not in this pillar). After the SEAS driver
emits `projected_velocity_<RES>m<SUFFIX>/Cycle000000/proc000000.vtu`:

```bash
# 3x3 sidecar | baseline | lcfar5000 | lcfar3000 | zgraded comparison.
# Writes results/pngs/compare_3x3_{Vp,Vs,density}_z<DEPTH>km.png
# for DEPTH in {0, 1, 2, 3, 4, 5, 7.5, 10, 15, 20} km.
python plot_sidecar_vs_baseline_multi_z.py

# Legacy diagnostics (only useful if the corresponding preview/
# variants from the v3-plan G-1/G-2 experiments are still on disk):
python plot_baseline_vs_g1g2.py                # baseline vs G-1 size-field
python plot_interp_5panel_multi_z.py           # interp diagnostics
```

The remaining `plot_*.py` scripts cover historical comparisons
(`plot_comparison_with_dg0.py`, `plot_old_vs_new_baseline*.py`) and
are kept for reproducing the v3 retraction-supporting figures.

##### What the four mesh variants do (read this before staring at the PNGs)

The same `velocity_safs.h5` feeds every projection — variants differ
only in **what the gmsh mesher did**.  Every variant uses the
distance-to-fault size field (Field[2]): tets are forced to ~1.5 km
inside a 3 km corridor of the fault and grow to a far-field cap as
you move away.  The variants differ in *how* the cap is set:

| Variant     | far-field cap                          | depth dependence?                            |
|-------------|----------------------------------------|----------------------------------------------|
| `nwcut`     | 10 km (baseline)                       | none — uniform top-to-bottom                 |
| `lcfar3000` | 3 km                                   | none — uniform top-to-bottom                 |
| `lcfar5000` | 5 km                                   | none — uniform top-to-bottom                 |
| `zgraded`   | 10 km, **plus** a depth-step ceiling:  | depth-step cap via Field[Min]:               |
|             | applied via `Field[Min]`               | surface zone (z > -500 m): ≤ 500 m           |
|             |                                        | basin body (-3 km..-500 m): ≤ 1000 m         |
|             |                                        | transition (-10 km..-3 km): ≤ 2000 m         |
|             |                                        | basement (z < -10 km): ≤ 3000 m              |

A vertical-slab cartoon of `zgraded` (the other three look uniform
top-to-bottom by comparison):

```
   z = +100 m  ┌──────────────────────────────────────────────┐  free surface
              │ ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒  │  ← surface ≤ 500 m
   z = -500 m │ ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓  │
              │ ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓  │  ← basin   ≤ 1 km
   z = -3 km  │ ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░  │
              │ ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░  │  ← transition ≤ 2 km
   z = -10 km │ ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒  │
              │ ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒  │  ← basement ≤ 3 km
   z = -41 km └──────────────────────────────────────────────┘  box bottom
```

Both rules (`fault-distance` and `depth-step`) are combined via
`Field[Min]`, so at every point the *smaller* target tet size wins:
the fault corridor always stays at `lc_near = 1500 m` regardless of
depth; outside the corridor, the depth cap kicks in.

**Why `zgraded` matters.**  The sidecar samples Vp/Vs/density at
**100 m vertical spacing inside the top 500 m** and at 250 m between
-1500 m and -500 m, then progressively coarsens with depth (1 km at
-3 km, 10 km below -30 km — see `velocity_safs.h5` `grid/z`).  A
non-zgraded mesh has ~1.5 km tall surface-zone tets that **average
over 15 sidecar slices** before the projector sees them: source
fidelity is real, but lost in averaging.  `zgraded` shrinks
surface-zone tets to ~500 m tall so the mesh can actually resolve
what the sidecar provides.

### Stage 4 — Project velocity onto the meshes (MFEM `seas_project_velocity_to_mesh`)

> *This stage produces the `preview/` .vtu files consumed by Stage 3.*

```bash
# In a fresh shell — the MFEM driver needs `mfem-dev`, not `pythonenv`:
conda activate mfem-dev
cd miniapps/seas
make seas_project_velocity_to_mesh    # if not already built

# MFEM's gmsh reader rejects msh4 with "vertices indices are not
# unique".  Convert each msh4 file to msh2 once (pythonenv has gmsh):
conda activate pythonenv
mkdir -p /tmp/safs_msh2
for stem in \
    safs_fault_box_nwcut_500m safs_fault_box_nwcut_1000m safs_fault_box_nwcut_2000m \
    safs_fault_box_nwcut_500m_lcfar3000 safs_fault_box_nwcut_1000m_lcfar3000 safs_fault_box_nwcut_2000m_lcfar3000 \
    safs_fault_box_nwcut_500m_lcfar5000 safs_fault_box_nwcut_1000m_lcfar5000 safs_fault_box_nwcut_2000m_lcfar5000 \
    safs_fault_box_nwcut_500m_zgraded   safs_fault_box_nwcut_1000m_zgraded   safs_fault_box_nwcut_2000m_zgraded; do
    gmsh ../../meshing/results/msh/${stem}.msh -format msh2 -save \
         -o /tmp/safs_msh2/${stem}.msh -v 0
done

# Then project each onto the regenerated sidecar (run from
# preview/).  Swap the ``cvmh`` segment to project against a
# different CVM version.
cd velocity/results/cvmh/preview
for res in 500 1000 2000; do
    for suffix in "" "_lcfar3000" "_lcfar5000" "_zgraded"; do
        stem="safs_fault_box_nwcut_${res}m${suffix}"
        ../../../../../seas_project_velocity_to_mesh \
            --mesh /tmp/safs_msh2/${stem}.msh \
            --sidecar ../velocity_safs.h5 \
            --out projected_velocity_${res}m${suffix}
    done
done
```

What this produces (under `velocity/results/cvmh/preview/`):

- `projected_velocity_<RES>m<SUFFIX>/projected_velocity_<RES>m<SUFFIX>.pvd`
- `projected_velocity_<RES>m<SUFFIX>/Cycle000000/proc000000.vtu` —
  H1-P2 (default `--order 2`) PointData fields Vp, Vs, density,
  lambda, mu.

Notes:

- **All 12 projections read the same `velocity_safs.h5`.**  The
  sidecar is mesh-agnostic; different variants differ only in the
  tets they sample it onto.
- **MFEM driver = `mfem-dev` env**; **gmsh-format conversion = `pythonenv` env.**
- The converted msh2 files in `/tmp/safs_msh2/` are reproducible from
  the canonical msh4 files in `meshing/results/msh/` and are NOT
  tracked.

### Tests

```bash
cd velocity/code
pytest -q test_data_projection.py test_build_size_field.py    # 27 tests, pure numpy
pytest -q test_fault_zone_metric.py                            # needs parse_mfem_vtu (vtk)
```

`parse_mfem_vtu.py` is an in-tree helper that wraps
`vtkXMLUnstructuredGridReader` — meshio does not parse MFEM's binary
VTU 2.2 dialect.  The plot scripts and `fault_zone_metric.py` import
it via the `sys.path.insert(HERE, …)` shim at the top of each file.
