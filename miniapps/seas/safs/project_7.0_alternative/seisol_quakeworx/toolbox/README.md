# Toolbox — SAFS SeisSol (QuakeWorx) post/pre-processing utilities

Self-contained Python tools. Requirements: `python3` with `numpy` and
`h5py` only (no vtk/meshio/netCDF4 — the `.nc` material file is netCDF4,
which IS HDF5, so h5py reads it directly). The one exception is
`generate_velocity_nc_from_raw`, which additionally needs `scipy`,
`pyproj`, and `netCDF4`.

## Setting up the conda environment ("pythonenv")

One command, using the spec shipped next to this README:

```bash
conda env create -f toolbox/environment.yml     # creates "pythonenv"
conda activate pythonenv
python -c "import numpy, scipy, h5py, pyproj, netCDF4; print('env OK')"
```

Manual alternative (same packages, latest compatible versions):

```bash
conda create -n pythonenv -c conda-forge \
    python=3.13 numpy h5py scipy pyproj netcdf4
conda activate pythonenv
```

Notes:

- `environment.yml` pins the exact versions the shipped artifacts were
  produced/verified with. The pins only truly matter for
  `generate_velocity_nc_from_raw`: a BIT-identical regeneration of
  `safs_material_cvm.nc` requires the same `scipy` (Qhull/Delaunay
  triangulation) and `pyproj` (EPSG transform) versions. With newer
  versions the output is still correct — just not guaranteed byte-equal.
  If the solver rejects a pin on your platform, relax it to
  major.minor (e.g. `scipy=1.16`).
- Optional, for independently verifying generated VTUs with the same
  reader ParaView uses: `conda install -n pythonenv -c conda-forge
  vtk meshio` (large download; none of the tools need it to run).
- No pip packages required; everything is on conda-forge.

All paths are command-line arguments or are resolved relative to the
script's own location, so this folder tree is relocatable — keep the
`toolbox/` subfolder structure and the case folder
(`safs_seisol_v2_0_0_RSSRW/`) side by side as uploaded. Example commands
below are run from the folder that contains `toolbox/`.

## h5_to_vtu/puml_h5_to_vtu.py

Convert a SeisSol PUML mesh (`.puml.h5`, pumgen output) to ParaView VTU
files: bulk tetrahedra + deduplicated dynamic-rupture (fault, BC 3)
surface; `--all-bcs` also writes free-surface / absorbing surfaces.

```bash
python3 toolbox/h5_to_vtu/puml_h5_to_vtu.py \
    safs_seisol_v2_0_0_RSSRW/safs_mesh.puml.h5 \
    --out-dir toolbox/h5_to_vtu --all-bcs
```

## on_fault_stress_projection/project_stress_to_vtu.py

Project the easi `!ConstantMap` initial-stress tensor (SeisSol
convention: compression-negative Pa, effective) onto the PUML fault and
write fault/bulk stress VTUs + a summary JSON. Outputs are in SEAS
convention (compression-positive MPa; Tandem fault basis s = up x n,
d = s x n; tau_strike + = right-lateral).

```bash
python3 toolbox/on_fault_stress_projection/project_stress_to_vtu.py \
    safs_seisol_v2_0_0_RSSRW/safs_initial_stress.yaml \
    safs_seisol_v2_0_0_RSSRW/safs_mesh.puml.h5 \
    --out-dir toolbox/on_fault_stress_projection \
    --P-p-MPa 20 --hypocenter 606971 3707270 -4965.62
```

## nc_to_vtu/nc_to_vtu.py

Convert the ASAGI material netCDF (rectilinear CVM grid) to a ParaView
VTU of hexahedral cells with point data rho_kg_m3, mu_GPa, lambda_GPa and
derived Vs_m_s, Vp_m_s, poisson_ratio. The full SAFS grid gives a
~675 MB file; `--stride 2` subsamples every 2nd node per axis (~1/8 the
size), `--no-derived` drops the derived fields.

```bash
python3 toolbox/nc_to_vtu/nc_to_vtu.py \
    safs_seisol_v2_0_0_RSSRW/safs_material_cvm.nc \
    --out-dir toolbox/nc_to_vtu
```

## generate_velocity_nc_from_raw/generate_velocity_nc_from_raw.py

Full reproduction of the CVM material pipeline in ONE file: raw
`velocity_raw_*.bp` ASCII slices (`raw_data/multiscale_statewise_cvm/`)
-> EPSG:4326 -> UTM 11N reprojection + linear resampling onto the
1500 m rectilinear grid (+100 m surface clone) -> elastic moduli at the
nodes (mu = rho Vs^2, lambda = rho (Vp^2 - 2 Vs^2)) -> uniform-z
resample (dz 250 m, z in [-45000, 0]) -> ASAGI NetCDF with compound
data(z,y,x) {rho; mu; lambda}. Defaults reproduce the shipped
`safs_seisol_v2_0_0_RSSRW/safs_material_cvm.nc` BIT-IDENTICALLY
(verified via `--compare-to`). Needs scipy + pyproj + netCDF4.

```bash
python3 toolbox/generate_velocity_nc_from_raw/generate_velocity_nc_from_raw.py \
    --mesh-puml safs_seisol_v2_0_0_RSSRW/safs_mesh.puml.h5 \
    --compare-to safs_seisol_v2_0_0_RSSRW/safs_material_cvm.nc
# output: toolbox/generate_velocity_nc_from_raw/safs_material_cvm.nc
# options: --dz 50 (keep every source z-level), --write-sidecar X.h5,
#          --grid-dx / --z-min / --z-max / --dtype
```

## hypocenter_facts/hypocenter_facts.py

Print the "hypocenter facts" sheet (sigma_N_eff, |tau|, mu_apparent from
the stress projection; Vs, rho, shear modulus, Poisson ratio from the
CVM material netCDF) at the nucleation location. Imports the projection
module from the sibling folder — keep the toolbox structure intact.

```bash
python3 toolbox/hypocenter_facts/hypocenter_facts.py \
    --json toolbox/hypocenter_facts/hypocenter_facts.json
# material sampling: --sampling trilinear (default; = easi linear) | nearest
```
