# NW-Cut Volume Mesh — ALT6 Fault inside an Adjusted Box

**Date generated:** 2026-05-08
**Mesher:** gmsh 4.15.0 (in `pythonenv`)
**Driver:** `code_meshing/run_nwcut_meshing.py`
**Template:** `code_meshing/safs_fault_box_nwcut.geo`
**Inputs:** `data_cutnwfault/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_{500,1000,2000}m_clean_clip_nwcut.stl`
**Outputs:** `code_meshing/safs_fault_box_nwcut_{500,1000,2000}m.{msh, _bulk.vtu, _fault.vtu}`

## Adjusted bounding box (per resolution, fitted to the cut trace)

The bounding box is now derived from the *cut* fault's bbox — not the
uncut ALT6 strip's bbox — so the lateral extent shrunk from ~443 km × ~296 km
to ~258 km × ~147 km. Padding around the fault matches
`safs_fault_box_freesurface_clip.geo` exactly:

| Pad parameter | Value | Purpose |
|---|---|---|
| `PAD_XY`     | 50 km   | horizontal margin (absorbing far-field) |
| `PAD_TOP`    | 100 m   | thin overburden — keeps fault top edge strictly interior to box top face |
| `PAD_BOTTOM` | 25 km   | rock below deepest fault vertex |

**Cut-fault bbox per resolution** (input to the box):

| Resolution | x [m] | y [m] | z [m] | Width × Height × Depth |
|---|---|---|---|---|
| 500 m  | `[364865.58, 622329.50]` | `[3692278.00, 3838985.00]` | `[-16607.61, 0.00]` | `257.46 × 146.71 × 16.61 km` |
| 1000 m | `[364830.83, 622329.50]` | `[3692278.00, 3838944.94]` | `[-16607.43, 0.00]` | `257.50 × 146.67 × 16.61 km` |
| 2000 m | `[364748.63, 622329.50]` | `[3692278.00, 3838955.11]` | `[-16607.45, 0.00]` | `257.58 × 146.68 × 16.61 km` |

**Final padded box** (each input bbox + the pads above): all three meshes
land in a domain of roughly `357 × 247 × 41.7 km`. Compared to the
existing `safs_fault_box_freesurface_clip.geo` at uncut `543 × 396 × 41.7 km`,
the new box is **~57 % the volume** of the old one.

## Mesh statistics

Mesh-size controls (identical to existing `safs_fault_box_freesurface_clip.geo`):

| Field | Value |
|---|---|
| `LC_NEAR`     | 1500 m |
| `LC_FAR`      | 10000 m |
| `DIST_INNER`  | 3000 m |
| `DIST_OUTER`  | 40000 m |
| `Mesh.Algorithm`   | 6 (Frontal-Delaunay 2-D) |
| `Mesh.Algorithm3D` | 1 (Delaunay 3-D) |

### 2000 m fault resolution

- Fault triangles in the .msh: **2 697**
- Bulk tetrahedra: **126 422**
- gmsh wall-clock: **2.6 s** (-nt 7)
- .msh size: **5.36 MB**  (compare: uncut clipped 2000m .msh = 10.09 MB)
- Bulk edge length [m]: **min = 77.47**, median = 2 821.31, mean = 4 450.82, max = 21 214.22
- Fault edge length [m]: **min = 77.47**, median = 2 000.80, mean = 1 990.09, max = 3 168.66
- Tet quality (isoperimetric, 1 = regular tet): **min = 0.0978**, mean = 0.7128, max = 0.7937
- Quality bins: `q < 0.1`: **1 tet (0.00 %)** ; `0.1 ≤ q < 0.3`: 638 tets (0.50 %) ; `q ≥ 0.3`: 125 783 tets (**99.50 %**)

### 1000 m fault resolution

- Fault triangles in the .msh: **10 620**
- Bulk tetrahedra: **172 109**
- gmsh wall-clock: **4.1 s** (-nt 7)
- .msh size: **7.38 MB**
- Bulk edge length [m]: **min = 8.64**, median = 2 226.44, mean = 3 616.40, max = 21 107.50
- Fault edge length [m]: **min = 8.64**, median = 1 000.41, mean = 1 001.09, max = 1 491.35
- Tet quality: **min = 0.0665**, mean = 0.7060, max = 0.7937
- Quality bins: `q < 0.1`: **5 tets (0.00 %)** ; `0.1 ≤ q < 0.3`: 442 tets (0.26 %) ; `q ≥ 0.3`: 171 662 tets (**99.74 %**)

### 500 m fault resolution

- Fault triangles in the .msh: **42 389**
- Bulk tetrahedra: **317 914**
- gmsh wall-clock: **11.5 s** (-nt 7)
- .msh size: **13.99 MB**
- Bulk edge length [m]: **min = 9.85**, median = 1 212.50, mean = 2 308.72, max = 20 500.69
- Fault edge length [m]: **min = 9.85**, median = 500.34, mean = 500.94, max = 1 091.34
- Tet quality: **min = 0.0643**, mean = 0.6946, max = 0.7937
- Quality bins: `q < 0.1`: **13 tets (0.00 %)** ; `0.1 ≤ q < 0.3`: 519 tets (0.16 %) ; `q ≥ 0.3`: 317 382 tets (**99.83 %**)

## Cross-resolution observations

- **All three meshes are valid:** every resolution produces a closed
  manifold, the fault is correctly embedded as a 2-D constraint
  (Physical Surface tag 101), and the box top / bottom / sides receive
  Physical Surface tags 102 / 103 / 104. `Physical Volume("rock", 1)`
  contains every tet.
- **Mean tet quality is consistent (~0.70)** across resolutions, which
  is what gmsh's Delaunay 3-D mesher typically delivers with an
  embedded discrete surface — this matches what the existing
  `safs_fault_box_freesurface_clip.{msh, vtu}` shows.
- **Min tet edge size shrinks at finer resolution** (77.5 m → 8.6 m →
  9.9 m). The 8–10 m mins at 1000 m and 500 m come from sliver fault
  triangles created by the NW cut: the linear-interpolation clipper
  produces clipped vertices very close to surviving keep-side vertices
  whenever the drop vertex was barely on the drop side. The `% of bad
  tets (q < 0.1)` is still essentially zero (≤ 0.005 % across all
  resolutions), so the slivers are isolated and do not propagate. If
  these slivers become a CFL bottleneck, run a one-shot pymeshlab
  edge-collapse pass on the cut STL to weld vertices closer than ~10 m
  before remeshing the volume.
- **Tet count scales sub-linearly with fault triangle count** (2 697 /
  126k vs 10 620 / 172k vs 42 389 / 318k tets), because the bulk mesh
  size is set by the same `LC_NEAR=1500 m`/`LC_FAR=10000 m` controls
  for every resolution.

## Reproducing

```bash
conda activate pythonenv          # gmsh 4.15.0, meshio 5.3.5, pymeshlab
cd /Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_alternative/code_meshing
python run_nwcut_meshing.py                     # all three resolutions
python run_nwcut_meshing.py --res 2000          # smoke run only
```

The wrapper looks up `gmsh` next to the running `python` (so an
explicit `conda activate pythonenv` from your shell is enough — no
shell-side `PATH` munging is required).

## Files written

```
code_meshing/safs_fault_box_nwcut.geo               (template, 6.3 KB)
code_meshing/run_nwcut_meshing.py                   (driver, 7.3 KB)
code_meshing/safs_fault_box_nwcut_2000m.msh         (5.36 MB)
code_meshing/safs_fault_box_nwcut_2000m_bulk.vtu    (3.69 MB, has `quality` cell-data)
code_meshing/safs_fault_box_nwcut_2000m_fault.vtu   (0.62 MB)
code_meshing/safs_fault_box_nwcut_1000m.msh         (7.38 MB)
code_meshing/safs_fault_box_nwcut_1000m_bulk.vtu    (5.01 MB)
code_meshing/safs_fault_box_nwcut_1000m_fault.vtu   (0.85 MB)
code_meshing/safs_fault_box_nwcut_500m.msh          (13.99 MB)
code_meshing/safs_fault_box_nwcut_500m_bulk.vtu     (9.29 MB)
code_meshing/safs_fault_box_nwcut_500m_fault.vtu    (1.64 MB)
```

The `_bulk.vtu` carries the per-tet `quality` field (isoperimetric tet
quality, computed by `code_preprocess/msh_to_vtu.py`). Open it in
ParaView and threshold by `quality` to localise any sliver tet visually.
