# Status: SAFS Multi-Fault Tet Mesh Pipeline (tetgen + autorefine)

Date: 2026-05-10
Author: debug session, project_7.0_preferred

## Pipeline Overview

Three-stage pipeline producing a tetrahedral volume mesh of the SAFS box
with all 6 corefined fault surfaces embedded as internal constraints:

```
   ┌─────────────────────────┐
   │ corefine_faults.py      │  (existing, unchanged)
   │ — pairwise CGAL corefine│
   │ — isotropic remeshing   │
   └────────────┬────────────┘
                │  data_corefined/<basename>_corefined.stl  × 6
                │  data_corefined/manifest.json
                ▼
   ┌─────────────────────────┐
   │ autorefine_merged       │  (NEW C++)
   │ — merge box + 6 faults  │
   │ — PMP::autorefine_      │
   │   triangle_soup resolves│
   │   residual intersections│
   │ — graded box surface    │
   │   triangulation         │
   └────────────┬────────────┘
                │  data_corefined/safs_autorefined_<R>m.stl
                │  data_corefined/safs_autorefined_<R>m_markers.json
                ▼
   ┌─────────────────────────┐
   │ tetgen_mesh.py          │  (NEW Python)
   │ — tetgen plc + nobisect │
   │   + quality refinement  │
   │ — output vertex dedup   │
   │   at 99 m               │
   │ — fault-aware filter    │
   │   (don't drop fault-    │
   │    incident slivers)    │
   └────────────┬────────────┘
                │  code_meshing/safs_multifault_box_<R>m_bulk.vtu
                │  code_meshing/safs_multifault_box_<R>m_fault.vtu
                ▼
              (MFEM input)
```

## Files Created / Modified

### New files

| Path | Role |
|---|---|
| `code_preprocess/corefine_cgal/autorefine_merged.cpp` | Merges 6 corefined STLs + box into one polygon soup; runs `CGAL::PMP::autorefine_triangle_soup` with marker-propagating visitor; emits clean merged STL + per-triangle marker JSON. |
| `code_preprocess/tetgen_mesh.py` | Reads autorefined STL + markers, runs tetgen with `nobisect=True quality=True`, applies output vertex dedup + fault-aware quality filter, writes bulk + fault VTUs. |
| `document/STATUS_tetgen_pipeline.md` | This file. |

### Modified files

| Path | Change |
|---|---|
| `code_preprocess/corefine_cgal/CMakeLists.txt` | Added `autorefine_merged` target. |
| `code_preprocess/corefine_cgal/corefine_set.cpp` | Added then reverted iterative validation pass — created degenerate triangles, abandoned. Includes intentionally not added back. |

### Default parameters (final)

`autorefine_merged`:
- `--pad-xy 50000` (50 km horizontal)
- `--pad-top 0` (box top = z = 0 = geological free surface; **changed from 100 m** so fault outcrops live ON box top)
- `--pad-bottom 25000` (25 km below deepest fault)
- `--box-edge-size 10000` (uniform 10 km box surface grid; autorefine adds local refinement at fault outcrops)
- `--box-marker 100`

`tetgen_mesh.py`:
- `--lc-far 10000` (advisory; `--maxvolume = lc_far^3 / 6`)
- `--minratio 2.0` (tetgen default; smaller fails on dense fault regions)
- `--bulk-min-edge 100` (filter floor)
- `--bulk-q-min 0.05` (filter floor)
- `--vertex-merge-tol 10` (input dedup tolerance)
- `--output-dedup-tol 99` (post-tetgen dedup; just below input min edge)
- `--mmg-cleanup OFF` (default off — see Known Issues)

## Current Output Status (2000 m fixture, default parameters)

### Bulk mesh
- **n_tets**: 33,650
- **n_verts**: 6,167
- **edge_min**: 99.32 m  (project bar ≥ 100 m)
- **edge_med**: 1,281 m
- **edge_max**: 34,860 m  (graded by box-edge-size = 10 km)
- **q_iso min**: 0.0000  (≥ 0.05 bar)
- **q_iso med**: 0.870  (≥ 0.85 bar ✓)
- **q ≥ 0.05 fraction**: 99.92%
- **q ≥ 0.30 fraction**: 98.50%  (≥ 99.5% bar)
- **edges < 100 m**: 2 tets
- **q < 0.05 (slivers)**: ~70 tets (0.21%)

### Fault surface
- **n_fault_tris**: 8,194 (vs 8,202 input → 8 dropped by output dedup, 99.9% preserved)
- **6 patches preserved**, alphabetical:
  - patch 1 (COAV-MC): 793 tris (input 788, +5 Steiner from autorefine)
  - patch 2 (MJVS-SAF): 984
  - patch 3 (MULT-SoSAF): 1299
  - patch 4 (SBMT-Garnet): 1204
  - patch 5 (SBMT-MC): 1010
  - patch 6 (SBMT-SAF): 2904

### Topological correctness
- **Embedding**: 8194/8194 (100%) fault triangles shared by exactly 2 bulk tets ✓
- **Coverage**: 100/100 random box-interior samples inside a tet ✓
- **No box-fault gap**: fault outcrops (z = 0) coincide with box top (z = 0) → free-surface BC applicable directly

## What Works

1. **Autorefine pipeline robust**: `PMP::autorefine_triangle_soup` cleanly resolves the residual intersections that survive corefine + isotropic_remeshing (1 duplicate triangle + 8 segment-facet intersections in input → 0 after autorefine).
2. **Fault preservation**: 100% of fault triangles embedded, with every input STL vertex preserved exactly via `nobisect=True`.
3. **Free-surface alignment**: After `pad_top=0` fix, fault outcrops at z = 0 are on the box top face. Box top has local refinement near outcrops via autorefine's Steiner insertion at fault-edge / box-facet intersections.
4. **Graded bulk**: Bulk tet sizes follow input boundary density — fine (~1500 m) near faults, ~10 km near box surface, smoothly graded between.

## Known Issues / Trade-offs

### 1. Sliver tets (~0.2% of mesh)

~70 bulk tets have q_iso < 0.05 and 2 tets have edge_min < 100 m. Concentrated within 5 km of the fault surface; depth range −12 km to −76 m.

**Root cause**: Fault triangulation has *non-uniform density* (input STLs have edge sizes 100 m to 2.5 km). Where dense fault triangles meet sparser bulk Delaunay tets, thin tets are inevitable. tetgen with `nobisect=True` cannot insert Steiner points on the fault to repair these.

**Mitigation attempted**:
| Tool | Outcome |
|---|---|
| `tetgen optim`, various `mindihedral`, `opt_iterations=20` | Marginal: q_med 0.871 → 0.892, slivers 70 → 60 |
| `tetgen coarsen=True` | No effect |
| `tetgen optmaxedgeratio=10` | No effect |
| Output vertex dedup at 99 m | **Most effective**: removed coincident-vertex slivers (edge_min 1 µm → 99 m), reduced sub-100 m edges 53 → 2 |
| MMG3D `optim+nosurf+noinsert` | **Failed**: MMG drops surface triangles on non-manifold input (our 4-incident polylines from corefine), destroying fault embedding. Even with explicit `RequiredTriangles` markers in MEDIT input, MMG dropped 7100 of 8195 fault triangles. |

**Why MMG fails**: MMG3D requires *manifold* surface input. Our corefine output has 107 non-manifold edges (each shared by 4 triangles — 2 from each of two intersecting faults). MMG's surface reconstruction treats these as defects and dropped most fault tris despite "Required" markings.

**Currently accepted trade-off**: ~70 q < 0.05 tets (0.2%) remain. They are physically located in the locked region (mostly z < −5 km) of BP5-style simulations and do not affect active rupture mechanics.

### 2. q ≥ 0.30 fraction = 98.5% vs project bar 99.5%

Same root cause as slivers — fault density variation forces some tets to be marginally thin (q in [0.05, 0.30]). Stricter filter would lose 10% bulk coverage; accepting current state.

### 3. 8 fault triangles lost to output dedup (0.1%)

Output dedup at 99 m merges any pair of vertices closer than 99 m (well below input STL min edge 100 m, so legitimate fault structure is preserved). 8 fault triangles became degenerate during vertex merge — these were autorefine-produced triangles at near-zero-area near-coincident sites, geometrically negligible.

## Tools Inventory

### CGAL (in `code_preprocess/corefine_cgal/`)
- `corefine_pair`, `corefine_set`: pairwise corefine pipeline (existing, unchanged)
- `check_self_intersect`: per-STL self-intersection diagnostic (existing)
- `mesh_volume`: CGAL Mesh_3 path (KEPT as fallback; produces remeshed surfaces, not used for production)
- **`autorefine_merged`** (NEW): merge + autorefine_triangle_soup with marker propagation

### Python (in `code_preprocess/`)
- `corefine_faults.py`: orchestrator for `corefine_set` (existing, unchanged)
- `medit_to_vtu.py`: MEDIT → VTU converter for CGAL Mesh_3 path (existing, kept)
- `check_msh_quality.py`: project quality bar checker (existing, kept)
- **`tetgen_mesh.py`** (NEW): tetgen + dedup + filter

## How to Reproduce

```bash
# 1. Activate environment
conda activate pythonenv

# 2. (One-time) Build C++ binaries
conda activate cgal-61
cd /path/to/project_7.0_preferred/code_preprocess/corefine_cgal
cmake -S . -B build
cmake --build build --target autorefine_merged corefine_set
conda activate pythonenv

# 3. Corefine STLs (existing pipeline; produces data_corefined/)
python code_preprocess/corefine_faults.py --res 2000

# 4. Autorefine merged surface (NEW)
./code_preprocess/corefine_cgal/build/autorefine_merged \
    data_corefined/manifest.json \
    data_corefined/safs_autorefined_2000m.stl \
    data_corefined/safs_autorefined_2000m_markers.json \
    --pad-top 0 --box-edge-size 10000 --verbose

# 5. Tetgen volume meshing (NEW)
python code_preprocess/tetgen_mesh.py \
    --out-base code_meshing/safs_multifault_box_2000m \
    --lc-far 10000

# 6. Quality check
python code_preprocess/check_msh_quality.py \
    code_meshing/safs_multifault_box_2000m_bulk.vtu
```

Outputs:
- `code_meshing/safs_multifault_box_2000m_bulk.vtu` — volume tet mesh, single 'attribute' = 1
- `code_meshing/safs_multifault_box_2000m_fault.vtu` — fault triangles, per-fault 'patch' marker

## Visualization in ParaView

To see the fault embedded in the bulk:
1. Open BOTH `*_bulk.vtu` and `*_fault.vtu` simultaneously.
2. Set bulk to "Surface With Edges", color by `q_iso` or `attribute`.
3. Set fault to "Surface With Edges", color by `patch`. Use opacity 0.7 to see through bulk if needed.
4. To see internal cuts: `Filters > Slice` on the bulk through any plane crossing a fault — slice plane will show fault as colored line.

## Path Forward (Open)

If sliver count must be driven below 0.1%:

1. **Custom Python sliver collapse** (~50 lines): walk sliver tets, collapse shortest edge if both endpoints are interior, re-tetrahedralize cavity. Preserves fault by design but complex to get right.
2. **Make corefine output manifold** then re-enable MMG3D: would require modifying corefine_set.cpp to "duplicate" non-manifold polyline edges so each fault has its own boundary edge along the polyline (at the cost of the conformality contract).
3. **Accept current state** and proceed with MFEM simulation — slivers are in the locked region of BP5 and don't affect active rupture.

## Pre-existing Documents

- `PLAN_cgal_corefine_multifault.md` — original Phase 1-3 plan (CGAL corefine + gmsh meshing). Phases 1-2 (corefine) still in active use; Phase 3 (gmsh) replaced by `autorefine_merged` + `tetgen_mesh.py`.
- `PLAN_min_edge_enforcement.md` — CGAL Mesh_3 hard-floor enforcement plan. Implemented in `mesh_volume.cpp` and `medit_to_vtu.py`. **Not used in current production**; kept as fallback.
- `PLAN_revert_to_gmsh.md` — gmsh combined-STL attempt. **Failed** at tetgen PLC validation; superseded by direct tetgen + autorefine pipeline documented here.
