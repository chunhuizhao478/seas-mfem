# project_7.0_alternative — three-pillar layout

The SAFS alternative pipeline (NW-cut single-fault representation) is
organised around three independent data streams, each laid out the same
way: **raw** inputs → **code** that processes them → **results** the
downstream consumers read.

```
project_7.0_alternative/
├── meshing/        fault + box mesh                        (CFM .ts → STL → Gmsh .msh)
├── velocity/       CVM-H velocity → projected sidecar      (.bp → .h5/.vtr)
├── stress/         community stress model → on-fault       (.csv → stress_safs.h5)
└── document/       cross-cutting plans and reviews

each pillar:
  raw/       inputs from external sources (CFM, CVM-H, CSM)
  code/      Python scripts that turn raw → results
  docs/      pillar-specific plans, explore notes, review docs
  results/   artifacts consumed downstream (by gmsh, by the SEAS C++ runtime)
```

## Data flow

```
   meshing/raw/*.ts
         │  ts_to_stl.py
         ▼
   meshing/raw/*_unclipped.stl
         │  clean_freesurface_mesh.py
         ▼
   meshing/results/stl_cleaned/*_clean_clip.stl
         │  nw_cut_strip.py        (anchor from project_7.0_preferred mesh)
         ▼
   meshing/results/stl_nwcut/*_clean_clip_nwcut.stl
         │  run_nwcut_meshing.py   (gmsh + size_field.pos from velocity/)
         ▼
   meshing/results/{msh,vtu}/safs_fault_box_nwcut_{500,1000,2000}m.{msh,*.vtu}
                │                                              │
                │                                              │
   velocity/raw/velocity_raw_*.bp                              │
         │  build_velocity_cvmh.py  ◄──── mesh-bbox guard ─────┘
         ▼
   velocity/results/velocity_safs.{h5,vtr}
         │  build_size_field.py ──► meshing/results/msh/*.size_field.pos
         │
   stress/raw/CSM_Johnson_Hearn_raw_*km.csv
         │  build_stress_safs.py     (uses meshing/results/vtu/*_bulk.vtu)
         ▼
   stress/results/stress_safs.h5  +  per-resolution {500m,1000m,2000m}/
```

## Where to run scripts

- Meshing pipeline (Gmsh):     `conda activate pythonenv; cd meshing/code`
- Velocity projection:         `conda activate pythonenv; cd velocity/code`
- Stress projection + tests:   `conda activate pythonenv; cd stress/code`

Run tests for each pillar from inside that pillar's `code/` directory:
`pytest -q`.

## Python imports

The `code/` directories are flat — there is no `data_projection`
package any more. Scripts import sibling modules directly
(`from crs import geographic_to_utm11n`). Out-of-pillar imports add the
target `code/` dir to `sys.path` first (see `meshing/code/run_nwcut_meshing.py`
for the pattern that pulls `velocity/code/{crs,build_size_field}`).
