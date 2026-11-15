# cmame_revision: granite disc and cylinder, solid and hydromech

One TOML per MOOSE case of farms branch `cdms` in `configs/`. The header of each
file names the MOOSE input it reproduces, and all values are copied from that
input (and from the job script's overrides where a production run used them).
Run from `miniapps/pff`: `pff_driver -c cmame_revision/configs/<case>.toml`.

| File | MOOSE input | Notes |
| --- | --- | --- |
| `configs/2d_solid_benchmark.toml` | `pulsepower/cmame_revision/2d_solid/benchmark/elasticity.i` | two seeded slots (d_lb = 0.9) |
| `configs/3d_solid_em_5J.toml` | `pulsepower/cmame_revision/3d_solid/em_5J/elasticity.i` | dashpots on outer, top, bottom |
| `configs/2d_hydromech_benchmark.toml` | `pulsepower/cmame_revision/2d_hydromech/benchmark/elasticity_E1d25.i` | undrained hole, Darcy–Poiseuille permeability |
| `configs/3d_hydromech_em_5J.toml` | `pulsepower/cmame_revision/3d_hydromech/em_5J/elasticity_mesh2x.i` | normal-strain (Heider 2021) permeability, p0 = 0.0965 MPa |

Planned subfolders, as in seas benchmark folders: `mesh/` (.geo scripts),
`gold/<run>_jobid<N>_git<hash>/` (reference runs), `benchmark_data/` (MOOSE
CSVs used by verification tests), `benchmark_document/` (write-ups).

Boundaries are given by Gmsh physical name (spe_hftc, as MOOSE uses them) or
by number (cmame_revision). MFEM's Gmsh reader needs format 2.2: generate with
`msh_version="2.2"` (spe_hftc mesh scripts) or `gmsh -format msh22`. Mesh
files stay out of git, as in farms.

Still to add: the 2D/3D variants (Weibull and Gaussian Gc, free faces, the
`permeability_formula/` and drained/undrained studies, the other `em_*J`
energies). They differ from these files in a few keys.
