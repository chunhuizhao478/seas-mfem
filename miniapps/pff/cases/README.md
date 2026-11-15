# Cases

One TOML per MOOSE case of farms branch `cdms`. The header of each file names
the MOOSE input it reproduces, and all values are copied from that input (and
from the job script's overrides where a production run used them).

| File | MOOSE input | Notes |
| --- | --- | --- |
| `spe_hftc/cube_homog_eta0p45_absorb.toml` | `spe_hftc/code_3d/elasticity_absorbing.i` | Vista/LS6 3D production: η = 0.45, ζ_pulse = 1e-8, dt_quiet = 5e-7 |
| `spe_hftc/section_homog_eta0p45_absorb.toml` | `spe_hftc/code_2d/elasticity_absorbing.i` | Vista 2D homogeneous η = 0.45 |
| `cmame_revision/2d_solid_benchmark.toml` | `cmame_revision/2d_solid/benchmark/elasticity.i` | two seeded slots (d_lb = 0.9) |
| `cmame_revision/3d_solid_em_5J.toml` | `cmame_revision/3d_solid/em_5J/elasticity.i` | dashpots on outer, top, bottom |
| `cmame_revision/2d_hydromech_benchmark.toml` | `cmame_revision/2d_hydromech/benchmark/elasticity_E1d25.i` | undrained hole, Darcy–Poiseuille permeability |
| `cmame_revision/3d_hydromech_em_5J.toml` | `cmame_revision/3d_hydromech/em_5J/elasticity_mesh2x.i` | normal-strain (Heider 2021) permeability, p0 = 0.0965 MPa |

Boundaries are given by Gmsh physical name (spe_hftc, as MOOSE uses them) or
by number (cmame_revision). MFEM's Gmsh reader needs format 2.2: generate with
`msh_version="2.2"` (spe_hftc mesh scripts) or `gmsh -format msh22`. Mesh
files stay out of git, as in farms.

Still to add: the 2D/3D variants (Weibull and Gaussian Gc, free faces, the
`permeability_formula/` and drained/undrained studies, the other `em_*J`
energies). They differ from these files in a few keys.
