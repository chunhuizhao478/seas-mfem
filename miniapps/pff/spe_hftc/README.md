# spe_hftc: dry acrylic cube under a 20-pulse train

One TOML per MOOSE case of farms branch `cdms` in `configs/`. The header of each
file names the MOOSE input it reproduces, and all values are copied from that
input (and from the job script's overrides where a production run used them).
Run from `miniapps/pff`: `pff_driver -c spe_hftc/configs/<case>.toml`.

| File | MOOSE input | Notes |
| --- | --- | --- |
| `configs/cube_homog_eta0p45_absorb.toml` | `pulsepower/spe_hftc/code_3d/elasticity_absorbing.i` | Vista/LS6 3D production: η = 0.45, ζ_pulse = 1e-8, dt_quiet = 5e-7 |
| `configs/section_homog_eta0p45_absorb.toml` | `pulsepower/spe_hftc/code_2d/elasticity_absorbing.i` | Vista 2D homogeneous η = 0.45 |

Planned subfolders, as in seas benchmark folders: `mesh/` (.geo scripts),
`gold/<run>_jobid<N>_git<hash>/` (reference runs), `benchmark_data/` (MOOSE
CSVs used by verification tests), `benchmark_document/` (write-ups).

Boundaries are given by Gmsh physical name (spe_hftc, as MOOSE uses them) or
by number (cmame_revision). MFEM's Gmsh reader needs format 2.2: generate with
`msh_version="2.2"` (spe_hftc mesh scripts) or `gmsh -format msh22`. Mesh
files stay out of git, as in farms.

