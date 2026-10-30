# UCERF 2-D map view: San Andreas and San Jacinto

Input files for `config/fd_cycle/ucerf2d_saf_sjf_fd.toml`, a plane-strain
map-view cycle (QD -> dynamic -> QD) on the UCERF3 FM3.1 traces cut at
10 km depth.  Setup and choices: `/mnt/project-files/ucerf-2d/` (report
`ucerf-faults-2d-rev2.pdf` and later revisions).

## Files

| file | what |
|---|---|
| `mesh/ucerf2d_saf_sjf_500m_tilted_D415km_f10.msh` | gmsh 2.2, 163,524 nodes, 326,130 triangles; tags 3 faults, 5 the two long walls, 1 the two ends |
| `mesh/ucerf2d_saf_sjf_500m_tilted_D415km_f10.json` | the builder's record of that mesh: 18 sections, 822 km of fault, 1,644 fault edges, 10 tips, box 1393 x 893 km tilted to azimuth 123.0 deg, corner fillets R = 10 km (49 corners, 4 with a smaller R, traces moved at most 0.44 km) |
| `petsc_ts_rk45_ucerf.cfg` | QD time stepping: rk 5dp, atol 1e-7 (infinity norm), dt_max 10 yr |

sha256 of the mesh:
`9f6dd56245f3a625948e0ababa59dcbe5480f34415cbff763f7e1a7041e6f4ed`

## How the mesh was made

`/mnt/project-files/ucerf-2d/rounded-mesh-2026-10-29/build_mesh.py`
(sha256 `eedcd7aba642a0508e885ec6e05dc47a6ef20112c171ea96ebbe41b011daa9c7`),
set `typeA`, 500 m on the faults, a box tilted along the San Andreas with its
walls 415 km from every fault point, and `--fillet 10e3` (Chunhui Zhao,
2026-10-29: round the trace corners).  It reads the 10 km slice
`slice_FM3_1.json` and `analyze.py` from `/mnt/project-files/ucerf-2d/`.

## Runs

Long runs go on the Mac, from a pushed commit, with the output on its own
disk.  From `miniapps/seas`:

    mpirun -np 4 ./seas_fd_cycle_driver --deck config/fd_cycle/ucerf2d_saf_sjf_fd.toml
