# UCERF 2-D map view: San Andreas and San Jacinto

Input files for `config/fd_cycle/ucerf2d_saf_sjf_fd.toml` (run8, R 10 km
fillets) and `config/fd_cycle/ucerf2d_saf_sjf_500m_R70_avg_tip_fd.toml`
(R 70 km smoothed traces), plane-strain map-view cycles (QD -> dynamic ->
QD) on the UCERF3 FM3.1 traces cut at 10 km depth.  Setup and choices: `/mnt/project-files/ucerf-2d/` (report
`ucerf-faults-2d-rev2.pdf` and later revisions).

## Files

| file | what |
|---|---|
| `mesh/ucerf2d_saf_sjf_500m_tilted_D415km_f10.msh` | gmsh 2.2, 163,524 nodes, 326,130 triangles; tags 3 faults, 5 the two long walls, 1 the two ends |
| `mesh/ucerf2d_saf_sjf_500m_tilted_D415km_f10.json` | the builder's record of that mesh: 18 sections, 822 km of fault, 1,644 fault edges, 10 tips, box 1393 x 893 km tilted to azimuth 123.0 deg, corner fillets R = 10 km (49 corners, 4 with a smaller R, traces moved at most 0.44 km) |
| `mesh/ucerf2d_saf_sjf_500m_tilted_D415km_R70.msh` | gmsh 2.2, 163,170 nodes, 325,422 triangles, same tags; the traces smoothed to a minimum radius R = 70 km, 1,645 fault edges, h_min 194 m, at most 0.43 deg per fault vertex |
| `mesh/ucerf2d_saf_sjf_500m_tilted_D415km_R70.json` | the builder's record of that mesh: 18 sections in 5 chains, 822.6 km of fault, 10 tips, box 1394 x 894 km tilted to azimuth 123.0 deg, traces moved at most 6.24 km (San Gorgonio Pass; saf1 at most 2.74 km) |
| `petsc_ts_rk45_ucerf.cfg` | QD time stepping: rk 5dp, atol 1e-7 (infinity norm), dt_max 10 yr |

sha256 of the meshes:
- `_f10.msh` `9f6dd56245f3a625948e0ababa59dcbe5480f34415cbff763f7e1a7041e6f4ed`
- `_R70.msh` `d3df54b3379f5042d79897a593a91af419647c6d5526c000f57cb67b2e4a7dbd`

## How the mesh was made

`/mnt/project-files/ucerf-2d/rounded-mesh-2026-10-29/build_mesh.py`
(sha256 `eedcd7aba642a0508e885ec6e05dc47a6ef20112c171ea96ebbe41b011daa9c7`),
set `typeA`, 500 m on the faults, a box tilted along the San Andreas with its
walls 415 km from every fault point, and `--fillet 10e3` (Chunhui Zhao,
2026-10-29: round the trace corners).  It reads the 10 km slice
`slice_FM3_1.json` and `analyze.py` from `/mnt/project-files/ucerf-2d/`.

The R 70 km mesh (Chunhui Zhao, 2026-11-01: "we can smooth the trace
bends if needed"; the radius from the simple cases D6f-D10f, where the bend
stress at an arc end is about mu' D / R):
`/mnt/project-files/ucerf-2d/smoothed-mesh-2026-11-01/build_mesh_smooth.py`
(sha256 `6864180017a8f68142fdaf6db13d532c28715871bd57d3ec761dd2b272ad8b2d`),
run in that folder, which holds its inputs (`slice_FM3_1.json`, `analyze.py`,
`parse.py`) and a README with the method:

    python3 build_mesh_smooth.py --set typeA --drop "Garlock,Elsinore,Lytle Creek" \
      --lcf 500 --grade 0.1 --d0 10e3 --hmax 5e3 --tilt "San Andreas" \
      --pad 415e3 --pad_along 415e3 --fillet 0 --smooth_R 70e3 \
      --out ucerf2d_saf_sjf_500m_tilted_D415km_R70.msh

`--smooth_R` replaces the tangent angle of each joined chain by the mean of
its lower and upper 1/R-Lipschitz envelopes, so it turns at most 1/R per unit
length, then resamples every 500 m.  gmsh meshes on 4 threads, so a rebuild
gives the same traces but not the same triangles; the committed file is the
record.

## Runs

Long runs go on the Mac, from a pushed commit, with the output on its own
disk.  From `miniapps/seas`:

    mpirun -np 4 ./seas_fd_cycle_driver --deck config/fd_cycle/ucerf2d_saf_sjf_fd.toml
    mpirun -np 8 ./seas_fd_cycle_driver --deck config/fd_cycle/ucerf2d_saf_sjf_500m_R70_avg_tip_fd.toml
