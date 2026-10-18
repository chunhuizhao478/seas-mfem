# TPV104 gold: RK45 mixed flux p2, 200 m (Frontera job 7765380)

- Code: seas-mfem commit 04674fc6d5ee2d14de89237d24e90a10bc8b99d3 (04674fc, on GitHub in origin/nitsche-benchmarks);
  driver seas_spatial_dyn_driver via ibrun. Built, as the sbatch says, with `make seas_spatial_dyn_driver` in
  miniapps/seas (the binary actually used is not recorded).
- Deck: jobs/tpv104_p2_rk45_mixedflux_12N_600r_normal.sbatch (same file as in the repo at 04674fc) with
  tpv104/configs/tpv104_spatial_p2_rk45_mixedflux.toml. Rate-and-state slip law with strong rate weakening (a 0.01 / b
  0.014 and V_w 0.1 m/s in the core, a 0.02 and V_w 1.0 m/s outside, f_w 0.2, Dc 0.4 m, f0 0.6); tau 40 MPa, sigma_n
  120 MPa; 45 MPa gradual overstress over 1 s at 7.5 km depth. Numerics: order 2, interior_flux "scalar", mixed_flux
  "adjacent", time_integrator "rk45", cfl 0.25 (cfl_safety "dg"), auto dt. The sbatch adds --deriv-cache, --tfinal 12
  (TOML: 0.2 s) and --mesh (TOML default: the 1000 m mesh).
- Mesh: tpv104/mesh/tpv104_200m.msh, not in the repo (the sbatch runs the gmsh command below if it is missing);
  tpv104/mesh/tpv104_200m.geo at 04674fc: a copy of tpv102_200m.geo per its header, lc_fault 200 m, msh 2.2. Element
  count, sha256 and gmsh version: not recorded.
- Run: Frontera (named in the sbatch), partition normal, account EAR20006, 12 nodes x 600 ranks, wall limit 48:00:00,
  job 7765380. Run date, wall time used and restart segments: not recorded (run commit dated 2026-06-03).
- results/: the 9 SCEC station files tpv104_station_x2_<strike km>_x3_<depth km>.dat, 26,505 rows, 4.528e-4 s apart, t
  = 0 to 12 s. No units in the header; values are SI (h-shear-stress 4.0e7 Pa at t = 0). First lines of
  tpv104_station_x2_-12_x3_12.dat:
  ```text
  # TPV104 station: x2_-12_x3_12 (along_strike=-12000, down_dip=12000)
  # Columns (SCEC TPV104 trace layout): t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate  v-shear-stress  n-stress  psi
  # BP5 convention: h = strike = component 2, v = dip = component 1.  TPV104 is pure strike-slip so v-slip / v-slip-rate / v-shear-stress ≈ 0.
  ```
- plots/: tpv104_overview.png and tpv104_x2_<..>_x3_<..>.png for the 9 stations (plot script not recorded).

## Reproduce

```bash
git checkout 04674fc6d5ee2d14de89237d24e90a10bc8b99d3
# MFEM itself must be built with PETSc + MUMPS (config/config.mk; see build_frontera.sh)
(cd miniapps/seas && make seas_spatial_dyn_driver)
(cd miniapps/seas && gmsh -format msh22 -3 tpv104/mesh/tpv104_200m.geo -o tpv104/mesh/tpv104_200m.msh)
sbatch miniapps/seas/jobs/tpv104_spatial/tpv104_p2_rk45_mixedflux_12N_600r_normal.sbatch
```

Output: miniapps/seas/tpv104/out_p2_rk45_mf_normal_<jobid>/tpv104_station_*.dat. TPV104_CONFIG / _MESH / _OUT /
_TFINAL override the defaults; TPV104_RESTART=<previous OUT>/cp continues past the wall limit (same -n 600). On
another cluster adjust the module lines, `LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:...`, ibrun, -p normal and -A EAR20006;
the sbatch has no absolute paths.

Gold folder added in f127ef9 (2026-06-04); README added 2026-10-17 by Chunhui Zhao's request (trace and reproduce every gold from its commit).
