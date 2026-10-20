# TPV205 gold: RK45 mixed flux p1, 200 m (Frontera jobs 7762879 + 7763019)

- Code: seas-mfem commit 04674fc6d5ee2d14de89237d24e90a10bc8b99d3 (04674fc, on GitHub in origin/nitsche-benchmarks);
  driver seas_spatial_dyn_driver via ibrun. Built, as the sbatch says, with `make seas_spatial_dyn_driver` in
  miniapps/seas (the binary actually used is not recorded).
- Deck: jobs/tpv205_p1_rk45_mixedflux_8N_400r_dev.sbatch (same file as in the repo at 04674fc) with
  tpv205/configs/tpv205_spatial_rk45_mixedflux.toml. Linear slip weakening (mu_s 0.677, mu_d 0.525, d_c 0.40 m); tau
  70 MPa, sigma_n 120 MPa, 3 x 3 km patches at 7.5 km depth of 81.6 MPa (x = 0), 78.0 MPa (x = -7.5 km) and 62.0 MPa
  (x = +7.5 km); barrier outside |x| = 15 km and below 15 km. Numerics: order 1, interior_flux "scalar", mixed_flux
  "adjacent", time_integrator "rk45", cfl 0.25 (cfl_safety "dg"), auto dt. The sbatch adds --deriv-cache, --tfinal 12
  (TOML: 0.2 s) and --mesh (same 200 m mesh as the TOML).
- Mesh: tpv205/mesh/tpv2053d_200m.msh, not in the repo (the sbatch runs the gmsh command below if it is missing);
  tpv205/mesh/tpv2053d_200m.geo at 04674fc: lc_fault and lc_nucl 200 m, refined around the three patches, msh 2.2.
  Element count, sha256 and gmsh version: not recorded.
- Run: Frontera (named in the sbatch), partition flex, account EAR20006, 16 nodes x 800 ranks, wall limit 01:00:00;
  the file name (8N_400r_dev) and its header ("development queue, 8N/400r/2 h") do not match these #SBATCH lines.
  Jobs: the folder name gives 7762879 and 7763019 ("combined"). The sbatch runs one segment per submission and resumes
  from a checkpoint with TPV205_RESTART, which fits a first job plus a restart, but which job wrote which time range
  and how the station files were joined is not recorded; the files show no repeated or backward time stamps. Run dates
  and wall time used: not recorded (run commit dated 2026-06-03).
- results/: the 16 SCEC station files tpv205_station_x2_<strike km>_x3_<depth km>.dat, 15,977 rows, 7.512e-4 s apart,
  t = 0 to 12 s. No units in the header; values are SI (h-shear-stress 7.0e7 Pa at t = 0). First lines of
  tpv205_station_x2_-12_x3_0.dat:
  ```text
  # TPV205 station: x2_-12_x3_0 (along_strike=-12000, down_dip=0)
  # Columns (SCEC TPV5 trace layout, LSW): t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate  v-shear-stress  n-stress  mu_eff
  # BP5 convention: h = strike = component 2, v = dip = component 1.  TPV205 is right-lateral strike-slip; v-* channels remain near 0.  mu_eff is the LSW μ(δ) at the current slip magnitude δ = sqrt(slip1² + slip2²).
  ```
- plots/: tpv205_overview.png and tpv205_x2_<..>_x3_<..>.png for the 16 stations (plot script not recorded).

## Reproduce

```bash
git checkout 04674fc6d5ee2d14de89237d24e90a10bc8b99d3
# MFEM itself must be built with PETSc + MUMPS (config/config.mk; see build_frontera.sh)
(cd miniapps/seas && make seas_spatial_dyn_driver)
(cd miniapps/seas && gmsh -format msh22 -3 tpv205/mesh/tpv2053d_200m.geo -o tpv205/mesh/tpv2053d_200m.msh)
sbatch miniapps/seas/jobs/tpv205_spatial/tpv205_p1_rk45_mixedflux_8N_400r_dev.sbatch
TPV205_RESTART=tpv205/out_rk45_mf_<first jobid>/cp sbatch miniapps/seas/jobs/tpv205_spatial/tpv205_p1_rk45_mixedflux_8N_400r_dev.sbatch
```

Output: miniapps/seas/tpv205/out_rk45_mf_<jobid>/tpv205_station_*.dat; the second command continues a run cut by the 1
h wall into a new out_rk45_mf_<jobid> (same -n 800). TPV205_CONFIG / _MESH / _OUT / _TFINAL override the defaults. On
another cluster adjust the module lines, `LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:...`, ibrun, -p flex and -A EAR20006; the
sbatch has no absolute paths.

Gold folder added in f127ef9 (2026-06-04); README added 2026-10-17 by Chunhui Zhao's request (trace and reproduce every gold from its commit).
