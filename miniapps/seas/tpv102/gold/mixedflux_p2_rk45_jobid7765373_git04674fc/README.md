# TPV102 gold: RK45 mixed flux p2, 200 m (Frontera job 7765373)

- Code: seas-mfem commit 04674fc6d5ee2d14de89237d24e90a10bc8b99d3 (04674fc, on GitHub in origin/nitsche-benchmarks);
  driver seas_spatial_dyn_driver via ibrun. Built, as the sbatch says, with `make seas_spatial_dyn_driver` in
  miniapps/seas (the binary actually used is not recorded).
- Deck: jobs/tpv102_p2_rk45_mixedflux_12N_600r_normal.sbatch (same file as in the repo at 04674fc) with
  tpv102/configs/tpv102_spatial_p2_rk45_mixedflux.toml. Rate-and-state aging law (a 0.008 / b 0.012 in the core, a
  0.016 outside 30 x 15 km, Dc 0.02 m, f0 0.6); tau 75 MPa, sigma_n 120 MPa; 25 MPa gradual overstress over 1 s at 7.5
  km depth. Numerics: order 2, interior_flux "scalar", mixed_flux "adjacent", time_integrator "rk45", cfl 0.25
  (cfl_safety "dg"), auto dt. The sbatch adds --deriv-cache, --tfinal 12 (TOML: 0.2 s) and --mesh (TOML default: the
  1000 m mesh).
- Mesh: tpv102/mesh/tpv102_200m.msh, not in the repo (the sbatch runs the gmsh command below if it is missing);
  tpv102/mesh/tpv102_200m.geo at 04674fc: lc_fault 200 m, far field 5 km, msh 2.2. Element count, sha256 and gmsh
  version: not recorded.
- Run: Frontera (named in the sbatch), partition normal, account EAR20006, 12 nodes x 600 ranks, wall limit 48:00:00,
  job 7765373. Run date, wall time used and restart segments: not recorded (run commit dated 2026-06-03).
- results/: the 9 SCEC station files tpv102_station_flt_<strike km>_<depth km>.dat (n = negative), 26,505 rows,
  4.528e-4 s apart, t = 0 to 12 s. No units in the header; values are SI (tau2 7.5e7 Pa at t = 0). First lines of
  tpv102_station_flt_0_12.dat:
  ```text
  # TPV102 station: flt_0_12 (along_strike=0, down_dip=12000)
  # Columns: time slip1 slip2 V1 V2 tau1 tau2 sigma_n log10_theta
  # BP5 convention: 1 = dip (tangent1), 2 = strike (tangent2). TPV102 is pure strike-slip so V1/slip1/tau1 ~ 0 and the interesting physics is in column 2.
  ```
- plots/: tpv102_overview.png and tpv102_flt_<station>.png for the 9 stations (plot script not recorded).

## Reproduce

```bash
git checkout 04674fc6d5ee2d14de89237d24e90a10bc8b99d3
# MFEM itself must be built with PETSc + MUMPS (config/config.mk; see build_frontera.sh)
(cd miniapps/seas && make seas_spatial_dyn_driver)
(cd miniapps/seas && gmsh -format msh22 -3 tpv102/mesh/tpv102_200m.geo -o tpv102/mesh/tpv102_200m.msh)
sbatch miniapps/seas/jobs/tpv102_spatial/tpv102_p2_rk45_mixedflux_12N_600r_normal.sbatch
```

Output: miniapps/seas/tpv102/out_p2_rk45_mf_normal_<jobid>/tpv102_station_*.dat. TPV102_CONFIG / _MESH / _OUT /
_TFINAL override the defaults; TPV102_RESTART=<previous OUT>/cp continues past the wall limit (same -n 600). On
another cluster adjust the module lines, `LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:...`, ibrun, -p normal and -A EAR20006;
the sbatch has no absolute paths.

Gold folder added in f127ef9 (2026-06-04); README added 2026-10-17 by Chunhui Zhao's request (trace and reproduce every gold from its commit).
