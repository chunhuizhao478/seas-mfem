# TPV102 gold: ADER O2 upwind p1, symmetric mesh, 200 m (Frontera job 7774565)

- Code: seas-mfem commit b09dbfcc14e5000d666fa0be59b6c13621c889c1 (b09dbfc, on GitHub in origin/nitsche-benchmarks);
  driver seas_spatial_dyn_driver via ibrun. Built, as the sbatch says, by jobs/spatial_dyn_build.sbatch clean-rebuild
  (deletes the driver's objects, then `make seas_spatial_dyn_driver -j48`; the binary actually used is not recorded).
- Deck: jobs/tpv102_p1_pureupwind_200m_hybrid_normal.sbatch (same file as in the repo at b09dbfc) with
  tpv102/configs/tpv102_spatial.toml. Rate-and-state aging law (a 0.008 / b 0.012 in the core, a 0.016 outside 30 x 15
  km, Dc 0.02 m, f0 0.6); tau 75 MPa, sigma_n 120 MPa; 25 MPa gradual overstress over 1 s at 7.5 km depth. Numerics:
  order 1, ader_order 2, interior_flux "scalar", fault_iterator "substep", cfl 0.25 (cfl_safety "dg"). The TOML says
  mixed_flux "adjacent"; the sbatch overrides it with --mixed-flux none and adds --ader-order 2 --time-integrator ader
  --fault-overint 0 --deriv-cache --shared-ck-recursion --tfinal 12 --mesh (no over-integration, no resample).
- Mesh: tpv102/mesh/tpv102_200m_hybrid.msh, not in the repo (the sbatch runs the gmsh command below if it is missing);
  tpv102/mesh/tpv102_200m_hybrid.geo at b09dbfc: y-mirror structured prism strip each side of the fault (5 layers over
  1 km, lc_fault 200 m) in an unstructured box (lc 5 km), ~170 MB per the sbatch. Element count, sha256 and gmsh
  version: not recorded.
- Run: Frontera (the sbatch does not name the machine; the commit that added these golds, 4eeec87, calls them Frontera references), partition normal, account EAR20006,
  10 nodes x 500 ranks, wall limit 48:00:00, job 7774565. Run date, wall time used and restart segments: not recorded
  (run commit dated 2026-06-05).
- results/: the 9 SCEC station files tpv102_station_flt_<strike km>_<depth km>.dat (n = negative), 96,893 rows,
  1.239e-4 s apart, t = 0 to 12 s. No units in the header; values are SI (tau2 7.5e7 Pa at t = 0). First lines of
  tpv102_station_flt_0_12.dat:
  ```text
  # TPV102 station: flt_0_12 (along_strike=0, down_dip=12000)
  # Columns: time slip1 slip2 V1 V2 tau1 tau2 sigma_n log10_theta
  # BP5 convention: 1 = dip (tangent1), 2 = strike (tangent2). TPV102 is pure strike-slip so V1/slip1/tau1 ~ 0 and the interesting physics is in column 2.
  ```
- plots/: tpv102_overview.png and tpv102_flt_<station>.png for the 9 stations (plot script not recorded).

## Reproduce

```bash
git checkout b09dbfcc14e5000d666fa0be59b6c13621c889c1
# MFEM itself must be built with PETSc + MUMPS (config/config.mk; see build_frontera.sh)
sbatch miniapps/seas/jobs/spatial_dyn_build.sbatch clean-rebuild
(cd miniapps/seas && gmsh -format msh22 -3 tpv102/mesh/tpv102_200m_hybrid.geo -o tpv102/mesh/tpv102_200m_hybrid.msh)
sbatch miniapps/seas/jobs/tpv102_spatial/upwind_ader/tpv102_p1_pureupwind_200m_hybrid_normal.sbatch
```

Output: miniapps/seas/tpv102/out_p1_aderO2_pu_200m_hybrid_<jobid>/tpv102_station_*.dat. Mesh, config and tfinal are
set at the top of the sbatch; TPV102_RESTART=<previous OUT>/cp continues past the wall limit (same -n 500). On another
cluster adjust the module lines, `LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:...`, ibrun, -p normal and -A EAR20006; the sbatch
has no absolute paths.

Gold folder added in 4eeec87 (2026-06-12); README added 2026-10-17 by Chunhui Zhao's request (trace and reproduce every gold from its commit).
