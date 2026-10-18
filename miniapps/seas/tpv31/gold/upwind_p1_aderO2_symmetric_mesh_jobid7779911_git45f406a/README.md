# TPV31 gold: ADER O2 upwind p1, matrix Riemann, symmetric mesh, 50 m (Frontera job 7779911)

- Code: seas-mfem commit 45f406a65b6d03b87f4be61a9b9c1c4f9d5f938a (45f406a, on GitHub in origin/nitsche-benchmarks);
  driver seas_spatial_dyn_driver via ibrun. Built, as the sbatch says, by jobs/spatial_dyn_build.sbatch clean-rebuild
  (deletes the driver's objects, then `make seas_spatial_dyn_driver -j48`; the binary actually used is not recorded).
- Deck: jobs/tpv31_p1_pureupwind_50m_hybrid_normal.sbatch (same file as in the repo at 45f406a) with
  tpv31/configs/tpv31.toml. Linear slip weakening (mu_s 0.58, mu_d 0.45, d_c 0.18 m), cohesion taper in the top 2.4
  km; depth_profile_1d material (steps at 2.4, 5, 10 km); depth_proportional prestress; instantaneous overstress 4.95
  MPa, radius 1.4 km + 0.6 km taper at 7.5 km depth. Numerics: order 1, ader_order 2, mixed_flux "none", interior_flux
  "matrix", fault_iterator "substep", cfl 0.5 (cfl_safety "dg"), tfinal 15 s. The sbatch adds --ader-order 2
  --time-integrator ader --mixed-flux none --fault-overint 0 --deriv-cache --shared-ck-recursion --tfinal 15 --mesh.
- Mesh: tpv31/mesh/tpv31_50m_hybrid.msh, not in the repo (the sbatch runs the gmsh command below if it is missing);
  tpv31/mesh/tpv31_50m_hybrid.geo at 45f406a: y-mirror structured prism strip each side of the fault (5 layers of 50
  m) in an unstructured box (lc 5 km), ~350 MB per the sbatch. Element count, sha256 and gmsh version: not recorded.
- Run: Frontera (named in the sbatch), partition normal, account EAR20006, 10 nodes x 500 ranks, wall limit 48:00:00,
  job 7779911. Run date, wall time used and restart segments: not recorded (run commit dated 2026-06-08).
- results/: the 30 station files tpv31_station_faultst<sss>dp<ddd>.dat (strike 0/6/12 km x ten depths, in 100 m:
  faultst060dp024 = 6000 m along strike, 2400 m down dip), 238,423 rows, 6.291e-5 s apart, t = 0 to 15 s; the header
  gives Pa, n-stress compression positive. First lines of tpv31_station_faultst000dp000.dat:
  ```text
  # TPV31 station: faultst000dp000 (along_strike=0, down_dip=0)
  # Columns (SCEC TPV31 trace layout, LSW): t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate  v-shear-stress  n-stress  mu_eff
  # Canonical frame: h = strike = component 2, v = dip = component 1.  TPV31 is right-lateral strike-slip; v-* channels remain near 0.  Stresses in Pa, n-stress compression-positive; mu_eff is the LSW μ(δ) at δ = sqrt(slip1² + slip2²).
  ```
- plots/: tpv31_overview.png and tpv31_faultst<sss>dp<ddd>.png for the 30 stations (plot script not recorded).

## Reproduce

```bash
git checkout 45f406a65b6d03b87f4be61a9b9c1c4f9d5f938a
# MFEM itself must be built with PETSc + MUMPS (config/config.mk; see build_frontera.sh)
sbatch miniapps/seas/jobs/spatial_dyn_build.sbatch clean-rebuild
(cd miniapps/seas && gmsh -format msh22 -3 tpv31/mesh/tpv31_50m_hybrid.geo -o tpv31/mesh/tpv31_50m_hybrid.msh)
sbatch miniapps/seas/jobs/tpv31_spatial/upwind_ader/tpv31_p1_pureupwind_50m_hybrid_normal.sbatch
```

Output: miniapps/seas/tpv31/out_p1_aderO2_pu_50m_hybrid_<jobid>/. Mesh and config are set at the top of the sbatch;
TPV31_TFINAL overrides 15 s and TPV31_RESTART=<previous OUT>/cp continues past the wall limit (same -n 500). On
another cluster adjust the module lines, `LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:...`, ibrun, -p normal and -A EAR20006;
the sbatch has no absolute paths.

Gold folder added in 4eeec87 (2026-06-12); README added 2026-10-17 by Chunhui Zhao's request (trace and reproduce every gold from its commit).
