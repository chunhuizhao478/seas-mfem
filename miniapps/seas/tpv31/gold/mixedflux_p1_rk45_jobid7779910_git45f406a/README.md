# TPV31 gold: RK45 mixed flux p1, matrix Riemann, contrast guard 0.05, 50 m (Frontera job 7779910)

- Code: seas-mfem commit 45f406a65b6d03b87f4be61a9b9c1c4f9d5f938a (45f406a, on GitHub in origin/nitsche-benchmarks);
  driver seas_spatial_dyn_driver via ibrun. Built, as the sbatch says, with `make seas_spatial_dyn_driver` in
  miniapps/seas (the binary actually used is not recorded).
- Deck: jobs/tpv31_p1_rk45_mixedflux_50m_normal_contrastguard.sbatch (same file as in the repo at 45f406a) with
  tpv31/configs/tpv31_rk_mixedflux.toml. Linear slip weakening (mu_s 0.58, mu_d 0.45, d_c 0.18 m), cohesion taper in
  the top 2.4 km; depth_profile_1d material (steps at 2.4, 5, 10 km); depth_proportional prestress; instantaneous
  overstress 4.95 MPa, radius 1.4 km + 0.6 km taper at 7.5 km depth. Numerics: order 1, interior_flux "matrix",
  mixed_flux "adjacent", cfl 0.5 (cfl_safety "dg"), auto dt, tfinal 15 s, no time_integrator key, seam_continuous =
  true (a no-op per the sbatch). The sbatch adds --time-integrator rk45 --mixed-flux-contrast-tol 0.05 --tfinal 15
  --mesh, no --deriv-cache.
- Mesh: tpv31/mesh/tpv31_50m.msh, not in the repo (the sbatch runs the gmsh command below if it is missing);
  tpv31/mesh/tpv31_50m.geo at 45f406a: lc_fault 50 m, lc_far 5 km, msh 2.2, fault reaching the free surface. Element
  count, sha256 and gmsh version: not recorded.
- Run: Frontera (named in the sbatch), partition normal, account EAR20006, 10 nodes x 500 ranks, wall limit 48:00:00,
  job 7779910. Run date, wall time used and restart segments: not recorded (run commit dated 2026-06-08).
- results/: the 30 station files tpv31_station_faultst<sss>dp<ddd>.dat (strike 0/6/12 km x ten depths, in 100 m:
  faultst060dp024 = 6000 m along strike, 2400 m down dip), 44,580 rows, 3.365e-4 s apart, t = 0 to 15 s; the header
  gives Pa, n-stress compression positive. First lines of tpv31_station_faultst000dp000.dat:
  ```text
  # TPV31 station: faultst000dp000 (along_strike=0, down_dip=0)
  # Columns (SCEC TPV31 trace layout, LSW): t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate  v-shear-stress  n-stress  mu_eff
  # Canonical frame: h = strike = component 2, v = dip = component 1.  TPV31 is right-lateral strike-slip; v-* channels remain near 0.  Stresses in Pa, n-stress compression-positive; mu_eff is the LSW μ(δ) at δ = sqrt(slip1² + slip2²).
  ```
- plots/: tpv31_overview.png and tpv31_faultst<sss>dp<ddd>.png for the 30 stations (plot script not recorded).
- Name: the folder says mixedflux_p1_rk45; the job is the contrast-guard variant (--mixed-flux-contrast-tol 0.05).

## Reproduce

```bash
git checkout 45f406a65b6d03b87f4be61a9b9c1c4f9d5f938a
# MFEM itself must be built with PETSc + MUMPS (config/config.mk; see build_frontera.sh)
(cd miniapps/seas && make seas_spatial_dyn_driver)
(cd miniapps/seas && gmsh -format msh22 -3 tpv31/mesh/tpv31_50m.geo -o tpv31/mesh/tpv31_50m.msh)
sbatch miniapps/seas/jobs/tpv31_spatial/mixedflux_rk45/tpv31_p1_rk45_mixedflux_50m_normal_contrastguard.sbatch
```

Output: miniapps/seas/tpv31/out_p1_rk45_mf_normal_contrastguard_<jobid>/. TPV31_CONFIG / _MESH / _OUT / _TFINAL
override the defaults; TPV31_RESTART=<previous OUT>/cp continues past the wall limit (same -n 500). On another cluster
adjust the module lines, `LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:...`, ibrun, -p normal and -A EAR20006; the sbatch has no
absolute paths.

Gold folder added in 4eeec87 (2026-06-12); README added 2026-10-17 by Chunhui Zhao's request (trace and reproduce every gold from its commit).
