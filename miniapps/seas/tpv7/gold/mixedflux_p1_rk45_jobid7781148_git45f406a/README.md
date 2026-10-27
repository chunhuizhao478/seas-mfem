# TPV7 gold: RK45 mixed flux p1, matrix Riemann, contrast guard 0.05, standard mesh, 200 m (Frontera job 7781148)

- Code: seas-mfem commit 45f406a65b6d03b87f4be61a9b9c1c4f9d5f938a (45f406a, on GitHub in origin/nitsche-benchmarks);
  driver seas_spatial_dyn_driver via ibrun. Built, as the sbatch says, by jobs/spatial_dyn_build.sbatch clean-rebuild
  (deletes the driver's objects, then `make seas_spatial_dyn_driver -j48`; the binary actually used is not recorded).
- Deck: jobs/tpv7_p1_rk45_mixedflux_200m.sbatch (same file as in the repo at 45f406a) with tpv7/configs/tpv7.toml. It
  is tpv6.toml with a low-contrast far side (vp 5000 / vs 2887 / rho 2670; near side vp 6000 / vs 3464 / rho 2670),
  tfinal 8.0 s and [problem] tag "tpv6": linear slip weakening (mu_s 0.677, mu_d 0.525, d_c 0.40 m), tau 70 MPa,
  sigma_n 120 MPa, 81.6 MPa 3 x 3 km patch at 7.5 km depth. TOML numerics: order 1, time_integrator "ader", ader_order
  2, mixed_flux "none", interior_flux "matrix", cfl 0.25 (cfl_safety "dg"); the sbatch overrides with
  --time-integrator rk45 --mixed-flux adjacent --mixed-flux-contrast-tol 0.05 --tfinal 12.0 (TOML: 8.0 s) --mesh
  tpv6/mesh/tpv6_200m.msh (TOML: the symmetric mesh), no --deriv-cache.
- Mesh: tpv6/mesh/tpv6_200m.msh, not in the repo (the sbatch runs the gmsh command below if it is missing);
  tpv6/mesh/tpv6_200m.geo at 45f406a: standard (asymmetric) mesh, lc_fault 200 m (first-pass per its header; spec 100
  m), lc_far 5 km, no material tags, shared by TPV6 and TPV7. Element count, sha256 and gmsh version: not recorded.
- Run: Frontera (the sbatch does not name the machine; the commit that added these golds, 4eeec87, calls them Frontera references), partition normal, account EAR20006,
  4 nodes x 224 ranks, wall limit 12:00:00, job 7781148. Run date, wall time used and restart segments: not recorded
  (run commit dated 2026-06-08).
- results/: 10 files tpv6_<nearside|farside>_x2_<strike km>_x3_<depth km>.dat (named tpv6_ because tpv7.toml uses the
  TPV6 writer; header "TPV6/7 station"): 5 stations x 2 sides (nearside = strong, larger-impedance side), 16,083 rows,
  7.462e-4 s apart, t = 0 to 12 s; header units: m/s, m, MPa, n-stress compression positive. First lines of
  tpv6_farside_x2_-12_x3_0.dat:
  ```text
  # TPV6/7 station x2_-12_x3_0 farside (along_strike=-12000, down_dip=0)
  # nearside = STRONG (larger Zp) side; farside = WEAK side.
  # Columns: t h-disp h-vel h-stress v-disp v-vel v-stress n-disp n-vel n-stress
  # h=strike, v=dip, n=normal; vel m/s, disp m, stress MPa, n-stress COMPRESSION-POSITIVE (drdg3d ref is negative — the viz flips it).  Velocity = imposed split-node Godunov velocity; disp = trapezoidal time-integral.
  ```
- plots/: tpv7_x2_<..>_x3_<..>_<side>.png for the 10 files, no overview (plot script not recorded).

## Reproduce

```bash
git checkout 45f406a65b6d03b87f4be61a9b9c1c4f9d5f938a
# MFEM itself must be built with PETSc + MUMPS (config/config.mk; see build_frontera.sh)
sbatch miniapps/seas/jobs/spatial_dyn_build.sbatch clean-rebuild
(cd miniapps/seas && gmsh -format msh22 -3 tpv6/mesh/tpv6_200m.geo -o tpv6/mesh/tpv6_200m.msh)
sbatch miniapps/seas/jobs/tpv7_spatial/mixedflux_rk45/tpv7_p1_rk45_mixedflux_200m.sbatch
```

Output: miniapps/seas/tpv7/out_p1_rk45_mixedflux_200m_<jobid>/. Mesh and config are set at the top of the sbatch;
TPV7_TFINAL overrides 12 s and TPV7_RESTART=<previous OUT>/cp continues past the wall limit (same -n 224). On another
cluster adjust the module lines, `LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:...`, ibrun, -p normal and -A EAR20006; the sbatch
has no absolute paths.

Gold folder added in 4eeec87 (2026-06-12); README added 2026-10-17 by Chunhui Zhao's request (trace and reproduce every gold from its commit).
