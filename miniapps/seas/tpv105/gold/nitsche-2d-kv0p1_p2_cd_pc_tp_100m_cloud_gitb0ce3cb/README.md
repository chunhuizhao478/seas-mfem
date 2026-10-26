# TPV105-2D gold: Nitsche dynamic rupture on a 2-D mesh with thermal pressurization, Kelvin-Voigt damping, p2, 100 m (cloud run at b0ce3cb6)

- Code: seas-mfem branch 2d-inplane, commit b0ce3cb6d0b0a437ac773d2f61f27241b1771838 (b0ce3cb6, on GitHub),
  built and run in a Claude Code cloud container by jobs/run_gold_cloud.sh: the CLAUDE.md guard (git fetch,
  clean tracked tree, HEAD on GitHub), every `*.o` deleted, seas_nitsche_dyn_driver built at HEAD with
  jobs/mk.sh (the container's make wrapper around cycle/cycle.mk; MFEM built from the same checkout, whose
  files outside miniapps/seas are those of the branch point 3926e9c8), the tree checked unchanged after
  the build. Binary sha256 e7fd5b6f7cdf248bf9064d3d895b88858d31dff7e423b008d4b0f20564d171bf (the same
  binary as the TPV205-2D gold of this commit).
- Deck: jobs/tpv105_2d_nitsche.toml (= tpv105/configs/tpv105_2d_nitsche.toml at b0ce3cb6). Plane strain in the
  (x, y) plane, fault on y = 0; rate and state, slip law with strong rate weakening (a 0.01 inside
  |x| <= 15 km rising to 0.02 over w = 3 km, b 0.014, L 0.4 m, V0 1e-6 m/s, f0 0.6, f_w 0.2, V_w 0.1 m/s
  inside rising to 1.0 m/s), V_init 1e-16 m/s, sigma_n 120 MPa, tau 21 MPa; nucleation overstress 60 MPa
  in a 1.5 km radius about x = -10.5 km, ramped over 1 s; thermal pressurization (TPV105 eqs 11-12:
  alpha_th 1e-6 and alpha_hy 1e-4 m^2/s, rho c 2.7e6 J/(m^3 K), Lambda 0.1 MPa/K, h 0.01 m, T_ini 483.15 K,
  p_ini 80 MPa; grid dz0 0.5 mm = h/20, growth 1.05, z_max 0.4 m, 77 unknowns per fault node, backward
  Euler); Prakash-Clifton v* 1 m/s, L 0.01 m; [numerics] kv_eta 0.1, kv_scale "node", kv_length
  "inscribed"; cfl 0.5 of the 2-D C(2) = 0.365; absorbing outer boundary; `up = [0, 0, 1]` puts the
  in-plane slip in the strike slot (slip2 = SCEC h-slip).
- Mesh: tpv105/mesh/tpv105_2d_100m.msh (in the repo, committed in phase B2), 56,652 triangles, 604 fault
  edges, sha256 184a4d1b0f3ff570f2ac149a8be4e5199f59301fa97b95bf0f45c0108bd92565.
- Run: cloud container (Intel Xeon @ 2.10 GHz, 4 cores), `mpirun -np 4 --bind-to none`, OMP_NUM_THREADS=1,
  2026-10-26 16:42-16:59 UTC; `--order 2 --log-every 2000`; 0-15 s, 47,913 steps, dt 3.130720707e-4 s,
  central difference ("cd"), driver wall 1027 s (20.9 ms/step on the slowest rank; the fault work with
  thermal pressurization is 3.0 ms/step), max|V| 12.6 m/s, max T 2348 K, max p 159.1 MPa, no NaN;
  sigma_n_eff stays in [112.5, 125.6] MPa and sigma_n_eff - (p - p_ini) >= 40.9 MPa.
- results/: the 7 on-fault station files (t slip1 slip2 V1 V2 tau1 tau2 sigma_n psi T p; 1 = dip,
  2 = strike; SI, T in K; sigma_n = the effective normal stress the law sees, -T.n - (p - p_ini),
  compression positive; tpv105_2d_station_conventions.txt), the station snap list (every station within
  6.6 mm of a fault node), driver.log and RUN_INFO.txt; stations.sha256 holds the sha256 of each station
  file. The `*_copies.dat` files and the initial-state CSVs the run also wrote are not kept.
- Reproducibility: the station data (all lines but the commit-stamped headers) are byte for byte those of
  the review 8 preview run (the same code before it was committed) in the same container.
- plots/: tpv105_2d_stations.pdf and its page PNGs (misfit overview, then per station 3 x 3: slip rate,
  front zoom, slip, shear stress, effective normal stress, shear / normal, temperature, pore pressure,
  psi; MDSBI black, SPECFEM2D grey, this run red), compare_vs_ref.txt and summary.json, made by
  scripts/plot_tpv105_2d_stations.py (command below).
- vs MDSBI 50 m (tpv105/benchmark_data/MDSBI2D_50m), window 0-12.5 s (SPECFEM2D ends at 12.5 s); the
  SPECFEM2D 50 m values against MDSBI show the spread between two reference codes:

  | measure | this run vs MDSBI | SPECFEM2D vs MDSBI |
  |---|---|---|
  | arrival (1 mm/s), -15 to +15 km | -0.012 to +0.015 s | 0 to +0.018 s |
  | arrival at -18 / +18 km | +0.011 / +0.037 s | +0.020 / +0.030 s |
  | final slip, -15 to +15 km | -2.19 to -0.53 % | -2.34 to +0.72 % |
  | final slip at -18 / +18 km | -7.60 / -6.08 % | -12.35 / +1.79 % |
  | peak slip rate | -3.6 to +5.3 % | +2.0 to +10.8 % |
  | peak temperature | -2.0 to -0.3 % | -1.1 to +2.6 % |
  | pore pressure at 12.5 s | -0.91 to -0.15 MPa | -0.56 to +1.18 MPa |

  Final slip lies between the two references at the five stations from -18 to +9 km (at -18 km they
  differ by 12 %: MDSBI 1.255 m, SPECFEM 1.100 m, ours 1.160 m). At +15 and +18 km, the two ends of the
  taper into the velocity-strengthening barrier on the far side from the nucleation, ours is below both
  (+15 km: 4.527 m against 4.628 and 4.661 m; +18 km: 2.401 m against 2.556 and 2.602 m).
- The p1 100 m run of the preview (same code) was 0.05 s early with peaks 10-30 % low; p2 is the gold.
- Status: TPV105-2D is completed for the 2-D in-plane design (gate G2; Chunhui Zhao 2026-10-26: "gold
  tpv205-2d, tpv105-2d results p=2, and mark as completed"); see document/2d_inplane/BENCHMARK_STATUS.md.

## Reproduce

In a checkout of b0ce3cb6 with MFEM built (MPI, MUMPS, METIS) and the build variables of jobs/mk.sh adapted
to the machine, from miniapps/seas:

    make -f cycle/cycle.mk seas_nitsche_dyn_driver        # plus the MFEM_BUILD_DIR etc. of jobs/mk.sh
    mpirun -np 4 ./seas_nitsche_dyn_driver --config tpv105/configs/tpv105_2d_nitsche.toml \
        --order 2 --log-every 2000 --output-dir <out>

jobs/run_gold_cloud.sh is the script as run (container paths). Plots (from miniapps/seas, D = this folder):

    PLOT_NOTE="gold p2 100 m, commit b0ce3cb6, cloud 4 ranks|grey = SPECFEM2D 50 m minus MDSBI 50 m" \
    python3 -I scripts/plot_tpv105_2d_stations.py $PWD $D/plots 15 \
        MDSBI2D_50m=tpv105/benchmark_data/MDSBI2D_50m SPECFEM2D_50m=tpv105/benchmark_data/SPECFEM2D_50m \
        -- gold_p2_100m=$D/results

- Requested by Chunhui Zhao (2-D in-plane design, phase D gate G2, 2026-10-26).
