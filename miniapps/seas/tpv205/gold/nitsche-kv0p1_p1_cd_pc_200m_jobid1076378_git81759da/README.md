# TPV205 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p1, 200 m (Vista job 1076378)

- Code: seas-mfem branch nitsche-benchmarks, commit 81759da2fb36f180fbf3bdd087318b6d680b732f (81759da, on GitHub),
  built on TACC Vista by
  jobs/vista_nitsche/build_verify_vista.sbatch (clean checkout, BUILD_COMMIT_VISTA = 81759da).
  Driver seas_nitsche_dyn_driver.
- Deck and job script: jobs/ (= tpv205/configs/tpv205_nitsche.toml and jobs/vista_nitsche/{tpv205_200m_p1.sbatch,
  run_nitsche_vista.sh} at 81759da). Linear slip weakening (mu_s 0.677, mu_d 0.525, d_c 0.40 m), three
  prestress patches at 7.5 km depth; Prakash-Clifton v* 1 m/s, L 0.01 m; [numerics] kv_eta 0.1,
  kv_scale "node", kv_length "inscribed"; cfl 0.5.
- Mesh: tpv205/mesh/tpv2053d_200m.msh, made on Vista by build_verify_vista.sbatch (gmsh 4.15.0, msh 2.2):
  1,894,583 tets, 26,228 fault faces, sha256 30c63b923781177830330a36aebcc4159cfc8043374ea4dcd883a9b2eec13538
  (not in the repo).
- Run: Vista gg, 4 nodes / 576 ranks, 2026-10-16 3:46-4:05 pm PDT (resubmission of 1075918, which stopped in
  the job's git fetch), 0-12 s, 41,614 steps, dt 2.88365e-4 s, central difference ("cd"), 1102 s wall
  (292 s of it fault output), max|V| 8.12 m/s, no NaN.
- results/: the 16 SCEC station files as the job wrote them (t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path;
  1 = dip, 2 = strike; SI; sigma_n compression positive; tpv205_station_conventions.txt), the station snap
  list, driver.log and RUN_INFO.txt. Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  fault.vtkhdf is not here.
- plots/: tpv205_200m_p1_vista1076378_stations.pdf and its page PNGs (slip-rate overview of the 16 stations,
  then per station slip rate, front zoom, slip, shear stress, normal stress, dip slip rate; this run vs
  DRDG3D 100 m and 200 m and the mixed-flux p2 gold); compare.txt / compare_summary.json from
  compare_nitsche_fault_stations.py vs DRDG3D 100 m, DRDG3D 200 m and both mixed-flux golds.
- vs DRDG3D 100 m (median |misfit| of 16 stations): arrival 0.043 s (early at 15 of 16, equal at the hypocentre), final slip 1.46 %,
  peak h-slip-rate 19.3 % (low at all 16); noise RMS sigma_n / tau_h / tau_v 0.028 / 0.017 / 0.025 MPa
  (DRDG3D 100 m tau_h / tau_v 0.017 / 0.021). Fails the pass test (arrival > 0.1 s at 2 stations; slip > 2 %
  at 5 stations: the four at +-12 km, 2.8 to 5.4 % low, and x2 = -7.5 km at the surface, 2.5 %; peak within
  10 % at 2 of 16), as expected of the coarsest run.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-benchmarks \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-benchmarks
    cd seas-mfem-benchmarks && git checkout 81759da2fb36f180fbf3bdd087318b6d680b732f
    cd miniapps/seas/jobs/vista_nitsche && sbatch build_verify_vista.sbatch   # build, tpv2053d_200m.msh, gate tests
    sbatch tpv205_200m_p1.sbatch

build_verify_vista.sbatch generates tpv205/mesh/tpv2053d_200m.msh from the committed tpv205/mesh/tpv2053d_200m.geo with Vista's gmsh 4.15.0;
gmsh output depends on the platform, so compare its sha256 with the one above. The job refuses to start unless
the checkout is clean, HEAD is on GitHub and the driver was built at HEAD. Compare (from miniapps/seas, with
D = this folder; rerunning it on these results/ reproduces compare.txt):

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv205/stations/tpv205_fault_stations.csv \
        --prefix tpv205 --run vista_p1_200m=$D/results \
        --ref drdg3d_100m=tpv205/benchmark_data/DRDG3D_100m_O4 --ref drdg3d_200m=tpv205/benchmark_data/DRDG3D_200m_O4 \
        --gold mixedflux_p1=tpv205/gold/mixedflux_p1_rk45_jobidcombined_7762879_7763019_git04674fc/results \
        --gold mixedflux_p2=tpv205/gold/mixedflux_p2_rk45_jobid7762880_git04674fc/results \
        --peak-min-stations 14 --json compare_summary.json > compare.txt

- Requested by Chunhui Zhao (TPV205/TPV104 benchmarks, 2026-10-16).
