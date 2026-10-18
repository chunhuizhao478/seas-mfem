# TPV104 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p1, 200 m (Vista job 1075919)

- Code: seas-mfem branch nitsche-benchmarks, commit 81759da2fb36f180fbf3bdd087318b6d680b732f (81759da, on GitHub),
  built on TACC Vista by
  jobs/vista_nitsche/build_verify_vista.sbatch (clean checkout, BUILD_COMMIT_VISTA = 81759da).
  Driver seas_nitsche_dyn_driver.
- Deck and job script: jobs/ (= tpv104/configs/tpv104_nitsche.toml and jobs/vista_nitsche/{tpv104_200m_p1.sbatch,
  run_nitsche_vista.sh} at 81759da). Rate and state, slip law with strong rate weakening (f_w 0.2);
  Prakash-Clifton v* 1 m/s, L 0.01 m; [numerics] kv_eta 0.1, kv_scale "node", kv_length "inscribed"; cfl 0.5.
- Mesh: tpv104/mesh/tpv104_200m.msh, made on Vista by build_verify_vista.sbatch (gmsh 4.15.0, msh 2.2):
  2,466,549 tets, 37,928 fault faces, sha256 10cecfd73d0c7cb93cfa188aac08e4a99a6ec81d439a3363bd7175e3a8ae4aa0
  (not in the repo).
- Run: Vista gg, 4 nodes / 512 ranks, 2026-10-16 2:51-3:19 pm PDT, 0-12 s, 48,961 steps, dt 2.45096e-4 s,
  central difference ("cd"), 1641 s wall (430 s of it fault output), max|V| 17.14 m/s, no NaN.
- results/: the 9 SCEC station files as the job wrote them (t slip1 slip2 V1 V2 tau1 tau2 sigma_n psi;
  1 = dip, 2 = strike; SI; sigma_n compression positive; tpv104_station_conventions.txt), the station snap
  list, driver.log and RUN_INFO.txt. Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  fault.vtkhdf (7.7 GB) is not here.
- plots/: tpv104_200m_p1_vista1075919_stations.pdf and its page PNGs (slip-rate overview of the 9 stations, then
  per station slip rate, front zoom, slip, shear stress, normal stress, psi; this run vs DRDG3D, SeisSol and
  the mixed-flux p2 gold); compare.txt / compare_summary.json from compare_nitsche_fault_stations.py vs
  DRDG3D, SeisSol and both mixed-flux golds.
- vs DRDG3D (median |misfit| of 9 stations): arrival 0.079 s (early at all 9), final slip 0.49 %,
  peak h-slip-rate 20.1 % (low at 8 of 9; within 10 % only at the hypocentre); noise RMS
  sigma_n / tau_h / tau_v 0.071 / 0.024 / 0.043 MPa (DRDG3D 0.019 / 0.030 / 0.025). Fails the pass test
  (arrival > 0.1 s at 3 stations, slip > 2 % at the 4 stations at +-12 km), as expected of the coarsest run.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-benchmarks \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-benchmarks
    cd seas-mfem-benchmarks && git checkout 81759da2fb36f180fbf3bdd087318b6d680b732f
    cd miniapps/seas/jobs/vista_nitsche && sbatch build_verify_vista.sbatch   # build, tpv104_200m.msh, gate tests
    sbatch tpv104_200m_p1.sbatch

build_verify_vista.sbatch generates tpv104/mesh/tpv104_200m.msh from the committed tpv104/mesh/tpv104_200m.geo with Vista's gmsh 4.15.0;
gmsh output depends on the platform, so compare its sha256 with the one above. The job refuses to start unless
the checkout is clean, HEAD is on GitHub and the driver was built at HEAD. Compare (from miniapps/seas, with
D = this folder; rerunning it on these results/ reproduces compare.txt):

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv104/stations/tpv104_fault_stations.csv \
        --prefix tpv104 --run vista_p1_200m=$D/results \
        --ref drdg3d=tpv104/benchmark_data/DRDG3D --ref seissol=tpv104/benchmark_data/seisol \
        --gold mixedflux_p1=tpv104/gold/mixedflux_p1_rk45_jobid7765845_git04674fc/results \
        --gold mixedflux_p2=tpv104/gold/mixedflux_p2_rk45_jobid7765380_git04674fc/results \
        --json compare_summary.json > compare.txt

- Requested by Chunhui Zhao (TPV205/TPV104 benchmarks, 2026-10-16).
