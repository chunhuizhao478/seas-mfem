# TPV104 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p2, 200 m (Vista job 1075917)

- Code: seas-mfem branch nitsche-benchmarks, commit 81759da2fb36f180fbf3bdd087318b6d680b732f (81759da, on GitHub),
  built on TACC Vista by
  jobs/vista_nitsche/build_verify_vista.sbatch (clean checkout, BUILD_COMMIT_VISTA = 81759da).
  Driver seas_nitsche_dyn_driver. Same commit, deck and mesh as the 200 m p1 gold
  (nitsche-kv0p1_p1_cd_pc_200m_jobid1075919_git81759da); only the order differs.
- Deck and job script: jobs/ (= tpv104/configs/tpv104_nitsche.toml and jobs/vista_nitsche/{tpv104_200m_p2.sbatch,
  run_nitsche_vista.sh} at 81759da; the sbatch and deck are byte-identical to the copies the job saved). Rate and
  state, slip law with strong rate weakening (f_w 0.2); fault ends at x = +-18 km (the old TPV102 geometry);
  Prakash-Clifton v* 1 m/s, L 0.01 m; [numerics] kv_eta 0.1, kv_scale "node", kv_length "inscribed"; cfl 0.5.
- Mesh: tpv104/mesh/tpv104_200m.msh, made on Vista by build_verify_vista.sbatch (gmsh 4.15.0, msh 2.2):
  2,466,549 tets, 37,928 fault faces, sha256 10cecfd73d0c7cb93cfa188aac08e4a99a6ec81d439a3363bd7175e3a8ae4aa0
  (not in the repo). The mesh is not mirror-symmetric about the fault (see the normal-stress note below).
- Run: Vista gg, 16 nodes / 2048 ranks (128 per node), account EAR26024, 2026-10-17 3:36-4:40 am PDT
  (1:03:42 job wall, about 5.6 SU), 0-12 s, 104,388 steps, dt 1.149563888e-4 s, central difference ("cd"),
  3785 s driver wall (418 s of it fault output), max|V| 27.62 m/s, no NaN.
- results/: the 9 SCEC station files (t slip1 slip2 V1 V2 tau1 tau2 sigma_n psi; 1 = dip, 2 = strike; SI;
  sigma_n compression positive; tpv104_station_conventions.txt), the station snap list, driver.log and
  RUN_INFO.txt. stations.sha256 holds the sha256 of each station file as written on Vista (the files were
  moved in two parts each and checked against it after joining). Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  fault.vtkhdf (Vista scratch) is not here.
- plots/: tpv104_200m_p2_vista1075917_stations.pdf and its page PNGs (slip-rate overview of the 9 stations, then
  per station slip rate, front zoom, slip, shear stress, normal stress, psi; this run vs DRDG3D, SeisSol, the
  mixed-flux p2 gold, the 200 m p1 gold and the 100 m p1 run 1076418, which stopped at t = 8.08 s);
  compare.txt / compare_summary.json from compare_nitsche_fault_stations.py vs DRDG3D, SeisSol, both mixed-flux
  golds and the 200 m p1 gold (command below).
- vs DRDG3D (median |misfit| of 9 stations): arrival 0.012 s (late at 6 of 9), final slip 0.49 % (low at 8 of 9),
  peak h-slip-rate 4.8 %; noise RMS sigma_n / tau_h / tau_v 0.048 / 0.030 / 0.033 MPa (DRDG3D 0.019 / 0.030 /
  0.025). Arrival within 0.1 s everywhere PASS (worst +0.056 s at x2 = 12 km, 12 km depth); slip within 2 %
  everywhere PASS (worst -1.45 % at x2 = -12 km, 3 km depth); peak within 10 % at 7 of 9. The two misses:
  x2 = 0, 12 km depth, 12.3 % high, where the references disagree (SeisSol 34.0 % above DRDG3D, mixed-flux p2
  23.6 % above); x2 = 12 km, 12 km depth, 13.5 % low (SeisSol 6.9 % low). Reference spread SeisSol vs DRDG3D:
  peak within 10 % at 8 of 9. Every median is below the 200 m p1 gold's (arrival 0.079 s, slip 0.49 %,
  peak 20.1 %).
- Fault-wide normal stress (driver.log): 5.1-196 MPa, against 120 MPa at rest; the minimum fell from 45 MPa to
  about 5 MPa between t = 8.3 and 8.6 s, when the rupture reached the fault ends, and stayed there to 12 s; the
  stations see 0.048 MPa RMS. The log gives only the range; that the minimum sits at the fault end is inferred
  from the 100 m p1 run 1076418, which went tensile at x = 18 km at t = 8.08 s. Cause: the tets on the two
  sides of each fault face differ in shape on this gmsh mesh; a mesh mirrored across the fault keeps 120 MPa
  exactly (1000 m tests, 2026-10-17). Mirror-mesh reruns (22 km fault with a velocity-strengthening border) use
  the scripts of commit 41879e8.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-benchmarks \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-benchmarks
    cd seas-mfem-benchmarks && git checkout 81759da2fb36f180fbf3bdd087318b6d680b732f
    cd miniapps/seas/jobs/vista_nitsche && sbatch build_verify_vista.sbatch   # build, tpv104_200m.msh, gate tests
    sbatch tpv104_200m_p2.sbatch

build_verify_vista.sbatch generates tpv104/mesh/tpv104_200m.msh from the committed tpv104/mesh/tpv104_200m.geo with Vista's gmsh 4.15.0;
gmsh output depends on the platform, so compare its sha256 with the one above. The job refuses to start unless
the checkout is clean, HEAD is on GitHub and the driver was built at HEAD. Compare (from miniapps/seas, with
D = this folder; rerunning it on these results/ reproduces compare.txt apart from its last line, the path it
wrote the JSON to, and compare_summary.json to round-off):

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv104/stations/tpv104_fault_stations.csv \
        --prefix tpv104 --run vista_p2_200m=$D/results \
        --ref drdg3d=tpv104/benchmark_data/DRDG3D --ref seissol=tpv104/benchmark_data/seisol \
        --gold mixedflux_p1=tpv104/gold/mixedflux_p1_rk45_jobid7765845_git04674fc/results \
        --gold mixedflux_p2=tpv104/gold/mixedflux_p2_rk45_jobid7765380_git04674fc/results \
        --gold vista_p1_200m=tpv104/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid1075919_git81759da/results:nitsche \
        --json compare_summary.json > compare.txt

- Requested by Chunhui Zhao (TPV205/TPV104 benchmarks, 2026-10-17).
