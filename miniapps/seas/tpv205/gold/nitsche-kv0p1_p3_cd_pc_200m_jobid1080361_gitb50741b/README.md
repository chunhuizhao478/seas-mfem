# TPV205 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p3, 200 m (Vista job 1080361)

- Code: seas-mfem branch nitsche-benchmarks, commit b50741be112a90e98a792abcc69202263cba5bdf (b50741b, on GitHub),
  built on TACC Vista by jobs/vista_nitsche/build_verify_vista.sbatch (job 1080351, clean checkout,
  BUILD_COMMIT_VISTA = b50741b, every gate suite 0 failed). b50741b is trunk a4c9dee, which has the Nitsche
  speedups on by default (the run logs "fast fault on, fused K/K_eta off, dense blocks on"), plus the p3 job
  scripts. Driver seas_nitsche_dyn_driver. Same deck, mesh and 16 x 144 rank layout as the 200 m p2 gold
  (nitsche-kv0p1_p2_cd_pc_200m_jobid1075916_git81759da); the order and the speedups differ, and the speedups
  change results only at round-off (profiling job 1079640: speedups off is bitwise the p2 gold).
- Deck and job script: jobs/ (= tpv205/configs/tpv205_nitsche.toml and jobs/vista_nitsche/{tpv205_200m_p3.sbatch,
  run_nitsche_vista.sh} at b50741b; the sbatch and deck are byte-identical to the copies the job saved, and the
  deck is unchanged since 81759da). Linear slip weakening (mu_s 0.677, mu_d 0.525, d_c 0.40 m), three prestress
  patches at 7.5 km depth; Prakash-Clifton v* 1 m/s, L 0.01 m; [numerics] kv_eta 0.1, kv_scale "node",
  kv_length "inscribed"; cfl 0.5.
- Mesh: tpv205/mesh/tpv2053d_200m.msh, made on Vista by build_verify_vista.sbatch (gmsh 4.15.0, msh 2.2):
  1,894,583 tets, 26,228 fault faces, sha256 30c63b923781177830330a36aebcc4159cfc8043374ea4dcd883a9b2eec13538
  (not in the repo), the same file as the p1 and p2 golds. The mesh is not mirror-symmetric about the fault.
- Run: Vista gg, 16 nodes / 2304 ranks (144 per node), account EAR24034, 2026-10-18 11:34 pm to 2026-10-19
  12:52 am PDT (1:18:03 job wall, 6.9 SU), 0-12 s, 150,076 steps, dt 7.995960894e-5 s, central difference
  ("cd"), 4663 s driver wall (119 s of it fault output), max|V| 10.70 m/s, no NaN.
- results/: the 16 SCEC station files (t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip, 2 = strike;
  SI; sigma_n compression positive; tpv205_station_conventions.txt), the station snap list, driver.log,
  RUN_INFO.txt, the job's stdout tpv205_200m_p3_1080361.out (and .err, a module notice only) and meminfo_kB.log
  (first node, every 2 s: epoch s, MemTotal kB, MemAvailable kB), the sources of the speed numbers below.
  stations.sha256 holds the sha256 of each station file as written on Vista (the files were moved in three parts
  each and checked against it after joining). Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  fault.vtkhdf (600 frames, Vista scratch) is not here.
- plots/: tpv205_200m_p3_vista1080361_stations.pdf and its page PNGs (slip-rate overview of the 16 stations,
  then per station slip rate, front zoom, slip, shear stress, normal stress, dip slip rate; this run in red vs
  DRDG3D 100 m and 200 m, the mixed-flux p2 gold, the 200 m and 100 m Nitsche p1 golds and the 200 m p2 gold in
  green); compare.txt / compare_summary.json from compare_nitsche_fault_stations.py vs DRDG3D 100 m, DRDG3D
  200 m, both mixed-flux golds and both Nitsche p1 golds, with the p2 gold as the run before this one (command
  below). DRDG3D 100 m ends at 11.9999 s, so every source is measured to there.
- vs DRDG3D 100 m (median |misfit| of 16 stations): arrival 0.014 s (late at 14 of 16), final slip 0.25 %
  (low at all 16), peak h-slip-rate 1.7 %; noise RMS sigma_n / tau_h / tau_v 0.015 / 0.013 / 0.018 MPa
  (DRDG3D 100 m tau_h / tau_v 0.017 / 0.021). PASSES all four pass tests: arrival within 0.1 s everywhere
  (worst +0.041 s at x2 = 7.5 km, 7.5 km depth); slip within 2 % everywhere (worst -0.61 % at x2 = 12 km,
  7.5 km depth); peak within 10 % at 15 of 16 (need 14); every median below the p2 gold's (arrival 0.023 s,
  slip 0.47 %, peak 1.8 %). The one miss is x2 = -12 km on the surface, 18.6 % low (p2 gold 28.1 % low,
  mixed-flux p2 16.4 % low, DRDG3D 200 m 2.3 % low). The p2 gold's other miss, x2 = 12 km at 7.5 km depth
  (11.4 % low), is 1.4 % low here.
- Fault-wide normal stress (driver.log): 53.4-185.5 MPa over the run, reached by t = 8 s and held to 12 s,
  against 120 MPa at rest (p2 gold 58.8-183 MPa); the stations see 0.015 MPa RMS. Cause as in the p2 gold
  README: the tets on the two sides of each fault face differ in shape on this gmsh mesh.

## Speed

From results/tpv205_200m_p3_1080361.out, driver.log and meminfo_kB.log; job wall from sacct; Vista gg charges
0.33 SU per node-hour. The p2 column is the p2 gold, run before the speedups, at the same ranks.

| | p3, this run (b50741b) | p2 gold 1075916 (81759da) |
|---|---|---|
| steps to 12 s, dt | 150,076, 7.996e-5 s | 88,725, 1.353e-4 s |
| job wall (sacct) | 1:18:03 | 44:54 |
| driver wall | 4663 s (31.1 ms per step) | 2670 s (30.1 ms per step) |
| time in ode.Step per rank, slowest | 4094 s = 27.28 ms per step | 27.28 ms per step (speedups off, 1 s rerun 1079640) |
| same, with the speedups | (this run) | 9.36 ms per step (1 s rerun 1079640) |
| time in ode.Step per rank, min / mean | 22.16 / 25.04 ms per step (max / mean 1.09) | not logged |
| step 1 (with operator assembly) | 6.6 s | not logged |
| setup plus the first 1000 steps | 129 s (driver.log wall at step 1000) | |
| fault work, busiest rank | 3.53 ms per step | |
| fault output (rank 0) | 600 frames, 119 s (0.20 s per frame) | 599 frames, 293 s (0.49 s per frame) |
| node-hours, SU | 20.8, 6.9 | 12.0, 4.0 |
| first-node memory peak | 145.1 of 237.6 GB | not logged |
| output folder on Vista | 4.9 GB | |

- A p3 step now costs what a p2 step cost before the speedups (27.28 ms at the slowest rank, 31 ms of driver
  wall per step), so the p3 run takes 1.74 x the p2 gold's job wall, mostly from its 1.69 x more steps
  (dt shrinks by C(3)/9 against C(2)/4 = 0.591).
- With the speedups on, p3 costs 2.91 x p2 per step at this rank count (27.28 against 9.36 ms). The job script
  estimated 39.7 ms from the 4.24 x ratio measured in the cloud on a 5000 m mesh; Vista came in 31 % lower.
- The driver's "time in ode.Step" excludes setup; the slowest rank sets the pace.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-benchmarks \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-p3
    cd seas-mfem-p3 && git checkout b50741be112a90e98a792abcc69202263cba5bdf
    cd miniapps/seas/jobs/vista_nitsche && sbatch build_verify_vista.sbatch   # build, tpv2053d_200m.msh, gate tests
    sbatch tpv205_200m_p3.sbatch

build_verify_vista.sbatch generates tpv205/mesh/tpv2053d_200m.msh from the committed tpv205/mesh/tpv2053d_200m.geo
with Vista's gmsh 4.15.0; gmsh output depends on the platform, so compare its sha256 with the one above (the job
refuses a mesh whose sha256 is not listed). The job refuses to start unless the checkout is clean, HEAD is on
GitHub and the driver was built at HEAD. Compare (from miniapps/seas, with D = this folder; rerunning it on these
results/ reproduces compare.txt apart from its last line, the path it wrote the JSON to, and compare_summary.json
to round-off; checked 2026-10-19):

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv205/stations/tpv205_fault_stations.csv \
        --prefix tpv205 --run vista_p2_200m=tpv205/gold/nitsche-kv0p1_p2_cd_pc_200m_jobid1075916_git81759da/results \
        --run vista_p3_200m=$D/results \
        --ref drdg3d_100m=tpv205/benchmark_data/DRDG3D_100m_O4 --ref drdg3d_200m=tpv205/benchmark_data/DRDG3D_200m_O4 \
        --gold mixedflux_p1=tpv205/gold/mixedflux_p1_rk45_jobidcombined_7762879_7763019_git04674fc/results \
        --gold mixedflux_p2=tpv205/gold/mixedflux_p2_rk45_jobid7762880_git04674fc/results \
        --gold vista_p1_200m=tpv205/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid1076378_git81759da/results:nitsche \
        --gold vista_p1_100m=tpv205/gold/nitsche-kv0p1_p1_cd_pc_100m_jobid1076419_git4803a90/results:nitsche \
        --peak-min-stations 14 --json compare_summary.json > compare.txt

- Requested by Chunhui Zhao (TPV205/TPV104 benchmarks, 2026-10-19).
