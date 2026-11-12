# TPV104 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p3, 200 m (Vista job 1080362)

- Code: seas-mfem branch nitsche-benchmarks, commit b50741be112a90e98a792abcc69202263cba5bdf (b50741b, on GitHub),
  built on TACC Vista by jobs/vista_nitsche/build_verify_vista.sbatch (job 1080351, clean checkout,
  BUILD_COMMIT_VISTA = b50741b, every gate suite 0 failed). b50741b is trunk a4c9dee, which has the Nitsche
  speedups on by default (the run logs "fast fault on, fused K/K_eta off, dense blocks on"), plus the p3 job
  scripts. Driver seas_nitsche_dyn_driver. Same deck, mesh and 16 x 128 rank layout as the 200 m p2 gold
  (nitsche-kv0p1_p2_cd_pc_200m_jobid1075917_git81759da); the order and the speedups differ, and the speedups
  change results only at round-off (checked on TPV205, profiling job 1079640).
- Deck and job script: jobs/ (= tpv104/configs/tpv104_nitsche.toml and jobs/vista_nitsche/{tpv104_200m_p3.sbatch,
  run_nitsche_vista.sh} at b50741b; the sbatch and deck are byte-identical to the copies the job saved, and the
  deck is unchanged since 81759da). Rate and state, slip law with strong rate weakening (f_w 0.2); fault ends at
  x = +-18 km (the old TPV102 geometry); Prakash-Clifton v* 1 m/s, L 0.01 m; [numerics] kv_eta 0.1, kv_scale
  "node", kv_length "inscribed"; cfl 0.5.
- Mesh: tpv104/mesh/tpv104_200m.msh, made on Vista by build_verify_vista.sbatch (gmsh 4.15.0, msh 2.2):
  2,466,549 tets, 37,928 fault faces, sha256 10cecfd73d0c7cb93cfa188aac08e4a99a6ec81d439a3363bd7175e3a8ae4aa0
  (not in the repo), the same file as the p1 and p2 golds. The mesh is not mirror-symmetric about the fault.
- Run: Vista gg, 16 nodes / 2048 ranks (128 per node), account EAR24034, 2026-10-18 11:34 pm to 2026-10-19
  1:28 am PDT (1:54:18 job wall, 10.1 SU), 0-12 s, 176,571 steps, dt 6.796163869e-5 s, central difference
  ("cd"), 6839 s driver wall (192 s of it fault output), max|V| 37.48 m/s, no NaN.
- results/: the 9 SCEC station files (t slip1 slip2 V1 V2 tau1 tau2 sigma_n psi; 1 = dip, 2 = strike; SI;
  sigma_n compression positive; tpv104_station_conventions.txt), the station snap list, driver.log, RUN_INFO.txt,
  the job's stdout tpv104_200m_p3_1080362.out (and .err, a module notice only) and meminfo_kB.log (first node,
  every 2 s: epoch s, MemTotal kB, MemAvailable kB), the sources of the speed numbers below. stations.sha256
  holds the sha256 of each station file as written on Vista (the files were moved in three parts each and
  checked against it after joining). Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  fault.vtkhdf (600 frames, Vista scratch) is not here.
- plots/: tpv104_200m_p3_vista1080362_stations.pdf and its page PNGs (slip-rate overview of the 9 stations, then
  per station slip rate, front zoom, slip, shear stress, normal stress, psi; this run in red vs DRDG3D, SeisSol,
  the mixed-flux p2 gold, the 200 m p1 gold, the 100 m p1 run 1076418 (stopped at t = 8.08 s) and the 200 m p2
  gold in green); compare.txt / compare_summary.json from compare_nitsche_fault_stations.py vs DRDG3D, SeisSol,
  both mixed-flux golds and the 200 m p1 gold, with the p2 gold as the run before this one (command below).
  SeisSol ends at 11.9987 s, so every source is measured to there.
- vs DRDG3D (median |misfit| of 9 stations): arrival 0.028 s (late at 7 of 9), final slip 0.27 % (low at 7 of 9),
  peak h-slip-rate 6.9 %; noise RMS sigma_n / tau_h / tau_v 0.035 / 0.023 / 0.022 MPa (DRDG3D 0.019 / 0.030 /
  0.025). Arrival within 0.1 s everywhere PASS (worst +0.052 s at x2 = 0, 12 km depth); slip within 2 %
  everywhere PASS (worst -0.55 % at x2 = -12 km and 12 km, 3 km depth); peak within 10 % at 7 of 9. The two
  misses: x2 = -12 km, 3 km depth, 13.8 % high (SeisSol 9.3 % high, mixed-flux p2 17.6 % high); x2 = 12 km,
  3 km depth, 13.2 % high (SeisSol 3.9 % high, mixed-flux p2 20.9 % high). The p2 gold's two misses are inside
  10 % here: x2 = 0, 12 km depth, 6.9 % high (p2 12.3 %) and x2 = 12 km, 12 km depth, 3.5 % low (p2 13.5 %).
  Medians along p2 -> p3 FAIL to shrink: slip 0.49 -> 0.27 % shrinks, arrival 0.012 -> 0.028 s and peak
  4.8 -> 6.9 % grow. For scale, SeisSol vs DRDG3D: max |dt| 0.106 s, max |slip| 2.29 %, peak within 10 % at
  8 of 9.
- Fault-wide normal stress (driver.log): 11.8-211.8 MPa over the run, against 120 MPa at rest (p2 gold
  5.1-196 MPa); the minimum fell from 74.9 MPa at t = 8.0 s to 12.6 MPa at t = 10.0 s and stayed there to 12 s;
  the stations see 0.035 MPa RMS. Cause as in the p2 gold README: the tets on the two sides of each fault face
  differ in shape on this gmsh mesh.

## Speed

From results/tpv104_200m_p3_1080362.out, driver.log and meminfo_kB.log; job wall from sacct; Vista gg charges
0.33 SU per node-hour. The p2 column is the p2 gold, run before the speedups, at the same ranks.

| | p3, this run (b50741b) | p2 gold 1075917 (81759da) |
|---|---|---|
| steps to 12 s, dt | 176,571, 6.796e-5 s | 104,388, 1.150e-4 s |
| job wall (sacct) | 1:54:18 | 1:03:42 |
| driver wall | 6839 s (38.7 ms per step) | 3785 s (36.3 ms per step) |
| time in ode.Step per rank, slowest | 6313 s = 35.75 ms per step | not logged |
| time in ode.Step per rank, min / mean | 29.65 / 33.30 ms per step (max / mean 1.07) | not logged |
| step 1 (with operator assembly) | 37.7 s | not logged |
| setup plus the first 1000 steps | 190 s (driver.log wall at step 1000) | |
| fault work, busiest rank | 6.08 ms per step | |
| fault output (rank 0) | 600 frames, 192 s (0.32 s per frame) | 599 frames, 418 s (0.70 s per frame) |
| node-hours, SU | 30.5, 10.1 | 17.0, 5.6 |
| first-node memory peak | 177.7 of 237.6 GB | not logged |
| output folder on Vista | 7.8 GB | |

- As for TPV205, a p3 step now costs about what a p2 step cost before the speedups (38.7 against 36.3 ms of
  driver wall per step), so the p3 run takes 1.79 x the p2 gold's job wall, mostly from its 1.69 x more steps.
  No TPV104 p2 run with the speedups exists to compare per step.
- The job script estimated 47 to 58 ms per step from the TPV205 p2 measurement scaled to this mesh; Vista came
  in at 35.75 ms.
- Step 1 took 37.7 s here and 8.5 s in the TPV102 p3 run on the same mesh, ranks and order (job 1080363);
  the difference is not investigated.
- The driver's "time in ode.Step" excludes setup; the slowest rank sets the pace.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-benchmarks \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-p3
    cd seas-mfem-p3 && git checkout b50741be112a90e98a792abcc69202263cba5bdf
    cd miniapps/seas/jobs/vista_nitsche && sbatch build_verify_vista.sbatch   # build, tpv104_200m.msh, gate tests
    sbatch tpv104_200m_p3.sbatch

build_verify_vista.sbatch generates tpv104/mesh/tpv104_200m.msh from the committed tpv104/mesh/tpv104_200m.geo with
Vista's gmsh 4.15.0; gmsh output depends on the platform, so compare its sha256 with the one above (the job refuses
a mesh whose sha256 is not listed). The job refuses to start unless the checkout is clean, HEAD is on GitHub and
the driver was built at HEAD. Compare (from miniapps/seas, with D = this folder; rerunning it on these results/
reproduces compare.txt apart from its last line, the path it wrote the JSON to, and compare_summary.json to
round-off; checked 2026-10-19):

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv104/stations/tpv104_fault_stations.csv \
        --prefix tpv104 --run vista_p2_200m=tpv104/gold/nitsche-kv0p1_p2_cd_pc_200m_jobid1075917_git81759da/results \
        --run vista_p3_200m=$D/results \
        --ref drdg3d=tpv104/benchmark_data/DRDG3D --ref seissol=tpv104/benchmark_data/seisol \
        --gold mixedflux_p1=tpv104/gold/mixedflux_p1_rk45_jobid7765845_git04674fc/results \
        --gold mixedflux_p2=tpv104/gold/mixedflux_p2_rk45_jobid7765380_git04674fc/results \
        --gold vista_p1_200m=tpv104/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid1075919_git81759da/results:nitsche \
        --json compare_summary.json > compare.txt

- Requested by Chunhui Zhao (TPV205/TPV104 benchmarks, 2026-10-19).
