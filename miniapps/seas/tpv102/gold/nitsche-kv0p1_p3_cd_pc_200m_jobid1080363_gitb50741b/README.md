# TPV102 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p3, 200 m (Vista job 1080363)

- Code: seas-mfem branch nitsche-benchmarks, commit b50741be112a90e98a792abcc69202263cba5bdf (b50741b, on GitHub),
  built on TACC Vista by jobs/vista_nitsche/build_verify_vista.sbatch (job 1080351, clean checkout,
  BUILD_COMMIT_VISTA = b50741b, every gate suite 0 failed). b50741b is trunk a4c9dee, which has the Nitsche
  speedups on by default (the run logs "fast fault on, fused K/K_eta off, dense blocks on"), plus the p3 job
  scripts. Driver seas_nitsche_dyn_driver.
- Deck and job script: jobs/ (= tpv102/configs/tpv102_nitsche_pc_p2_kv0p1.toml and
  jobs/vista_nitsche/{tpv102_200m_p3.sbatch, run_nitsche_vista.sh} at b50741b; the sbatch and deck are
  byte-identical to the copies the job saved). The deck is the one of the damped p2 gold
  (nitsche-kv0p1_p2_cd_pc_200m_jobid7995488_git2e0c04a, Frontera, stopped at 9 s) plus the [output.stations]
  block added since (last changed in c906a8a); the job sets the order (3) and the 200 m mesh, and runs to the
  deck's own 12 s. Rate and state, aging law (a 0.008 in
  the core, 0.016 outside, b 0.012, Dc 0.02 m), 75 MPa shear and 120 MPa normal prestress, 25 MPa compact-bell
  nucleation; Prakash-Clifton v* 1 m/s, L 0.01 m; [numerics] kv_eta 0.1, kv_scale "node", kv_length
  "inscribed"; cfl 0.5; fault output every 0.1 s.
- Mesh: tpv104/mesh/tpv104_200m.msh, made on Vista by build_verify_vista.sbatch (gmsh 4.15.0, msh 2.2):
  2,466,549 tets, 37,928 fault faces, sha256 10cecfd73d0c7cb93cfa188aac08e4a99a6ec81d439a3363bd7175e3a8ae4aa0
  (not in the repo). tpv104_200m.geo equals tpv102_200m.geo apart from comments, so this is the Vista mesh of
  tpv102_200m.geo, the same file as the TPV104 p1, p2 and p3 golds. The earlier TPV102 golds used the Frontera
  mesh of the same .geo (2,466,377 tets, sha256 eea25a27...). The mesh is not mirror-symmetric about the fault.
- Run: Vista gg, 16 nodes / 2048 ranks (128 per node), account EAR24034, 2026-10-18 11:34 pm to 2026-10-19
  1:26 am PDT (1:51:27 job wall, 9.8 SU), 0-12 s, 176,571 steps, dt 6.79616e-5 s, central difference ("cd"),
  6668 s driver wall (38 s of it fault output), max|V| 21.03 m/s, max|V| dt / Dc 0.071, no NaN.
- results/: the 9 SCEC station files (t slip1 slip2 V1 V2 tau1 tau2 sigma_n log10_theta; 1 = dip, 2 = strike;
  SI; sigma_n compression positive; tpv102_station_conventions.txt), the station snap list, driver.log,
  RUN_INFO.txt, the job's stdout tpv102_200m_p3_1080363.out (and .err, a module notice only) and meminfo_kB.log
  (first node, every 2 s: epoch s, MemTotal kB, MemAvailable kB), the sources of the speed numbers below.
  stations.sha256 holds the sha256 of each station file as written on Vista (the files were moved in three parts
  each and checked against it after joining). Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  fault.vtkhdf (120 frames, Vista scratch) is not here.
- plots/: tpv102_200m_p3_vista1080363_stations.pdf and its page PNGs (slip-rate overview of the 9 stations, then
  per station slip rate, front zoom, slip, shear stress, normal stress, state column; this run in red vs DRDG3D,
  PyLith, the mixed-flux p2 gold, the damped p1 200 m gold (Frontera 7995743) and the damped p2 200 m gold
  (Frontera 7995488, to t = 9 s) in green); compare.txt / compare_summary.json from
  compare_nitsche_fault_stations.py vs DRDG3D, PyLith, both mixed-flux golds and the damped p1 200 m gold over
  the full 0-12 s (command below). The p2 gold ends at 9 s, so it is plotted only.
- vs DRDG3D (median |misfit| of 9 stations): arrival 0.027 s (late at 7 of 9), final slip 0.16 % (low at 6 of 9),
  peak h-slip-rate 2.7 %; noise RMS sigma_n / tau_h / tau_v 0.023 / 0.025 / 0.033 MPa (DRDG3D 0.016 / 0.035 /
  0.050). Arrival within 0.1 s everywhere PASS (worst +0.053 s at x2 = 0, 12 km depth). Slip within 2 %
  everywhere FAILS only at x2 = 12 km and -12 km, 3 km depth, both 2.05 % low (PyLith 1.82 % and 1.85 % low,
  mixed-flux p2 1.89 % and 1.90 % low); every other station is within 0.5 %. Peak within 10 % at 8 of 9; the
  miss is x2 = 12 km, 3 km depth, 11.9 % high (mixed-flux p2 16.2 % high, PyLith 4.9 % low). For scale, PyLith
  vs DRDG3D: max |dt| 0.172 s, max |slip| 1.85 %, peak within 10 % at 6 of 9. The damped p1 200 m gold has
  medians 0.141 s, 0.40 % and 31.8 % (peaks 17 to 38 % low at 8 of 9 stations).
- Fault-wide normal stress (driver.log): 55.5-189.7 MPa over the run, against 120 MPa at rest; the range was
  97-138 MPa at t = 8 s and 57-182 MPa at t = 10 s, and widened slowly to 12 s; the stations see 0.023 MPa RMS.
  Cause as in the TPV104 p2 gold README: the tets on the two sides of each fault face differ in shape on this
  gmsh mesh.

## Speed

From results/tpv102_200m_p3_1080363.out, driver.log and meminfo_kB.log; job wall from sacct; Vista gg charges
0.33 SU per node-hour. No earlier TPV102 run on Vista exists (the TPV102 golds ran on Frontera), so the second
column is the TPV104 p3 gold, the same mesh, ranks, order and code.

| | TPV102 p3, this run (b50741b) | TPV104 p3 gold 1080362 (b50741b) |
|---|---|---|
| steps to 12 s, dt | 176,571, 6.796e-5 s | 176,571, 6.796e-5 s |
| job wall (sacct) | 1:51:27 | 1:54:18 |
| driver wall | 6668 s (37.8 ms per step) | 6839 s (38.7 ms per step) |
| time in ode.Step per rank, slowest | 6186 s = 35.03 ms per step | 6313 s = 35.75 ms per step |
| time in ode.Step per rank, min / mean | 28.39 / 32.38 ms per step (max / mean 1.08) | 29.65 / 33.30 ms per step |
| step 1 (with operator assembly) | 8.5 s | 37.7 s |
| setup plus the first 1000 steps | 146 s (driver.log wall at step 1000) | 190 s |
| fault work, busiest rank | 6.00 ms per step | 6.08 ms per step |
| fault output (rank 0) | 120 frames, 38 s (0.32 s per frame) | 600 frames, 192 s (0.32 s per frame) |
| node-hours, SU | 29.7, 9.8 | 30.5, 10.1 |
| first-node memory peak | 171.6 of 237.6 GB | 177.7 of 237.6 GB |
| output folder on Vista | 2.1 GB | 7.8 GB |

- Most of the 171 s of driver wall this run saves against TPV104 is fault output (120 frames at 0.1 s against
  600 at 0.02 s); the per-step cost is the same to 2 %.
- The job script estimated 47 to 58 ms per step (as for TPV104); Vista came in at 35.03 ms.
- The step-1 gap to TPV104 (8.5 against 37.7 s on the same mesh and ranks) is not investigated.
- The driver's "time in ode.Step" excludes setup; the slowest rank sets the pace.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-benchmarks \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-p3
    cd seas-mfem-p3 && git checkout b50741be112a90e98a792abcc69202263cba5bdf
    cd miniapps/seas/jobs/vista_nitsche && sbatch build_verify_vista.sbatch   # build, tpv104_200m.msh, gate tests
    sbatch tpv102_200m_p3.sbatch

build_verify_vista.sbatch generates tpv104/mesh/tpv104_200m.msh from the committed tpv104/mesh/tpv104_200m.geo with
Vista's gmsh 4.15.0; gmsh output depends on the platform, so compare its sha256 with the one above (the job refuses
a mesh whose sha256 is not listed). The job refuses to start unless the checkout is clean, HEAD is on GitHub and
the driver was built at HEAD. Compare (from miniapps/seas, with D = this folder; rerunning it on these results/
reproduces compare.txt apart from its last line, the path it wrote the JSON to, and compare_summary.json to
round-off; checked 2026-10-19):

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv102/stations/tpv102_fault_stations.csv \
        --prefix tpv102 --run vista_p3_200m=$D/results \
        --ref drdg3d=tpv102/benchmark_data/scec_drdg3d --ref pylith=tpv102/benchmark_data/scec_pylith \
        --gold mixedflux_p1=tpv102/gold/mixedflux_p1_rk45_jobid7765843_git04674fc/results:nitsche \
        --gold mixedflux_p2=tpv102/gold/mixedflux_p2_rk45_jobid7765373_git04674fc/results:nitsche \
        --gold nitsche_p1_200m=tpv102/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid7995743_git2e0c04a/results:nitsche \
        --json compare_summary.json > compare.txt

- Requested by Chunhui Zhao (TPV205/TPV104 benchmarks, 2026-10-19).
