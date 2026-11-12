# TPV31 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p2, 200 m (Vista job 1096105)

- Code: seas-mfem branch nitsche-benchmarks, commit 8eef3333058893a47e875c2c9249661f83cfe771 (8eef333, on GitHub),
  built on TACC Vista by jobs/vista_nitsche/build_verify_vista.sbatch (job 1096104, 2026-10-24, clean checkout,
  BUILD_COMMIT_VISTA = 8eef333, every gate suite and smoke run passed). 8eef333 is the p1 gold commit 9607712 plus
  the p1 gold itself (cf0a868c), the trunk merge 746b9c4f (rank 0 mesh setup: only rank 0 reads the serial mesh and
  sends each rank its part) and the committed-memory log (8eef333). The Nitsche speedups are on by default (the run
  logs "fast fault on, fused K/K_eta off, dense blocks on"). Driver seas_nitsche_dyn_driver.
- Deck and job script: jobs/ (= tpv31/configs/tpv31_nitsche.toml and jobs/vista_nitsche/{tpv31_200m_p2.sbatch,
  run_nitsche_vista.sh} at 8eef333; the deck is unchanged since 9607712). Same deck and mesh as the p1 gold
  (nitsche-kv0p1_p1_cd_pc_200m_jobid1085893_git9607712, which lists the TPV31 setup); only the order and the rank
  layout differ.
- Mesh: tpv31/mesh/tpv31_nitsche_200m.msh, made on Vista by build_verify_vista.sbatch from the committed .geo
  (gmsh 4.15.0, built-in kernel, msh 2.2): 2,029,790 tets, 26,248 fault faces, sha256
  5bd0b52ddd5080c531e2fc1e96b7b5d86edee0b8957c55f11e7d5294c02b9a71 (not in the repo), the same file as the p1 gold.
- Run: Vista gg, 16 nodes / 2304 ranks (144 per node, 126.9 k tets per node), account EAR24034, 2026-10-24
  1:37 to 2:01 am PDT (24:15 job wall, 2.1 SU), 0-15 s, 124,835 steps, dt 1.201591584e-4 s, central difference
  ("cd"), 1433.6 s driver wall (62.6 s of it fault output), max|V| 2.838 m/s, no NaN. The first try at this layout
  (Vista 1087801, 9607712) failed in setup when every rank read the serial mesh; the rank 0 mesh setup fixed that.
- results/: the 30 SCEC station files (strike 0, 6, 12 km x depths 0, 0.2, 0.5, 1, 2.4, 3, 5, 7.5, 10, 12 km;
  t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip, 2 = strike; SI; sigma_n compression positive;
  tpv31_station_conventions.txt), the station snap list (nearest fault node, up to 56 m from the requested point),
  driver.log, RUN_INFO.txt, the job's stdout tpv31_200m_p2_1096105.out (and .err, a module notice only) and
  meminfo_kB.log (first node, every 2 s: epoch s, MemTotal, MemAvailable, CommitLimit, Committed_AS, all kB), the
  sources of the speed numbers below. stations.sha256 holds the sha256 of each station file as written on Vista
  (the files came off Vista through the Mac in 15 MB parts and were checked against it after joining). Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  fault.vtkhdf (300 frames, Vista scratch) is not here.
- plots/: tpv31_200m_p2_vista1096105_stations.pdf and its page PNGs (misfit maps of arrival, final slip and peak
  slip rate over the 30 stations for this run and the p1 gold, then per station slip rate, front zoom, slip, shear
  stress, normal stress and shear/normal ratio; this run in red, drawn last, over the p1 gold in orange, vs SeisSol
  in black and EQdyna in gray); compare.txt / compare_summary.json from compare_nitsche_fault_stations.py vs SeisSol
  and EQdyna with the p1 gold as the run before this one (command below). SeisSol ends at 14.9964 s, so every source
  is measured to there. compare.txt equals the job's comparison on Vista except for the pass-3 line, which now gives
  the 24-station requirement, and the last line (the JSON path).
- vs SeisSol (median |misfit| of 30 stations): arrival 0.010 s (early at 23 of 30), final slip 0.91 % (low at 28),
  peak h-slip-rate 5.28 % (low at 21); noise RMS sigma_n / tau_h / tau_v 0.005 / 0.004 / 0.005 MPa (SeisSol tau_h /
  tau_v 0.004 / 0.005). Passes pass tests 1, 3 and 4, fails 2:
  1. arrival within 0.1 s everywhere (worst -0.096 s at x2 = 12 km, 2.4 km depth);
  2. slip within 2 % at 28 of 30: x2 = 6 km at 12 km depth and x2 = 12 km at 10 km depth are both 2.06 % low;
  3. peak within 10 % at 24 of 30 (need 24); the misses are 5 of the 10 stations at 12 km along strike (worst
     -21.1 % at 3 km depth) and x2 = 0 at 0.2 km depth (+10.7 %);
  4. every median below the p1 gold's (0.113 s, 1.73 %, 21.8 %).
  The p1 gold's along-strike growth is mostly gone: median arrival -0.005 / -0.006 / -0.010 s and peak +1.3 / -2.7 /
  -10.2 % at 0 / 6 / 12 km. For scale, EQdyna vs SeisSol: max |dt| 0.030 s, max |slip| 2.02 %, peak within 10 % at
  25 of 30.
- Fault-wide normal stress (driver.log): 24.5-81.1 MPa at rest (the mu/mu_0 scaling), 12.8-103.1 MPa over the run
  (p1 gold 14.5-100.1 MPa), reached within 2 % by t = 8.2 s and held to 15 s.

## Speed

From results/tpv31_200m_p2_1096105.out, driver.log and meminfo_kB.log; job wall from sacct; Vista gg charges 0.33 SU
per node-hour. The p1 column is the p1 gold on 4 nodes (507 k tets per node).

| | p2, this run (8eef333) | p1 gold 1085893 (9607712) |
|---|---|---|
| nodes x ranks | 16 x 144 | 4 x 144 |
| steps to 15 s, dt | 124,835, 1.202e-4 s | 58,551, 2.562e-4 s |
| job wall (sacct) | 24:15 | 9:41 |
| driver wall | 1433.6 s (11.48 ms per step) | 559.55 s (9.56 ms per step) |
| time in ode.Step per rank, slowest | 1397.3 s = 11.19 ms per step | 7.32 ms per step |
| time in ode.Step per rank, min / mean | 9.58 / 11.18 ms per step | 5.80 / 6.80 ms per step |
| step 1 (with operator assembly) | 2.47 s | 0.92 s |
| setup (slowest rank) | 38.9 s (rank 0 mesh read 17.1 s, ParMesh 12.5 s, partition 6.1 s) | |
| setup plus the first 1000 steps | 49.5 s (driver.log wall at step 1000) | 65.5 s |
| fault work, busiest rank | 1.75 ms per step | 1.06 ms per step |
| fault output (rank 0) | 300 frames, 62.6 s (0.209 s per frame) | 300 frames, 62.0 s |
| dense blocks, all ranks | 79.1 GB | |
| node-hours, SU | 6.47, 2.1 | 0.65, 0.21 |
| first-node memory peak in use | 62.3 of 237.6 GiB | 141.6 of 237.6 GiB |
| first-node memory peak committed | 51.8 of 118.8 GiB (step 1, operator assembly) | not logged |
| output folder on Vista | 3.3 GB | 2.7 GB |

- p2 takes 2.13 x the p1 steps (dt shrinks by C(2)/4 against C(1)/1 = 0.469) and, at 4 x the nodes, 1.53 x the p1
  cost per step at the slowest rank.
- The committed-memory peak is the first step's operator assembly; Vista compute nodes refuse allocations past the
  118.8 GiB commit limit, so this sets how many tets a node can take.
- The driver's "time in ode.Step" excludes setup; the slowest rank sets the pace.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-benchmarks \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-tpv6731
    cd seas-mfem-tpv6731 && git checkout 8eef3333058893a47e875c2c9249661f83cfe771
    cd miniapps/seas/jobs/vista_nitsche && sbatch build_verify_vista.sbatch   # build, tpv31_nitsche_200m.msh, gate tests
    sbatch tpv31_200m_p2.sbatch

build_verify_vista.sbatch generates tpv31/mesh/tpv31_nitsche_200m.msh from the committed .geo with Vista's gmsh
4.15.0; gmsh output depends on the platform, so compare its sha256 with the one above (the run script refuses a
mesh whose sha256 is not in MESH_SHA256_VISTA). The job refuses to start unless the checkout is clean, HEAD is on
GitHub and the driver was built at HEAD. Compare (from miniapps/seas, with D = this folder and P1 = the p1 gold
folder; rerunning it on these results/ reproduces compare.txt and compare_summary.json; checked 2026-10-25):

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv31/stations/tpv31_fault_stations.csv \
        --prefix tpv31 --run p1_200m=$P1/results --run p2_200m=$D/results \
        --ref seisol=tpv31/benchmark_data/scec_seisol --ref eqdyna=tpv31/benchmark_data/scec_eqdyna \
        --peak-min-stations 24 --t-end 15 --json $D/plots/compare_summary.json > $D/plots/compare.txt

- Requested by Chunhui Zhao (TPV6/TPV7/TPV31 benchmarks, 2026-10-25).
