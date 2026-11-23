# TPV31 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p3, 200 m (Vista job 1096106)

- Code: seas-mfem branch nitsche-benchmarks, commit 8eef3333058893a47e875c2c9249661f83cfe771 (8eef333, on GitHub),
  built on TACC Vista by jobs/vista_nitsche/build_verify_vista.sbatch (job 1096104, 2026-10-24, clean checkout,
  BUILD_COMMIT_VISTA = 8eef333, every gate suite and smoke run passed), the same build as the p2 gold
  (nitsche-kv0p1_p2_cd_pc_200m_jobid1096105_git8eef333, which says what 8eef333 adds to the p1 gold commit). The
  Nitsche speedups are on by default (the run logs "fast fault on, fused K/K_eta off, dense blocks on"). Driver
  seas_nitsche_dyn_driver.
- Deck and job script: jobs/ (= tpv31/configs/tpv31_nitsche.toml and jobs/vista_nitsche/{tpv31_200m_p3.sbatch,
  run_nitsche_vista.sh} at 8eef333; the deck is unchanged since 9607712). Same deck and mesh as the p1 gold
  (nitsche-kv0p1_p1_cd_pc_200m_jobid1085893_git9607712, which lists the TPV31 setup) and the same deck, mesh, build
  and 16 x 144 layout as the p2 gold; only the order differs.
- Mesh: tpv31/mesh/tpv31_nitsche_200m.msh, made on Vista by build_verify_vista.sbatch from the committed .geo
  (gmsh 4.15.0, built-in kernel, msh 2.2): 2,029,790 tets, 26,248 fault faces, sha256
  5bd0b52ddd5080c531e2fc1e96b7b5d86edee0b8957c55f11e7d5294c02b9a71 (not in the repo), the same file as the p1 and
  p2 golds.
- Run: Vista gg, 16 nodes / 2304 ranks (144 per node, 126.9 k tets per node), account EAR24034, 2026-10-24
  5:00 to 6:50 am PDT (1:50:08 job wall, 9.7 SU), 0-15 s, 211,157 steps, dt 7.10374899e-5 s, central difference
  ("cd"), 6584.6 s driver wall (61.5 s of it fault output), max|V| 3.051 m/s, no NaN. The first try at this layout
  (Vista 1087802, 9607712) failed in setup when every rank read the serial mesh; the rank 0 mesh setup fixed that.
- results/: the 30 SCEC station files (strike 0, 6, 12 km x depths 0, 0.2, 0.5, 1, 2.4, 3, 5, 7.5, 10, 12 km;
  t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip, 2 = strike; SI; sigma_n compression positive;
  tpv31_station_conventions.txt), the station snap list (nearest fault node, up to 45 m from the requested point),
  driver.log, RUN_INFO.txt, the job's stdout tpv31_200m_p3_1096106.out (and .err, a module notice only) and
  meminfo_kB.log (first node, every 2 s: epoch s, MemTotal, MemAvailable, CommitLimit, Committed_AS, all kB), the
  sources of the speed numbers below. stations.sha256 holds the sha256 of each station file as written on Vista
  (the files came off Vista through the Mac in 15 MB parts and were checked against it after joining). Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  fault.vtkhdf (299 frames, Vista scratch) is not here.
- plots/: tpv31_200m_p3_vista1096106_stations.pdf and its page PNGs (misfit maps of arrival, final slip and peak
  slip rate over the 30 stations for this run, the p2 gold and the p1 gold, then per station slip rate, front zoom,
  slip, shear stress, normal stress and shear/normal ratio; this run in red, drawn last, over the p2 gold in green
  and the p1 gold in orange, vs SeisSol in black and EQdyna in gray); compare.txt / compare_summary.json from
  compare_nitsche_fault_stations.py vs SeisSol and EQdyna with the p1 and p2 golds as the runs before this one
  (command below). SeisSol ends at 14.9964 s, so every source is measured to there. compare.txt equals the job's
  comparison on Vista except for the pass-3 line, which now gives the 24-station requirement, and the last line
  (the JSON path).
- vs SeisSol (median |misfit| of 30 stations): arrival 0.0083 s (early at 25 of 30), final slip 0.42 % (low at 24),
  peak h-slip-rate 2.45 % (low at 19); noise RMS sigma_n / tau_h / tau_v 0.003 / 0.003 / 0.005 MPa (SeisSol tau_h /
  tau_v 0.004 / 0.005). PASSES all four pass tests:
  1. arrival within 0.1 s everywhere (worst -0.048 s at x2 = 12 km, 1 km depth);
  2. slip within 2 % everywhere (worst -1.67 % at x2 = 0, 7.5 km depth);
  3. peak within 10 % at 29 of 30 (need 24); the one miss is x2 = 0 at 0.2 km depth, 11.5 % low (p2 gold 10.7 %
     high, p1 gold 5.2 % high, EQdyna 1.9 % high);
  4. every median below the p2 gold's (0.010 s, 0.91 %, 5.28 %).
  The misfit no longer grows along strike: median arrival -0.007 s at 0, 6 and 12 km, peak +0.0 / -0.9 / -1.7 %.
  The p2 gold's two slip misses (2.06 % low) are 1.02 % and 1.56 % low here, and its five peak misses at 12 km
  along strike are all within 7.6 %. For scale, EQdyna vs SeisSol: max |dt| 0.030 s, max |slip| 2.02 %, peak within 10 % at
  25 of 30; this run's arrival is further than 0.030 s from SeisSol at 6 of 30 stations.
- Fault-wide normal stress (driver.log): 24.5-81.1 MPa at rest (the mu/mu_0 scaling), 12.4-103.1 MPa over the run
  (p2 gold 12.8-103.1 MPa), reached within 2 % by t = 8.2 s and held to 15 s.

## Speed

From results/tpv31_200m_p3_1096106.out, driver.log and meminfo_kB.log; job wall from sacct; Vista gg charges 0.33 SU
per node-hour. The p2 column is the p2 gold, same build and layout.

| | p3, this run (8eef333) | p2 gold 1096105 (8eef333) |
|---|---|---|
| nodes x ranks | 16 x 144 | 16 x 144 |
| steps to 15 s, dt | 211,157, 7.104e-5 s | 124,835, 1.202e-4 s |
| job wall (sacct) | 1:50:08 | 24:15 |
| driver wall | 6584.6 s (31.18 ms per step) | 1433.6 s (11.48 ms per step) |
| time in ode.Step per rank, slowest | 6542.9 s = 30.99 ms per step | 11.19 ms per step |
| time in ode.Step per rank, min / mean | 26.64 / 30.96 ms per step | 9.58 / 11.18 ms per step |
| step 1 (with operator assembly) | 8.03 s | 2.47 s |
| setup (slowest rank) | 38.8 s (rank 0 mesh read 17.0 s, ParMesh 12.5 s, partition 6.1 s) | 38.9 s |
| setup plus the first 1000 steps | 92.5 s (driver.log wall at step 1000) | 49.5 s |
| fault work, busiest rank | 6.28 ms per step | 1.75 ms per step |
| fault output (rank 0) | 299 frames, 61.5 s (0.206 s per frame) | 300 frames, 62.6 s |
| dense blocks, all ranks | 314.8 GB | 79.1 GB |
| node-hours, SU | 29.4, 9.7 | 6.47, 2.1 |
| first-node memory peak in use | 132.3 of 237.6 GiB | 62.3 of 237.6 GiB |
| first-node memory peak committed | 115.2 of 118.8 GiB (step 1, operator assembly) | 51.8 of 118.8 GiB |
| output folder on Vista | 4.1 GB | 3.3 GB |

- p3 takes 1.69 x the p2 steps (dt shrinks by C(3)/9 against C(2)/4 = 0.591) at 2.77 x the p2 cost per step
  (TPV205 p3 against p2 with the speedups: 2.91 x), so 4.5 x the p2 job wall.
- The committed-memory peak (the first step's operator assembly; steady state about 68 GiB) came within 3.6 GiB of
  Vista's 118.8 GiB commit limit, past which compute nodes refuse allocations. Only the first node is logged, and its
  share of the mesh need not be the largest, so this layout has no margin at p3: fewer than 16 nodes, or a finer
  mesh per node, will not fit.
- The driver's "time in ode.Step" excludes setup; the slowest rank sets the pace.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-benchmarks \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-tpv6731
    cd seas-mfem-tpv6731 && git checkout 8eef3333058893a47e875c2c9249661f83cfe771
    cd miniapps/seas/jobs/vista_nitsche && sbatch build_verify_vista.sbatch   # build, tpv31_nitsche_200m.msh, gate tests
    sbatch tpv31_200m_p3.sbatch

build_verify_vista.sbatch generates tpv31/mesh/tpv31_nitsche_200m.msh from the committed .geo with Vista's gmsh
4.15.0; gmsh output depends on the platform, so compare its sha256 with the one above (the run script refuses a
mesh whose sha256 is not in MESH_SHA256_VISTA). The job refuses to start unless the checkout is clean, HEAD is on
GitHub and the driver was built at HEAD. Compare (from miniapps/seas, with D = this folder and P1 / P2 = the p1 /
p2 gold folders; rerunning it on these results/ reproduces compare.txt and compare_summary.json; checked
2026-10-25):

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv31/stations/tpv31_fault_stations.csv \
        --prefix tpv31 --run p1_200m=$P1/results --run p2_200m=$P2/results --run p3_200m=$D/results \
        --ref seisol=tpv31/benchmark_data/scec_seisol --ref eqdyna=tpv31/benchmark_data/scec_eqdyna \
        --peak-min-stations 24 --t-end 15 --json $D/plots/compare_summary.json > $D/plots/compare.txt

- Requested by Chunhui Zhao (TPV6/TPV7/TPV31 benchmarks, 2026-10-25).
