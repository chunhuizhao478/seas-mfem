# TPV27 gold: Nitsche dynamic rupture with Drucker-Prager viscoplasticity, Kelvin-Voigt damping, p3, 200 m (Vista job 1154788)

- Code: seas-mfem branch nitsche-dynamic-rupture-nonlinearity, commit 72121c997abec6cfe1cb09a9c7f746c47ad313a6
  (72121c9, on GitHub; also contained in scec-benchmarks since 9f562eab), built on TACC Vista by
  jobs/vista_rheology/build_tpv2627_vista.sbatch (job 1154787, 2026-11-19 10:32 to 10:45 PST, clean checkout,
  "build OK: commit 72121c9..., tracked changes: 0", "all checks passed"; gcc/14.2.0 + openmpi/5.0.11, the module
  stack after Vista removed HCOLL on 2026-11-11; its informational 200 m load check on one node was killed out of
  memory, as in every build since 1099060, and does not gate). Driver seas_nitsche_dyn_driver; the run logs "fast
  fault on, fused K/K_eta off, dense blocks on".
- Deck and job scripts: jobs/ (= tpv27/configs/tpv27_nitsche_kv0p1.toml, jobs/vista_rheology/{tpv27_200m_p3.sbatch,
  build_tpv2627_vista.sbatch} and jobs/vista_nitsche/run_nitsche_vista.sh at 72121c9; the deck is byte-identical to
  the copy the run wrote). TPV27 = the TPV26 deck (forced rupture from (-5, -10) km, depth-dependent effective stress
  with fluid pressure, linear slip weakening on slip path, Prakash-Clifton regularization v* = 1 m/s, L = 0.01 m,
  Kelvin-Voigt eta_bar = 0.1 over the whole domain) with SCEC Drucker-Prager viscoplasticity (exact exponential
  Duvaut-Lions update), the plastic state at the 10 nodes of degree 2 in each tet (representation "nodal_pm1").
- Mesh: tpv26/mesh/tpv2627_mirror_200m.msh, made on Vista by the build from the committed tpv26/mesh/tpv2627_half.geo
  (gmsh 4.15.0, the y >= 0 half mirrored across the fault): 2,997,834 tets, sha256
  ed8ce90b30a55999bd34a5ae8c7fa90da3063e8d1377f270e2c988559de4f599 (not in the repo), the same mesh as the TPV26 gold
  1099061 and TPV27 p1 1114572 / p2 1114571. The 12 on-fault stations are mesh vertices (snap 0 m).
- Run: Vista gg, 32 nodes / 4096 ranks (128 per node; 32 is the gg QOS per-job limit), account EAR24034,
  2026-11-20 04:01 to 05:55 PST (1 h 54 min job wall, about 20 SU), 0 to 13 s, 203,364 steps, dt 6.392508e-5 s,
  central difference ("cd"), 6804.4 s driver wall, exit 0, max|V| 4.806 m/s, sigma_n_eff min -1.16 MPa (near the
  free surface), no NaN. At 13 s: 241,567 yielding nodes, 436,513 elements with plastic strain, max eta 0.00688.
- results/: the 12 SCEC on-fault station files (t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip,
  2 = strike; SI; sigma_n compression positive; tpv27_station_conventions.txt), the 6 off-fault station files
  (tpv27_offfault_body*.dat), the station snap list, driver.log and RUN_INFO.txt. files.sha256 is the sha256 of each
  file as written on Vista; the files came off Vista through the Mac (38 of 38 matched there) and were checked again
  here. Not here: the per-station _copies.dat files (12 x 39 MB; on the Mac), the job's stdout and memory log, and
  fault.vtkhdf (649 frames) and rheology.vtkhdf (plastic strain and eta every 0.5 s); they stay in
  $SCRATCH/nitsche_benchmarks/tpv27_mirror200m_p3_1154788 on Vista (12 GB).
- plots/: tpv27_200m_p123_vista_stations.pdf (sha256 9dbb5c39aad5828a18806c8ce5221164218be888ee938f8751287593083949b1)
  and its 14 page PNGs: misfit maps of arrival, final slip and peak slip rate, then per station slip rate, front
  zoom, slip, shear stress, normal stress and shear/normal ratio; p1 (1114572), p2 (1114571) and this run p3 (drawn
  last) vs DRDG3D 200 m O4 in black. compare.txt and compare_summary.json from scripts/compare_nitsche_fault_stations.py
  at 6e35a58 (runs p1, p2, p3 in that order); sigma_n_late.txt, the mean sigma_n over 10.5 to 11.9997 s per station
  for p1, p2, p3 and DRDG3D (DRDG3D sign flipped to compression positive). DRDG3D ends at 11.9997 s, so every measure
  is to t = 11.9997 s. compare.txt's last line is the JSON path on the Mac where it was written.
- vs DRDG3D (median |misfit| over the 10 slipping stations; the two surface stations never slip in any run):

  | order | job | arrival | final slip | peak h-slip-rate |
  |---|---|---|---|---|
  | p1 | 1114572 | 0.510 s | 1.17 % | 16.6 % |
  | p2 | 1114571 | 0.034 s | 0.28 % | 2.9 % |
  | p3 | 1154788 (this gold) | 0.019 s | 0.31 % | 2.3 % |

  1. arrival within 0.1 s at every station (worst +0.043 s). PASS.
  2. slip within 2 % everywhere (worst -0.57 %). PASS.
  3. peak within 10 % at 11 of 12; the miss is x2 = 15 km, 15 km depth at +15.5 % (4.452 against 3.854 m/s; p2
     +5.1 %, p1 -21.3 %). The other stations are within -0.3 to +8.2 %.
  4. medians shrink p1 -> p2 -> p3 for arrival and peak; the slip median goes 0.28 % -> 0.31 %.
- Normal stress after rupture (sigma_n_late.txt): p3 is within 0.11 MPa of DRDG3D at 11 of 12 stations and +0.42 MPa
  at x2 = 15 km, 15 km depth. p2 sat up to 1.54 MPa above DRDG3D there and p1 changed sign at several stations, so
  those offsets were the resolution of the plastic zone. (TPV26, elastic on the mirror mesh, has no change in
  sigma_n.)

## Commands

Build and run (Vista login node, clean checkout of 72121c9 in $WORK/seas-mfem-nonlinearity):

```
cd $WORK/seas-mfem-nonlinearity/miniapps/seas/jobs/vista_rheology
B=$(sbatch --parsable build_tpv2627_vista.sbatch | tail -n 1)
sbatch --dependency=afterok:$B tpv27_200m_p3.sbatch
```

Comparison (Mac, pythonenv), from a checkout of miniapps/seas at 6e35a58, with P1, P2, P3 the results folders of
1114572, 1114571 and 1154788; plot_complex_stations_p123.py (sha256
2c316ee118d4adc7d8e380888beac828f41c0042122f862905c504a112d52fab) is the complex-geometries thread's multi-order
station page script, not in the repo:

```
python3 scripts/compare_nitsche_fault_stations.py --stations tpv26/stations/tpv26_fault_stations.csv \
    --prefix tpv27 --run p1=P1 --run p2=P2 --run p3=P3 --ref DRDG3D=tpv27/benchmark_data/scec_drdg3d \
    --t-end 11.9997 --json OUT/compare_summary.json > OUT/compare.txt
```

Then plot_complex_stations_p123.py on P1, P2, P3, the DRDG3D folder and OUT/compare.txt, in the same form as the
p1/p2 PDF of 2026-11-02 (Mac: tpv27_200m/p12_6e35a58); its exact argument list was not recorded with this gold.

## Speed

From results/driver.log; Vista gg charges 0.33 SU per node-hour.

| | this run (72121c9) |
|---|---|
| nodes x ranks | 32 x 128 |
| steps to 13 s, dt | 203,364, 6.393e-5 s |
| driver wall | 6804.4 s (33.5 ms per step) |
| time in ode.Step per rank, slowest | 6670.6 s = 32.80 ms per step |
| time in ode.Step per rank, min / mean | 29.19 / 32.79 ms per step |
| step 1 (with operator assembly) | 5.34 s |
| setup (slowest rank) | 60.2 s (rank 0 mesh read 23.3 s, partition 11.1 s, ParMesh 21.7 s) |
| fault work, busiest rank | 3.03 ms per step |
| fault output (rank 0) | 649 frames, 209.7 s (0.323 s per frame) |
| dense blocks, all ranks | 467.7 GB |
| first-node memory peak committed | 68.8 of 118.8 GiB (in use 99.1 GiB; the job's memory log on Vista) |

TPV27 p2 1114571 (24 x 128) ran 17.1 ms per step; per step p3 costs 2.56 x p2 per node-equivalent
(32.8 x 32 / (17.1 x 24)).
