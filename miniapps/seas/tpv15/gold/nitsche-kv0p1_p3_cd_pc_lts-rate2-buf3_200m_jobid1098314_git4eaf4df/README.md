
# TPV15 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p3, LTS, 200 m (Vista job 1098314)

- Code: seas-mfem branch nitsche-dynamic-rupture-complex-geometries, commit
  4eaf4df769dfb4eb8279deb02f4efa75dc1fe804 (4eaf4df, on GitHub), built on TACC Vista by
  jobs/vista_complex/build_verify_complex_vista.sbatch (job 1086877, 2026-10-21, clean checkout, every check
  passed). Driver seas_nitsche_dyn_driver with --lts rate2 (4 clusters, 3 buffer layers). The branch's later
  commits up to 311658f5 do not touch the driver or this deck; 94dede26 adds the opt-in [output.stations] dt
  (and dt = "1e-4s" in the deck), which this run predates.
- Deck and job script: jobs/ (= tpv15/configs/tpv15_nitsche.toml and jobs/vista_complex/{tpv15_200m_p3.sbatch,
  run_complex_vista.sh} at 4eaf4df; the deck copy the job wrote is byte-identical). Main fault on Y = 0,
  X in [-16, 12] km; branch at 30 degrees from 100 m to 12 km along it from the junction (the spec's gap; 50 m
  off the main plane); both 15 km deep. Per-fault prestress by [[stress_region]] boxes, slip-path weakening O3,
  zero-strength tension floor O2, Kelvin-Voigt eta_bar 0.1 and Prakash-Clifton (v_star 1 m/s, L 0.01 m); no
  plasticity.
- Mesh: tpv15/mesh/tpv15_200m.msh, made on Vista by build_verify_complex_vista.sbatch from the committed .geo
  (gmsh 4.15.0, msh 2.2): 2,529,921 tets, 55,185 fault faces, graded to 50 m cells at the junction gap, sha256
  5b9f7f0a03ae996291e2829fdd765314332c20e9ec0cfa81c0feb4292125ecec (not in the repo).
- Run: Vista gg, 24 nodes / 2304 ranks (96 per node), account EAR24034, started 2026-10-24 4:51 pm PDT
  (4:04:38 job wall, 32.3 SU), 0-12 s. LTS: dt0 = 7.99268e-06 s (the global CFL step); the fault (110,370 tets
  plus 210,616 buffer tets) sits in cluster 0, so it steps at dt0: 1,501,376 steps. Driver wall 14,658.6 s
  (221.8 s of it fault output, 600 frames). max|V| 12.66 m/s, no NaN. The tension floor acted on 970 M
  DOF-steps (sigma_n_eff down to -199.8 MPa, which the friction law saw as 0); where on the fault that tension
  sits was not checked (it is at none of the stations).
- results/: the 14 SCEC station files (8 on the main fault, 6 on the branch; tpv15_station_<main|branch>_x2_<a>_
  x3_<b>.dat, x2 = km from the junction along each fault, x3 = depth km; t slip1 slip2 V1 V2 tau1 tau2 sigma_n
  slip_path; 1 = dip, 2 = strike; SI; sigma_n compression positive; tpv15_station_conventions.txt), the snap list,
  driver.log, RUN_INFO.txt, the job's .out/.err and meminfo_kB.log.
  - THINNED: the run wrote a row every fault step, 1,501,377 rows per station (3.1 GB in all). The files here
    keep every 13th row plus the header and last rows (115,492 rows, about 0.104 ms apart, 240 MB), thinned on
    Vista with the awk command in THIN_INFO.txt. THIN_INFO.txt holds the sha256 of every original and thinned
    file; stations.sha256 holds the thinned files' sha256 as written on Vista. The originals stay on Vista
    scratch, which can be purged. SeisSol's output is 0.32 ms apart, so the comparison is unaffected.
  - Left out: the per-station DG copy envelopes (*_copies.dat) and fault.vtkhdf (both on Vista scratch, 15 GB).
- reference_seissol_gap/: the SeisSol TPV15 run that keeps the junction gap (CVWS user pelties, SeisSol rev 810,
  300 m, 18-Dec-2012), 13 station files exactly as downloaded from SCEC CVWS on 2026-10-26, with MANIFEST.txt
  (source, original names, sha256). SeisSol has no file for branch x2 = 9 km, x3 = 0.
- plots/: tpv15_200m_p3_vista1098314_vs_seissolgap_stations.pdf (17 pages) and its PNGs: misfit maps, an overview
  per fault, then per station (a) strike slip rate, (b) front zoom, (c) strike slip, (d) strike shear stress,
  (e) normal stress and (f) dip slip rate; this run in red against SeisSol with gap in black.
  plot_complex_stations_p123.py v4 made them. compare.txt / compare_summary.json: compare_nitsche_fault_stations.py
  at 6e35a58, this run against SeisSol with gap. For the record: compare_p123_vs_seissol_nogap.* (p1, p2, p3
  against the SeisSol no-gap run used before 2026-10-26), refs_vs_faultmod50.txt (each reference against
  FaultMod 50 m) and cvws_tpv15_inventory.txt (the four CVWS TPV15 uploads).
- vs SeisSol with gap (median |misfit| over the stations where both rupture): arrival 0.0095 s, final slip
  2.14 %, peak strike-slip rate 6.07 %.
  0. signs of final slip and t = 0 stresses agree: PASS;
  1. arrival within 0.1 s: every measured station is within 0.02 s; the script reports FAIL only because SeisSol
     lacks branch x2 = 9 km, x3 = 0;
  2. slip within 2 %: FAIL. Branch 1.2 to 3.6 % low, main fault before the junction 1.2 to 1.3 % low, and main
     x2 = 2 km, x3 = 0 +18.4 % (0.78 m against 0.66 m, same arrival and peak);
  3. peak within 10 % at 12 of 13 stations (miss: branch x2 = 2 km, x3 = 0, -13.5 %).
  The references disagree with each other by similar amounts (refs_vs_faultmod50.txt): FaultMod 100 m and
  SeisSol with gap differ from FaultMod 50 m by about 2 to 6 % in slip, and at main x2 = 2 km the slip ranges
  over 0.015 m (FaultMod 50 m), 0.26 m (FaultMod 100 m), 0.66 m (SeisSol with gap) and 0.78 m (this run).
- driver.log prints a "Y=0 Face Classification Audit" warning (325 y = 0 faces with no BC): the main fault covers
  only part of the plane y = 0, so those are ordinary interior faces.

## Speed

From results/tpv15_200m_p3_1098314.out, driver.log and meminfo_kB.log; job wall from sacct; Vista gg charges
0.33 SU per node-hour. The p2 column is run 1098313, same build, deck, mesh, LTS and layout (not here).

| | p3, this run | p2 run 1098313 |
|---|---|---|
| nodes x ranks | 24 x 96 | 24 x 96 |
| dt0, steps to 12 s | 7.993e-6 s, 1,501,376 | 1.352e-5 s, 887,608 |
| job wall (sacct) | 4:04:38 | 58:11 |
| driver wall | 14,658.6 s | 3,471.9 s |
| time in ode.Step per rank, slowest | 11,989.2 s = 7.99 ms per step | 2,353.1 s = 2.65 ms per step |
| time in ode.Step, min / mean | 4.91 / 7.39 ms per step | 1.62 / 2.36 ms per step |
| fault work, busiest rank | 0.59 ms per step | |
| fault output (rank 0) | 600 frames, 221.8 s | 600 frames, 226.5 s |
| node-hours, SU | 97.9, 32.3 | 23.3, 7.7 |
| first-node memory peak in use | 100.3 of 237.5 GiB (114 s after start) | 104.7 GiB |
| output folder on Vista | 15 GB | 13 GB |

- p3 takes 1.69 x the p2 steps at 3.0 x the p2 cost per step. The step count is set by the 50 m cells at the
  junction, which put the fault in the finest LTS cluster.
- The run script logs MemTotal and MemAvailable only, not Committed_AS.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-dynamic-rupture-complex-geometries \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-complex-geometries
    cd seas-mfem-complex-geometries && git checkout 4eaf4df769dfb4eb8279deb02f4efa75dc1fe804
    cd miniapps/seas/jobs/vista_complex && sbatch build_verify_complex_vista.sbatch
    sbatch tpv15_200m_p3.sbatch

Compare the generated mesh's sha256 with the one above (the run script refuses others). Thin the station files
with the awk command in results/THIN_INFO.txt. The comparison and plots were made on Chunhui Zhao's Mac on
2026-10-26 from miniapps/seas at 6e35a58 (D = this folder):

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv15/stations/tpv15_fault_stations.csv \
        --prefix tpv15 --run p3=$D/results --ref SeisSolGap=$D/reference_seissol_gap --by-name --signed \
        --t-end 12 --json $D/plots/compare_summary.json > $D/plots/compare.txt
    python3 $D/plots/plot_complex_stations_p123.py --seas . --prob tpv15 --comp h --t-end 12 \
        --tag tpv15_200m_p3_vista1098314_vs_seissolgap --out $D/plots \
        --ref "SeisSolGap|$D/reference_seissol_gap" \
        --run "p3|Nitsche p=3 200 m LTS (Vista 1098314)|$D/results" --cmp $D/plots/compare.txt

- Requested by Chunhui Zhao (complex fault geometries, 2026-11-01): "gold the seissol with gap and mfem results".
