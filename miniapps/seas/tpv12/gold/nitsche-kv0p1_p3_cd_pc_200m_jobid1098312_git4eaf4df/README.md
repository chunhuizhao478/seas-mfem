# TPV12 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p3, 200 m (Vista job 1098312)

- Code: seas-mfem branch nitsche-dynamic-rupture-complex-geometries, commit
  4eaf4df769dfb4eb8279deb02f4efa75dc1fe804 (4eaf4df, on GitHub), built on TACC Vista by
  jobs/vista_complex/build_verify_complex_vista.sbatch (job 1086877, 2026-10-21, clean checkout,
  BUILD_COMMIT_VISTA = VERIFIED_COMMIT_COMPLEX_VISTA = 4eaf4df, every check passed). The Nitsche speedups are on
  by default (the run logs "fast fault on, fused K/K_eta off, dense blocks on"). Driver seas_nitsche_dyn_driver.
  The branch's next commit, 6e35a58, changes only scripts/compare_nitsche_fault_stations.py.
- Deck and job script: jobs/ (= tpv12/configs/tpv12_nitsche.toml and jobs/vista_complex/{tpv12_200m_p3.sbatch,
  run_complex_vista.sh} at 4eaf4df). The deck sets the TPV12 depth-dependent effective prestress through option O1
  [prestress_depth_table], the zero-strength tension floor O2 ([fault_tension] mode = "zero_strength_floor"),
  slip-path weakening O3 ([slip_weakening_distance] measure = "slip_path"), Kelvin-Voigt eta_bar 0.1 (scale node,
  inscribed length) and the Prakash-Clifton regularization (v_star 1 m/s, L 0.01 m); no plasticity.
- Mesh: tpv12/mesh/tpv12_200m.msh, made on Vista by build_verify_complex_vista.sbatch from the committed .geo
  (gmsh 4.15.0, msh 2.2): 1,880,734 tets, 25,528 fault faces, sha256
  de6ed89eb495194ba1b0c3d88b7b78ea9871bf0a6d38173bc4137b3c57376319 (not in the repo).
- Run: Vista gg, 16 nodes / 2304 ranks (144 per node, 117.5 k tets per node, every rank reads the serial mesh on
  this branch), account EAR24034, 2026-10-24 3:21 to 4:13 pm PDT (51:45 job wall, 4.6 SU), 0-8 s, 94,738 steps,
  dt 8.4444001e-5 s, central difference ("cd"), 3065.9 s driver wall (87.8 s of it fault output), max|V|
  49.50 m/s, no NaN. The .err holds "cannot lock ref" lines from git fetch: the six p2/p3 jobs started together
  in one checkout, and run_complex_vista.sh retries the fetch for that reason; the commit guard then passed.
- results/: the 10 SCEC station files (x2 = 0, 4.5, 12 km at x3 = 0 and 7.5 km, x2 = 0 at x3 = 1.5, 3, 4.5,
  12 km; x3 is distance down dip; t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip, 2 = strike; SI;
  sigma_n compression positive; tpv12_station_conventions.txt), the station snap list (nearest fault node, up to
  39.5 m from the requested point), driver.log, RUN_INFO.txt, the job's stdout tpv12_200m_p3_1098312.out (and
  .err) and meminfo_kB.log (first node, every 2 s: epoch s, MemTotal, MemAvailable, kB). stations.sha256 holds the
  sha256 of each station file as written on Vista (the files came off Vista through the Mac and were checked
  against it). Read with `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  Left out: the per-station DG copy envelopes (*_copies.dat, 10 x 18 MB) and fault.vtkhdf (399 frames); both
  are on Vista scratch, the envelopes also on Chunhui Zhao's Mac.
- plots/: tpv12_200m_p3_vista1098312_stations.pdf (12 pages) and its page PNGs: misfit maps of arrival, final
  slip and peak dip-slip rate at the 10 stations for this run, the p2 run (Vista 1098311, 200 m) and the p1 run
  (Vista 1086893, 100 m), then an overview and per station (a) dip slip rate, (b) front zoom, (c) dip slip,
  (d) dip shear stress, (e) normal stress and (f) strike slip rate; this run in red, drawn last, over p2 in green
  and p1 in orange dashed, vs DRDG3D in black. plot_complex_stations_p123.py made them (the benchmarks thread's
  TPV31 p1/p2/p3 station pages with input changes, listed in its docstring). compare.txt / compare_summary.json
  from compare_nitsche_fault_stations.py at 6e35a58 vs DRDG3D with p1 and p2 as the runs before this one
  (command below). DRDG3D ends at 7.9999 s, so every source is measured to there. The p1 and p2 runs are not
  golds and their station files are not in the repo.
- vs DRDG3D (median |misfit| of 10 stations): arrival 0.0053 s, final slip 0.04 %, peak dip-slip rate 2.07 %
  (p2 0.0060 s, 0.18 %, 6.64 %; p1 at 100 m 0.0254 s, 0.79 %, 18.60 %).
  0. signs of final slip and the t = 0 stresses agree everywhere: PASS;
  1. arrival within 0.1 s everywhere: PASS (worst -0.090 s at x2 = 12 km, x3 = 7.5 km);
  2. slip within 2 % where DRDG3D slips more than 0.1 m: PASS (worst -0.36 %, same station);
  3. peak within 10 % at 9 of 10 stations; the miss is x2 = 12 km, x3 = 7.5 km, 13.2 % low (p2 38.4 %, p1 54.1 %
     low), 3 km from the fault's side edge, where DRDG3D's peak is a 39.1 m/s spike;
  4. every median below p2's and p1's: PASS.
  No second TPV12 reference is in hand to give DRDG3D's own spread.
- Fault-wide (driver.log): max|V| 49.50 m/s (p2 36.4, p1 at 100 m 25.9) and sigma_n_eff down to -368 MPa (p2 -282,
  p1 -180) grow with the order. The tension floor held the law's sigma_n at 0 on 78.9 M DOF-steps. A 1000 m p1
  run (cloud, 2026-10-22) placed the tension at the corner where the fault's side edge meets the surface trace
  (X = 15 km), a trace-end singularity; this run was not checked for where its extremes sit. The stations,
  which are away from that corner, converge (medians above).
- driver.log prints a "Y=0 Face Classification Audit" warning (345 y = 0 faces with no BC). The audit is written
  for a vertical fault in the plane y = 0; the TPV12 fault meets y = 0 only along its surface trace, so these are
  ordinary interior faces of that plane.

## Speed

From results/tpv12_200m_p3_1098312.out, driver.log and meminfo_kB.log; job wall from sacct; Vista gg charges
0.33 SU per node-hour. The p2 column is the p2 run, same build, deck, mesh and layout (from its .out, not here).

| | p3, this run (4eaf4df) | p2 run 1098311 (4eaf4df) |
|---|---|---|
| nodes x ranks | 16 x 144 | 16 x 144 |
| steps to 8 s, dt | 94,738, 8.444e-5 s | 56,009, 1.43e-4 s |
| job wall (sacct) | 51:45 | 14:19 |
| driver wall | 3065.9 s (32.36 ms per step) | |
| time in ode.Step per rank, slowest | 2573.1 s = 27.16 ms per step | 10.2 ms per step |
| time in ode.Step per rank, min / mean | 21.79 / 25.39 ms per step | |
| step 1 (with operator assembly) | 8.53 s | |
| setup plus the first 1000 steps | 123.7 s (driver.log wall at step 1000) | |
| fault work, busiest rank | 3.57 ms per step | |
| fault output (rank 0) | 399 frames, 87.8 s (0.220 s per frame) | |
| dense blocks, all ranks | 292.3 GB | |
| node-hours, SU | 13.8, 4.6 | 3.8, 1.3 |
| first-node memory peak in use | 162.4 of 237.6 GiB (89 s after start, step 1) | 141.8 GiB |
| output folder on Vista | 3.5 GB | |

- p3 takes 1.69 x the p2 steps at 2.66 x the p2 cost per step, so 3.6 x the p2 job wall. The job script's
  estimate was 4.0 SU.
- This branch's run script logs MemTotal and MemAvailable only, not Committed_AS, so the margin to Vista's
  118.8 GiB per-node commit limit is not known; the run did not reach it.
- The driver's "time in ode.Step" excludes setup; the slowest rank sets the pace.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-dynamic-rupture-complex-geometries \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-complex-geometries
    cd seas-mfem-complex-geometries && git checkout 4eaf4df769dfb4eb8279deb02f4efa75dc1fe804
    cd miniapps/seas/jobs/vista_complex && sbatch build_verify_complex_vista.sbatch   # build, meshes, checks
    sbatch tpv12_200m_p3.sbatch

build_verify_complex_vista.sbatch generates tpv12/mesh/tpv12_200m.msh from the committed .geo with Vista's gmsh
4.15.0; gmsh output depends on the platform, so compare its sha256 with the one above (the run script refuses a
mesh whose sha256 is not in MESH_SHA256_COMPLEX_VISTA). The job refuses to start unless the checkout is clean,
HEAD is on GitHub and the driver was built and verified at HEAD.

The DRDG3D reference is not in the repo (tpv12/benchmark_data/README.md). The comparison and the plots were made
on Chunhui Zhao's Mac on 2026-10-25 (Python with numpy and matplotlib) from miniapps/seas at 6e35a58, with the
commands below; there, D = this folder, REF = the DRDG3D station folder, and P1 / P2 = the p1 / p2 station
folders stand in for the Mac paths that compare.txt prints. Its "# wrote" line names the JSON file before it was
renamed to compare_summary.json.

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv12/stations/tpv12_fault_stations.csv \
        --prefix tpv12 --run p1=$P1 --run p2=$P2 --run p3=$D/results --ref DRDG3D=$REF --by-name --signed \
        --t-end 8 --json $D/plots/compare_summary.json > $D/plots/compare.txt
    python3 $D/plots/plot_complex_stations_p123.py --seas . --prob tpv12 --comp v --t-end 8 \
        --tag tpv12_200m_p3_vista1098312 --out $D/plots --ref "DRDG3D|$REF" \
        --run "p1|Nitsche p=1 100 m (Vista 1086893)|$P1" --run "p2|Nitsche p=2 200 m (Vista 1098311)|$P2" \
        --run "p3|Nitsche p=3 200 m (Vista 1098312)|$D/results" --cmp $D/plots/compare.txt

- Requested by Chunhui Zhao (complex fault geometries, 2026-10-25).
