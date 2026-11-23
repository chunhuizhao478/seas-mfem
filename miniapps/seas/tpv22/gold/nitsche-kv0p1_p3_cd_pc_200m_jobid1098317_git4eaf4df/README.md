
# TPV22 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p3, 200 m (Vista job 1098317)

- Code: seas-mfem branch nitsche-dynamic-rupture-complex-geometries, commit
  4eaf4df769dfb4eb8279deb02f4efa75dc1fe804 (4eaf4df, on GitHub), built on TACC Vista by
  jobs/vista_complex/build_verify_complex_vista.sbatch (job 1086877, 2026-10-21, clean checkout,
  BUILD_COMMIT_VISTA = VERIFIED_COMMIT_COMPLEX_VISTA = 4eaf4df, every check passed). The Nitsche speedups are on
  by default (the run logs "fast fault on, fused K/K_eta off, dense blocks on"). Driver seas_nitsche_dyn_driver.
  The branch's later commits (6e35a58 compare script, 1421eb9a TPV12 gold) do not touch the driver or this deck.
- Deck and job script: jobs/ (= tpv22/configs/tpv22_nitsche.toml and jobs/vista_complex/{tpv22_200m_p3.sbatch,
  run_complex_vista.sh} at 4eaf4df; the deck copy the job wrote next to its output is byte-identical). Two vertical
  right-lateral faults with an extensional stepover: fault 1 on Y = 0, X in [-25, 5] km, fault 2 on Y = -1.6 km,
  X in [-5, 25] km, both 20 km deep. The deck sets the TPV22 prestress through option O1 [prestress_depth_table]
  (total stress: sigma_yy -60 MPa, sigma_xy 29.38 MPa to 15 km depth, then linear to 14.69 MPa at 20 km), the
  zero-strength tension floor O2 ([fault_tension] mode = "zero_strength_floor"), slip-path weakening O3
  ([slip_weakening_distance] measure = "slip_path"), the forced-rupture nucleation O4 ([forced_rupture]: r_crit
  3 km around (-10, 0, 10 km depth) on fault 1, 2424.8 m/s, smoothing 0.081, t0 0.5 s), locked fault borders O5
  ([fault_tip_lock], the free-surface traces stay free), Kelvin-Voigt eta_bar 0.1 (scale node, inscribed length)
  and the Prakash-Clifton regularization (v_star 1 m/s, L 0.01 m); no plasticity.
- Mesh: tpv22/mesh/tpv22_200m.msh, made on Vista by build_verify_complex_vista.sbatch from the committed .geo
  (gmsh 4.15.0, msh 2.2): 3,628,170 tets, 69,922 fault faces (610 shared between ranks), sha256
  4348fa1993fc11c4b0cf8b3a6b5703d99af376d0a88d741f119ea0e5d6c4a51c (not in the repo).
- Run: Vista gg, 32 nodes / 2304 ranks (72 per node, 113.4 k tets per node, every rank reads the serial mesh on
  this branch), account EAR24034, 2026-10-24 4:11 to 6:05 pm PDT (1:54:21 job wall, 20.1 SU), 0-15 s, 199,924
  steps, dt 7.502857152e-05 s, central difference ("cd"), 6842.5 s driver wall (322.2 s of it fault output),
  max|V| 8.63 m/s, no NaN. 3,409 of 699,220 fault DOFs locked on the borders; the forced rupture covers 16,475
  fault DOFs. sigma_n_eff stays in [9.72, 106.7] MPa, so the fault never went into tension and the tension floor
  never acted (0 DOF-steps).
- results/: the 18 SCEC station files (7 on fault 1, 11 on fault 2; tpv22_station_fault<k>_x2_<a>_x3_<b>.dat with
  x2 = global X km and x3 = depth km; t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip, 2 = strike, so
  h-slip is slip2, positive right-lateral on both faults; SI; sigma_n compression positive;
  tpv22_station_conventions.txt), the station snap list (nearest fault node, up to 35.6 m from the requested
  point), driver.log, RUN_INFO.txt, the job's stdout tpv22_200m_p3_1098317.out (and .err) and meminfo_kB.log
  (first node, every 2 s: epoch s, MemTotal, MemAvailable, kB). stations.sha256 holds the sha256 of each station
  file as written on Vista (sha256sum in the run folder on Vista; the files came off Vista through the Mac and
  match it). Read with `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  Left out: the per-station DG copy envelopes (*_copies.dat, 18 files) and fault.vtkhdf (748 frames); both are
  on Vista scratch (13 GB output folder), the envelopes also on Chunhui Zhao's Mac.
- plots/: tpv22_200m_p3_vista1098317_stations.pdf (21 pages) and its page PNGs: misfit maps of arrival, final
  slip and peak strike-slip rate at the 18 stations for this run, the p2 run (Vista 1098316, 200 m) and the p1
  run (Vista 1086895, 200 m), then an overview per fault and per station (a) strike slip rate, (b) front zoom,
  (c) strike slip, (d) strike shear stress, (e) normal stress and (f) dip slip rate; this run in red, drawn last,
  over p2 in green and p1 in orange dashed, vs DRDG3D in black. plot_complex_stations_p123.py made them (v3, the
  same file as in the TPV12 p3 gold; the benchmarks thread's TPV31 p1/p2/p3 station pages with input changes,
  listed in its docstring). The script ran with --tag tpv22_200m_p123_vista; the PDF and PNGs were renamed to
  this run's tag, their content unchanged (the tag sets only the file names). compare.txt / compare_summary.json
  from compare_nitsche_fault_stations.py at 6e35a58 vs DRDG3D with p1 and p2 as the runs before this one
  (command below). DRDG3D ends at 14.9996 s, so every source is measured to there. The p1 and p2 runs are not
  golds and their station files are not in the repo.
- vs DRDG3D (median |misfit| of 18 stations): arrival 0.0342 s, final slip 0.23 %, peak strike-slip rate 2.82 %
  (p2 0.0460 s, 0.47 %, 2.90 %; p1 0.1993 s, 4.93 %, 13.06 %).
  0. signs of final slip and the t = 0 stresses agree everywhere: PASS;
  1. arrival within 0.1 s everywhere: PASS (worst +0.098 s at fault 2, x2 = 5 km, x3 = 5 km);
  2. slip within 2 % where DRDG3D slips more than 0.1 m: PASS (worst -0.63 % at fault 2, x2 = 5 km, x3 = 15 km);
  3. peak within 10 % at 14 of 18 stations. The four misses are all above DRDG3D, and three of them grow from p2
     to p3: fault 1 (0, 0) +17.6 % (p2 +10.7 %, p1 +13.5 %), a narrow spike at the front; fault 1 (-5, 10)
     +14.7 % (p2 +8.7 %, p1 -10.0 %); fault 2 (0, 0) +10.2 % (p2 +7.9 %, p1 +24.6 %); fault 2 (20, 10) +14.0 %
     (p2 +8.4 %, p1 -12.6 %). Positions are (x2, x3) in km;
  4. every median below p2's and p1's: PASS.
  No second TPV22 reference is in hand to give DRDG3D's own spread, so whether those peaks are this run's
  overshoot or DRDG3D's under-resolution is not known.
- Fault-wide (driver.log): max|V| 8.63 m/s (p2 7.63) and the lowest sigma_n_eff 9.72 MPa (p2 21.2 MPa) both move
  with the order.
- driver.log prints a "Y=0 Face Classification Audit" warning (758 y = 0 faces with no BC). The audit is written
  for a fault covering the plane y = 0; fault 1 covers only X in [-25, 5] km of that plane, so the other y = 0
  faces are ordinary interior faces.

## Speed

From results/tpv22_200m_p3_1098317.out, driver.log and meminfo_kB.log; job wall from sacct; Vista gg charges
0.33 SU per node-hour. The p2 column is the p2 run, same build, deck, mesh and layout (from its .out and
driver.log, not here).

| | p3, this run (4eaf4df) | p2 run 1098316 (4eaf4df) |
|---|---|---|
| nodes x ranks | 32 x 72 | 32 x 72 |
| steps to 15 s, dt | 199,924, 7.503e-5 s | 118,194, 1.269e-4 s |
| job wall (sacct) | 1:54:21 | 31:52 |
| driver wall | 6842.5 s (34.23 ms per step) | 1893.7 s (16.02 ms per step) |
| time in ode.Step per rank, slowest | 6033.1 s = 30.18 ms per step | 1428.9 s = 12.09 ms per step |
| time in ode.Step per rank, min / mean | 24.61 / 28.27 ms per step | 8.70 / 11.48 ms per step |
| step 1 (with operator assembly) | 10.15 s | 2.72 s |
| setup plus the first 1000 steps | 190.1 s (driver.log wall at step 1000) | 155.0 s |
| fault work, busiest rank | 4.85 ms per step | 1.43 ms per step |
| fault output (rank 0) | 748 frames, 322.2 s (0.431 s per frame) | 748 frames, 323.9 s (0.433 s per frame) |
| dense blocks, all ranks | 556.7 GB | 139.9 GB |
| node-hours, SU | 61.0, 20.1 | 17.0, 5.6 |
| first-node memory peak in use | 124.5 of 237.6 GiB (148 s after start, in setup) | 124.6 GiB |
| output folder on Vista | 13 GB | 13 GB |

- p3 takes 1.69 x the p2 steps at 2.50 x the p2 cost per step, so 4.2 x the p2 stepping time; setup and the
  fault output cost about the same at both orders, so the job wall is 3.6 x p2's. The job script's estimate was
  35.4 SU plus the fault output.
- The memory peak comes in setup, before step 1000, and is the same at p2 and p3. While stepping, the first node
  holds about 47 to 54 GiB at p3 and 33 GiB at p2. So the peak does not come from the order; that it comes from
  every rank reading the serial mesh is inferred, not measured.
- This branch's run script logs MemTotal and MemAvailable only, not Committed_AS, so the margin to Vista's
  118.8 GiB per-node commit limit is not known; the run did not reach it.
- The driver's "time in ode.Step" excludes setup; the slowest rank sets the pace.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-dynamic-rupture-complex-geometries \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-complex-geometries
    cd seas-mfem-complex-geometries && git checkout 4eaf4df769dfb4eb8279deb02f4efa75dc1fe804
    cd miniapps/seas/jobs/vista_complex && sbatch build_verify_complex_vista.sbatch   # build, meshes, checks
    sbatch tpv22_200m_p3.sbatch

build_verify_complex_vista.sbatch generates tpv22/mesh/tpv22_200m.msh from the committed .geo with Vista's gmsh
4.15.0; gmsh output depends on the platform, so compare its sha256 with the one above (the run script refuses a
mesh whose sha256 is not in MESH_SHA256_COMPLEX_VISTA). The job refuses to start unless the checkout is clean,
HEAD is on GitHub and the driver was built and verified at HEAD.

The DRDG3D reference is not in the repo (tpv22/benchmark_data/README.md). The comparison and the plots were made
on Chunhui Zhao's Mac on 2026-10-25 (Python with numpy and matplotlib) from miniapps/seas at 6e35a58, with the
commands below; there, D = this folder, REF = the DRDG3D station folder, and P1 / P2 = the p1 / p2 station
folders stand in for the Mac paths that compare.txt prints. Its "# wrote" line names the JSON file before it was
renamed to compare_summary.json; the plot tag was tpv22_200m_p123_vista before the files were renamed.

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv22/stations/tpv22_fault_stations.csv \
        --prefix tpv22 --run p1=$P1 --run p2=$P2 --run p3=$D/results --ref DRDG3D=$REF --by-name --signed \
        --t-end 15 --json $D/plots/compare_summary.json > $D/plots/compare.txt
    python3 $D/plots/plot_complex_stations_p123.py --seas . --prob tpv22 --comp h --t-end 15 \
        --tag tpv22_200m_p3_vista1098317 --out $D/plots --ref "DRDG3D|$REF" \
        --run "p1|Nitsche p=1 200 m (Vista 1086895)|$P1" --run "p2|Nitsche p=2 200 m (Vista 1098316)|$P2" \
        --run "p3|Nitsche p=3 200 m (Vista 1098317)|$D/results" --cmp $D/plots/compare.txt

- Requested by Chunhui Zhao (complex fault geometries, 2026-10-25): "tpv22 p3 results look good to me, gold the
  results, commit and push".
