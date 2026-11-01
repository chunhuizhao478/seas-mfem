# TPV26 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p2, 200 m (Vista job 1099061)

- Code: seas-mfem branch nitsche-dynamic-rupture-nonlinearity, commit 3ab9e59b730f69cfd6e14dc69206c42d301fe8d5
  (3ab9e59, on GitHub), built on TACC Vista by jobs/vista_rheology/build_tpv2627_vista.sbatch (job 1099060,
  2026-10-24 22:22 to 22:51 PDT, clean checkout, BUILD_COMMIT_VISTA = 3ab9e59, "all checks passed"; its
  informational 200 m load check on one node was killed at 230 GiB committed, as in build 1086883, and does not
  gate). Driver seas_nitsche_dyn_driver; the run logs "fast fault on, fused K/K_eta off, dense blocks on".
- Deck and job scripts: jobs/ (= tpv26/configs/tpv26_nitsche_kv0p1.toml, jobs/vista_rheology/{tpv26_200m_p2.sbatch,
  build_tpv2627_vista.sbatch} and jobs/vista_nitsche/run_nitsche_vista.sh at 3ab9e59). TPV26 elastic: forced
  rupture from (-5, -10) km, depth-dependent effective stress with fluid pressure, linear slip weakening on slip
  path, Prakash-Clifton regularization (v* = 1 m/s, L = 0.01 m), Kelvin-Voigt eta_bar = 0.1 over the whole domain,
  rheology table present with no model enabled.
- Mesh: tpv26/mesh/tpv2627_mirror_200m.msh, made on Vista by build 1099060 from the committed
  tpv26/mesh/tpv2627_half.geo with tpv26/mesh/make_tpv2627_meshes.sh (gmsh 4.15.0, built-in kernel, msh 2.2; the
  y >= 0 half mirrored across the fault): 2,997,834 tets, 46,338 fault faces, sha256
  ed8ce90b30a55999bd34a5ae8c7fa90da3063e8d1377f270e2c988559de4f599 (not in the repo). The 12 on-fault stations are
  mesh vertices (snap distance 0 m at all 12, results/tpv26_station_snap.txt). The earlier run 1094998 (b618bbb5,
  stations not vertices) snapped up to 57 m away and read sigma_n 0.3 to 1 % high at five stations; this gold
  replaces it.
- Run: Vista gg, 24 nodes / 3072 ranks (128 per node), account EAR26024, 2026-10-24 22:52 to 23:18 PDT (about
  26 min job wall, about 3.4 SU), 0 to 13 s, 120,228 steps, dt 1.081286e-4 s, central difference ("cd"),
  1527.5 s driver wall, max|V| 6.577 m/s, no NaN.
- results/: the 12 SCEC on-fault station files (t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip,
  2 = strike; SI; sigma_n compression positive; tpv26_station_conventions.txt), the 6 off-fault station files
  (tpv26_offfault_body*.dat), the station snap list, driver.log and RUN_INFO.txt. files.sha256 is the sha256 of each
  file as written on Vista; the files came off Vista through the Mac and were checked against it there and again
  in the cloud. Not here: the per-station _copies.dat files (12 x 23 MB, the DG copies the station rows average),
  the job's stdout and memory log, and fault.vtkhdf (649 frames); they stay in
  $SCRATCH/nitsche_benchmarks/tpv26_mirror200m_p2_1099061 on Vista (copies files also on the Mac).
- plots/: tpv26_200m_p2_vista1099061_stations.pdf and its page PNGs (misfit maps of arrival, final slip and peak
  slip rate over the 12 stations, then per station slip rate, front zoom, slip, shear stress, normal stress and
  shear/normal ratio; this run in red, drawn last, vs DRDG3D 200 m O4 in black); compare.txt and
  compare_summary.json from scripts/compare_nitsche_fault_stations.py at 6e35a58; t0_normal_stress.txt, the first
  row of each station file against the first row of each DRDG3D file (DRDG3D n-stress sign flipped to compression
  positive). DRDG3D ends at 11.9997 s, so every measure is to t = 11.9997 s. compare.txt's last line is the JSON
  path on the Mac where it was written.
- vs DRDG3D (median |misfit| of 12 stations): arrival 0.034 s, final slip 0.35 %, peak h-slip-rate 4.19 %.
  1. arrival within 0.1 s at 11 of 12; the miss is x2 = 10 km, 10 km depth at -0.313 s by the 1 mm/s rule, where
     a small rise crosses 1 mm/s before the front; the main rise on its page is about 0.05 s ahead of DRDG3D.
  2. slip within 2 % everywhere (worst -0.95 % at x2 = 15 km, 15 km depth). PASS.
  3. peak within 10 % at 7 of 12; 12 to 23 % low at x2 = -15, 5, 10, 15 km (10 km depth) and x2 = 15 km
     (15 km depth), the sharp-front stations, as on the old mesh.
- Stress at t = 0 (t0_normal_stress.txt): sigma_n and shear equal DRDG3D to 1e-4 MPa at 0, 5 and 10 km depth
  (87.8204 / 13.8316 MPa at 5 km, 175.6409 / 27.6633 MPa at 10 km). At 15 km, where the taper starts, sigma_n is
  263.468 and 263.461 MPa against DRDG3D 263.389 (+0.03 %) and shear 41.51 and 41.49 against 41.33 MPa (+0.4 %);
  these sit on the deck formula, and tpv26/benchmark_data notes DRDG3D is about 0.4 % below the formula in shear
  there.

## Commands

Build and run (Vista login node, clean checkout of 3ab9e59 in $WORK/seas-mfem-nonlinearity):

```
cd $WORK/seas-mfem-nonlinearity/miniapps/seas/jobs/vista_rheology
B=$(sbatch --parsable build_tpv2627_vista.sbatch | tail -n 1)
sbatch --dependency=afterok:$B tpv26_200m_p2.sbatch
```

Comparison (Mac, pythonenv, matplotlib 3.10.6), from a checkout of miniapps/seas at 6e35a58 with this folder's
results/ as RUN; plot_complex_stations.py (sha256 44893758de364a97f59139584a09dacb4edd350f1a208550f305f9252ac34732)
is the complex-geometries thread's station page script, not in the repo:

```
python3 scripts/compare_nitsche_fault_stations.py --stations tpv26/stations/tpv26_fault_stations.csv \
    --prefix tpv26 --run p2_200m=RUN --ref DRDG3D=tpv26/benchmark_data/scec_drdg3d --t-end 11.9997 \
    --json OUT/compare_summary.json > OUT/compare.txt
python3 plot_complex_stations.py . RUN OUT "Nitsche p=2 200 m (Vista 1099061)" 11.9997 tpv26 \
    tpv26/benchmark_data/scec_drdg3d DRDG3D h tpv26_200m_p2_vista1099061 OUT/compare.txt
```

## Speed

From results/driver.log; Vista gg charges 0.33 SU per node-hour.

| | this run (3ab9e59) |
|---|---|
| nodes x ranks | 24 x 128 |
| steps to 13 s, dt | 120,228, 1.081e-4 s |
| driver wall | 1527.5 s (12.70 ms per step) |
| time in ode.Step per rank, slowest | 1475.9 s = 12.28 ms per step |
| time in ode.Step per rank, min / mean | 10.49 / 12.27 ms per step |
| step 1 (with operator assembly) | 1.73 s |
| setup (slowest rank) | 58.3 s (rank 0 mesh read 23.4 s, partition 9.6 s, ParMesh 21.7 s) |
| fault work, busiest rank | 1.60 ms per step |
| fault output (rank 0) | 649 frames, 215.3 s (0.332 s per frame) |
| dense blocks, all ranks | 116.8 GB |
| first-node memory peak committed | 40.1 of 118.8 GiB (the job's memory log on Vista) |

The old-mesh run 1094998 (same layout and deck, 2,995,976 tets) took 97,755 steps at 12.8 ms per step: Vista's
gmsh made a smaller smallest element with the station vertices, so dt fell 19 % and the step count rose by the same.
