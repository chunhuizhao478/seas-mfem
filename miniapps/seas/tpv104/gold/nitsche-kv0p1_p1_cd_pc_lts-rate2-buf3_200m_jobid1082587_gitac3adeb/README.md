# TPV104 gold: Nitsche dynamic rupture with local time stepping, p1, 200 m (Vista job 1082587)

- Code: seas-mfem branch claude/nitsche-profiling-w336yq, commit ac3adeb9b9e5781b9110139eac2ba2e39bbfdd74 (ac3adeb, on
  GitHub; the code change is da56f4b, "local time stepping, fault buffer layers against the 2-step ripple"), built on
  TACC Vista by jobs/vista_profiling/build_lts_vista.sbatch (job 1082383: clean checkout, BUILD_COMMIT_VISTA =
  ac3adeb; the six gate test sets and the LTS suite, 24 of 24, passed; both 200 m gold meshes matched meshes.sha256).
  Driver seas_nitsche_dyn_driver.
- Time stepping: `--lts rate2` (multirate leapfrog, rate 2, at most 4 clusters) with the default `--lts-fault-buffer
  3`: the 75,856 fault tets and the 157,010 tets of 3 face-neighbour layers around them are in cluster 0. Tets per
  cluster 239,340 / 219,781 / 1,031,158 / 976,270; ideal speedup 3.38x. dt0 = 2.450956970e-04 s, equal to global
  stepping's dt; 48,968 ticks to 12.0018 s. Everything else as the global-stepping gold this one sits beside,
  nitsche-kv0p1_p1_cd_pc_200m_jobid1075919_git81759da: same deck, mesh, order and rank layout.
- Deck and job script: jobs/ (= tpv104/configs/tpv104_nitsche.toml, jobs/vista_profiling/lts_tpv104_200m_p1.sbatch and
  run_lts_vista.sh at ac3adeb; the deck is byte-identical to the copy the job saved) and the job's stdout and stderr,
  lts_tpv104_200m_p1_1082587.out / .err.
- Mesh: tpv104/mesh/tpv104_200m.msh, the global-stepping gold's mesh (2,466,549 tets, 37,928 fault faces), sha256
  10cecfd73d0c7cb93cfa188aac08e4a99a6ec81d439a3363bd7175e3a8ae4aa0 (not in the repo; build_lts_vista.sbatch copies it
  from $WORK/seas-mfem-benchmarks and checks it against jobs/vista_profiling/meshes.sha256).
- Run: Vista gg, 4 nodes / 512 ranks (128 per node), account EAR26024, 2026-10-19 4:39-4:51 pm PDT (11:25 job wall,
  about 0.25 SU). The job first runs global stepping (the defaults) 0-1.0 s at the same commit as the speed and
  round-off reference, then LTS 0-12 s. Resubmission of 1082387 with node i617-131 excluded (`sbatch
  --exclude=i617-131`, Chunhui Zhao 2026-10-19). 1082387 ran on i617-[031,053,131-132]; both its runs stopped on
  i617-131 (METIS "Memory allocation failed" in the global-stepping run, std::bad_alloc right after the LTS setup),
  with no station output; its stdout is in jobs/vista_profiling/results/.
- Speed (busiest rank, time in ode.Step after step 1): LTS 8.14 ms per step against global stepping 11.56 ms per step
  at the same dt, 1.42x. Driver wall: LTS 518 s for 0-12 s; global stepping 143 s for 0-1 s (4,081 steps). The patch
  kicks cost 0.45 of global stepping per tick and the halo wait takes up to 0.69 of a rank's stepping, so the speedup
  stays below the ideal 3.38x. First-node memory peak 126.3 GB (LTS), 145.3 GB (global stepping) of 237.6 GB. max|V|
  17.1446 m/s, no NaN.
- results/: the 9 SCEC station files of the LTS run (t slip1 slip2 V1 V2 tau1 tau2 sigma_n psi; 1 = dip, 2 = strike;
  SI; sigma_n compression positive; tpv104_station_conventions.txt), the station snap list, the LTS driver.log,
  gts_driver.log (the 1 s global-stepping run) and RUN_INFO.txt. stations.sha256 holds the sha256 of each station file
  as written on Vista (checked after the transfer). fault.vtkhdf (Vista scratch) is not here. Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
- Against the global-stepping gold 1075919, 0-12 s (plots/lts_vs_gts_gold.txt, from
  jobs/vista_profiling/judge_lts_gold.py; criteria of the LTS buffer-fix review, with the 10 Pa floor of criterion 4
  that Chunhui Zhao accepted on 2026-10-20): ALL PASS. 1 SCEC pass and fail the same as the gold's; 2 arrival within
  0.98 ms, peak slip rate within 0.114 %, final slip within 0.092 % at every station; 3 noise RMS within 3.3 % of the
  gold's; 4 the period-2 part of the worst station within 2x the gold's or below 10 Pa in every 1 s window: tau_h
  above 2x in 2 windows (9-10 s, 11-12 s), at most 2.8 Pa against the gold's 0.9-1.1 Pa; sigma_n above 2x in 3 windows
  (9-12 s), at most 3.4 Pa against the gold's 0.8-1.0 Pa. All are late windows (from 7 s on), where the gold's own
  level is a few Pa; inferred to be what the 3 buffer layers leave of the seam ripple.
- Whole traces (plots/compare_gold.txt, max |LTS - gold| / max |gold| per station, groups as compare_stations.py
  --groups, 48,961 rows): worst slip 2.061e-03, V 1.749e-02, tau 1.164e-02, sigma_n 1.976e-03. The LTS run before the
  buffer fix (job 1081280, 2d09c44) gave 1.840e-02, 1.888e-01, 8.595e-02, 1.371e-02. The job's own checks: global
  stepping against the gold over 0-1 s 1.000e-10 (round-off); LTS against global stepping over 0-1 s: slip 1.444e-04,
  V 1.886e-04, tau 4.056e-05, sigma_n 6.593e-06.
- vs DRDG3D (plots/compare.txt): median |misfit| arrival 0.0775 s, final slip 0.56 %, peak h-slip-rate 20.15 %; pass
  tests 1 FAIL, 2 FAIL, 3 peak within 10% at 1 of 9 stations (pass and fail the same as the global-stepping gold's,
  criterion 1 above).
- plots/: tpv104_200m_p1_lts_vista1082587_stations.pdf and its page PNGs (slip-rate overview of the 9 stations, then
  per station (a)-(f) slip rate, front zoom, slip, shear stress, normal stress and dip slip rate against the
  global-stepping gold and DRDG3D, (g)-(j) LTS - gold of slip rate, slip, shear stress and normal stress, (k)-(l) the
  period-2 part per 0.25 s; jobs/vista_profiling/plot_lts_gold.py); compare.txt / compare_summary.json from
  compare_nitsche_fault_stations.py against the references, both mixed-flux golds and the global-stepping gold;
  compare_gold.txt; lts_vs_gts_gold.txt.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b claude/nitsche-profiling-w336yq \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-profiling
    cd seas-mfem-profiling && git checkout ac3adeb9b9e5781b9110139eac2ba2e39bbfdd74
    cd miniapps/seas/jobs/vista_profiling && sbatch build_lts_vista.sbatch   # build, gold meshes, gate + LTS tests
    sbatch --dependency=afterok:<build job> lts_tpv104_200m_p1.sbatch

The jobs refuse to start unless the checkout is clean, HEAD is on GitHub, the driver was built at HEAD and the mesh
sha256 is in meshes.sha256. Checks and plots (from miniapps/seas, at the commit that adds this folder, with
D = this folder and G = tpv104/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid1075919_git81759da):

    python3 -I jobs/vista_profiling/compare_stations.py $G/results $D/results --tmax 12.0 --time-tol 2e-9 --groups
    python3 -I jobs/vista_profiling/judge_lts_gold.py tpv104 $G/results $D/results /tmp/judge > $D/plots/lts_vs_gts_gold.txt
    python3 -I jobs/vista_profiling/plot_lts_gold.py tpv104 $D/results $G/results tpv104/benchmark_data/DRDG3D \
        $D/plots/tpv104_200m_p1_lts_vista1082587 'LTS p=1 200 m (Vista 1082587)' 'gold (Vista 1075919)' 'DRDG3D'
    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv104/stations/tpv104_fault_stations.csv \
        --prefix tpv104 --run lts_p1_200m=$D/results \
        --ref drdg3d=tpv104/benchmark_data/DRDG3D --ref seissol=tpv104/benchmark_data/seisol \
        --gold mixedflux_p1=tpv104/gold/mixedflux_p1_rk45_jobid7765845_git04674fc/results \
        --gold mixedflux_p2=tpv104/gold/mixedflux_p2_rk45_jobid7765380_git04674fc/results \
        --gold gts_p1_200m=$G/results:nitsche \
        --json $D/plots/compare_summary.json > $D/plots/compare.txt

Rerunning these reproduces lts_vs_gts_gold.txt and compare.txt (apart from compare.txt's last line, the path it
wrote the JSON to) and the plots. compare_gold.txt holds the per-station numbers of the first command for this run
and for round 1 (job 1081280, whose station files stayed off the repo).

- Requested by Chunhui Zhao (Nitsche profiling, plan item 7: local time stepping golds, 2026-10-20).
