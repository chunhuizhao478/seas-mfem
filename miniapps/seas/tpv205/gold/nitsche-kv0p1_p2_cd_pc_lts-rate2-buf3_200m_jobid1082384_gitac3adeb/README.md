# TPV205 gold: Nitsche dynamic rupture with local time stepping, p2, 200 m (Vista job 1082384)

- Code: seas-mfem branch claude/nitsche-profiling-w336yq, commit ac3adeb9b9e5781b9110139eac2ba2e39bbfdd74 (ac3adeb, on
  GitHub; the code change is da56f4b, "local time stepping, fault buffer layers against the 2-step ripple"), built on
  TACC Vista by jobs/vista_profiling/build_lts_vista.sbatch (job 1082383: clean checkout, BUILD_COMMIT_VISTA =
  ac3adeb; the six gate test sets and the LTS suite, 24 of 24, passed; both 200 m gold meshes matched meshes.sha256).
  Driver seas_nitsche_dyn_driver.
- Time stepping: `--lts rate2` (multirate leapfrog, rate 2, at most 4 clusters) with the default `--lts-fault-buffer
  3`: the 52,456 fault tets and the 106,114 tets of 3 face-neighbour layers around them are in cluster 0. Tets per
  cluster 169,867 / 313,803 / 712,947 / 697,966; ideal speedup 3.20x. dt0 = 1.352508279e-04 s, equal to global
  stepping's dt; 88,728 ticks to 12.0005 s. Everything else as the global-stepping gold this one sits beside,
  nitsche-kv0p1_p2_cd_pc_200m_jobid1075916_git81759da: same deck, mesh, order and rank layout.
- Deck and job script: jobs/ (= tpv205/configs/tpv205_nitsche.toml, jobs/vista_profiling/lts_tpv205_200m_p2.sbatch and
  run_lts_vista.sh at ac3adeb; the deck is byte-identical to the copy the job saved) and the job's stdout and stderr,
  lts_tpv205_200m_p2_1082384.out / .err.
- Mesh: tpv205/mesh/tpv2053d_200m.msh, the global-stepping gold's mesh (1,894,583 tets, 26,228 fault faces), sha256
  30c63b923781177830330a36aebcc4159cfc8043374ea4dcd883a9b2eec13538 (not in the repo; build_lts_vista.sbatch copies it
  from $WORK/seas-mfem-benchmarks and checks it against jobs/vista_profiling/meshes.sha256).
- Run: Vista gg, 16 nodes / 2304 ranks (144 per node), account EAR26024, 2026-10-20 4:29-4:45 am PDT (16:20 job wall,
  about 1.44 SU). The job first runs global stepping (the defaults) 0-1.0 s at the same commit as the speed and
  round-off reference, then LTS 0-12 s. Run with node i617-131 excluded (`scontrol update ExcNodeList=i617-131` while
  it waited, Chunhui Zhao 2026-10-19), after that node failed the first p1 jobs out of memory.
- Speed (busiest rank, time in ode.Step after step 1): LTS 5.24 ms per step against global stepping 9.42 ms per step
  at the same dt, 1.80x. Driver wall: LTS 739 s for 0-12 s; global stepping 199 s for 0-1 s (7,394 steps). The patch
  kicks cost 0.48 of global stepping per tick and the halo wait takes up to 0.57 of a rank's stepping, so the speedup
  stays below the ideal 3.20x. First-node memory peak 116.3 GB (LTS), 138.4 GB (global stepping) of 237.6 GB. max|V|
  10.2045 m/s, no NaN.
- results/: the 16 SCEC station files of the LTS run (t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip, 2 =
  strike; SI; sigma_n compression positive; tpv205_station_conventions.txt), the station snap list, the LTS
  driver.log, gts_driver.log (the 1 s global-stepping run) and RUN_INFO.txt. stations.sha256 holds the sha256 of each
  station file as written on Vista (checked after the transfer). fault.vtkhdf (Vista scratch) is not here. Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
- Against the global-stepping gold 1075916, 0-12 s (plots/lts_vs_gts_gold.txt, from
  jobs/vista_profiling/judge_lts_gold.py; criteria of the LTS buffer-fix review, with the 10 Pa floor of criterion 4
  that Chunhui Zhao accepted on 2026-10-20): ALL PASS. 1 SCEC pass and fail the same as the gold's; 2 arrival within
  0.41 ms, peak slip rate within 0.051 %, final slip within 0.020 % at every station; 3 noise RMS within 0.2 % of the
  gold's; 4 the period-2 part of the worst station within 2x the gold's or below 10 Pa in every 1 s window: tau_h
  above 2x in 1 window (7-8 s), at most 6.0 Pa against the gold's 2.5 Pa; sigma_n above 2x in 5 windows (7-12 s), at
  most 4.3 Pa against the gold's 0.6-1.4 Pa. All are late windows (from 7 s on), where the gold's own level is a few
  Pa; inferred to be what the 3 buffer layers leave of the seam ripple.
- Whole traces (plots/compare_gold.txt, max |LTS - gold| / max |gold| per station, groups as compare_stations.py
  --groups, 88,725 rows): worst slip 6.488e-04, V 1.104e-02, tau 3.460e-03, sigma_n 3.024e-04. The LTS run before the
  buffer fix (job 1081277, 2d09c44) gave 3.168e-03, 1.689e-02, 1.757e-02, 1.599e-03. The job's own checks: global
  stepping against the gold over 0-1 s 2.691e-11 (round-off); LTS against global stepping over 0-1 s: slip 1.404e-04,
  V 1.613e-04, tau 1.671e-05, sigma_n 5.296e-06.
- vs DRDG3D 100 m (plots/compare.txt): median |misfit| arrival 0.0232 s, final slip 0.49 %, peak h-slip-rate 1.78 %;
  pass tests 1 PASS, 2 PASS, 3 peak within 10% at 14 of 16 stations PASS (pass and fail the same as the
  global-stepping gold's, criterion 1 above).
- plots/: tpv205_200m_p2_lts_vista1082384_stations.pdf and its page PNGs (slip-rate overview of the 16 stations, then
  per station (a)-(f) slip rate, front zoom, slip, shear stress, normal stress and dip slip rate against the
  global-stepping gold and DRDG3D 100 m, (g)-(j) LTS - gold of slip rate, slip, shear stress and normal stress,
  (k)-(l) the period-2 part per 0.25 s; jobs/vista_profiling/plot_lts_gold.py); compare.txt / compare_summary.json
  from compare_nitsche_fault_stations.py against the references, both mixed-flux golds and the global-stepping gold;
  compare_gold.txt; lts_vs_gts_gold.txt.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b claude/nitsche-profiling-w336yq \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-profiling
    cd seas-mfem-profiling && git checkout ac3adeb9b9e5781b9110139eac2ba2e39bbfdd74
    cd miniapps/seas/jobs/vista_profiling && sbatch build_lts_vista.sbatch   # build, gold meshes, gate + LTS tests
    sbatch --dependency=afterok:<build job> lts_tpv205_200m_p2.sbatch

The jobs refuse to start unless the checkout is clean, HEAD is on GitHub, the driver was built at HEAD and the mesh
sha256 is in meshes.sha256. Checks and plots (from miniapps/seas, at the commit that adds this folder, with
D = this folder and G = tpv205/gold/nitsche-kv0p1_p2_cd_pc_200m_jobid1075916_git81759da):

    python3 -I jobs/vista_profiling/compare_stations.py $G/results $D/results --tmax 12.0 --time-tol 2e-9 --groups
    python3 -I jobs/vista_profiling/judge_lts_gold.py tpv205 $G/results $D/results /tmp/judge > $D/plots/lts_vs_gts_gold.txt
    python3 -I jobs/vista_profiling/plot_lts_gold.py tpv205 $D/results $G/results tpv205/benchmark_data/DRDG3D_100m_O4 \
        $D/plots/tpv205_200m_p2_lts_vista1082384 'LTS p=2 200 m (Vista 1082384)' 'gold (Vista 1075916)' 'DRDG3D 100 m O4'
    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv205/stations/tpv205_fault_stations.csv \
        --prefix tpv205 --run lts_p2_200m=$D/results \
        --ref drdg3d_100m=tpv205/benchmark_data/DRDG3D_100m_O4 --ref drdg3d_200m=tpv205/benchmark_data/DRDG3D_200m_O4 \
        --gold mixedflux_p1=tpv205/gold/mixedflux_p1_rk45_jobidcombined_7762879_7763019_git04674fc/results \
        --gold mixedflux_p2=tpv205/gold/mixedflux_p2_rk45_jobid7762880_git04674fc/results \
        --gold gts_p2_200m=$G/results:nitsche --peak-min-stations 14 \
        --json $D/plots/compare_summary.json > $D/plots/compare.txt

Rerunning these reproduces lts_vs_gts_gold.txt and compare.txt (apart from compare.txt's last line, the path it
wrote the JSON to) and the plots. compare_gold.txt holds the per-station numbers of the first command for this run
and for round 1 (job 1081277, whose station files stayed off the repo).

- Requested by Chunhui Zhao (Nitsche profiling, plan item 7: local time stepping golds, 2026-10-20).
