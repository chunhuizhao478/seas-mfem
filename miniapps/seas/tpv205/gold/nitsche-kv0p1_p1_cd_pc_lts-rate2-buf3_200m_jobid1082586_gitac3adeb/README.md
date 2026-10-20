# TPV205 gold: Nitsche dynamic rupture with local time stepping, p1, 200 m (Vista job 1082586)

- Code: seas-mfem branch claude/nitsche-profiling-w336yq, commit ac3adeb9b9e5781b9110139eac2ba2e39bbfdd74 (ac3adeb, on
  GitHub; the code change is da56f4b, "local time stepping, fault buffer layers against the 2-step ripple"), built on
  TACC Vista by jobs/vista_profiling/build_lts_vista.sbatch (job 1082383: clean checkout, BUILD_COMMIT_VISTA =
  ac3adeb; the six gate test sets and the LTS suite, 24 of 24, passed; both 200 m gold meshes matched meshes.sha256).
  Driver seas_nitsche_dyn_driver.
- Time stepping: `--lts rate2` (multirate leapfrog, rate 2, at most 4 clusters) with the default `--lts-fault-buffer
  3`: the 52,456 fault tets and the 106,114 tets of 3 face-neighbour layers around them are in cluster 0. Tets per
  cluster 169,867 / 313,803 / 712,947 / 697,966; ideal speedup 3.20x. dt0 = 2.883649727e-04 s, equal to global
  stepping's dt; 41,616 ticks to 12.0006 s. Everything else as the global-stepping gold this one sits beside,
  nitsche-kv0p1_p1_cd_pc_200m_jobid1076378_git81759da: same deck, mesh, order and rank layout.
- Deck and job script: jobs/ (= tpv205/configs/tpv205_nitsche.toml, jobs/vista_profiling/lts_tpv205_200m_p1.sbatch and
  run_lts_vista.sh at ac3adeb; the deck is byte-identical to the copy the job saved) and the job's stdout and stderr,
  lts_tpv205_200m_p1_1082586.out / .err.
- Mesh: tpv205/mesh/tpv2053d_200m.msh, the global-stepping gold's mesh (1,894,583 tets, 26,228 fault faces), sha256
  30c63b923781177830330a36aebcc4159cfc8043374ea4dcd883a9b2eec13538 (not in the repo; build_lts_vista.sbatch copies it
  from $WORK/seas-mfem-benchmarks and checks it against jobs/vista_profiling/meshes.sha256).
- Run: Vista gg, 4 nodes / 576 ranks (144 per node), account EAR26024, 2026-10-19 4:31-4:39 pm PDT (8:07 job wall,
  about 0.18 SU). The job first runs global stepping (the defaults) 0-1.0 s at the same commit as the speed and
  round-off reference, then LTS 0-12 s. Resubmission of 1082386 with node i617-131 excluded (`sbatch
  --exclude=i617-131`, Chunhui Zhao 2026-10-19). 1082386 ran on i617-[031,053,131-132]; its 1 s global-stepping run
  stopped at startup on i617-131 (METIS "Memory allocation failed"), while its LTS run finished with station files
  bitwise equal to this run's; its stdout is in jobs/vista_profiling/results/.
- Speed (busiest rank, time in ode.Step after step 1): LTS 6.07 ms per step against global stepping 7.87 ms per step
  at the same dt, 1.30x. Driver wall: LTS 359 s for 0-12 s; global stepping 100 s for 0-1 s (3,468 steps). The patch
  kicks cost 0.48 of global stepping per tick and the halo wait takes up to 0.61 of a rank's stepping, so the speedup
  stays below the ideal 3.20x. First-node memory peak 112.7 GB (LTS), 136.3 GB (global stepping) of 237.6 GB. max|V|
  8.12351 m/s, no NaN.
- results/: the 16 SCEC station files of the LTS run (t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip, 2 =
  strike; SI; sigma_n compression positive; tpv205_station_conventions.txt), the station snap list, the LTS
  driver.log, gts_driver.log (the 1 s global-stepping run) and RUN_INFO.txt. stations.sha256 holds the sha256 of each
  station file as written on Vista (checked after the transfer). fault.vtkhdf (Vista scratch) is not here. Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
- Against the global-stepping gold 1076378, 0-12 s (plots/lts_vs_gts_gold.txt, from
  jobs/vista_profiling/judge_lts_gold.py; criteria of the LTS buffer-fix review, with the 10 Pa floor of criterion 4
  that Chunhui Zhao accepted on 2026-10-20): ALL PASS. 1 SCEC pass and fail the same as the gold's; 2 arrival within
  1.15 ms, peak slip rate within 0.098 %, final slip within 0.047 % at every station; 3 noise RMS within 0.6 % of the
  gold's; 4 the period-2 part of the worst station within 2x the gold's or below 10 Pa in every 1 s window: tau_h
  above 2x in 1 window (11-12 s), at most 5.1 Pa against the gold's 1.6 Pa; sigma_n above 2x in 5 windows (7-12 s), at
  most 4.3 Pa against the gold's 0.4-1.9 Pa. All are late windows (from 7 s on), where the gold's own level is a few
  Pa; inferred to be what the 3 buffer layers leave of the seam ripple.
- Whole traces (plots/compare_gold.txt, max |LTS - gold| / max |gold| per station, groups as compare_stations.py
  --groups, 41,614 rows): worst slip 1.293e-03, V 1.569e-02, tau 7.180e-03, sigma_n 8.850e-04. The LTS run before the
  buffer fix (job 1081279, 2d09c44) gave 1.031e-02, 2.222e-02, 1.052e-02, 1.378e-03. The job's own checks: global
  stepping against the gold over 0-1 s 2.965e-12 (round-off); LTS against global stepping over 0-1 s: slip 3.334e-04,
  V 3.263e-04, tau 4.127e-05, sigma_n 2.338e-05.
- vs DRDG3D 100 m (plots/compare.txt): median |misfit| arrival 0.0420 s, final slip 1.50 %, peak h-slip-rate 19.36 %;
  pass tests 1 FAIL, 2 FAIL, 3 peak within 10% at 2 of 16 stations FAIL (pass and fail the same as the global-stepping
  gold's, criterion 1 above).
- plots/: tpv205_200m_p1_lts_vista1082586_stations.pdf and its page PNGs (slip-rate overview of the 16 stations, then
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
    sbatch --dependency=afterok:<build job> lts_tpv205_200m_p1.sbatch

The jobs refuse to start unless the checkout is clean, HEAD is on GitHub, the driver was built at HEAD and the mesh
sha256 is in meshes.sha256. Checks and plots (from miniapps/seas, at the commit that adds this folder, with
D = this folder and G = tpv205/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid1076378_git81759da):

    python3 -I jobs/vista_profiling/compare_stations.py $G/results $D/results --tmax 12.0 --time-tol 2e-9 --groups
    python3 -I jobs/vista_profiling/judge_lts_gold.py tpv205 $G/results $D/results /tmp/judge > $D/plots/lts_vs_gts_gold.txt
    python3 -I jobs/vista_profiling/plot_lts_gold.py tpv205 $D/results $G/results tpv205/benchmark_data/DRDG3D_100m_O4 \
        $D/plots/tpv205_200m_p1_lts_vista1082586 'LTS p=1 200 m (Vista 1082586)' 'gold (Vista 1076378)' 'DRDG3D 100 m O4'
    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv205/stations/tpv205_fault_stations.csv \
        --prefix tpv205 --run lts_p1_200m=$D/results \
        --ref drdg3d_100m=tpv205/benchmark_data/DRDG3D_100m_O4 --ref drdg3d_200m=tpv205/benchmark_data/DRDG3D_200m_O4 \
        --gold mixedflux_p1=tpv205/gold/mixedflux_p1_rk45_jobidcombined_7762879_7763019_git04674fc/results \
        --gold mixedflux_p2=tpv205/gold/mixedflux_p2_rk45_jobid7762880_git04674fc/results \
        --gold gts_p1_200m=$G/results:nitsche --peak-min-stations 14 \
        --json $D/plots/compare_summary.json > $D/plots/compare.txt

Rerunning these reproduces lts_vs_gts_gold.txt and compare.txt (apart from compare.txt's last line, the path it
wrote the JSON to) and the plots. compare_gold.txt holds the per-station numbers of the first command for this run
and for round 1 (job 1081279, whose station files stayed off the repo).

- Requested by Chunhui Zhao (Nitsche profiling, plan item 7: local time stepping golds, 2026-10-20).
