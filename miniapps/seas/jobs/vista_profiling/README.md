# Vista jobs for the Nitsche profiling speedups

Branch `claude/nitsche-profiling-w336yq` (Chunhui Zhao, 2026-10-17). The jobs check
and time plan items 1, 3, 4, 5, 6 and 8a of the Nitsche profiling plan on Vista
(account EAR26024, partition `gg`, 144 cores per node). The plan and the reports
are in `miniapps/seas/document/nitsche_profiling_dev/`.

## Files

| File | What it does |
|---|---|
| `build_profiling_vista.sbatch` | 1 node. The CLAUDE.md guard; `build_vista.sh` (MFEM, `seas_nitsche_dyn_driver`, the six gate tests), then `seas_bench_nitsche_kernels`; the six gate tests; copies the 200 m gold mesh from `$WORK/seas-mfem-benchmarks` and checks its sha256 against `meshes.sha256`; makes the 1000 m bench mesh; runs the bench on 16 ranks (each new call against the call it replaces) and the fault writers at the 200 m fault size on 1 rank (`[H]`). |
| `speed_compare_tpv205_200m_p2.sbatch` | 16 nodes x 144 ranks, as gold 1075916. TPV205, 200 m gold mesh, p = 2, 0 to 1.0 s, deck `tpv205/configs/tpv205_nitsche.toml`. Runs `off` (`--fast-fault 0 --fused-spmv 0 --fast-fault-output 0`), `on` (the defaults) and `weighted` (`--fault-weight 0.25`) in one allocation, then compares the stations: `off` against the first 1.0 s of gold 1075916, `on` and `weighted` against `off`. |
| `speed_dense_tpv205_200m_p2.sbatch` | Plan item 8a (2026-10-18), the same 16 x 144 ranks, mesh, deck and 1.0 s. Runs `off` (`--fast-fault 0 --fused-spmv 0 --fast-fault-output 0 --dense-blocks 0`), `item4` (`--dense-blocks 0`, the defaults before item 8a) and `dense` (the defaults), then compares the stations: `off` against the first 1.0 s of gold 1075916, `item4` and `dense` against `off`, and `dense` against `item4`. |
| `build_lts_vista.sbatch` | Plan item 7, local time stepping (2026-10-19), 1 node. The CLAUDE.md guard; `build_vista.sh`, then `seas_test_nitsche_lts`; the six gate tests and the LTS suite (L-1 to L-10); copies the TPV205 and TPV104 200 m gold meshes from `$WORK/seas-mfem-benchmarks` and checks their sha256 against `meshes.sha256`. |
| `lts_tpv205_200m_p2.sbatch`, `lts_tpv104_200m_p2.sbatch`, `lts_tpv205_200m_p1.sbatch`, `lts_tpv104_200m_p1.sbatch` | Plan item 7, P5: the four new LTS golds, each with the deck, mesh, order and rank layout of the global-stepping gold it replaces (1075916, 1075917, 1076378, 1075919). Each runs `gts` (the defaults, 0 to 1.0 s, the speed reference at this commit) and `lts` (`--lts rate2`, 0 to 12.0 s, with the driver's default `--lts-fault-buffer 3` from the commit after 2d09c44) through `run_lts_vista.sh`, then compares the stations: `gts` against the first 1.0 s of the gold, `lts` against `gts` over 1.0 s and against the gold over 12 s. |
| `run_lts_vista.sh` | Shared body of the four LTS jobs. The LTS station comparisons use `--time-tol 2e-9`: the stations print times as `%.10e`, whose quantum at t >= 10 s is 1e-9 s, so the 1e-9 of 2d09c44 stopped the 12 s comparison ("times differ"). |
| `compare_stations.py` | Station comparison used by the speed jobs (standard library only). `--time-tol` and `--groups` for the LTS jobs. |
| `judge_lts_gold.py` | Item 7 (2026-10-20), run off Vista from `miniapps/seas`: checks an LTS run against the global-stepping gold it sits beside over 12 s with the criteria of the LTS buffer-fix review (SCEC pass and fail as the gold; arrival within 5 ms, peak slip rate within 1 %, final slip within 0.5 %; noise RMS within 10 %; period-2 part within 2x the gold's or below 10 Pa in every 1 s window). Output `plots/lts_vs_gts_gold.txt` of each LTS gold folder. |
| `plot_lts_gold.py` | Item 7 (2026-10-20): station plots of an LTS run against its global-stepping gold and DRDG3D (one PDF and its page PNGs in each LTS gold folder's `plots/`). |
| `m0_lts_tpv205_200m_p2.sbatch` | M0 of the speed-next plan (rev 1.1, §6; 2026-10-20): 16 x 144 ranks as LTS gold 1082384. TPV205, 200 m gold mesh, p = 2, `--lts rate2`, 0 to 1.0 s, four variants with existing flags only: `base` (the defaults), `every10` (`--station-every 10`), `first` (`--station-copies first`) and `deflate0` (deck `tpv205/configs/tpv205_nitsche_m0_deflate0.toml`, fault deflate level 0). Measures where the time outside the step goes (stations or not) and the frame cost; `base` and `deflate0` stations are compared with the first 1.0 s of 1082384 (meant to be bitwise equal). |
| `m0_lts_tpv205_200m_p1.sbatch` | M0, 4 x 144 ranks as LTS gold 1082586: TPV205, 200 m, p = 1, `--lts rate2`, 0 to 12.0 s with the deflate-level-0 deck copy: the whole-run frame cost; stations compared with 1082586 (meant to be bitwise equal). |
| `run_m0_vista.sh` | Shared body of the two M0 jobs and the two V1 jobs (`TAG` names the output folder). |
| `build_map_vista.sbatch` | Profiling (2026-10-20): builds `seas_nitsche_dyn_driver_g` (the driver with `-g`, Makefile target of the same name) after `build_lts_vista.sbatch`, and checks that its object without the debug sections is byte-identical to the production driver's. |
| `map_lts_tpv205_200m_p1.sbatch` | Profiling: TPV205, 200 m gold mesh, p = 1, `--lts rate2`, 0 to 3.0 s, 4 x 144 ranks as LTS gold 1082586, under Linaro MAP (`map --profile`): time per function and part of the code over the whole run, setup included. |
| `map_lts_tpv205_200m_p2.sbatch` | Profiling: the same at p = 2, 0 to 2.0 s, 16 x 144 ranks as LTS gold 1082384. |
| `run_map_vista.sh` | Shared body of the two profiling jobs: the guard, the MAP run (`.map` file, `map_functions.csv`, `profile.json`), the station check against the gold, and a perf fallback on four ranks if MAP fails. |
| `perf_wrap.sh` | The perf fallback's per-rank wrapper. |
| `build_speed_next_vista.sbatch` | Speed-next O1, O2 and the setup timers (2026-10-20), 1 node: `build_lts_vista.sbatch` plus `seas_test_station_output` and its test (O2-1 to O2-4 on 1, 2, 3 and 4 ranks). |
| `v1_lts_tpv205_200m_p1.sbatch` | Speed-next V1: TPV205, 200 m gold mesh, p = 1, `--lts rate2`, 0 to 12.0 s, 4 x 144 ranks as LTS gold 1082586. Variants `o2` (the defaults: copy average over the ranks holding copies, rows flushed every 100) and `o1o2` (the same with `--paraview-fault-deflate-level none`); stations compared with 1082586 over the 12 s (meant to be bitwise equal); wall times against its 359 s. |
| `v1_lts_tpv205_200m_p2.sbatch` | V1 at p = 2: 16 x 144 ranks as LTS gold 1082384, against its 739 s. |
| `meshes.sha256` | sha256 of the TPV205 and TPV104 200 m gold meshes. |
| `results/` | Job stdout of the runs below, as relayed from Vista. |

## Submit

```
cd $WORK && git clone --recurse-submodules -b claude/nitsche-profiling-w336yq \
    https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-profiling
cd seas-mfem-profiling/miniapps/seas/jobs/vista_profiling
sbatch build_profiling_vista.sbatch
# after it passes, at the same commit:
sbatch speed_compare_tpv205_200m_p2.sbatch   # items 1-6
sbatch speed_dense_tpv205_200m_p2.sbatch     # item 8a
# item 7 (local time stepping), at the LTS commit:
sbatch build_lts_vista.sbatch
# after it prints "all checks passed", at the same commit:
sbatch lts_tpv205_200m_p2.sbatch; sbatch lts_tpv104_200m_p2.sbatch
sbatch lts_tpv205_200m_p1.sbatch; sbatch lts_tpv104_200m_p1.sbatch
# (Vista's sbatch prints its banner on stdout before the job id, hence the tail -n 1)
# M0 of the speed-next plan, at the M0 commit: the build, then both runs after it
B=$(sbatch --parsable build_lts_vista.sbatch | tail -n 1)
sbatch --dependency=afterok:$B m0_lts_tpv205_200m_p2.sbatch
sbatch --dependency=afterok:$B m0_lts_tpv205_200m_p1.sbatch
# profiling, at the same commit: the -g driver after the build, then both MAP runs after it
G=$(sbatch --parsable --dependency=afterok:$B build_map_vista.sbatch | tail -n 1)
sbatch --dependency=afterok:$G map_lts_tpv205_200m_p1.sbatch
sbatch --dependency=afterok:$G map_lts_tpv205_200m_p2.sbatch
# speed-next V1, at the implementation commit, once the M0 and profiling jobs have finished
# (one clone, so one commit at a time): the build and tests, then both runs after it
V=$(sbatch --parsable build_speed_next_vista.sbatch | tail -n 1)
sbatch --dependency=afterok:$V v1_lts_tpv205_200m_p1.sbatch
sbatch --dependency=afterok:$V v1_lts_tpv205_200m_p2.sbatch
# profiling again at the MAP route commit (MAP starts the job itself, through ibrun or else srun;
# run_map_vista.sh), after the speed-next build V above at the same commit: the -g driver, then both runs
G=$(sbatch --parsable --dependency=afterok:$V build_map_vista.sbatch | tail -n 1)
sbatch --dependency=afterok:$G map_lts_tpv205_200m_p1.sbatch
sbatch --dependency=afterok:$G map_lts_tpv205_200m_p2.sbatch
```

Each speed and LTS job stops if the driver was built at another commit (`BUILD_COMMIT_VISTA`)
or if the mesh sha256 differs from `meshes.sha256`.

## Results (18 to 21 Oct 2026)

Clone `$WORK/seas-mfem-profiling`, modules gcc/14.2.0, cuda/12.6, openmpi/5.0.5,
hypre/2.31.0, metis/5.1.0.3, phdf5/1.14.6.

| Job | Commit | Script | Result |
|---|---|---|---|
| 1078824 | 6d41ba7 | `build_profiling_vista.sbatch` | Build and the six gate test sets pass (74, 22, 49, 13, 9, 42 passed; 0 failed). The 16-rank bench crashed after its `[V]` section: MFEM partial assembly deletes the mesh's face-neighbour data, which the traction checks need on ranks with shared fault faces. Bench only; fixed in 0c83fb5. |
| 1078825 | 6d41ba7 | `speed_compare_tpv205_200m_p2.sbatch` | `results/tpv205_200m_p2_speed_1078825.out`. `off` equals gold 1075916 bitwise at all 16 stations; `on` and `weighted` differ from `off` by at most 6.6e-11 relative. Time in `ode.Step` on the slowest rank: 27.00 ms per step `off`, 16.38 ms `on` (1.65x), 16.26 ms `weighted` (1.66x). Wall time of the 1 s run: 321 s `off`, 244 s `on`. First-node memory peak: 147, 146, 161 GB of 238 GB. |
| 1078847 | 0c83fb5 | `build_profiling_vista.sbatch` | Build and gate tests pass. `results/bench_np16_1078847.log`: 9 of 9 checks pass on 16 ranks (10 shared fault faces); `ComputeTractionFast` 38.5x and `AddSlipLoadFast` 47x faster than the calls they replace. `results/bench_h_1078847.txt`: fault writers at the 200 m fault size. |
| 1079599 | 7cab085 | `build_profiling_vista.sbatch` | Build and the six gate test sets pass (74, 22, 49, 13, 9, 42 passed; 0 failed). Bench on 16 ranks (TPV205 1000 m, p = 1): 14 of 14 checks pass; `[D]` dense blocks give `K u + K_eta v` 1.52x faster than the two hypre products and 1.57x faster than the item 4 sweep (minimum of 60 pairs). The bench log stays on Vista (`miniapps/seas/bench_np16_1079599.log` in the clone). |
| 1079640 | 7cab085 | `speed_dense_tpv205_200m_p2.sbatch` | `results/tpv205_200m_p2_dense_1079640.out`. `off` equals gold 1075916 bitwise at all 16 stations; `item4` and `dense` differ from `off` by at most 6.6e-11 and 2.7e-11 relative. Time in `ode.Step` on the slowest rank: 27.28 ms per step `off`, 16.91 ms `item4`, 9.36 ms `dense` (1.81x against `item4`, 2.91x against `off`). Wall time of the 1 s run: 322, 257, 197 s. Dense blocks 68.9 GiB over all ranks; first-node memory peak 138.6, 138.6, 138.2 GiB of 237.6 GiB. |
| 1081276 | 2d09c44 | `build_lts_vista.sbatch` | Build, the six gate test sets (74, 22, 49, 13, 9, 42 passed) and the LTS suite (24 of 24) pass; both 200 m meshes match `meshes.sha256`. 7:30 on 1 node. |
| 1081277, 1081278, 1081279, 1081280 | 2d09c44 | `lts_tpv205_200m_p2`, `lts_tpv104_200m_p2`, `lts_tpv205_200m_p1`, `lts_tpv104_200m_p1` | `results/lts_<case>_<job>.out` and `.err`. `gts` against the old gold over 1 s: 2.7e-11, 1.2e-11, 3.0e-12, 1.0e-10 relative; `lts` against `gts` over 1 s: 1.1e-3, 2.4e-3, 1.9e-3, 7.0e-3 (max V). The job's 12 s comparison stopped at "times differ" (the tolerance above) and was redone off Vista with 2e-9. Arrival, peak slip rate and final slip match the old golds (arrival within 3.5 ms but one station, peak within 2 %, slip within 1.1 %), with the same SCEC pass and fail. But the LTS stations carry a ripple of period 2 dt0 that global stepping does not have (tau_h up to 0.9 MPa at the surface stations of TPV205 p2, stick-slip chatter that moved the arrival at x2 = 0 on the surface 55 ms early), so these runs are not golds. The cause and the fix (`--lts-fault-buffer`, default 3) are in the comment of `solver/nitsche_lts_run.hpp`; the four jobs are to be rerun with it. Time in `ode.Step` on the slowest rank, `gts` then `lts`: 9.89 and 5.27 ms per step (TPV205 p2), 11.94 and 6.02 (TPV104 p2), 7.87 and 5.80 (TPV205 p1), 15.39 and 7.91 (TPV104 p1); the maxima include waiting on rank 0's fault frame writes (124-197 s of each 12 s `lts` run). |
| 1082383 | ac3adeb | `build_lts_vista.sbatch` | Build, the six gate test sets (74, 22, 49, 13, 9, 42 passed) and the LTS suite (24 of 24) pass; both 200 m meshes match `meshes.sha256`. |
| 1082386, 1082387 | ac3adeb | `lts_tpv205_200m_p1`, `lts_tpv104_200m_p1` | Stopped out of memory at startup on node i617-131 (METIS partitioning in the `gts` runs, `std::bad_alloc` after the LTS setup in 1082387's `lts` run); `results/lts_tpv205_200m_p1_1082386.out` and `.err`, `results/lts_tpv104_200m_p1_1082387.out` (its `.err` stayed on Vista). The `lts` run of 1082386 finished, with station files bitwise equal to those of 1082586. Chunhui Zhao had both rerun with `sbatch --exclude=i617-131` (1082586, 1082587) and the waiting p2 jobs kept off that node (`scontrol update ExcNodeList=i617-131`). |
| 1082384, 1082385, 1082586, 1082587 | ac3adeb | `lts_tpv205_200m_p2`, `lts_tpv104_200m_p2`, `lts_tpv205_200m_p1`, `lts_tpv104_200m_p1` | The four LTS golds, `tpv205/gold/` and `tpv104/gold/` `nitsche-kv0p1_p<order>_cd_pc_lts-rate2-buf3_200m_jobid<job>_gitac3adeb/` (stations, job stdout, checks, plots and a README each). `gts` against the old gold over 1 s: 2.7e-11, 1.2e-11, 3.0e-12, 1.0e-10 relative. Against the global-stepping golds over 12 s (`judge_lts_gold.py`): arrival within 1.15 ms, peak slip rate within 0.22 %, final slip within 0.09 %, noise RMS within 3.3 %, the same SCEC pass and fail; the period-2 part matches the gold's during the rupture and stays below 10 Pa after it (up to 7 Pa against the gold's 0.4-2.5 Pa). Largest station difference from the gold (max V): 1.1e-2 to 2.3e-2, against 1.7e-2 to 1.9e-1 before the buffer fix. Time in `ode.Step` on the slowest rank, `gts` then `lts`: 9.42 and 5.24 ms per step (TPV205 p2, 1.80x), 11.96 and 6.56 (TPV104 p2, 1.82x), 7.87 and 6.07 (TPV205 p1, 1.30x), 11.56 and 8.14 (TPV104 p1, 1.42x). |
| 1085105 | eb6b16a | `build_lts_vista.sbatch` | Build, the six gate test sets (74, 22, 49, 13, 9, 42 passed) and the LTS suite (24 of 24) pass; both 200 m meshes match `meshes.sha256`. 7:35 on 1 node. |
| 1085110 | eb6b16a | `m0_lts_tpv205_200m_p1.sbatch` | M0 at p = 1 with the fault frames at deflate level 0, 0 to 12.0 s (41,616 ticks): wall 252.6 s against the 359.2 s of LTS gold 1082586 (level 1), 1.42x. Time in `ode.Step` on the busiest rank 146.2 s (3.51 ms per tick) against the gold's 252.8 s; halo wait at most 66.1 s against 136.6 s. `fault.vtkhdf` 15.2 GB. Stations bitwise equal to 1082586 (16 of 16). First-node memory peak 120.3 of 237.6 GB. |
| 1085111 | eb6b16a | `build_map_vista.sbatch` | `seas_nitsche_dyn_driver_g` built, 8 debug sections. Without its debug sections the -g object differed in bytes from the production object, with every section the same size; the job compares the disassembly from the next commit on. |
| 1085112 | eb6b16a | `map_lts_tpv205_200m_p1.sbatch` | No profile. Every rank stopped at startup (`/home1/apps/forge/24.1.1/lib/libstdc++.so.6: version 'GLIBCXX_3.4.32' not found`): module map/24.1.1 puts Forge's libstdc++ ahead of gcc 14's. The perf fallback stopped the same way, and the job still exited 0. Both fixed in `run_map_vista.sh` at the next commit. |
| 1085109 | eb6b16a | `m0_lts_tpv205_200m_p2.sbatch` | M0 at p = 2, 0 to 1.0 s (7,400 ticks), four variants. Wall: `base` 160.2 s, `every10` 144.4 s, `first` 145.5 s, `deflate0` 156.1 s. The 50 fault frames took 6.34 s at deflate level 1 (0.127 s each) and 1.96 s at level 0 (0.039 s each); `fault.vtkhdf` 133.6 MB and 1,280.1 MB. Time in `ode.Step` on the busiest rank: 41.79, 45.61, 52.59 and 37.69 s. The stations of `base` and `deflate0` are bitwise equal to the first 1.0 s of 1082384 (16 of 16), and all four end states are the same. |
| 1085113 | eb6b16a | `map_lts_tpv205_200m_p2.sbatch` | Cancelled before it started (Chunhui Zhao, 2026-10-21), since it would have stopped as 1085112 did. |
| 1085934 | f81f374 | `build_lts_vista.sbatch` | Build, the six gate test sets (74, 22, 49, 13, 9, 42 passed) and the LTS suite (24 of 24) pass; both 200 m meshes match `meshes.sha256`. |
| 1085935 | f81f374 | `build_map_vista.sbatch` | `seas_nitsche_dyn_driver_g` built. Without its debug sections the -g object differs in bytes from the production object, with every section the same size and the same disassembly. |
| 1085936 | f81f374 | `map_lts_tpv205_200m_p1.sbatch` | The driver now loads gcc 14's libstdc++, but MAP exited 8 s after "Waiting for licence" without a profile, and its message was lost (both addressed at the next commit). The perf fallback ran 0 to 3.0 s (wall 138.2 s; `ode.Step` at most 56.2 s), with stations bitwise equal to the first 3.0 s of 1082586 and perf profiles of ranks 0, 1, 288 and 575 (`perf/`, on Vista). The job exited 1. |
| 1085937 | f81f374 | `map_lts_tpv205_200m_p2.sbatch` | As 1085936 at p = 2 (16 x 144 ranks, 0 to 2.0 s): MAP exited 9 s after "Waiting for licence". The perf fallback ran (wall 219.3 s, 14,792 ticks; `ode.Step` at most 79.5 s), with stations bitwise equal to the first 2.0 s of 1082384 (16 of 16) and perf profiles of ranks 0, 1, 1152 and 2303. The job exited 1. |
| 1086943 | 488052f | `build_lts_vista.sbatch` | Build and all gate tests pass (0 failed); both 200 m meshes match `meshes.sha256`. 7:28 on 1 node. |
| 1086944 | 488052f | `build_map_vista.sbatch` | `seas_nitsche_dyn_driver_g` built, 8 debug sections; without them every section has the size and the disassembly of the production object. |
| 1086945 | 488052f | `map_lts_tpv205_200m_p1.sbatch` | MAP checked out its licence and started its express launch with `--mpi=openmpi`, then stopped: "pindirect exited before reaching MPIR_Breakpoint" (express launch waits for mpirun's MPIR interface, which Open MPI 5's mpirun no longer has; addressed at the next commit). The perf fallback ran 0 to 3.0 s (wall 144.5 s; `ode.Step` at most 51.9 s), with stations bitwise equal to the first 3.0 s of 1082586 (16 of 16). The job exited 1. |
| 1086946 | 488052f | `map_lts_tpv205_200m_p2.sbatch` | As 1086945 at p = 2: MAP stopped after 9 s with the same "pindirect exited before reaching MPIR_Breakpoint". The perf fallback ran 0 to 2.0 s (wall 223.9 s; `ode.Step` at most 87.7 s), with stations bitwise equal to the first 2.0 s of 1082384 (16 of 16) and perf profiles of ranks 0, 1, 1152 and 2303. The job exited 1. |
| 1087793 | 08ec423 | `build_speed_next_vista.sbatch` | Build, the six gate test sets (74, 22, 49, 13, 9, 42 passed), the LTS suite (24 of 24) and the station output test (168 checks) pass; both 200 m meshes match `meshes.sha256`. 7:32 on 1 node. |
| 1087795 | 08ec423 | `build_map_vista.sbatch` | `seas_nitsche_dyn_driver_g` built, 8 debug sections; without them every section has the size and the disassembly of the production object. |
| 1087796 | 08ec423 | `map_lts_tpv205_200m_p1.sbatch` | No profile: MAP refused both routes at its argument check, "--nodes can only be specified with --queue" (`--profile` runs without the queue; fixed at the next commit). The perf fallback ran 0 to 3.0 s (wall 132.3 s; `ode.Step` at most 60.3 s), with stations bitwise equal to the first 3.0 s of 1082586 (16 of 16). Setup, max over ranks: 76.4 s, of which the mesh read 54.7 s, the partition 7.6 s and the LTS clusters 5.5 s. The job exited 1. |

The run outputs (fault files, stations, driver logs) stayed on Vista in
`$SCRATCH/nitsche_profiling/tpv205_200m_p2_speed_1078825/{off,on,weighted}` and
`$SCRATCH/nitsche_profiling/tpv205_200m_p2_dense_1079640/{off,item4,dense}`.
Scratch is purged, and these runs are not golds. The bench log in `results/`
has its repeated lines collapsed (libcuda "(ignored)" lines and duplicated
audit / `[K-DIAG]` blocks). `results/tpv205_200m_p2_dense_1079640.out` is the
job stdout as relayed from Vista, with the run-directory paths shortened to
`.../` and the `du` listing joined onto one line. The analysis is in §7 of
`miniapps/seas/document/nitsche_profiling_dev/nitsche-speedups-impl-rev1.1.pdf`
(items 1 to 6) and §6 of `nitsche-dense-blocks-8a-rev1.1.pdf` in the same
folder (item 8a).
