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
```

Each speed and LTS job stops if the driver was built at another commit (`BUILD_COMMIT_VISTA`)
or if the mesh sha256 differs from `meshes.sha256`.

## Results (18 and 19 Oct 2026)

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
