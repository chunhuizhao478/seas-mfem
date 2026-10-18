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
| `compare_stations.py` | Station comparison used by the speed jobs (standard library only). |
| `meshes.sha256` | sha256 of the TPV205 200 m gold mesh. |
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
```

Each speed job stops if the driver was built at another commit (`BUILD_COMMIT_VISTA`)
or if the mesh sha256 differs from `meshes.sha256`.

## Results (18 Oct 2026)

Clone `$WORK/seas-mfem-profiling`, modules gcc/14.2.0, cuda/12.6, openmpi/5.0.5,
hypre/2.31.0, metis/5.1.0.3, phdf5/1.14.6.

| Job | Commit | Script | Result |
|---|---|---|---|
| 1078824 | 6d41ba7 | `build_profiling_vista.sbatch` | Build and the six gate test sets pass (74, 22, 49, 13, 9, 42 passed; 0 failed). The 16-rank bench crashed after its `[V]` section: MFEM partial assembly deletes the mesh's face-neighbour data, which the traction checks need on ranks with shared fault faces. Bench only; fixed in 0c83fb5. |
| 1078825 | 6d41ba7 | `speed_compare_tpv205_200m_p2.sbatch` | `results/tpv205_200m_p2_speed_1078825.out`. `off` equals gold 1075916 bitwise at all 16 stations; `on` and `weighted` differ from `off` by at most 6.6e-11 relative. Time in `ode.Step` on the slowest rank: 27.00 ms per step `off`, 16.38 ms `on` (1.65x), 16.26 ms `weighted` (1.66x). Wall time of the 1 s run: 321 s `off`, 244 s `on`. First-node memory peak: 147, 146, 161 GB of 238 GB. |
| 1078847 | 0c83fb5 | `build_profiling_vista.sbatch` | Build and gate tests pass. `results/bench_np16_1078847.log`: 9 of 9 checks pass on 16 ranks (10 shared fault faces); `ComputeTractionFast` 38.5x and `AddSlipLoadFast` 47x faster than the calls they replace. `results/bench_h_1078847.txt`: fault writers at the 200 m fault size. |

The run outputs (fault files, stations, driver logs) stayed on Vista in
`$SCRATCH/nitsche_profiling/tpv205_200m_p2_speed_1078825/{off,on,weighted}`.
Scratch is purged, and these runs are not golds. The bench log in `results/`
has its repeated lines collapsed (libcuda "(ignored)" lines and duplicated
audit / `[K-DIAG]` blocks). The analysis is in §7 of
`miniapps/seas/document/nitsche_profiling_dev/nitsche-speedups-impl-rev1.1.pdf`.
