# BP5-FD p = 2, 0 to 1800 yr, Vista 576 ranks (jobs 1094927 + 1094928)

Produced from commit `3926e9c81c126214d9c12d4d9136d0ebec2f7d34` (branch `fully-dynamic-cycle`) of https://github.com/chunhuizhao478/seas-mfem:
https://github.com/chunhuizhao478/seas-mfem/commit/3926e9c81c126214d9c12d4d9136d0ebec2f7d34

The Vista checkout had no uncommitted changes in tracked files; the only untracked files were the build's
`miniapps/seas/seas_fd_cycle_driver` and `miniapps/seas/BUILD_COMMIT_VISTA` (= 3926e9c8). The job script refuses
to run otherwise (clean tree, HEAD on GitHub, binaries built at HEAD). Every output file header records the commit.

This README was written on 2026-10-26, after the run, from the job logs, sacct and the commit.

## Contents

- `scec_upload/`: the SCEC upload files of the BP5 description, section 4 (30 files): `fltst_strk*dp*.dat` (10 on-fault
  stations), `blkst_strk*.dat` (9 off-fault stations), `global.dat` (domain A source), `catalog.dat` (7 events),
  `slip_2_depth.dat`, `slip_3_depth.dat`, `stress_2_depth.dat`, `stress_3_depth.dat` and the `_strike` four (profiles),
  `rupture.dat` (event 1, domain A). Headers: problem SEAS Benchmark BP5-FD, code seas-mfem (seas_fd_cycle_driver),
  version 3926e9c8..., modeler Chunhui Zhao, element_size 1000 m, polynomial_degree 2.
- `jobs/`: the job script, submit script, build script and deck exactly as at the commit, plus the Slurm out/err of the
  build (1094923), tests (1094925), check (1094926) and the two run links (1094927, 1094928).
- `logs/`: the driver logs of both links (`run_<job>.log`), per-node memory logs and the regime switch file.
- `compare/`: station comparison against BICyclE and SBEM to 1800 yr (PDF, table, plotting script).
- `MANIFEST.sha256`: sha256 of every file here.

Not included (too large; on Vista in the run folder): `bp5_fd_p2_1800yr_fault_snapshots.vtkhdf` (83.2 GB),
`bp5_fd_p2_1800yr_surface_snapshots.vtkhdf` (12.6 GB), the raw `*.txt` outputs (308 MB) and the checkpoints.

## Jobs

- Cluster: Vista (TACC), partition gg, account ear20006, 4 nodes x 144 = 576 MPI ranks (ibrun).
- Build 1094923 (`vista_build.sbatch`, FDC_PART=build), tests 1094925 (FDC_PART=tests; `test_cycle_bench_run: PASSED`,
  0 failing steps), check 1094926 (`vista_p2_check.sbatch`; `p2 check: 0 failed`).
- Link 1: job 1094927, nodes i618-[044,051-053], 2026-10-24T08:51:23 to 2026-10-26T05:51:40 (CDT), t = 0 to 1575 yr,
  stopped at the 45 h wall limit with a checkpoint (exit 8, the designed hand-off), COMPLETED.
- Link 2: job 1094928, nodes i615-[021-024], 2026-10-26T05:52:13 to 13:04:14 (CDT), 1575 to 1800 yr, restarted from
  link 1's checkpoint. The driver finished with `status: done`, final checkpoint at t = 56803680000.00 s = t_final.
  Slurm reports this job FAILED (exit 9) and the spare link 1094929 was cancelled: a false alarm in the job script's
  end check, which compares the summary's 3-digit `final time: 5.68e+10` against `t_final: 5.68037e+10`. The run itself
  reached t_final; the DONE marker and the script's automatic SCEC export were skipped because of it, so the export was
  run by hand (below). The end check is fixed in the commit "fd_cycle: the Vista p = 2 chain ends on the driver's
  "status: done"" (the job script only; the driver is unchanged).
- Submitted with `bash miniapps/seas/jobs/fd_cycle/submit_vista_p2.sh` on the login node (build, tests, check and
  three links, chained with afterok).

## Inputs

- Deck: `miniapps/seas/config/fd_cycle/bp5_fd_p2_1800yr.toml` (copy in `jobs/`): order 2, t_final 5.680368e10 s,
  IC scec_fd + eta_V_init, QD regime MUMPS + PETSc TS rk5dp (atol 1e-8) with the Green's-function H-matrix stage
  operator, dynamic regime global central-difference stepping (LTS off), switch at V_on 1e-2 / V_off 1e-3 m/s.
- Mesh (path relative to `miniapps/seas`, scaled x 1000), sha256:

```
f413b957cf1c6d02d2c84b65fc05d1f99232301ed8083cb7d8533c9fa442c3b4  bp5/mesh/reference/bp5_tandem_exact.msh
```

- H-matrix cache (read, not built, by both links): `$SCRATCH/seas-qd-green/bp5_fd_p2_h1e-10.bin`, 40697028608 bytes
  (Vista: /scratch/10024/zhaochun/seas-qd-green/), built by the check job's code at this commit.

## Build

- Binary `miniapps/seas/seas_fd_cycle_driver`, sha256 275416441822d68d01c353fc9b4e32ff74127e08eae42cba16f3d50e53d11a44.
- Modules: gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5. Deps in `$WORK/seas-mfem-fdcycle-deps`.
- Built by `sbatch vista_build.sbatch` (FDC_PART=build) from `miniapps/seas/jobs/fd_cycle`; its log is
  `jobs/fdc_build_1094923.out` (`[build_vista] done: deps mfem seas at commit 3926e9c8...`).

## Run

From `miniapps/seas` (as in `vista_p2_1800yr.sbatch`):

```
ibrun ./seas_fd_cycle_driver --deck config/fd_cycle/bp5_fd_p2_1800yr.toml \
   --output-dir $WORK/seas-mfem-fdcycle-runs/p2_1800_1094927 \
   --qd-green-hcache $SCRATCH/seas-qd-green/bp5_fd_p2_h1e-10.bin --wall-limit 162000 --log-every 1000
# link 2 adds: --restart $WORK/seas-mfem-fdcycle-runs/p2_1800_1094927/bp5_fd_p2_1800yr
```

Run summary: 705046 TS steps, 30320 rejections, 7 dynamic phases / events, fault and surface snapshots 23250 each.

## SCEC export

From `miniapps/seas` at the same commit, on Vista:

```
python3 -I cycle/scec_bp5_export.py $WORK/seas-mfem-fdcycle-runs/p2_1800_1094927 bp5_fd_p2_1800yr \
   $WORK/seas-mfem-fdcycle-runs/p2_1800_1094927/scec_upload --modeler "Chunhui Zhao" --element-size "1000 m" --order 2
```

(the job script's own call adds `--version <HEAD>`, the same commit the headers already give). 44723 rows per time
series, 0 repeated times dropped; profiles 3171 rows; rupture.dat 6160 positions from the event-1 file.

## Results against BICyclE and SBEM (see `compare/`)

- Event starts (yr): 0, 291.81, 583.50, 875.15, 1166.78, 1458.42, 1750.06. Recurrence settles at 291.6 yr against
  297.4 (BICyclE) and 298.3 (SBEM), about 2 % short.
- Event slip at the rupture stations within about 1-3 %; total slip at 1800 yr 1.1-1.8 % low at every station.
- Known open items: edge stations (strk+-36dp+00, strk+00dp+22) creep more after each event than the references;
  a spurious sigma_n drift (it should stay at 25 MPa), smaller than at p = 1 but still present.
- The dynamic phases have no absorbing boundary: the side walls (x = +-200 km, y = +-100 km) hold the plate-loading
  displacement and the bottom (100 km depth) is traction-free, so waves reflect back to the fault, about 30 to 90 s
  after leaving it (estimated from the wave speeds), within dynamic phases of 326 to 330 s (`logs/*_switch.txt`).
  Chunhui Zhao, 2026-10-26: dashpots on all outer walls are the next change; this run is from before it.
