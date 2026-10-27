#!/bin/bash
# G3 step 1: seas-mfem BP3-QD d60 thrust, p = 1, 1500 yr, 4 ranks, at 9503fe25 (QD stages
# through the 2D discrete Green's function).  Lock taken and released by locked_exec.sh
# (released on any exit); peak RSS and wall go to g3-bp3/run_stats.log.
source /Users/chunhuizhao/miniforge/etc/profile.d/conda.sh; conda activate mumps-petsc
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
G=/Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/g3-bp3
R=$G/seas_p1_d60_thrust
cd /Users/chunhuizhao/projects/seas-mfem-2d-inplane/miniapps/seas
CMD="mpirun -np 4 ./seas_fd_cycle_driver --deck config/fd_cycle/bp3_qd_d60_thrust.toml --output-dir $R/out"
{ echo "start $(date)"; echo "commit $(git rev-parse HEAD) dirty $(git status --porcelain --untracked-files=no | wc -l)"
  echo "binary sha256 $(shasum -a 256 seas_fd_cycle_driver | cut -d' ' -f1)"; echo "env mumps-petsc, OMP_NUM_THREADS=1"; echo "cmd $CMD"; } > $R/run_info.txt
t0=$(date +%s); echo $t0 > $R/t0.txt
$G/locked_exec.sh $CMD > $R/run.log 2>&1; rc=$?
echo "end $(date) rc $rc wall $(( $(date +%s) - t0 )) s" >> $R/run_info.txt
echo "rc $rc"
