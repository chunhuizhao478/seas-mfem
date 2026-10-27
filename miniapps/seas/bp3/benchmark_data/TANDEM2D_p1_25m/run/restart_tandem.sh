#!/bin/bash
# Restart of the G3 Tandem run with ulimit -Ss 65520 (the hard maximum; 65536 was refused),
# keeping the Mac run lock: take over the lock owner first, stop the old run by its PIDs
# (prterun 55633 and tandem ranks 55634-55637, children of script 55600), move the partial
# output aside, start again.
source /Users/chunhuizhao/miniforge/etc/profile.d/conda.sh; conda activate tandem-petsc321
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
G=/Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/g3-bp3/tandem_p1_d60_thrust
TD=/Users/chunhuizhao/projects/tandem
LOCK=/Users/chunhuizhao/projects/.seas-mac-run.lock
grep -q "pid 55600 " $LOCK/owner || { echo "ERROR: lock not held by our run 55600: $(cat $LOCK/owner)"; exit 1; }
echo "pid $$ thread 2d-inplane (G3 Tandem p1 BP3 d60 thrust 1500 yr, restart with stack 65520 KB) since $(date)" > $LOCK/owner
have=1
release() { [ "$have" = 1 ] && grep -q "pid $$ " "$LOCK/owner" 2>/dev/null && rm -f "$LOCK/owner" && rmdir "$LOCK" && echo "lock released $(date)"; }
trap release EXIT
kill -TERM 55633 55634 55635 55636 55637 2>/dev/null
for i in $(seq 1 30); do kill -0 55600 2>/dev/null || break; sleep 1; done
kill -0 55634 2>/dev/null && { kill -KILL 55633 55634 55635 55636 55637 2>/dev/null; sleep 2; }
kill -0 55600 2>/dev/null && { echo "ERROR: old script 55600 still alive"; exit 1; }
echo "old run stopped $(date)"
mkdir -p $G/aborted_stack7MiB
mv $G/out $G/run.log $G/run_info.txt $G/t0.txt $G/run_tandem.out $G/aborted_stack7MiB/ 2>/dev/null
[ -d $G/checkpoint ] && mv $G/checkpoint $G/aborted_stack7MiB/
mkdir -p $G/out
ulimit -Ss 65520 || { echo "ERROR: ulimit -Ss 65520 failed"; exit 1; }
cd $G
CMD="mpirun -np 4 $TD/build_2d_p1/app/tandem bp3.toml --petsc -options_file $TD/examples/options/lu_mumps.cfg -options_file $TD/examples/options/rk45.cfg -ts_monitor"
{ echo "start $(date)"; echo "restart: the first start (22:22:10 PDT, ulimit -Ss 65536 refused, stack 8176 KB) was stopped after about 2 min and moved to aborted_stack7MiB/; this run uses ulimit -Ss 65520 (the hard maximum)"
  echo "tandem commit $(git -C $TD rev-parse HEAD) (status: $(git -C $TD status --short | grep -v '^??' | tr '\n' ' '))"
  echo "PETSc env $CONDA_PREFIX ($(grep -E 'define PETSC_VERSION_(MAJOR|MINOR|SUBMINOR) ' $CONDA_PREFIX/include/petscversion.h | awk '{printf $3"."}'))"
  echo "ulimit -s $(ulimit -s)"; echo "cmd $CMD"; } > run_info.txt
t0=$(date +%s); echo $t0 > t0.txt
$CMD > run.log 2>&1; rc=$?
echo "end $(date) rc $rc wall $(( $(date +%s) - t0 )) s" >> run_info.txt
echo "rc $rc wall $(( $(date +%s) - t0 )) s"
