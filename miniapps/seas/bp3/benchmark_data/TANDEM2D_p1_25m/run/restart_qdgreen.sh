#!/bin/bash
# Restart of the G3 Tandem run in QDGreen mode (precomputed discrete Green's function; BP3's
# boundary is linear in t and bp3.toml has boundary_linear = true, which QDGreen requires).
# Keeps the Mac run lock: take over the owner first, stop the default-mode run by its PIDs
# (prterun 55768, tandem 55769-55772, script 55732), move its output aside, start again with
# --mode QDGreen and [gf_checkpoint] prefix = "gf_checkpoint" (inside this folder).
source /Users/chunhuizhao/miniforge/etc/profile.d/conda.sh; conda activate tandem-petsc321
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
G=/Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/g3-bp3/tandem_p1_d60_thrust
TD=/Users/chunhuizhao/projects/tandem
LOCK=/Users/chunhuizhao/projects/.seas-mac-run.lock
grep -q "pid 55732 " $LOCK/owner || { echo "ERROR: lock not held by our run 55732: $(cat $LOCK/owner)"; exit 1; }
echo "pid $$ thread 2d-inplane (G3 Tandem p1 BP3 d60 thrust 1500 yr, QDGreen) since $(date)" > $LOCK/owner
have=1
release() { [ "$have" = 1 ] && grep -q "pid $$ " "$LOCK/owner" 2>/dev/null && rm -f "$LOCK/owner" && rmdir "$LOCK" && echo "lock released $(date)"; }
trap release EXIT
kill -TERM 55768 55769 55770 55771 55772 2>/dev/null
for i in $(seq 1 30); do kill -0 55732 2>/dev/null || break; sleep 1; done
kill -0 55769 2>/dev/null && { kill -KILL 55768 55769 55770 55771 55772 2>/dev/null; sleep 2; }
kill -0 55732 2>/dev/null && { echo "ERROR: old script 55732 still alive"; exit 1; }
echo "default-mode run stopped $(date)"
mkdir -p $G/aborted_default_mode
mv $G/out $G/run.log $G/run_info.txt $G/t0.txt $G/restart_tandem.out $G/aborted_default_mode/ 2>/dev/null
[ -d $G/checkpoint ] && mv $G/checkpoint $G/aborted_default_mode/
mkdir -p $G/out
grep -q '^\[gf_checkpoint\]' $G/bp3.toml || printf '\n# Added for the QDGreen restart (2026-10-27): Green'"'"'s-function files inside this run folder.\n[gf_checkpoint]\nprefix = "gf_checkpoint"\n' >> $G/bp3.toml
ulimit -Ss 65520 || { echo "ERROR: ulimit -Ss 65520 failed"; exit 1; }
cd $G
CMD="mpirun -np 4 $TD/build_2d_p1/app/tandem bp3.toml --mode QDGreen --petsc -options_file $TD/examples/options/lu_mumps.cfg -options_file $TD/examples/options/rk45.cfg -ts_monitor"
{ echo "start $(date)"
  echo "history: 1st start 22:22:10 PDT (stack 7 MiB) -> aborted_stack7MiB/; 2nd start 22:24:03 PDT (QD, stack 63 MiB, about 0.5 s/step, projected 30-40 h) -> aborted_default_mode/; this run: --mode QDGreen with [gf_checkpoint] prefix = gf_checkpoint, otherwise unchanged"
  echo "tandem commit $(git -C $TD rev-parse HEAD) (status: $(git -C $TD status --short | grep -v '^??' | tr '\n' ' '))"
  echo "PETSc env $CONDA_PREFIX ($(grep -E 'define PETSC_VERSION_(MAJOR|MINOR|SUBMINOR) ' $CONDA_PREFIX/include/petscversion.h | awk '{printf $3"."}'))"
  echo "ulimit -s $(ulimit -s)"; echo "cmd $CMD"; } > run_info.txt
t0=$(date +%s); echo $t0 > t0.txt
$CMD > run.log 2>&1; rc=$?
echo "end $(date) rc $rc wall $(( $(date +%s) - t0 )) s" >> run_info.txt
echo "rc $rc wall $(( $(date +%s) - t0 )) s"
