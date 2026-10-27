#!/bin/bash
# G3 step B: Tandem p = 1 BP3-QD d60 thrust (bp3_d60_reverse), 1500 yr, 4 ranks, on
# seas-mfem's bp3_d60_25m.msh, the example's options (lu_mumps.cfg + rk45.cfg).
# Holds the Mac run lock for the whole run.
source /Users/chunhuizhao/miniforge/etc/profile.d/conda.sh; conda activate tandem-petsc321
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ulimit -Ss 65536
G=/Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/g3-bp3/tandem_p1_d60_thrust
TD=/Users/chunhuizhao/projects/tandem
LOCK=/Users/chunhuizhao/projects/.seas-mac-run.lock
have=0
release() { [ "$have" = 1 ] && grep -q "pid $$ " "$LOCK/owner" 2>/dev/null && rm -f "$LOCK/owner" && rmdir "$LOCK" && echo "lock released $(date)"; }
trap release EXIT
until mkdir "$LOCK" 2>/dev/null; do echo "$(date) waiting for lock: $(cat $LOCK/owner 2>/dev/null)"; sleep 60; done
have=1; echo "pid $$ thread 2d-inplane (G3 Tandem p1 BP3 d60 thrust 1500 yr) since $(date)" > "$LOCK/owner"
busy() { ps -axo pid=,comm= | awk '{n=split($2,p,"/"); b=p[n]} b ~ /^seas_/ || b ~ /bp5/ || b=="prterun" || b=="mpirun" || b=="mpiexec" || b=="tandem" {print}'; }
while [ -n "$(busy)" ]; do echo "$(date) waiting for: $(busy)"; sleep 60; done
cd $G
CMD="mpirun -np 4 $TD/build_2d_p1/app/tandem bp3.toml --petsc -options_file $TD/examples/options/lu_mumps.cfg -options_file $TD/examples/options/rk45.cfg -ts_monitor"
{ echo "start $(date)"; echo "tandem commit $(git -C $TD rev-parse HEAD) (status: $(git -C $TD status --short | grep -v '^??' | tr '\n' ' '))"
  echo "PETSc env $CONDA_PREFIX ($(grep -E 'define PETSC_VERSION_(MAJOR|MINOR|SUBMINOR) ' $CONDA_PREFIX/include/petscversion.h | awk '{printf $3"."}'))"
  echo "ulimit -s $(ulimit -s)"; echo "cmd $CMD"; } > run_info.txt
t0=$(date +%s); echo $t0 > t0.txt
$CMD > run.log 2>&1; rc=$?
echo "end $(date) rc $rc wall $(( $(date +%s) - t0 )) s" >> run_info.txt
echo "rc $rc wall $(( $(date +%s) - t0 )) s"
