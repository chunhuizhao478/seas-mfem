#!/bin/bash
#SBATCH -J denselu
#SBATCH -o denselu_%j.out
#SBATCH -e denselu_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 01:30:00
#SBATCH -A EAR20006
# Timing of the shared dense LU of the implicit QD fault solve on Vista
# (Chunhui Zhao, 2026-11-10: "I would like to emphasize profiling and
# speedup for the current workflow"): [cycle.qd] dense_lu = "shared"
# (cycle/par_dense_lu.hpp, ScaLAPACK) against "serial" (every rank factors).
#   1. seas_test_par_dense_lu (built here, one file): the checks L1, L2 on
#      144 ranks, then the factor and solve times at M = 692 (sjf1) and
#      M = 3,290 (UCERF) on 144, 64, 36 and 16 ranks.
#   2. The UCERF implicit deck for 300 QD steps from t = 0, serial and
#      shared (the same steps; the driver's esdirk line gives Newton
#      iterations, LU factors and seconds per step).  Output in
#      <checkout>-runs/denselu_<job id>/{serial,shared}; nothing is kept as
#      a result.
# The binaries must be built (vista_build.sbatch) at this code: the guard
# of vista_2d_cycle.sbatch.  Cost: 1 node x <= 1.5 h on gg = <= 0.5 SU
# (EAR20006).  Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
DEPS_ROOT=${DEPS_ROOT:-$ROOT-deps}
PETSC_PREFIX=$DEPS_ROOT/petsc-${PETSC_VERSION:-3.24.6}
cd "$ROOT" || exit 2

# ---- GitHub-first guard (CLAUDE.md), and the binaries match this code ----
# No fetch: the compute nodes may have no internet; fetch on the login node.
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
HEAD=$(git rev-parse HEAD)
BUILT=$(cat miniapps/seas/BUILD_COMMIT_VISTA 2>/dev/null)
git cat-file -e "$BUILT^{commit}" 2>/dev/null \
   || { echo "ERROR: binaries built at '$BUILT', not a commit here; run vista_build.sbatch"; exit 2; }
CODE=$(git diff --name-only "$BUILT" HEAD \
       | grep -v -E '^miniapps/seas/(config/fd_cycle|jobs/fd_cycle|jobs/2d_inplane|simple_cases/configs|simple_cases/scripts)/')
[ -z "$CODE" ] || { echo "ERROR: code changed since the build at $BUILT:"; echo "$CODE"; exit 2; }
echo "commit $HEAD, binaries built at $BUILT"

module purge
for m in $MODULES; do module load $m || exit 2; done
cd miniapps/seas || exit 2
OUT=$ROOT-runs/denselu_$SLURM_JOB_ID
mkdir -p "$OUT/serial" "$OUT/shared" || exit 2
echo "host $(hostname)  ranks $SLURM_NTASKS  output $OUT"
fails=0

# ---- 1. the LU alone ----
make -f tests/2d_inplane/2d_inplane.mk SEAS_FS_LIB= HDF5_DIR=$PETSC_PREFIX seas_test_par_dense_lu \
   > "$OUT/build_test.log" 2>&1 || { echo "ERROR: seas_test_par_dense_lu did not build"; tail -20 "$OUT/build_test.log"; exit 2; }
ibrun ./seas_test_par_dense_lu -n 692 -reps 5 > "$OUT/lu_144.log" 2>&1; rc=$?
echo "seas_test_par_dense_lu (checks, 144 ranks) exit $rc"; cat "$OUT/lu_144.log"
[ $rc -eq 0 ] || fails=$((fails + 1))
for n in 144 64 36 16; do
   for M in 692 3290; do
      echo "--- $n ranks, M = $M"
      ibrun -n $n ./seas_test_par_dense_lu -only-time -n $M -reps 5 2>&1 | grep -A4 "L3 process grid"
   done
done

# ---- 2. the UCERF implicit deck, 300 QD steps, serial and shared ----
B=config/fd_cycle/ucerf2d_saf_sjf_500m_R70_avg_tip_contact_pml0p01_quench_retry_esdirk
for k in serial shared; do
   D=${B}_fd.toml; [ $k = shared ] && D=${B}_sharedlu_fd.toml
   s=$(date +%s)
   ibrun ./seas_fd_cycle_driver --deck "$D" --output-dir "$OUT/$k" --max-steps 300 --log-every 50 \
      > "$OUT/$k/run.log" 2>&1; rc=$?
   echo "--- UCERF esdirk, dense_lu $k: exit $rc after $(( $(date +%s) - s )) s"
   grep -E "commit:|integrator esdirk|MUMPS: analysis|QD-GREEN\] built|esdirk in the QD regime|esdirk \(option I\)|final time|TS steps|step rejections|status" "$OUT/$k/run.log"
   [ $rc -eq 0 ] || { fails=$((fails + 1)); tail -20 "$OUT/$k/run.log"; }
done
# the two runs take the same steps: compare the last accepted step
for k in serial shared; do echo "$k: $(grep -E '^ +300 ' "$OUT/$k/run.log" | tail -1)"; done
du -sh "$OUT"
echo "denselu timing: $fails failed step(s)"
exit $fails
