#!/bin/bash
#SBATCH --job-name=qdtiming
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=1 --ntasks-per-node=20 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=04:00:00
#SBATCH --output=qdtiming_%j.out
#SBATCH --error=qdtiming_%j.err
# Timing of the implicit QD step's speed-up keys on Discovery (USC CARC,
# SCEC condo), the condo twin of jobs/2d_inplane/vista_dense_lu_timing.sbatch
# (Chunhui Zhao, 2026-11-10: "I thought we will upload it and run on
# Discovery"; "let's proceed with the speedup options"):
#   1. seas_test_par_dense_lu (built here, one file): the checks L1, L2 on
#      20 ranks, then the factor and solve times at M = 692 (sjf1) and
#      M = 3,290 (UCERF) on 20, 16, 9 and 4 ranks.
#   2. The UCERF implicit deck for QD_STEPS (default 300) QD steps from
#      t = 0, three times: as it is ("serial"), with [cycle.qd] dense_lu =
#      "shared" and with [cycle.qd] stage_guess = "extrapolate".  The first
#      run writes the Green's function file and the others read it, so all
#      three step on the same G.  The driver's esdirk line gives Newton
#      iterations, LU factors and seconds per step; PETSc's -log_view is not
#      on (the decks' options file has none).  Output in
#      <checkout>-runs/qdtiming_<job id>/{serial,shared,guess}; nothing is
#      kept as a result.
# One node, 20 ranks (Vista ran the deck on 144): the serial factor costs
# the same on any rank count, the shared one has 20 cores to use.
# The binaries must be built (discovery_2d_build.sbatch) at this code: the
# guard of discovery_2d_cycle.sbatch.  Cost: 1 node x <= 4 h on the SCEC
# condo (scec_608).  Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
MODULES="gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 hdf5/1.14.3 cmake/3.29.4 python/3.11.9"
STEPS=${QD_STEPS:-300}
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd "$ROOT" || exit 2

# ---- GitHub-first guard (CLAUDE.md), and the binaries match this code ----
# No fetch: run as checked out (the login-node submit fetched).
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
HEAD=$(git rev-parse HEAD)
BUILT=$(cat miniapps/seas/BUILD_COMMIT_DISCOVERY 2>/dev/null)
git cat-file -e "$BUILT^{commit}" 2>/dev/null \
   || { echo "ERROR: binaries built at '$BUILT', not a commit here; run discovery_2d_build.sbatch"; exit 2; }
CODE=$(git diff --name-only "$BUILT" HEAD \
       | grep -v -E '^miniapps/seas/(config/fd_cycle|jobs/fd_cycle|jobs/2d_inplane|simple_cases/configs|simple_cases/scripts)/')
[ -z "$CODE" ] || { echo "ERROR: code changed since the build at $BUILT:"; echo "$CODE"; exit 2; }
echo "commit $HEAD, binaries built at $BUILT"

module purge
for m in $MODULES; do module load $m || exit 2; done
[ -n "${HDF5_ROOT:-}" ] && [ -f "$HDF5_ROOT/include/hdf5.h" ] \
   || { echo "ERROR: HDF5_ROOT not set by the hdf5 module"; exit 2; }
cd miniapps/seas || exit 2
OUT=$ROOT-runs/qdtiming_$SLURM_JOB_ID
mkdir -p "$OUT/serial" "$OUT/shared" "$OUT/guess" || exit 2
echo "host $(hostname)  cpu $(lscpu | grep -m1 'Model name' | xargs)  ranks $SLURM_NTASKS  output $OUT"
RUN="mpirun --bind-to core"
fails=0

# ---- 1. the LU alone ----
make -f tests/2d_inplane/2d_inplane.mk SEAS_FS_LIB= HDF5_DIR=$HDF5_ROOT seas_test_par_dense_lu \
   > "$OUT/build_test.log" 2>&1 || { echo "ERROR: seas_test_par_dense_lu did not build"; tail -20 "$OUT/build_test.log"; exit 2; }
$RUN -np "$SLURM_NTASKS" ./seas_test_par_dense_lu -n 692 -reps 5 > "$OUT/lu_checks.log" 2>&1; rc=$?
echo "seas_test_par_dense_lu (checks, $SLURM_NTASKS ranks) exit $rc"; cat "$OUT/lu_checks.log"
[ $rc -eq 0 ] || fails=$((fails + 1))
for n in 20 16 9 4; do
   [ "$n" -le "$SLURM_NTASKS" ] || continue
   for M in 692 3290; do
      echo "--- $n ranks, M = $M"
      $RUN -np $n ./seas_test_par_dense_lu -only-time -n $M -reps 5 2>&1 | grep -A4 "L3 process grid"
   done
done

# ---- 2. the UCERF implicit deck, QD steps only: serial, shared, guess ----
B=config/fd_cycle/ucerf2d_saf_sjf_500m_R70_avg_tip_contact_pml0p01_quench_retry_esdirk
G=$OUT/green_ucerf.bin
for k in serial shared guess; do
   case $k in
      serial) D=${B}_fd.toml ;;
      shared) D=${B}_sharedlu_fd.toml ;;
      guess)  D=${B}_sg_fd.toml ;;
   esac
   s=$(date +%s)
   $RUN -np "$SLURM_NTASKS" ./seas_fd_cycle_driver --deck "$D" --output-dir "$OUT/$k" \
      --max-steps "$STEPS" --log-every 50 --qd-green-cache "$G" > "$OUT/$k/run.log" 2>&1; rc=$?
   echo "--- UCERF esdirk, $k: exit $rc after $(( $(date +%s) - s )) s"
   grep -E "commit:|integrator esdirk|MUMPS: analysis|QD-GREEN\]|esdirk in the QD regime|esdirk \(option I\)|final time|TS steps|step rejections|status" "$OUT/$k/run.log"
   [ $rc -eq 0 ] || { fails=$((fails + 1)); tail -20 "$OUT/$k/run.log"; }
done
# serial and shared take the same steps; guess may differ within newton_tol:
# the last printed step of each
for k in serial shared guess; do echo "$k: $(grep -E "^ +$STEPS " "$OUT/$k/run.log" | tail -1)"; done
rm -f "$G"
du -sh "$OUT"
echo "qd timing: $fails failed step(s)"
exit $fails
