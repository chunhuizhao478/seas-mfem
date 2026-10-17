#!/bin/bash
#SBATCH --job-name=fdc_build
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=1 --ntasks-per-node=20 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=04:00:00
#SBATCH --output=fdc_build_%j.out
#SBATCH --error=fdc_build_%j.err
# Fully dynamic BP5 cycle on Discovery: build and tests on an SCEC condo node
# (Chunhui Zhao, 2026-10-16).  On the login node first:
#   mkdir -p /scratch1/$USER/fdcycle-discovery-20261016 && cd $_
#   git clone --recurse-submodules -b fully-dynamic-cycle \
#       https://github.com/chunhuizhao478/seas-mfem.git seas-mfem
#   cd seas-mfem && STAGE=fetch bash miniapps/seas/jobs/fd_cycle/build_discovery.sh
#   cd miniapps/seas/jobs/fd_cycle && sbatch discovery_build.sbatch
# Steps (logs in ../runs/build_<jobid>, next to the checkout):
#   1. build_discovery.sh deps mfem seas: PETSc 3.24.6 with MUMPS, MFEM,
#      cycle-all (driver, tests, 209 API pins, owned-paths check M1)
#   2. test-cycle-p3 on 1 and 2 ranks (bit for bit on 1 rank)
#   3. test-cycle-p2 (frame_map on 1, 2, 4, 8 ranks; the rest on 2)
#   4. test_cycle_restart.sh: the P1 gold deck, 25 + 25 steps vs 50, 1 rank
# Exits non-zero when any step fails, so discovery_switch.sbatch can wait on
# it with --dependency=afterok.  Cost: 1 node x <= 4 h on the SCEC condo.
# Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MODULES="gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 hdf5/1.14.3 cmake/3.29.4 python/3.11.9"
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
RUNS=$(dirname "$ROOT")/runs/build_$SLURM_JOB_ID
mkdir -p "$RUNS" || exit 2
echo "host $(hostname)  cpu $(grep -m1 'model name' /proc/cpuinfo)  logs $RUNS"
T0=$(date +%s)

# 1. build (the GitHub-first guard runs inside, before anything is built)
JOBS=20 STAGE="deps mfem seas" bash "$ROOT/miniapps/seas/jobs/fd_cycle/build_discovery.sh" \
   > "$RUNS/build.log" 2>&1
rc=$?; echo "build_discovery.sh exit $rc after $(( $(date +%s) - T0 )) s"; tail -25 "$RUNS/build.log"
[ $rc -eq 0 ] || exit $rc
echo "commit $(git -C "$ROOT" rev-parse HEAD)"

# The same modules for the tests.
module purge
for m in $MODULES; do module load $m || exit 2; done
cd "$ROOT/miniapps/seas"
MK="make -f cycle/cycle.mk SEAS_FS_LIB= HDF5_DIR=$HDF5_ROOT"
fails=0

# 2. P3 tests
s=$(date +%s)
$MK test-cycle-p3 CYCLE_P3_WORK="$RUNS/p3_work" > "$RUNS/test_p3.log" 2>&1; rc=$?
echo "test-cycle-p3 exit $rc in $(( $(date +%s) - s )) s"
grep -E "^  (ok  |FAIL)|control:|^  A:|test_cycle_p3: (PASSED|FAILED)" "$RUNS/test_p3.log"
[ $rc -eq 0 ] || fails=$((fails + 1))

# 3. P2 tests
s=$(date +%s)
$MK test-cycle-p2 > "$RUNS/test_p2.log" 2>&1; rc=$?
echo "test-cycle-p2 exit $rc in $(( $(date +%s) - s )) s"
grep -E "FAIL|PASSED|FAILED" "$RUNS/test_p2.log" | tail -20
[ $rc -eq 0 ] || fails=$((fails + 1))

# 4. P1 restart, bit for bit on 1 rank
s=$(date +%s)
MPIRUN=mpirun tests/cycle/test_cycle_restart.sh config/fd_cycle/bp5_gold.toml "$RUNS/restart" 25 1 \
   > "$RUNS/test_restart.log" 2>&1; rc=$?
echo "test_cycle_restart.sh exit $rc in $(( $(date +%s) - s )) s"
grep -E "PASSED|FAILED|differ" "$RUNS/test_restart.log" | tail -10
[ $rc -eq 0 ] || fails=$((fails + 1))

echo "total $(( $(date +%s) - T0 )) s, failing steps: $fails"
[ $fails -eq 0 ] || exit 1
