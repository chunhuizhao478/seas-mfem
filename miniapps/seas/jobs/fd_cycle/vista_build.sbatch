#!/bin/bash
#SBATCH -J fdc_build
#SBATCH -o fdc_build_%j.out
#SBATCH -e fdc_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 06:00:00
#SBATCH -A EAR20006
# Fully dynamic BP5 cycle on Vista: build and every cycle test on one gg node
# (Chunhui Zhao, 2026-10-21).  discovery_build.sbatch's steps, plus the
# parameter parity test and the trunk's T12 and hook spy tests, which the
# trunk merge (7baf0a4) was checked with.  On the login node first:
#   cd $WORK
#   git clone --recurse-submodules -b fully-dynamic-cycle \
#       https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-fdcycle
#   cd seas-mfem-fdcycle && STAGE=fetch bash miniapps/seas/jobs/fd_cycle/build_vista.sh
#   cd miniapps/seas/jobs/fd_cycle && sbatch vista_build.sbatch
# Steps (logs in <checkout>-runs/build_<jobid>, next to the checkout):
#   0. python3 with numpy, h5py and tomllib (CYCLE_PY), checked before the
#      build: the snapshot test reads its files back with it
#   1. build_vista.sh deps mfem seas: PETSc 3.24.6 with MUMPS and HDF5,
#      MFEM, cycle-all (driver, tests, 253 API pins, owned-paths check M1)
#   2. test-cycle-params-parity: the decks against BP5Params, 1 rank
#   3. test-cycle-p3 on 1 and 2 ranks (bit for bit on 1 rank)
#   4. test-cycle-p2 (frame_map on 1, 2, 4, 8 ranks; the rest on 2)
#   5. test_cycle_restart.sh: the P1 gold deck, 25 + 25 steps vs 50, 1 rank
#   6. test-cycle-p4 on 1 and 2 ranks (the QD -> dynamic hand-over, re-entry)
#   7. test-cycle-snapshots on 2 ranks: the fault and surface snapshot files
#      through a switch, restarts, a trim and a wall-limit stop, read back
#      (a run that cannot read them back fails this step)
#   8. trunk test T12 (test-nitsche-pf) on 1 and 2 ranks
#   9. trunk test L-10 (test-nitsche-hook-spy, hook H-B1) on 1, 2 and 4 ranks
#  10. test-cycle-bench on 1 and 2 ranks: the SCEC outputs of
#      cycle/bench_output.hpp on a two-face fault and a box (unit test)
#  11. test-cycle-bench-run on 4 ranks: the p = 2 deck's options in short
#      p = 1 runs (restart with the SCEC files, the QD velocity, the stage
#      operators, LTS against global stepping), MUMPS
# FDC_PART splits the job (submit_vista_p2.sh): "build" runs steps 0 and 1
# only, "tests" steps 0 and 2 to 11 on the binaries a build job made at
# this commit (it refuses to start if make would rebuild them, since
# vista_p2_check.sbatch may be running them), "all" (the default) both.
# Exits non-zero when any step fails.  Cost: 1 node x <= 6 h (<= 2 SU).
# Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
export MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
DEPS_ROOT=${DEPS_ROOT:-$ROOT-deps}
PETSC_PREFIX=$DEPS_ROOT/petsc-${PETSC_VERSION:-3.24.6}
PART=${FDC_PART:-all}
case "$PART" in all|build|tests) ;; *) echo "ERROR: FDC_PART '$PART' (all, build or tests)"; exit 2 ;; esac
RUNS=$ROOT-runs/build_$SLURM_JOB_ID
mkdir -p "$RUNS" || exit 2
echo "part: $PART"
echo "host $(hostname)  cpu $(lscpu | grep -m1 'Model name' | xargs)  logs $RUNS"
T0=$(date +%s)

# 0. The tests' python (Vista's default python3 has no numpy; the
# python3/3.11.8 module's interpreter, called by path as the QD jobs do).
# The build keeps the default python3 for PETSc's configure, as the QD build.
CYCLE_PY=${CYCLE_PY:-/opt/apps/gcc14/cuda12/python3/3.11.8/bin/python3}
"$CYCLE_PY" -c "import numpy, h5py, tomllib; print('python', __import__('sys').version.split()[0], 'numpy', numpy.__version__, 'h5py', h5py.__version__)" \
   || { echo "ERROR: $CYCLE_PY cannot import numpy, h5py and tomllib; set CYCLE_PY"; exit 2; }
mkdir -p "$RUNS/pybin" && ln -sf "$CYCLE_PY" "$RUNS/pybin/python3" || exit 2

# 1. build (the GitHub-first guard runs inside, before anything is built)
if [ "$PART" != tests ]; then
   JOBS=64 STAGE="deps mfem seas" bash "$ROOT/miniapps/seas/jobs/fd_cycle/build_vista.sh" \
      > "$RUNS/build.log" 2>&1
   rc=$?; echo "build_vista.sh exit $rc after $(( $(date +%s) - T0 )) s"; tail -25 "$RUNS/build.log"
   [ $rc -eq 0 ] || exit $rc
   echo "commit $(git -C "$ROOT" rev-parse HEAD)"
   [ "$PART" = build ] && { echo "build only (FDC_PART=build): the tests run in their own job"; exit 0; }
else
   # The tests only: the GitHub-first guard, and binaries built at HEAD
   [ -z "$(git -C "$ROOT" status --porcelain --untracked-files=no)" ] \
      || { echo "ERROR: uncommitted changes in $ROOT"; git -C "$ROOT" status --short; exit 2; }
   [ -n "$(git -C "$ROOT" branch -r --contains HEAD)" ] \
      || { echo "ERROR: HEAD $(git -C "$ROOT" rev-parse --short HEAD) is not on GitHub"; exit 2; }
   BUILT=$(cat "$ROOT/miniapps/seas/BUILD_COMMIT_VISTA" 2>/dev/null)
   [ "$BUILT" = "$(git -C "$ROOT" rev-parse HEAD)" ] \
      || { echo "ERROR: binaries built at '$BUILT', HEAD is $(git -C "$ROOT" rev-parse HEAD); run FDC_PART=build first"; exit 2; }
   echo "commit $BUILT, binaries from the build job"
fi

# The same modules for the tests, and CYCLE_PY as their python3.
for m in $MODULES; do module load $m || exit 2; done
export PATH="$RUNS/pybin:$PATH"
cd "$ROOT/miniapps/seas"
MK="make -f cycle/cycle.mk SEAS_FS_LIB= HDF5_DIR=$PETSC_PREFIX"
if [ "$PART" = tests ]; then
   $MK -q seas_fd_cycle_driver seas_test_cycle_params_parity seas_test_cycle_p2 \
          seas_test_cycle_p3 seas_test_cycle_p4 seas_test_cycle_bench \
      || { echo "ERROR: make would rebuild the build job's binaries"; exit 2; }
fi
fails=0

# 2. Parameter parity
s=$(date +%s)
$MK test-cycle-params-parity > "$RUNS/test_parity.log" 2>&1; rc=$?
echo "test-cycle-params-parity exit $rc in $(( $(date +%s) - s )) s"
grep -E "PASSED|FAIL" "$RUNS/test_parity.log" | tail -10
[ $rc -eq 0 ] || fails=$((fails + 1))

# 3. P3 tests
s=$(date +%s)
$MK test-cycle-p3 CYCLE_P3_WORK="$RUNS/p3_work" > "$RUNS/test_p3.log" 2>&1; rc=$?
echo "test-cycle-p3 exit $rc in $(( $(date +%s) - s )) s"
grep -E "^  (ok  |FAIL)|control:|^  A:|test_cycle_p3: (PASSED|FAILED)" "$RUNS/test_p3.log"
[ $rc -eq 0 ] || fails=$((fails + 1))

# 4. P2 tests
s=$(date +%s)
$MK test-cycle-p2 > "$RUNS/test_p2.log" 2>&1; rc=$?
echo "test-cycle-p2 exit $rc in $(( $(date +%s) - s )) s"
grep -E "FAIL|PASSED|FAILED" "$RUNS/test_p2.log" | tail -20
[ $rc -eq 0 ] || fails=$((fails + 1))

# 5. P1 restart, bit for bit on 1 rank
s=$(date +%s)
MPIRUN=mpirun tests/cycle/test_cycle_restart.sh config/fd_cycle/bp5_gold.toml "$RUNS/restart" 25 1 \
   > "$RUNS/test_restart.log" 2>&1; rc=$?
echo "test_cycle_restart.sh exit $rc in $(( $(date +%s) - s )) s"
grep -E "PASSED|FAILED|differ" "$RUNS/test_restart.log" | tail -10
[ $rc -eq 0 ] || fails=$((fails + 1))

# 6. P4 tests
s=$(date +%s)
$MK test-cycle-p4 > "$RUNS/test_p4.log" 2>&1; rc=$?
echo "test-cycle-p4 exit $rc in $(( $(date +%s) - s )) s"
grep -E "^  (ok  |FAIL)|test_cycle_p4: (PASSED|FAILED)" "$RUNS/test_p4.log"
[ $rc -eq 0 ] || fails=$((fails + 1))

# 7. Snapshot files (MUMPS as the production decks), read back
s=$(date +%s)
SOLVER=mumps $MK test-cycle-snapshots CYCLE_SNAP_WORK="$RUNS/snapshots" \
   > "$RUNS/test_snapshots.log" 2>&1; rc=$?
echo "test-cycle-snapshots exit $rc in $(( $(date +%s) - s )) s"
grep -E "NOTE|^[A-G]: exit|wall limit:|FAIL|test_cycle_snapshots: (PASSED|FAILED)" "$RUNS/test_snapshots.log"
if grep -q "files not read back" "$RUNS/test_snapshots.log"; then
   echo "FAILED: the snapshot files were not read back (python3 without numpy or h5py)"; rc=1
fi
[ $rc -eq 0 ] || fails=$((fails + 1))

# 8. Trunk test T12: FaultPointState as the cycle's checkpoints store it
s=$(date +%s)
$MK test-nitsche-pf > "$RUNS/test_nitsche_pf.log" 2>&1; rc=$?
echo "test-nitsche-pf exit $rc in $(( $(date +%s) - s )) s"
grep -E "Summary|FAIL" "$RUNS/test_nitsche_pf.log"
[ $rc -eq 0 ] || fails=$((fails + 1))

# 9. Trunk test L-10: the hook H-B1 spy, and its abort case
s=$(date +%s)
$MK test-nitsche-hook-spy > "$RUNS/test_nitsche_hook_spy.log" 2>&1; rc=$?
echo "test-nitsche-hook-spy exit $rc in $(( $(date +%s) - s )) s"
grep -E "Summary|aborted|FAIL" "$RUNS/test_nitsche_hook_spy.log"
[ $rc -eq 0 ] || fails=$((fails + 1))

# 10. The SCEC output writer (unit test)
s=$(date +%s)
$MK test-cycle-bench CYCLE_BENCH_WORK="$RUNS/bench" > "$RUNS/test_bench.log" 2>&1; rc=$?
echo "test-cycle-bench exit $rc in $(( $(date +%s) - s )) s"
grep -E "^  (ok  |FAIL)|test_cycle_bench: (PASSED|FAILED)" "$RUNS/test_bench.log"
[ $rc -eq 0 ] || fails=$((fails + 1))

# 11. The p = 2 deck's options in short runs
s=$(date +%s)
SOLVER=mumps $MK test-cycle-bench-run CYCLE_BENCH_RUN_WORK="$RUNS/bench_run" CYCLE_BENCH_RUN_NP=4 \
   > "$RUNS/test_bench_run.log" 2>&1; rc=$?
echo "test-cycle-bench-run exit $rc in $(( $(date +%s) - s )) s"
grep -E "^(---|PASSED|FAILED)|: exit|  FAIL|test_cycle_bench_run: " "$RUNS/test_bench_run.log"
[ $rc -eq 0 ] || fails=$((fails + 1))

echo "total $(( $(date +%s) - T0 )) s, failing steps: $fails"
[ $fails -eq 0 ] || exit 1
