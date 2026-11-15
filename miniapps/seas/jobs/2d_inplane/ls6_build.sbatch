#!/bin/bash
#SBATCH -J cycle2d_build
#SBATCH -o cycle2d_build_%j.out
#SBATCH -e cycle2d_build_%j.err
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 04:00:00
# 2D map-view cycle on Lonestar6 (Chunhui Zhao, 2026-11-15: "let's try on
# lonestar6", for the UCERF 2.5D relax5k 100-event run while Vista's queue
# waits): build seas_fd_cycle_driver on one normal node (2 x 64 EPYC 7763)
# and check it.  vista_build.sbatch with the Lonestar6 toolchain
# (jobs/fd_cycle/build_ls6.sh: gcc/13.2.0, impi/21.18, sequential MKL,
# -march=znver3).  No account line: submit with sbatch -A <allocation>
# (submit_ls6_2d_cycle.sh passes CYCLE_ACCOUNT).
# On the login node first (it has internet; the compute nodes may not):
#   cd $WORK
#   git clone --recurse-submodules -b fully-dynamic-fault-opening \
#       https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-fully-dynamic-fault-opening
#   cd seas-mfem-fully-dynamic-fault-opening && STAGE=fetch bash miniapps/seas/jobs/fd_cycle/build_ls6.sh
#   CYCLE_ACCOUNT=<allocation> CYCLE_DECKS=... bash miniapps/seas/jobs/2d_inplane/submit_ls6_2d_cycle.sh
# Steps (logs in <checkout>-runs/build_<jobid>, next to the checkout):
#   1. jobs/fd_cycle/build_ls6.sh, stages deps and mfem: PETSc 3.24.6 with
#      MUMPS and HDF5 into <checkout>-deps, MFEM in place (the GitHub-first
#      guard runs inside, before anything is built)
#   2. seas_fd_cycle_driver from fresh objects, then BUILD_COMMIT_LS6
#   3. the tests of vista_build.sbatch's step 3 (test-2d-deck,
#      test-prestress-bp3, test-2d-absorbing, test-absorbing-all,
#      test-sigma-floor, test-relaxation, test-tip-zone), and
#      tests/pml/pml.mk test-pml-switch: the rest state after a switch with
#      the band, global and LTS stepping, plane strain and 2.5D (the
#      substrate load fix 27ea7321), which the 2.5D + LTS run depends on
# Not here: vista_build.sbatch's steps 0 and 4 (the snapshot files read back
# with python).  No python on Lonestar6 has h5py (python/3.12.11 has numpy
# and tomllib only; checked 2026-11-15), and the snapshot writer is the code
# Vista's build 1145585 checked at this commit.
# Exits non-zero when any step fails.  Cost: 1 node x <= 4 h on normal
# = <= 4 SU.
# Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export MODULES="gcc/13.2.0 impi/21.18 mkl/26.1 cmake/3.31.5 python/3.12.11"
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
DEPS_ROOT=${DEPS_ROOT:-$ROOT-deps}
PETSC_PREFIX=$DEPS_ROOT/petsc-${PETSC_VERSION:-3.24.6}
RUNS=$ROOT-runs/build_$SLURM_JOB_ID
mkdir -p "$RUNS" || exit 2
echo "host $(hostname)  cpu $(lscpu | grep -m1 'Model name' | xargs)  logs $RUNS"
T0=$(date +%s)

# 1. PETSc and MFEM (the GitHub-first guard runs inside)
JOBS=64 STAGE="deps mfem" DEPS_ROOT="$DEPS_ROOT" \
   bash "$ROOT/miniapps/seas/jobs/fd_cycle/build_ls6.sh" > "$RUNS/build_deps_mfem.log" 2>&1
rc=$?; echo "build_ls6.sh deps mfem exit $rc after $(( $(date +%s) - T0 )) s"; tail -5 "$RUNS/build_deps_mfem.log"
[ $rc -eq 0 ] || exit $rc

# 2. The driver, from fresh objects (the shared Makefile has no header
# dependencies), with the toolchain of step 1
for m in $MODULES; do module load $m || exit 2; done
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
cd "$ROOT/miniapps/seas" || exit 2
case "$(git -C "$ROOT" submodule status miniapps/seas/extern/toml11)" in
   " "*) ;;
   *) echo "ERROR: toml11 submodule missing or not at its recorded commit"; exit 2 ;;
esac
find . -name '*.o' -not -path './extern/*' -delete
MK="make -f tests/pml/pml.mk SEAS_FS_LIB= HDF5_DIR=$PETSC_PREFIX"
s=$(date +%s)
$MK -j64 seas_fd_cycle_driver > "$RUNS/build_seas.log" 2>&1; rc=$?
echo "seas_fd_cycle_driver build exit $rc in $(( $(date +%s) - s )) s, $(grep -c warning "$RUNS/build_seas.log") warning lines"
[ $rc -eq 0 ] || { tail -30 "$RUNS/build_seas.log"; exit $rc; }
git -C "$ROOT" rev-parse HEAD > BUILD_COMMIT_LS6
echo "commit $(cat BUILD_COMMIT_LS6)"
fails=0

# 3. The tests
for t in test-2d-deck test-prestress-bp3 test-2d-absorbing test-absorbing-all \
         test-sigma-floor test-relaxation test-tip-zone test-pml-switch; do
   s=$(date +%s)
   $MK $t > "$RUNS/$t.log" 2>&1; rc=$?
   echo "$t exit $rc in $(( $(date +%s) - s )) s"
   grep -E "passed ===|[0-9]+ / [0-9]+ passed|: PASSED \(|: FAILED|^  FAIL|  FAILED" "$RUNS/$t.log" | tail -8
   [ $rc -eq 0 ] || fails=$((fails + 1))
done

echo "total $(( $(date +%s) - T0 )) s, failing steps: $fails"
[ $fails -eq 0 ] || exit 1
