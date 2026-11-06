#!/bin/bash
#SBATCH --job-name=fo_build
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=1 --ntasks-per-node=20 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=04:00:00
#SBATCH --output=fo_build_%j.out
#SBATCH --error=fo_build_%j.err
# 2D in-plane cycle build on Discovery (USC CARC, SCEC condo), the condo twin
# of jobs/2d_inplane/vista_build.sbatch (Chunhui Zhao, 2026-11-05: "also
# let's try discovery, which you should have access").  On the login node
# first (it has internet; the compute nodes may not):
#   mkdir -p /scratch1/$USER/seas-mfem-fully-dynamic-fault-opening && cd $_
#   git clone --recurse-submodules -b fully-dynamic-fault-opening \
#       https://github.com/chunhuizhao478/seas-mfem.git seas-mfem
#   cd seas-mfem && STAGE=fetch DEPS_ROOT=$PWD/../deps \
#       bash miniapps/seas/jobs/fd_cycle/build_discovery.sh
#   cd miniapps/seas/jobs/2d_inplane && sbatch discovery_2d_build.sbatch
# Steps (logs in <checkout>-runs/build_<jobid>, next to the checkout):
#   1. jobs/fd_cycle/build_discovery.sh, STAGE="deps mfem": PETSc 3.24.6
#      with MUMPS into DEPS_ROOT, MFEM in place (the GitHub-first guard runs
#      inside, before anything is built).  Its seas stage is not used:
#      cycle-all ends with fully-dynamic-cycle's rule M1 check
#      (cycle-owned-check), which refuses this branch's 2D files by design.
#   2. seas_fd_cycle_driver from fresh objects, then BUILD_COMMIT_DISCOVERY.
#   3. The Vista twin's 2D test targets (test-2d-deck, test-prestress-bp3,
#      test-2d-absorbing, test-absorbing-all, test-sigma-floor,
#      test-relaxation, test-tip-zone) plus this branch's fault-opening
#      suites (test-fault-opening, test-fault-contact, test-fault-contact-dyn,
#      test-fault-contact-qd).  test-handback-table and the snapshot
#      read-back need an h5py python, which Discovery has not staged: those
#      run in the cloud and in the Vista build job instead.
# Exits non-zero when any step fails, so discovery_2d_cycle.sbatch can wait
# on it with --dependency=afterok.  Cost: 1 node x <= 4 h on the SCEC condo
# (scec_608).
# Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export MODULES="gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 hdf5/1.14.3 cmake/3.29.4 python/3.11.9"
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
DEPS_ROOT=${DEPS_ROOT:-$(dirname "$ROOT")/deps}
RUNS=$ROOT-runs/build_$SLURM_JOB_ID
mkdir -p "$RUNS" || exit 2
echo "host $(hostname)  cpu $(lscpu | grep -m1 'Model name' | xargs)  logs $RUNS"
T0=$(date +%s)

# 1. PETSc and MFEM (the GitHub-first guard runs inside)
JOBS=20 STAGE="deps mfem" DEPS_ROOT="$DEPS_ROOT" \
   bash "$ROOT/miniapps/seas/jobs/fd_cycle/build_discovery.sh" > "$RUNS/build_deps_mfem.log" 2>&1
rc=$?; echo "build_discovery.sh deps mfem exit $rc after $(( $(date +%s) - T0 )) s"; tail -5 "$RUNS/build_deps_mfem.log"
[ $rc -eq 0 ] || exit $rc

# 2. The driver, from fresh objects (the shared Makefile has no header
# dependencies), with the toolchain of step 1
for m in $MODULES; do module load $m || exit 2; done
[ -n "${HDF5_ROOT:-}" ] && [ -f "$HDF5_ROOT/include/hdf5.h" ] \
   || { echo "ERROR: HDF5_ROOT not set by the hdf5 module"; exit 2; }
cd "$ROOT/miniapps/seas" || exit 2
case "$(git -C "$ROOT" submodule status miniapps/seas/extern/toml11)" in
   " "*) ;;
   *) echo "ERROR: toml11 submodule missing or not at its recorded commit"; exit 2 ;;
esac
find . -name '*.o' -not -path './extern/*' -delete
MK="make -f tests/2d_inplane/2d_inplane.mk SEAS_FS_LIB= HDF5_DIR=$HDF5_ROOT"
s=$(date +%s)
$MK -j20 seas_fd_cycle_driver > "$RUNS/build_seas.log" 2>&1; rc=$?
echo "seas_fd_cycle_driver build exit $rc in $(( $(date +%s) - s )) s, $(grep -c warning "$RUNS/build_seas.log") warning lines"
[ $rc -eq 0 ] || { tail -30 "$RUNS/build_seas.log"; exit $rc; }
git -C "$ROOT" rev-parse HEAD > BUILD_COMMIT_DISCOVERY
echo "commit $(cat BUILD_COMMIT_DISCOVERY)"
fails=0

# 3. The 2D test targets of the Vista twin, and this branch's fault-opening
# suites
for t in test-2d-deck test-prestress-bp3 test-2d-absorbing test-absorbing-all \
         test-sigma-floor test-relaxation test-tip-zone \
         test-fault-opening test-fault-contact test-fault-contact-dyn \
         test-fault-contact-qd; do
   s=$(date +%s)
   $MK $t > "$RUNS/$t.log" 2>&1; rc=$?
   echo "$t exit $rc in $(( $(date +%s) - s )) s"
   grep -E "passed ===|[0-9]+ / [0-9]+ passed|: PASSED \(|: FAILED|^  FAIL|  FAILED" "$RUNS/$t.log" | tail -8
   [ $rc -eq 0 ] || fails=$((fails + 1))
done

echo "total $(( $(date +%s) - T0 )) s, failing steps: $fails"
[ $fails -eq 0 ] || exit 1
