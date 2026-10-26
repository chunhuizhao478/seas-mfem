#!/bin/bash
#SBATCH -J qd_rock4_test
#SBATCH -o qd_rock4_test_%j.out
#SBATCH -e qd_rock4_test_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# QD speed S3: build tests/parallel/test_ts_rock4.cpp (the PETSc TS type
# rock4 of solver/ts_rock4.hpp; no mesh) in this clone and run it at 1, 2, 4
# and 8 ranks.  Each run must end with ALL PASSED, and the lines of its
# per-method check (y(1) of the nonlinear problem for every method of 5 to
# 124 stages, printed with 17 digits) must be the same at every rank count.
# Output in bp5/results_qd_speed/rock4_test_job<id>.  Needs the clone's
# MFEM build (jobs/vista/cpu_clone_build.sbatch); qd_speed_submit.sh
# STEPS=rock4 submits it after the build.
#
# GitHub is the version control: the checkout must be clean, its HEAD on
# GitHub (CLAUDE.md guard) and at the commit its upstream branch points to,
# or the job stops.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"

cd "${SEAS_MFEM_ROOT}"
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in ${SEAS_MFEM_ROOT}"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "commit $(git rev-parse HEAD)"
if [ "$(git rev-parse HEAD)" != "$(git rev-parse '@{u}' 2>/dev/null)" ]; then
   echo "ERROR: HEAD $(git rev-parse --short HEAD) is not its upstream" \
        "($(git rev-parse --short '@{u}' 2>/dev/null || echo none)); pull first." >&2
   exit 1
fi
echo "Checkout: ${SEAS_MFEM_ROOT} at $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD))"

cd "${SEAS_DIR}"
OUT="bp5/results_qd_speed/rock4_test_job${SLURM_JOB_ID}"
mkdir -p "${OUT}"
g0_provenance "${OUT}"

# A fresh object, so the test is compiled at this commit.
rm -f tests/parallel/test_ts_rock4.o seas_test_ts_rock4
if ! make seas_test_ts_rock4 > "${OUT}/build.log" 2>&1; then
   echo "ERROR: building seas_test_ts_rock4 failed:" >&2
   tail -n 30 "${OUT}/build.log" >&2
   exit 1
fi

fail=0
for np in 1 2 4 8; do
   ibrun -n "${np}" -o 0 ./seas_test_ts_rock4 > "${OUT}/test_np${np}.txt" 2>&1
   rc=$?
   echo "np ${np}: exit ${rc}, $(tail -n 1 "${OUT}/test_np${np}.txt")"
   if [ "${rc}" != 0 ] || ! grep -qx "ALL PASSED" "${OUT}/test_np${np}.txt"; then
      fail=1
   fi
done
grep '^  nonlinear ' "${OUT}/test_np1.txt" > "${OUT}/nonlinear_np1.txt"
if [ "$(wc -l < "${OUT}/nonlinear_np1.txt")" != 47 ]; then
   echo "ERROR: expected 47 nonlinear lines at 1 rank." >&2
   fail=1
fi
for np in 2 4 8; do
   if ! grep '^  nonlinear ' "${OUT}/test_np${np}.txt" \
        | cmp -s - "${OUT}/nonlinear_np1.txt"; then
      echo "ERROR: the nonlinear lines at ${np} ranks differ from 1 rank's." >&2
      fail=1
   fi
done
echo "ROCK4 test: $([ "${fail}" = 0 ] && echo PASS || echo FAIL)"
exit "${fail}"
