#!/bin/bash
#SBATCH -J tpv102_200m_400r_v2
#SBATCH -o tpv102_200m_400r_v2_%j.out
#SBATCH -e tpv102_200m_400r_v2_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV102 dev-queue v1-repro run with v2 fixes (v2 recommendation #3).
#
# This is the DISPOSITIVE test: repeat the failing configuration from
# results_200m_p1_3s_dev_job7664258 (the v1 symptom: rupture stops at
# the rank-3 / rank-293 partition seam) with the v1+v2 fixes in place
# and `--debug-qnorm` enabled.  If the fixes are correct, the fault-
# surface PVD at t=1.5 s should show the rupture front past |x| >= 1.5
# km — i.e., the front has crossed at least one partition seam — and
# the qnorm watch should show non-hypocenter ranks ramping up in step
# with the hypocenter rank.
#
# DO NOT submit this before the 4-rank (rec #1) and 50-rank (rec #2)
# sanity runs have passed.  Goal of rec #1/#2: catch an R-101
# consistency abort cheaply before burning 400-rank node-hours.
#
# WALL-TIME: 400 ranks × 200 m ≈ 0.1 s/step, tfinal=1.5 s is 5250
# steps → ~1.5 hr, fits in the 2 hr dev-queue window.
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch

export LC_ALL=C
export LANG=C

# R-403 fix: pin the submit dir before any `cd` so the post-run RESULT.txt
# check can locate the SBATCH -o file.  Slurm writes #SBATCH -o relative to
# the submit dir, but this script `cd`s to miniapps/seas before grepping —
# without this pin, the grep runs in the wrong directory and RESULT.txt
# reports "no [qnorm:watch] line" FAIL regardless of physics outcome.
SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
make seas_tpv102_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_200m_p1_1.5s_400r_v2_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "tfinal=1.5 s, dt=auto, order=1, 400 ranks, v1+v2 fixes, debug-qnorm"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 1.5 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv \
      --debug-qnorm

# R-306 fix: post-run machine-checkable success criterion for the
# dispositive 400-rank run.  PASS means (a) the R-101 shared-fault
# DOFData diagnostic did not trip, and (b) the final [qnorm:watch]
# line shows non-hypocenter ranks with non-trivial ||Q||_inf — i.e.,
# bulk wave energy actually crossed at least one partition seam.
# Without this check, `ibrun` exits 0 regardless of physics outcome
# and a regressed run hides behind a green CI signal.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
# R-403 fix: absolute path anchored at the submit dir (see SBATCH_LOG_DIR
# at top of script).  The prior relative path failed whenever the user
# submitted from anywhere other than miniapps/seas/.
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_200m_400r_v2_${SLURM_JOB_ID}.out"
if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}'). Hint: SLURM_SUBMIT_DIR='${SLURM_SUBMIT_DIR:-<unset>}'." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi
FINAL_WATCH=$(grep "\[qnorm:watch\]" "${JOB_LOG}" | tail -1)
# R-402 fix: dead-rank detector keyed to the driver's fixed %.3e format
# (`M.MMMe[+-]NN`).  Matches either an exact zero (`0.000e+00`) or a
# sufficiently tiny magnitude (exponent -20 through -99, i.e., values that
# are effectively dead for this simulation).  The prior `r[0-9]+=0(\b|...)`
# / `=1e-3[0-9]+` combo false-matched `rN=0.0001` and (after the %.3e
# format switch) missed scientific zeros like `rN=0.000e+00`.
if grep -q "R-101 shared-fault DOFData consistency FAILED" "${JOB_LOG}" \
   || grep -q "R-101 shared-fault DOFData: .* unpaired entries" "${JOB_LOG}"; then
   echo "FAIL: R-101 shared-fault consistency diagnostic tripped" > "${RESULT_FILE}"
elif [ -z "${FINAL_WATCH}" ]; then
   echo "FAIL: no [qnorm:watch] line in job output — driver never produced diagnostic" > "${RESULT_FILE}"
elif echo "${FINAL_WATCH}" | grep -qE "r[0-9]+=(0\.0{3}e\+[0-9]{2}|[0-9]\.[0-9]{3}e-[2-9][0-9])([[:space:]]|$)"; then
   echo "FAIL: final [qnorm:watch] shows dead ranks (suggests rupture did not cross partition seam):" > "${RESULT_FILE}"
   echo "      ${FINAL_WATCH}" >> "${RESULT_FILE}"
else
   echo "PASS: cross-rank wave propagation evidence present; no R-101 abort" > "${RESULT_FILE}"
   echo "      ${FINAL_WATCH}" >> "${RESULT_FILE}"
fi
cat "${RESULT_FILE}"
