#!/bin/bash
#SBATCH -J tpv102_200m_50r
#SBATCH -o tpv102_200m_50r_%j.out
#SBATCH -e tpv102_200m_50r_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 50
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV102 dev-queue intermediate-scale run (v2 recommendation #2).
#
# Purpose: exercise the full (v1+v2) fix on a mesh where Λ_dyn/h ≈ 0.8
# (borderline-resolved cohesive zone; rupture should at least start
# propagating past the nucleation patch).  50 ranks gives METIS enough
# headroom to create real cross-rank partition seams in both the bulk
# (always) and — depending on the partitioning — possibly across the
# fault as shared faces.
#
# Diagnostics expected:
#   - `--debug-qnorm` watch list shows ||Q||_inf growing on non-
#     hypocenter ranks within ~100 ms of nucleation (c_s × 0.1 s = 350 m,
#     comparable to partition neighborhoods at h=200 m).
#   - R-101 runtime diagnostic prints
#     `[R-101 check] shared-fault DOFData consistency OK: N pairs matched`
#     OR returns silently if no shared fault faces exist.  Both are
#     acceptable; an MFEM_ABORT from that diagnostic is a hard-fail
#     signal that the fix is incomplete.
#   - Fault-surface PVD every 0.1 s (volume PVD suppressed for disk).
#
# WALL-TIME BUDGET: ~2 hr dev queue; 200 m mesh at 50 ranks runs
# ~0.5 s/step — 1.5 s physical time ≈ 5250 steps → ~3.0 hr.  Set
# tfinal=1.0 s here so we complete within the 2 hr window; this
# still covers the 1 s nucleation ramp and ~0.3 s of post-nucleation
# propagation.  For 1.5 s or longer, use the normal queue.
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch

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

RESULT_DIR="tpv102/results_200m_p1_1.0s_50r_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "tfinal=1.0 s (dev-queue budget), dt=auto, order=1"
echo "debug-qnorm ON, fault-surface PVD only, output every 0.1 s"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 1.0 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv \
      --debug-qnorm

# R-306 fix: post-run machine-checkable success criterion.  PASS iff
# (a) no R-101 shared-fault consistency abort, and (b) the final
# [qnorm:watch] line shows non-hypocenter ranks carrying non-trivial
# ||Q||_inf (i.e., wave energy crossed at least one partition seam).
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
# R-403 fix: absolute path anchored at the submit dir (see SBATCH_LOG_DIR
# at top of script).  The prior relative path failed whenever the user
# submitted from anywhere other than miniapps/seas/.
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_200m_50r_${SLURM_JOB_ID}.out"
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
   echo "FAIL: final [qnorm:watch] shows dead ranks (rupture may not have crossed seam):" > "${RESULT_FILE}"
   echo "      ${FINAL_WATCH}" >> "${RESULT_FILE}"
else
   echo "PASS: cross-rank wave propagation evidence present; no R-101 abort" > "${RESULT_FILE}"
   echo "      ${FINAL_WATCH}" >> "${RESULT_FILE}"
fi
cat "${RESULT_FILE}"
