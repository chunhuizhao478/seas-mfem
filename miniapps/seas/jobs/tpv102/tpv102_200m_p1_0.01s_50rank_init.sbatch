#!/bin/bash
#SBATCH -J tpv102_init_50r
#SBATCH -o tpv102_init_50r_%j.out
#SBATCH -e tpv102_init_50r_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 50
#SBATCH -t 00:20:00
#SBATCH -A EAR20006

# TPV102 v8.0.0 debug plan — Phase 1 Escalation A-Frontera, step (1).
#
# Purpose: cheap init-only sanity on the 200 m mesh at 50 ranks.
# Confirms:
#   - driver builds and starts on Frontera with current HEAD
#   - [BUILD] banner writes to build_info.txt (R-501 / R-606)
#   - R-101 shared-fault DOFData consistency passes at 50-rank
#     scale (pairs matched > 0, max_rel_diff within tolerance)
#   - no MFEM_ABORT during ParMesh / face-neighbor setup
#
# NOT a physics test — `tfinal = 0.01 s` is pre-nucleation and
# cannot reproduce the v8.0.0 pinned-V_max signature.  Purely a
# plumbing check before the longer 100-rank / 400-rank jobs.
#
# WALL-CLOCK: ~5 min on 1 development node.
# SU COST: ~0.1 SU (1 node-hour allocated, < 1/6 of that used).
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_0.01s_50rank_init.sbatch

export LC_ALL=C
export LANG=C

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
# Ensure a clean rebuild so the new [BUILD] banner compiles in even if
# .o files were built before the banner patch.
rm -f drivers/tpv102_driver.o dynamic/wave_operator.o dynamic/fault_face_flux.o
make seas_tpv102_driver -j"${JOBS}" SEAS_EXTRA_CPPFLAGS="" \
   || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_init_50r_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase 1 Escalation A-Frontera step (1): init-only sanity (tfinal=0.01 s)"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 0.01 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.01 \
      --debug-qnorm

# Post-run machine-checkable success criterion for init-only.
# PASS iff:
#   (a) [BUILD] banner present (confirms new banner compiled in),
#   (b) [R-101 check] line shows pairs matched > 0 OR the SKIPPED
#       message ("no shared fault faces"),
#   (c) no MFEM_ABORT in the log.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_init_50r_${SLURM_JOB_ID}.out"
if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}'). Hint: SLURM_SUBMIT_DIR='${SLURM_SUBMIT_DIR:-<unset>}'." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

BUILD_LINE=$(grep "\[BUILD\]" "${JOB_LOG}" | head -1)
R101_LINE=$(grep -E "\[R-101 check\]|pairs matched" "${JOB_LOG}" | head -1)

if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData consistency FAILED" "${JOB_LOG}"; then
   echo "FAIL: MFEM abort during init" > "${RESULT_FILE}"
elif [ -z "${BUILD_LINE}" ]; then
   echo "FAIL: no [BUILD] banner — R-501 banner patch did not compile in" > "${RESULT_FILE}"
elif [ -z "${R101_LINE}" ]; then
   echo "FAIL: no R-101 diagnostic line (driver did not reach init or diagnostic suppressed)" > "${RESULT_FILE}"
else
   echo "PASS: init sanity complete" > "${RESULT_FILE}"
   echo "      ${BUILD_LINE}" >> "${RESULT_FILE}"
   echo "      ${R101_LINE}" >> "${RESULT_FILE}"
fi
cat "${RESULT_FILE}"
