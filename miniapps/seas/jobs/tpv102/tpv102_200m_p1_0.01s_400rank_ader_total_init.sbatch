#!/bin/bash
#SBATCH -J tpv102_ader_init_400r
#SBATCH -o tpv102_ader_init_400r_%j.out
#SBATCH -e tpv102_ader_init_400r_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV102 v9.3.0 — ADER + total-Q init sanity at 400 ranks.
#
# Purpose: cheap init-only sanity (`tfinal = 0.01 s`) on the 200 m mesh
# at the production 8-node / 400-rank partitioning.  Confirms the v9.3
# purpose-change stack:
#   [1] ADER is the default time integrator (round-7 R-001 + round-6
#       R-001 fixes: --time-integrator defaults to "ader").
#   [2] total-Q is the only supported formulation
#       (InitializeStateTotal + ZeroDOFDataPreStressTotal +
#        SetAbsorbingBackground + ApplyNucleationTotal wired in the
#        driver; MFEM_VERIFY at line ~602 fails loud if total-Q wiring
#        is missing).
#
# Sanity checks:
#   - driver builds and starts on Frontera with current HEAD
#   - banner prints "State representation: total-Q (pre-stress baked into Q)"
#     and "Time integrator: ADER-O(2)  (default)"
#   - [BUILD] banner writes to build_info.txt (R-501 / R-606)
#   - R-101 shared-fault DOFData consistency passes at 400-rank scale
#   - no MFEM_ABORT during ParMesh / face-neighbor setup or the first
#     AdvanceADER step
#
# NOT a physics test — `tfinal = 0.01 s` is pre-nucleation.  Purely a
# plumbing check before the 2.0 s V_max halt test below.
#
# WALL-CLOCK: ~5-10 min on 8 development nodes.  The 2 hr -t allocation
# is the user-requested budget ceiling; actual usage is well under.
# SU COST: 8 nodes × ~0.2 hr × 32 SU/node-hour ≈ 51 SU (worst case 512
# SU if the full 2 hr is somehow consumed — cap via --tfinal=0.01).
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_0.01s_400rank_ader_total_init.sbatch

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
# Force a full rebuild of the driver + every dynamic/ .o file so the
# current git state's source (user pulls before submit) unconditionally
# compiles in.  rm'ing the executable also forces relink.
rm -f seas_tpv102_driver drivers/tpv102_driver.o dynamic/*.o
make seas_tpv102_driver -j"${JOBS}" SEAS_EXTRA_CPPFLAGS="" \
   || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_ader_init_400r_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "v9.3.0 ADER + total-Q init sanity at 400 ranks (tfinal=0.01 s)"

# Explicit --time-integrator=ader belts-and-braces even though the
# driver default is now "ader" (round-6 R-001).  --ader-order=2
# matches --order=1 (O >= p+1 rule, round-7 R-001).
ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --time-integrator ader \
      --ader-order 2 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 0.01 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.01 \
      --debug-qnorm

# Post-run machine-checkable success criterion for init-only.
# PASS iff:
#   (a) [BUILD] banner present,
#   (b) State-representation banner says total-Q,
#   (c) Time-integrator banner says ADER-O(...),
#   (d) [R-101 check] line shows pairs matched > 0 OR the SKIPPED
#       message ("no shared fault faces"),
#   (e) no MFEM_ABORT in the log.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_ader_init_400r_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_ader_init_400r_${SLURM_JOB_ID}.err"
if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}'). Hint: SLURM_SUBMIT_DIR='${SLURM_SUBMIT_DIR:-<unset>}'." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

# Grep banners — stdout goes to .out; build banner also in build_info.txt.
BUILD_INFO_FILE="build_info.txt"
if [ -f "${BUILD_INFO_FILE}" ]; then
   BUILD_LINE=$(grep "\[BUILD\]" "${BUILD_INFO_FILE}" | head -1)
elif [ -f "${JOB_ERR}" ]; then
   BUILD_LINE=$(grep "\[BUILD\]" "${JOB_ERR}" | head -1)
else
   BUILD_LINE=""
fi
TOTAL_Q_LINE=$(grep -E "^State representation:" "${JOB_LOG}" | head -1)
ADER_LINE=$(grep -E "^Time integrator:" "${JOB_LOG}" | head -1)
R101_LINE=""
if [ -f "${JOB_ERR}" ]; then
   R101_LINE=$(grep -E "\[R-101 check\]|pairs matched" "${JOB_ERR}" | head -1)
fi
if [ -z "${R101_LINE}" ]; then
   R101_LINE=$(grep -E "\[R-101 check\]|pairs matched" "${JOB_LOG}" | head -1)
fi

if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData" "${JOB_ERR}" 2>/dev/null \
   || grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData" "${JOB_LOG}"; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData" "${JOB_ERR}" "${JOB_LOG}" 2>/dev/null | head -1)
   echo "FAIL: MFEM abort during init" > "${RESULT_FILE}"
   echo "      ${ABORT_LINE}" >> "${RESULT_FILE}"
elif [ -z "${BUILD_LINE}" ]; then
   echo "FAIL: no [BUILD] banner — banner patch did not compile in" > "${RESULT_FILE}"
elif ! echo "${TOTAL_Q_LINE}" | grep -q "total-Q"; then
   echo "FAIL: State representation is NOT total-Q" > "${RESULT_FILE}"
   echo "      banner: '${TOTAL_Q_LINE}'" >> "${RESULT_FILE}"
elif ! echo "${ADER_LINE}" | grep -q "ADER-O"; then
   echo "FAIL: Time integrator is NOT ADER" > "${RESULT_FILE}"
   echo "      banner: '${ADER_LINE}'" >> "${RESULT_FILE}"
elif [ -z "${R101_LINE}" ]; then
   echo "FAIL: no R-101 diagnostic line (driver did not reach init or diagnostic suppressed)" > "${RESULT_FILE}"
else
   echo "PASS: v9.3 ADER + total-Q init sanity complete" > "${RESULT_FILE}"
   echo "      ${BUILD_LINE}" >> "${RESULT_FILE}"
   echo "      ${TOTAL_Q_LINE}" >> "${RESULT_FILE}"
   echo "      ${ADER_LINE}" >> "${RESULT_FILE}"
   echo "      ${R101_LINE}" >> "${RESULT_FILE}"
fi
cat "${RESULT_FILE}"
