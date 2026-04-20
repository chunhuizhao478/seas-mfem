#!/bin/bash
#SBATCH -J tpv102_bisect_400r
#SBATCH -o tpv102_bisect_400r_%j.out
#SBATCH -e tpv102_bisect_400r_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV102 v8.0.0 debug plan — Phase 1 Tier 3 dispositive reproducer.
#
# Purpose: repeat the known-failing 400-rank configuration at
# tfinal=1.2 s (R-502) on the current v7.0.0 post-fix HEAD
# (commit c015048).  This is the bug-regime ground-truth run:
# it should reproduce the v8.0.0 signature — V_max pinned at
# 7.699 m/s at hypocenter, off-hypocenter stations at V_ini=1e-12,
# bulk `max‖Q‖_∞ ≈ 1.66e-8`.  If the bug does NOT reproduce here,
# the reproducer has drifted and the plan must be re-anchored.
#
# Pair with 100-rank Tier 2 bisection to identify scale-dependence
# per the decision table in the 100-rank sbatch header.
#
# Plan alignment:
#   - tfinal=1.2 s (was 1.5 s in the pre-existing 400rank_dev.sbatch;
#     R-502 says 1.2 s is sufficient for breakaway + saturation and
#     keeps us in the 2-hr dev-queue window).
#   - mesh tpv102_200m.msh, order=1, CFL=0.5 (driver auto dt).
#   - --debug-qnorm to diagnose cross-rank wave propagation.
#
# WALL-CLOCK: ~1.5 hr on 8 Cascade Lake nodes (5250-step
# equivalent at 400-rank).
# SU COST: ~400 SU (8 nodes × 1.7 hr × 32 SU/node-hour, budget).
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_1.2s_400rank_bisect.sbatch

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
rm -f drivers/tpv102_driver.o dynamic/wave_operator.o dynamic/fault_face_flux.o
make seas_tpv102_driver -j"${JOBS}" SEAS_EXTRA_CPPFLAGS="" \
   || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_200m_p1_1.2s_400r_bisect_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase 1 Tier 3 dispositive: 400 ranks / 8 nodes, tfinal=1.2 s"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 1.2 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv \
      --debug-qnorm

# Phase 1 success criterion [R-506] — same scoring as Tier 2.
# Expected: 0/9 stations slip (bug reproduces).  Any other result
# means the reproducer has changed.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_bisect_400r_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_bisect_400r_${SLURM_JOB_ID}.err"
if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}'). Hint: SLURM_SUBMIT_DIR='${SLURM_SUBMIT_DIR:-<unset>}'." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

# Early bailout for MFEM_ABORT (writes to stderr → .err).
if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData" "${JOB_ERR}" 2>/dev/null; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData" "${JOB_ERR}" | head -1)
   echo "FAIL: MFEM abort during Tier 3" > "${RESULT_FILE}"
   echo "      ${ABORT_LINE}" >> "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

SLIP_COUNT=0
SLIPPED_STATIONS=""
SILENT_STATIONS=""
for sta in flt_0_3 flt_0_7.5 flt_0_12 flt_9_7.5 flt_12_3 flt_12_12 \
           flt_n9_7.5 flt_n12_3 flt_n12_12; do
   STA_FILE="${RESULT_DIR}/${sta}.dat"
   if [ ! -f "${STA_FILE}" ]; then
      SILENT_STATIONS="${SILENT_STATIONS} ${sta}(missing)"
      continue
   fi
   V2=$(awk '!/^#/ && NF>=5 {v=$5} END {print (v=="" ? "0" : v)}' "${STA_FILE}")
   if awk "BEGIN {exit !( ($V2)*($V2) > 1e-12 )}"; then
      SLIP_COUNT=$((SLIP_COUNT + 1))
      SLIPPED_STATIONS="${SLIPPED_STATIONS} ${sta}"
   else
      SILENT_STATIONS="${SILENT_STATIONS} ${sta}"
   fi
done

FINAL_WATCH=$(grep "\[qnorm:watch\]" "${JOB_LOG}" | tail -1)
V_MAX_LINE=$(grep "V_max" "${JOB_LOG}" | tail -1)

if grep -q "R-101 shared-fault DOFData consistency FAILED\|R-101 shared-fault DOFData: .* unpaired entries" "${JOB_ERR}" 2>/dev/null \
   || grep -q "R-101 shared-fault DOFData consistency FAILED\|R-101 shared-fault DOFData: .* unpaired entries" "${JOB_LOG}"; then
   echo "FAIL: R-101 shared-fault consistency diagnostic tripped" > "${RESULT_FILE}"
elif grep -q "MFEM abort\|MFEM_ABORT" "${JOB_ERR}" 2>/dev/null \
   || grep -q "MFEM abort\|MFEM_ABORT" "${JOB_LOG}"; then
   echo "FAIL: MFEM abort during Tier 3" > "${RESULT_FILE}"
else
   echo "Tier 3 (400-rank, dispositive) propagation score: ${SLIP_COUNT} / 9 stations slipped" > "${RESULT_FILE}"
   echo "  Slipped:${SLIPPED_STATIONS}" >> "${RESULT_FILE}"
   echo "  Silent :${SILENT_STATIONS}" >> "${RESULT_FILE}"
   if [ -n "${V_MAX_LINE}" ]; then
      echo "  Final V_max line: ${V_MAX_LINE}" >> "${RESULT_FILE}"
   fi
   if [ -n "${FINAL_WATCH}" ]; then
      echo "  Final qnorm watch: ${FINAL_WATCH}" >> "${RESULT_FILE}"
   fi
   if [ "${SLIP_COUNT}" -eq 0 ]; then
      echo "  Expected: 0/9 (v8.0.0 bug reproduced — ground-truth anchor for diagnosis)" >> "${RESULT_FILE}"
   else
      echo "  UNEXPECTED: ${SLIP_COUNT}/9 slipped — v8.0.0 reproducer may have drifted" >> "${RESULT_FILE}"
   fi
fi
cat "${RESULT_FILE}"
