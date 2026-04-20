#!/bin/bash
#SBATCH -J tpv102_bisect_100r
#SBATCH -o tpv102_bisect_100r_%j.out
#SBATCH -e tpv102_bisect_100r_%j.err
#SBATCH -p development
#SBATCH -N 2
#SBATCH -n 100
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV102 v8.0.0 debug plan — Phase 1 Run 1b (Tier 2 bisection).
#
# Purpose: determine whether the v8.0.0 pinned-V_max bug manifests
# at intermediate rank count.  The bug is known to reproduce at
# 400-rank / tfinal=1.5 s (job 7665297).  Run 1a (np=1, separate
# sbatch) tests the interior-only path.  This run sits between:
# 100 ranks across 2 Cascade Lake nodes exercises BOTH the
# interior-fault branch (intra-node faces) AND the shared-fault
# branch (inter-node partition seams) without the full 8-node
# topology complexity.
#
# Decision mapping (combined with 1-rank Run 1a and 400-rank Tier 3):
#   - 1-rank 0/9, 100-rank 0/9, 400-rank 0/9  → Case A: interior
#     branch / Evaluate / accumulation bug (scale-invariant).
#   - 1-rank 4-9/9, 100-rank 0/9, 400-rank 0/9 → Case B: shared
#     branch bug (manifests whenever shared fault faces exist).
#   - 1-rank 4-9/9, 100-rank 4-9/9, 400-rank 0/9 → Case C: bug
#     requires 400-rank topology (inter-rank resonance / MPI
#     collective edge case).
#
# Per plan R-502: tfinal = 1.2 s is REQUIRED to cross breakaway
# (V_max > 0.1 m/s at t ≈ 1.0 s).  tfinal = 0.3 s does NOT
# reproduce the signature.
#
# Per plan R-505: 50-rank is intentionally SKIPPED for bisection;
# 50-rank × tfinal=1.2 s would exceed 2 hr wall-clock on dev queue.
# 100-rank / 2-nodes gives ~1.2 hr.
#
# WALL-CLOCK: ~1.2 hr on 2 Cascade Lake nodes.
# SU COST: ~100 SU (2 nodes × 1.5 hr × 32 SU/node-hour, budget).
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_1.2s_100rank_bisect.sbatch

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

RESULT_DIR="tpv102/results_200m_p1_1.2s_100r_bisect_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase 1 Run 1b (Tier 2 bisection): 100 ranks / 2 nodes, tfinal=1.2 s"

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

# Phase 1 success criterion [R-506]: station propagation score.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_bisect_100r_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_bisect_100r_${SLURM_JOB_ID}.err"
if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}'). Hint: SLURM_SUBMIT_DIR='${SLURM_SUBMIT_DIR:-<unset>}'." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

# Early bailout for MFEM_ABORT (writes to stderr → .err).
if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData" "${JOB_ERR}" 2>/dev/null; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData" "${JOB_ERR}" | head -1)
   echo "FAIL: MFEM abort during Run 1b (Tier 2)" > "${RESULT_FILE}"
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
   # Column 5 = V2 (strike component, BP5 convention).
   V2=$(awk '!/^#/ && NF>=5 {v=$5} END {print (v=="" ? "0" : v)}' "${STA_FILE}")
   if awk "BEGIN {exit !( ($V2)*($V2) > 1e-12 )}"; then
      SLIP_COUNT=$((SLIP_COUNT + 1))
      SLIPPED_STATIONS="${SLIPPED_STATIONS} ${sta}"
   else
      SILENT_STATIONS="${SILENT_STATIONS} ${sta}"
   fi
done

# Qnorm cross-rank tail watch.
FINAL_WATCH=$(grep "\[qnorm:watch\]" "${JOB_LOG}" | tail -1)

if grep -q "R-101 shared-fault DOFData consistency FAILED\|R-101 shared-fault DOFData: .* unpaired entries" "${JOB_ERR}" 2>/dev/null \
   || grep -q "R-101 shared-fault DOFData consistency FAILED\|R-101 shared-fault DOFData: .* unpaired entries" "${JOB_LOG}"; then
   echo "FAIL: R-101 shared-fault consistency diagnostic tripped" > "${RESULT_FILE}"
elif grep -q "MFEM abort\|MFEM_ABORT" "${JOB_ERR}" 2>/dev/null \
   || grep -q "MFEM abort\|MFEM_ABORT" "${JOB_LOG}"; then
   echo "FAIL: MFEM abort during Run 1b (Tier 2)" > "${RESULT_FILE}"
else
   echo "Run 1b Tier 2 propagation score: ${SLIP_COUNT} / 9 stations slipped" > "${RESULT_FILE}"
   echo "  Slipped:${SLIPPED_STATIONS}" >> "${RESULT_FILE}"
   echo "  Silent :${SILENT_STATIONS}" >> "${RESULT_FILE}"
   if [ -n "${FINAL_WATCH}" ]; then
      echo "  Final qnorm watch: ${FINAL_WATCH}" >> "${RESULT_FILE}"
   fi
   if [ "${SLIP_COUNT}" -eq 0 ]; then
      echo "  Decision: 0/9 — bug reproduces at 100-rank (NOT scale-exclusive)" >> "${RESULT_FILE}"
   elif [ "${SLIP_COUNT}" -lt 4 ]; then
      echo "  Decision: 1-3/9 partial — Case B-partial" >> "${RESULT_FILE}"
   else
      echo "  Decision: 4-9/9 healthy at 100-rank — bug is 400-rank-specific (Case C)" >> "${RESULT_FILE}"
   fi
fi
cat "${RESULT_FILE}"
