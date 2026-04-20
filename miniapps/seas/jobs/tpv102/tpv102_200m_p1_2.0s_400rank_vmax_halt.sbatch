#!/bin/bash
#SBATCH -J tpv102_vmax_halt_400r
#SBATCH -o tpv102_vmax_halt_400r_%j.out
#SBATCH -e tpv102_vmax_halt_400r_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV102 v8.0.0 Phase 1+2 confirmation — "does V_max halt or evolve?"
#
# Configuration: the original job 7665297 reproducer (400 ranks / 8 nodes
# / 200 m mesh) — but extended to tfinal = 2.0 s to observe V_max
# behaviour past breakaway.  Breakaway expected at t ~ 1 s; the remaining
# ~1 s captures whether V_max halts (pinned-V_max signature of v8.0.0)
# or evolves through the rupture phase (healthy).
#
# Post Phase 1+2 fixes:
#   - R-101 pairing now uses BP5 face-vertex-key (`shared_fault_key.hpp`).
#   - Phase 2 adds centroid-within-face tiebreaker so MFEM face-orientation
#     flips across ranks do not cause spurious pairing failures.
#   - Local np=14/30/50 runs at 1000 m / tfinal=0.05 s confirm R-101
#     clean (bit-exact), V_max growing monotonically, waves reaching all
#     ranks.  What remains to be tested on Frontera:
#       (a) does the Phase 1+2 verifier fix hold on the 200 m mesh at
#           400-rank METIS partition topology?
#       (b) does V_max grow through breakaway and evolve through the
#           rupture phase (i.e., does the v8.0.0 "pinned V_max" signature
#           resolve)?  Hence tfinal=2.0 s (longer than the 1.2 s
#           bisection runs) to capture post-breakaway behaviour.
#
# Decision rules:
#   - R-101 aborts → Phase 2 fix did not generalise to 200 m topology;
#     escalate diagnosis.
#   - R-101 clean, V_max stays at V_ini=1e-12 → classic v8.0.0 signature
#     reproduces; Phase 1+2 pairing fix was necessary but insufficient.
#     Open Phase 3 (bulk propagation / ghost Q / CFL / nucleation).
#   - R-101 clean, V_max breakaway at t ~ 1 s then evolves toward 1-10
#     m/s → v8.0.0 bug resolved by the Phase 1+2 verifier fixes (the
#     "drift" cascade we chased was entirely verifier-pairing artifact).
#
# WALL-CLOCK budget: ≤ 2 hr on 8 Cascade Lake nodes (8750-step
# equivalent at 400-rank with 200 m mesh dt_cfl ~ 0.23 ms).  User
# confirmed fits in dev-queue 2 hr window based on prior runs at
# this scale.
# SU COST: ~512 SU (8 nodes × 2 hr × 32 SU/node-hour, budget).
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_2.0s_400rank_vmax_halt.sbatch

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
# Force full rebuild: rm executable + ALL dynamic/ .o files so any
# source change (user pulls before submit) unconditionally picks up.
rm -f seas_tpv102_driver drivers/tpv102_driver.o dynamic/*.o
make seas_tpv102_driver -j"${JOBS}" SEAS_EXTRA_CPPFLAGS="" \
   || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_200m_p1_2.0s_400r_vmax_halt_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase 1+2 confirmation: 400 ranks / 8 nodes, tfinal=2.0 s (V_max halt test)"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 2.0 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv \
      --debug-qnorm

# Post-run classification.  Station scoring is the same 9-station test
# as the Tier 3 sbatch, but the decision table is oriented toward
# "did the rupture evolve?" rather than "did the v8.0.0 bug reproduce?".
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_vmax_halt_400r_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_vmax_halt_400r_${SLURM_JOB_ID}.err"
if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}'). Hint: SLURM_SUBMIT_DIR='${SLURM_SUBMIT_DIR:-<unset>}'." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

# Early bailout for MFEM_ABORT (writes to stderr → .err).
if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData" "${JOB_ERR}" 2>/dev/null; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData" "${JOB_ERR}" | head -1)
   echo "FAIL: MFEM abort during V_max halt test" > "${RESULT_FILE}"
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

# V_max timeline: look at every reported V_max line to see whether it
# pinned at V_ini (1e-12) or rose through breakaway.
V_MAX_INIT=$(grep "V_max" "${JOB_LOG}" | head -1)
V_MAX_FINAL=$(grep "V_max" "${JOB_LOG}" | tail -1)
V_MAX_PEAK=$(grep "V_max" "${JOB_LOG}" \
             | awk '{for(i=1;i<=NF;i++) if ($i ~ /^[0-9.e+-]+$/) v=$i; print v}' \
             | sort -g | tail -1)
FINAL_WATCH=$(grep "\[qnorm:watch\]" "${JOB_LOG}" | tail -1)

if grep -q "R-101 shared-fault DOFData consistency FAILED\|R-101 shared-fault DOFData: .* unpaired entries" "${JOB_ERR}" 2>/dev/null \
   || grep -q "R-101 shared-fault DOFData consistency FAILED\|R-101 shared-fault DOFData: .* unpaired entries" "${JOB_LOG}"; then
   echo "FAIL: R-101 shared-fault consistency diagnostic tripped" > "${RESULT_FILE}"
elif grep -q "MFEM abort\|MFEM_ABORT" "${JOB_ERR}" 2>/dev/null \
   || grep -q "MFEM abort\|MFEM_ABORT" "${JOB_LOG}"; then
   echo "FAIL: MFEM abort during V_max halt test" > "${RESULT_FILE}"
else
   echo "V_max halt test (400-rank / tfinal=2.0 s) propagation score: ${SLIP_COUNT} / 9 stations slipped" > "${RESULT_FILE}"
   echo "  Slipped:${SLIPPED_STATIONS}" >> "${RESULT_FILE}"
   echo "  Silent :${SILENT_STATIONS}" >> "${RESULT_FILE}"
   if [ -n "${V_MAX_INIT}" ]; then
      echo "  First  V_max line: ${V_MAX_INIT}" >> "${RESULT_FILE}"
   fi
   if [ -n "${V_MAX_PEAK}" ]; then
      echo "  Peak   V_max value: ${V_MAX_PEAK}" >> "${RESULT_FILE}"
   fi
   if [ -n "${V_MAX_FINAL}" ]; then
      echo "  Final  V_max line: ${V_MAX_FINAL}" >> "${RESULT_FILE}"
   fi
   if [ -n "${FINAL_WATCH}" ]; then
      echo "  Final qnorm watch: ${FINAL_WATCH}" >> "${RESULT_FILE}"
   fi
   # Decision:
   #   - 0/9 stations slipped AND V_max peak ~ 1e-12: v8.0.0 reproduces
   #     (Phase 1+2 fix was necessary but not sufficient; open Phase 3).
   #   - >=6/9 stations slipped AND V_max peak > 0.1 m/s: rupture evolved
   #     past breakaway; v8.0.0 resolved by Phase 1+2.
   #   - Intermediate: partial rupture — investigate station-by-station.
   if [ "${SLIP_COUNT}" -eq 0 ]; then
      echo "  Decision: 0/9 slipped — v8.0.0 signature persists; open Phase 3 (bulk propagation)" >> "${RESULT_FILE}"
   elif [ "${SLIP_COUNT}" -ge 6 ]; then
      echo "  Decision: ${SLIP_COUNT}/9 slipped — rupture evolved; v8.0.0 resolved by Phase 1+2" >> "${RESULT_FILE}"
   else
      echo "  Decision: ${SLIP_COUNT}/9 slipped — partial rupture; investigate per-station" >> "${RESULT_FILE}"
   fi
fi
cat "${RESULT_FILE}"
