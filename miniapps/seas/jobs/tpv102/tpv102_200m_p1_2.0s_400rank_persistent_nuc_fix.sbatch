#!/bin/bash
#SBATCH -J tpv102_persist_nuc_400r
#SBATCH -o tpv102_persist_nuc_400r_%j.out
#SBATCH -e tpv102_persist_nuc_400r_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV102 nucleation persistent-prestress channel fix — re-run of the
# 400-rank V_max halt test under the 2026-04-22 fix.
#
# What changed since
# `tpv102_200m_p1_2.0s_400rank_ader_total_vmax_halt.sbatch`:
#   - The bulk-Q add-delta nucleation path (`ApplyNucleationTotal(Q,
#     map, state, ...)`) is replaced by the persistent-prestress
#     channel (`ApplyNucleationTotalPrestress(dof_data, ..., t)` →
#     `DOFData.tau2_nuc` → additive read in
#     `FaultFaceFlux::EvaluateTotal`).  See
#       miniapps/seas/debug_document/tpv102_debug_document/
#         tpv102_nucleation_code_review_2026-04-22.md
#         tpv102_nucleation_code_fix_2026-04-22.md
#   - The legacy `FaultQPNodalMap` / `FaultQPNucleationState` /
#     `BuildFaultQPNodalMap` machinery remains in `tpv102_setup_total.hpp`
#     so the R-002 / R-009 / R-I06-003 unit tests stay green, but the
#     production driver no longer touches it (verified by `grep
#     ApplyNucleationTotal drivers/tpv102_driver.cpp` returning only
#     comment hits).
#
# Expected behaviour:
#   - V_max stays at 1e-12 m/s (V_ini) for t < ~0.5 s.
#   - As the nucleation ramp peaks at nuc_T = 1.0 s, V_max climbs by
#     ORDERS of magnitude — the fix's whole point is that the
#     friction solver now sees the FULL 25 MPa nucleation stress
#     (under the broken add-delta path it saw ~kPa, hence V_max
#     stalled at 2e-7 m/s in the prior 400-rank run).
#   - Breakaway around t ~ 0.5-1 s, V_max into the 1-10 m/s range.
#   - At least 6 / 9 SCEC fault stations should record nontrivial
#     slip rate by t = 2.0 s.
#
# WALL-CLOCK budget: ≤ 2 hr on 8 Cascade Lake dev-queue nodes.
# SU COST: ~512 SU (8 nodes × 2 hr × 32 SU/node-hour).
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_2.0s_400rank_persistent_nuc_fix.sbatch

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
# Force full rebuild: rm executable + ALL dynamic/ .o files so the
# 2026-04-22 fix unconditionally picks up.
rm -f seas_tpv102_driver drivers/tpv102_driver.o dynamic/*.o
make seas_tpv102_driver -j"${JOBS}" SEAS_EXTRA_CPPFLAGS="" \
   || { echo "BUILD FAILED"; exit 1; }

# Local fix-regression unit test in-batch (cheap, ~1 s, single-rank):
# guards against a stale build dir picking up the fix in the executable
# but not in the test object file.
make seas_test_persistent_nuc_prestress_channel -j"${JOBS}" \
   || { echo "FIX UNIT TEST BUILD FAILED"; exit 1; }
./seas_test_persistent_nuc_prestress_channel \
   || { echo "FIX UNIT TEST FAILED"; exit 1; }

RESULT_DIR="tpv102/results_persist_nuc_400r_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "TPV102 persistent-prestress fix: 400 ranks / 8 nodes, tfinal=2.0 s"

# CLI mirrors the working ADER+total-Q halt sbatch — only the source
# changed, not the runtime configuration.
ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --time-integrator ader \
      --ader-order 2 \
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

# Post-run classification: same shape as the prior halt sbatch so the
# RESULT.txt format stays comparable.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_persist_nuc_400r_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_persist_nuc_400r_${SLURM_JOB_ID}.err"
if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}'). Hint: SLURM_SUBMIT_DIR='${SLURM_SUBMIT_DIR:-<unset>}'." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

TOTAL_Q_LINE=$(grep -E "^State representation:" "${JOB_LOG}" | head -1)
ADER_LINE=$(grep -E "^Time integrator:" "${JOB_LOG}" | head -1)

if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData" "${JOB_ERR}" 2>/dev/null; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData" "${JOB_ERR}" | head -1)
   echo "FAIL: MFEM abort during persistent-prestress fix run" > "${RESULT_FILE}"
   echo "      ${ABORT_LINE}" >> "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

if ! echo "${TOTAL_Q_LINE}" | grep -q "total-Q"; then
   echo "FAIL: State representation banner shows NOT total-Q" > "${RESULT_FILE}"
   echo "      banner: '${TOTAL_Q_LINE}'" >> "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi
if ! echo "${ADER_LINE}" | grep -q "ADER-O"; then
   echo "FAIL: Time integrator banner shows NOT ADER" > "${RESULT_FILE}"
   echo "      banner: '${ADER_LINE}'" >> "${RESULT_FILE}"
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

V_MAX_INIT=$(grep "V_max" "${JOB_LOG}" | head -1)
V_MAX_FINAL=$(grep "V_max" "${JOB_LOG}" | tail -1)
V_MAX_PEAK=$(grep "V_max" "${JOB_LOG}" \
             | awk '{for(i=1;i<=NF;i++) if ($i ~ /^[0-9.e+-]+$/) v=$i; print v}' \
             | sort -g | tail -1)
FINAL_WATCH=$(grep "\[qnorm:watch\]" "${JOB_LOG}" | tail -1)

echo "TPV102 persistent-prestress fix (400-rank / tfinal=2.0 s)" > "${RESULT_FILE}"
echo "  Banners:" >> "${RESULT_FILE}"
echo "    ${TOTAL_Q_LINE}" >> "${RESULT_FILE}"
echo "    ${ADER_LINE}" >> "${RESULT_FILE}"
echo "  Station propagation: ${SLIP_COUNT} / 9 stations slipped" >> "${RESULT_FILE}"
echo "    Slipped:${SLIPPED_STATIONS}" >> "${RESULT_FILE}"
echo "    Silent :${SILENT_STATIONS}" >> "${RESULT_FILE}"
if [ -n "${V_MAX_INIT}"  ]; then echo "  First  V_max line: ${V_MAX_INIT}"  >> "${RESULT_FILE}"; fi
if [ -n "${V_MAX_PEAK}"  ]; then echo "  Peak   V_max value: ${V_MAX_PEAK}" >> "${RESULT_FILE}"; fi
if [ -n "${V_MAX_FINAL}" ]; then echo "  Final  V_max line: ${V_MAX_FINAL}" >> "${RESULT_FILE}"; fi
if [ -n "${FINAL_WATCH}" ]; then echo "  Final qnorm watch: ${FINAL_WATCH}" >> "${RESULT_FILE}"; fi

# Fix-specific decision branch: contrast against the broken-baseline
# 2e-7 m/s V_max the prior halt sbatch reported.
if [ "${SLIP_COUNT}" -eq 0 ]; then
   echo "  Decision: 0/9 slipped — fix did NOT recover rupture; investigate (possible second bug)" >> "${RESULT_FILE}"
elif [ "${SLIP_COUNT}" -ge 6 ]; then
   echo "  Decision: ${SLIP_COUNT}/9 slipped — persistent-prestress fix UNBLOCKED nucleation, rupture evolved" >> "${RESULT_FILE}"
else
   echo "  Decision: ${SLIP_COUNT}/9 slipped — partial rupture; compare to broken-baseline V_max stall (2e-7 m/s)" >> "${RESULT_FILE}"
fi
cat "${RESULT_FILE}"
