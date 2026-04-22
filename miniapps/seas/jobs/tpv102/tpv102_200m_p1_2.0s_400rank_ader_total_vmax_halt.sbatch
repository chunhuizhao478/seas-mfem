#!/bin/bash
#SBATCH -J tpv102_ader_vmax_400r
#SBATCH -o tpv102_ader_vmax_400r_%j.out
#SBATCH -e tpv102_ader_vmax_400r_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV102 v9.3.0 — ADER + total-Q production V_max halt test (400 ranks,
# 8 nodes, 200 m mesh, tfinal = 2.0 s).
#
# Purpose: mirror of `tpv102_200m_p1_2.0s_400rank_vmax_halt.sbatch`
# (the RK4 + fluctuation-Q reference from v8.0.0 Phase 1+2) but on
# the v9.3 purpose-changed stack:
#   [1] --time-integrator=ader (round-6 R-001 default flip + explicit
#       flag here for belt-and-braces against any future default change)
#   [2] total-Q formulation (wired automatically at driver init; the
#       driver's MFEM_VERIFY on wave.GetAbsorbingBackground() guarantees
#       no accidental fluctuation-Q fallback reaches wave.Mult or
#       AdvanceADER).
#
# Expected behaviour post-v9.3:
#   - Rupture initiates at t ~ 0.5-1 s (nucleation ramp peaks at
#     nuc_T = 1.0 s; the TPV102 hypocentre breaks at approximately
#     t ~ 0.5 s under the SCEC reference).
#   - V_max rises from 1e-12 m/s (V_ini) to 1-10 m/s during breakaway.
#   - Station output at the 9 SCEC fault stations shows non-trivial
#     slip rate (tau2 / V2 nonzero).
#   - ADER+total-Q should eliminate the v9.2.0 H-V92-K "_k4 − <field>"
#     drift by construction (plan Phase 7 §3) — the
#     `normal_stress_k4 - normal_stress` delta in the VTU must be zero
#     per cell under ADER.
#
# WALL-CLOCK budget: ≤ 2 hr on 8 Cascade Lake dev-queue nodes.  ADER
# at O=2 costs roughly the same per step as RK4; the extra predictor
# work offsets the single-stage advantage, so the 8750-step estimate
# from the reference sbatch holds.
# SU COST: ~512 SU (8 nodes × 2 hr × 32 SU/node-hour, budget).
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_2.0s_400rank_ader_total_vmax_halt.sbatch

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

RESULT_DIR="tpv102/results_ader_vmax_400r_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "v9.3.0 ADER + total-Q V_max halt test: 400 ranks / 8 nodes, tfinal=2.0 s"

# --time-integrator ader explicit, --ader-order 2 matches --order 1.
# --cfl 0.5 mirrors the RK4 reference; ADER's linear stability limit
# matches RK4's for this material / order (plan §Numerical constraints).
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

# Post-run classification: station-slip scoring + banner checks for
# ADER + total-Q.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_ader_vmax_400r_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_ader_vmax_400r_${SLURM_JOB_ID}.err"
if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}'). Hint: SLURM_SUBMIT_DIR='${SLURM_SUBMIT_DIR:-<unset>}'." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

# Banner sanity: ADER + total-Q must be the effective mode.
TOTAL_Q_LINE=$(grep -E "^State representation:" "${JOB_LOG}" | head -1)
ADER_LINE=$(grep -E "^Time integrator:" "${JOB_LOG}" | head -1)

# Early bailout for MFEM_ABORT (writes to stderr → .err).
if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData" "${JOB_ERR}" 2>/dev/null; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData" "${JOB_ERR}" | head -1)
   echo "FAIL: MFEM abort during V_max halt test" > "${RESULT_FILE}"
   echo "      ${ABORT_LINE}" >> "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

# Hard failure if ADER / total-Q is not the effective mode.
if ! echo "${TOTAL_Q_LINE}" | grep -q "total-Q"; then
   echo "FAIL: State representation banner shows NOT total-Q — production invariant violated" > "${RESULT_FILE}"
   echo "      banner: '${TOTAL_Q_LINE}'" >> "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi
if ! echo "${ADER_LINE}" | grep -q "ADER-O"; then
   echo "FAIL: Time integrator banner shows NOT ADER — production invariant violated" > "${RESULT_FILE}"
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

# V_max timeline from the stdout banner emitted every output-dt.
V_MAX_INIT=$(grep "V_max" "${JOB_LOG}" | head -1)
V_MAX_FINAL=$(grep "V_max" "${JOB_LOG}" | tail -1)
V_MAX_PEAK=$(grep "V_max" "${JOB_LOG}" \
             | awk '{for(i=1;i<=NF;i++) if ($i ~ /^[0-9.e+-]+$/) v=$i; print v}' \
             | sort -g | tail -1)
FINAL_WATCH=$(grep "\[qnorm:watch\]" "${JOB_LOG}" | tail -1)

echo "v9.3 ADER + total-Q V_max halt test (400-rank / tfinal=2.0 s)" > "${RESULT_FILE}"
echo "  Banners:" >> "${RESULT_FILE}"
echo "    ${TOTAL_Q_LINE}" >> "${RESULT_FILE}"
echo "    ${ADER_LINE}" >> "${RESULT_FILE}"
echo "  Station propagation: ${SLIP_COUNT} / 9 stations slipped" >> "${RESULT_FILE}"
echo "    Slipped:${SLIPPED_STATIONS}" >> "${RESULT_FILE}"
echo "    Silent :${SILENT_STATIONS}" >> "${RESULT_FILE}"
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
#   - 0/9 stations slipped AND V_max peak ~ 1e-12: no rupture (likely
#     a v9.3 regression, since the RK4+fluctuation baseline was
#     reported to evolve past breakaway).
#   - >=6/9 stations slipped AND V_max peak > 0.1 m/s: rupture evolved
#     under ADER + total-Q; v9.3 production stack validated.
#   - Intermediate: partial rupture — investigate station-by-station
#     and compare against the v8.0.0 / v9.2.0 RK4+fluctuation baseline.
if [ "${SLIP_COUNT}" -eq 0 ]; then
   echo "  Decision: 0/9 slipped under ADER + total-Q — investigate (possible v9.3 regression)" >> "${RESULT_FILE}"
elif [ "${SLIP_COUNT}" -ge 6 ]; then
   echo "  Decision: ${SLIP_COUNT}/9 slipped — rupture evolved under ADER + total-Q" >> "${RESULT_FILE}"
else
   echo "  Decision: ${SLIP_COUNT}/9 slipped — partial rupture; investigate per-station" >> "${RESULT_FILE}"
fi
cat "${RESULT_FILE}"
