#!/bin/bash
#SBATCH -J tpv102_v92_dt_qtr
#SBATCH -o tpv102_v92_dt_qtr_%j.out
#SBATCH -e tpv102_v92_dt_qtr_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 16:00:00
#SBATCH -A EAR20006

# TPV102 v9.2.0 plan Step 1 — dt-quartering companion to
# `tpv102_200m_p1_12.0s_400r_v92_dt_half.sbatch`.
#
# Identical to the v91 12 s baseline (job 7668434) EXCEPT `--cfl 0.125`
# (quarters dt vs the default cfl=0.5).  Paired with `_dt_half` to
# pin down the truncation order of the RK4/psi coupling:
#
#   | slip_dip ratio (v91 / dt_half) × (dt_half / dt_quarter) | Scaling |
#   |---------------------------------------------------------|---------|
#   | ≈ 16 × 16                                               | dt⁴ (pure RK4)       |
#   | ≈  4 ×  4                                               | dt² (psi split-step) |
#   | ≈  1 ×  1                                               | non-truncation       |
#
# Both runs compare the slip_dip peak at the 9 SCEC stations against
# the v91 baseline (1.83 m hypocentre peak, see §1.2 of plan).
#
# PREREQUISITE — source tree matches the one used for v91 job 7668434
# (same post-v9.1.0 head).  No v9.2.0 C0 source edits required.
#
# Plan references:
#   - debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_debug_plan.md
#     PRIORITY NOW Step 1; REVIEW round-5 R-V92-G01
#   - tpv102_200m_p1_12.0s_400r_v92_dt_half.sbatch (paired run)
#   - tpv102_200m_p1_12.0s_400rank_v91.sbatch (baseline)
#   - feedback_sbatch_modules.md
#   - feedback_frontera_approval.md — EXPLICIT APPROVAL REQUIRED before sbatch.
#
# WALL-CLOCK budget: ~16 hr (4× v91 baseline because nsteps quadruples).
# SU COST          : ~4096 SU (8 nodes × 16 hr × 32 SU/node-hr).
#                    If 16 hr exceeds the queue limit, reduce tfinal:
#                    `tfinal 4.0` → ~5 hr wall / ~1280 SU, still covers
#                    the rupture-tip growth window (§18 shows peaks
#                    already present at 4 s).
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_12.0s_400r_v92_dt_quarter.sbatch

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
rm -f seas_tpv102_driver drivers/tpv102_driver.o dynamic/*.o io/*.o fault/*.o
make seas_tpv102_driver -j"${JOBS}" SEAS_EXTRA_CPPFLAGS="" \
   || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_200m_p1_12.0s_400r_v92_dt_quarter_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase v9.2.0 Step 1 (dt/4): 400 ranks / 8 nodes / 200 m / P1 / tfinal=12.0 s / cfl=0.125"
echo "Purpose  : RK4-truncation order discrimination, paired with dt_half run"
echo "Baseline : v91 job 7668434 (cfl=0.5, slip_dip peak 1.83 m)"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.125 \
      --tfinal 12.0 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv \
      --debug-qnorm

# -------------------------------------------------------------------------
# Post-run classification for R-V92 Step 1 (dt/4 decision table).
# Same post-processing as _dt_half sbatch so results are directly
# comparable.  The combined output of both RESULT.txt files pins down
# the truncation order.
# -------------------------------------------------------------------------
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_v92_dt_qtr_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_v92_dt_qtr_${SLURM_JOB_ID}.err"

if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}')." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData\|ERROR: NaN" "${JOB_ERR}" 2>/dev/null; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData|ERROR: NaN" "${JOB_ERR}" | head -1)
   echo "FAIL/NAN: run aborted — classifies as CFL violation / bulk-operator issue" > "${RESULT_FILE}"
   echo "      ${ABORT_LINE}" >> "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

V91_PEAK=1.83
GLOBAL_PEAK=0
PER_STATION=""
for sta in flt_0_3 flt_0_7.5 flt_0_12 flt_9_7.5 flt_12_3 flt_12_12 \
           flt_n9_7.5 flt_n12_3 flt_n12_12; do
   STA_FILE="${RESULT_DIR}/tpv102_station_${sta}.dat"
   if [ ! -f "${STA_FILE}" ]; then
      STA_FILE="${RESULT_DIR}/${sta}.dat"
   fi
   if [ ! -f "${STA_FILE}" ]; then
      PER_STATION="${PER_STATION}  ${sta}: missing"$'\n'
      continue
   fi
   PK=$(awk '!/^#/ && NF>=2 {v=($2<0?-$2:$2); if (v>m) m=v} END{print m+0}' "${STA_FILE}")
   PER_STATION="${PER_STATION}  ${sta}: |slip_dip|_peak = ${PK} m"$'\n'
   if awk "BEGIN {exit !( ${PK} > ${GLOBAL_PEAK} )}"; then
      GLOBAL_PEAK="${PK}"
   fi
done

RATIO=$(awk "BEGIN {if (${GLOBAL_PEAK} > 0 && ${V91_PEAK} > 0) printf \"%.3f\", ${V91_PEAK} / ${GLOBAL_PEAK}; else print \"N/A\"}")

{
   echo "TPV102 v9.2.0 Step 1 — dt/4 Frontera result"
   echo "  cfl = 0.125 (v91 was 0.5)"
   echo "  v91 baseline slip_dip peak (hypo) : ${V91_PEAK} m"
   echo "  dt/4 slip_dip peak (global max)   : ${GLOBAL_PEAK} m"
   echo "  reduction factor v91 / dt_quarter : ${RATIO}"
   echo
   echo "Per-station |slip_dip| peaks:"
   printf "%s" "${PER_STATION}"
   echo
   echo "Decision — compare to sibling dt_half RESULT.txt:"
   echo "  (v91 / dt_half) × (dt_half / dt_quarter) ≈ 16×16=256 ==> RK4 truncation (dt^4)"
   echo "  (v91 / dt_half) × (dt_half / dt_quarter) ≈  4× 4= 16 ==> psi split (dt^2)"
   echo "  Both ratios ≈ 1                                     ==> non-truncation → Step 3"
   if awk "BEGIN {exit !( ${GLOBAL_PEAK} < 0.05 )}"; then
      echo "  THIS RUN: global peak < 0.05 m (>36× reduction) ==> strongly favours RK4 truncation (dt^4)"
   elif awk "BEGIN {exit !( ${GLOBAL_PEAK} < 0.2 )}"; then
      echo "  THIS RUN: 0.05 ≤ peak < 0.2 m   (9-36× reduction) ==> psi split (dt^2)"
   elif awk "BEGIN {exit !( ${GLOBAL_PEAK} < 1.5 )}"; then
      echo "  THIS RUN: peak ~ v91 (<2× reduction) ==> NOT time-integration ==> Step 3"
   else
      echo "  THIS RUN: peak ≥ v91 or NaN ==> CFL/bulk-operator issue"
   fi
   echo
   echo "See tpv102_debug_v9.2.0_debug_plan.md PRIORITY NOW Step 1 + Step 1 detail."
} > "${RESULT_FILE}"

cat "${RESULT_FILE}"
