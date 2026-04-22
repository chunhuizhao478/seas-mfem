#!/bin/bash
#SBATCH -J tpv102_v92_dt_half
#SBATCH -o tpv102_v92_dt_half_%j.out
#SBATCH -e tpv102_v92_dt_half_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 08:00:00
#SBATCH -A EAR20006

# TPV102 v9.2.0 plan Step 1 — dt-halving decisive test (HIGHEST PRIORITY)
# (tpv102_debug_v9.2.0_debug_plan.md Step 1, PRIORITY NOW rev-3h).
#
# Identical to the v91 12 s baseline (job 7668434) EXCEPT `--cfl 0.25`
# instead of the default `--cfl 0.5`.  Since the driver computes
# `cfl_internal = cfl_factor / (3*(2*order+1))` and `dt_cfl` scales
# linearly with cfl_factor (drivers/tpv102_driver.cpp:325-326), this
# halves dt relative to the v91 baseline.  All physics parameters,
# mesh, and output schedule are preserved so the comparison is 1:1.
#
# Purpose — classify R-V92 primary mechanism (plan §Step 1):
#
#   | slip_dip peak vs v91 (1.83 m)  | Classification                             | Next action       |
#   |--------------------------------|--------------------------------------------|-------------------|
#   | ~16× smaller (~0.11 m)         | Time-integration bug (4th-order truncate)  | Apply F01+F02     |
#   | ~8×  smaller (~0.23 m)         | psi splitting (F02 3rd-order)              | Integrate psi     |
#   | Unchanged (~1.8 m)             | NOT time-integration                       | → plan Step 3     |
#   | Larger / NaN                   | CFL violation / bulk-operator issue        | Audit & reduce CFL |
#
# Expected fixed-truncation scaling of RK4: slip_dip ∝ dt^4 if dominated
# by RK4 truncation; dt^2 if dominated by the psi operator split.  The
# classification table encodes both.
#
# PREREQUISITE — source tree matches the one used for v91 job 7668434
# (same post-v9.1.0 head).  No v9.2.0 C0 source edits required for
# this run; it is a pure time-step variation of the baseline build.
# If a v9.2.0 diagnostic patch (e.g. SEAS_DIAG_FAULT_SIGMA) has been
# merged to the branch, that is fine — this script does NOT compile
# with `-DSEAS_DIAG_FAULT_SIGMA`, so the diagnostic blocks stay dead.
#
# Plan references:
#   - debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_debug_plan.md
#     PRIORITY NOW Step 1; §Step 1 detail; REVIEW round-5 R-V92-G01
#   - tpv102_200m_p1_12.0s_400rank_v91.sbatch (baseline)
#   - feedback_sbatch_modules.md
#   - feedback_frontera_approval.md — EXPLICIT APPROVAL REQUIRED before
#     sbatch.
#
# WALL-CLOCK budget: ~8 hr (2× v91 baseline because nsteps doubles).
# SU COST          : ~2048 SU (8 nodes × 8 hr × 32 SU/node-hr).
#                    NOTE: plan header estimates 20 SU for the pair,
#                    but that is for the classifier only; production
#                    12 s / 200 m / 400 ranks at dt/2 costs ~2 ×
#                    v91's 1024 SU.  Override with a shorter --tfinal
#                    if budget constrained — `tfinal 4.0` still
#                    captures the rupture-tip growth through ~4 s,
#                    which §18 shows already exhibits the σ_n / slip
#                    anomaly at peak.  See Decision section at end.
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_12.0s_400r_v92_dt_half.sbatch

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
# Force full rebuild so the v91 post-build tree is definitely picked up.
rm -f seas_tpv102_driver drivers/tpv102_driver.o dynamic/*.o io/*.o fault/*.o
make seas_tpv102_driver -j"${JOBS}" SEAS_EXTRA_CPPFLAGS="" \
   || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_200m_p1_12.0s_400r_v92_dt_half_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase v9.2.0 Step 1 (dt/2): 400 ranks / 8 nodes / 200 m / P1 / tfinal=12.0 s / cfl=0.25"
echo "Purpose  : R-V92 primary-mechanism classification via RK4 truncation scaling"
echo "Baseline : v91 job 7668434 (cfl=0.5, slip_dip peak 1.83 m)"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.25 \
      --tfinal 12.0 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv \
      --debug-qnorm

# -------------------------------------------------------------------------
# Post-run classification for R-V92 Step 1 (dt-halving decision table).
# Extracts slip_dip (= slip1) peak at flt_0_7.5 and compares to v91's 1.83 m.
# -------------------------------------------------------------------------
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_v92_dt_half_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_v92_dt_half_${SLURM_JOB_ID}.err"

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

# Station .dat columns: 1=t, 2=slip1(dip), 3=slip2(strike), 4=V1(dip),
#                       5=V2(strike), 6=tau1, 7=tau2, 8=sigma_n, 9=psi
V91_PEAK=1.83
# Scan the nine SCEC stations for the maximum |slip_dip| — rupture may
# be asymmetric across stations post-cfl-change, so per-station peaks
# AND a global peak both go into the log.
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
   echo "TPV102 v9.2.0 Step 1 — dt/2 Frontera result"
   echo "  cfl = 0.25 (v91 was 0.5)"
   echo "  v91 baseline slip_dip peak (hypo) : ${V91_PEAK} m"
   echo "  dt/2 slip_dip peak (global max)   : ${GLOBAL_PEAK} m"
   echo "  reduction factor v91 / dt_half    : ${RATIO}"
   echo
   echo "Per-station |slip_dip| peaks:"
   printf "%s" "${PER_STATION}"
   echo
   echo "Decision (plan Step 1 table):"
   # dt/2 ⇒ RK4 truncation scales by 1/16, psi-split by 1/4
   if awk "BEGIN {exit !( ${GLOBAL_PEAK} < 0.2 )}"; then
      echo "  global peak < 0.2 m  (> 9× reduction) ==> RK4 truncation PRIMARY ==> Step 5 F01+F02"
   elif awk "BEGIN {exit !( ${GLOBAL_PEAK} < 0.5 )}"; then
      echo "  0.2 m ≤ peak < 0.5 m  (4-9× reduction) ==> psi splitting PRIMARY ==> F02 alone"
   elif awk "BEGIN {exit !( ${GLOBAL_PEAK} < 1.5 )}"; then
      echo "  peak ~ v91 (<2× reduction) ==> NOT time-integration ==> go to Step 3 (§5.1 bulk probe)"
   else
      echo "  peak ≥ v91 or NaN ==> CFL/bulk-operator issue ==> audit flux.Interior + reduce CFL"
   fi
   echo
   echo "See tpv102_debug_v9.2.0_debug_plan.md PRIORITY NOW Step 1 + Step 1 detail."
} > "${RESULT_FILE}"

cat "${RESULT_FILE}"
