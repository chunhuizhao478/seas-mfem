#!/bin/bash
#SBATCH -J tpv102_v92_postfix_12s_400r
#SBATCH -o tpv102_v92_postfix_12s_400r_%j.out
#SBATCH -e tpv102_v92_postfix_12s_400r_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 05:00:00
#SBATCH -A EAR20006

# TPV102 v9.2.0 Step 8 — 12 s post-fix confirmation run.
# Builds on v91 12 s baseline (job 7668434) with three driver-level
# fixes applied:
#
#   1. F01+F02 coupled RK4 on (Q, psi) — stage-wise dpsi/dt sampled via
#      AgingLawPsi::Rate, combined with classical (1,2,2,1)/6 Butcher
#      weights.  Replaces the pre-fix operator-split UpdateStateAnalytic
#      pattern (O(dt^2) coupling) with O(dt^4) coupling matching Q's
#      RK4 order.
#
#   2. R-V92-K01 endpoint re-evaluation — after RK4 combines Q(t+dt)
#      and psi(t+dt), a fresh wave.Mult(Q, k_endpoint) runs at the
#      endpoint state to overwrite dof_data.{V1, V2, slip_rate,
#      tau1_corr, tau2_corr, sigma_n_corr} with values instantaneously
#      self-consistent at t_{n+1}.  Replaces the pre-fix
#      Butcher-weighted mean of stage-i snapshots (which labelled the
#      sample time as t_{n+1} but held values at t_n + dt/2 — a
#      half-step phase lag).  Cost: 1 extra Mult per step (~25%).
#
#   3. R-V92-E02 `_k4` VTU CellData — fault-surface VTU now emits
#      `slip_rate_{dip,strike}_k4`, `traction_{dip,strike}_k4`,
#      `normal_stress_k4` alongside the averaged/endpoint fields.
#      Diagnostic only; `<field>_k4 − <field>` quantifies the
#      RK4-averaging residual per triangle in ParaView.
#
# Plan-closure gate (tpv102_debug_v9.2.0_debug_plan.md §20.5):
#   |σ_n − 120 MPa| ≤ 1 MPa  AND  |slip_dip| ≤ 10^-3 m at every one
#   of the 9 SCEC stations.
#
# v91 pre-fix baseline (job 7668434, tfinal = 12 s):
#   global |σ_n − 120 MPa|_peak ≈ 218 MPa  (VTU per-triangle max)
#   global |slip_dip|_peak      ≈ 1.83 m
#
# Classifier decision matrix (GLOBAL station peaks, absolute physics):
#
#   | Verdict          | σ_n threshold        | slip_dip threshold  |
#   |------------------|----------------------|---------------------|
#   | CLOSES (spec)    | ≤ 1 MPa every sta    | ≤ 10^-3 m every sta |
#   | CLOSES (loose)   | ≤ 45 MPa global      | ≤ 5·10^-3 m global  |
#   | PARTIAL          | ≤ 80 MPa global      | ≤ 5·10^-2 m global  |
#   | NO-CLOSE         | anything else        |                     |
#   | NAN/ABORT        | MFEM abort observed  |                     |
#
# CLOSES (spec) closes v9.2.0 per plan §20.5 → schedule
# tpv102_debug_v9.2.0_fix.md + tpv102_debug_v9.2.0_check.md.
# CLOSES (loose) means within Poisson bound + ~5× spec — F01+F02+K01
# eliminated the anomaly; continuing investigation is diminishing
# returns.  PARTIAL means residual ~v91-scale σ_n drift but O(50×)
# slip_dip reduction — psi-split + K01 were contributors but not the
# sole primary; continue Step 3 §5.1 bulk SYY.  NO-CLOSE means neither
# F01+F02 nor K01 was the primary mechanism; re-open with H-V92-G / V / U.
#
# PREREQUISITE:
#   Source tree must contain:
#     - drivers/tpv102_driver.cpp : F01+F02 coupled RK4 + R-V92-K01
#       endpoint re-eval blocks
#     - dynamic/fault_face_flux.cpp : psi-purity invariant
#     - io/paraview_output.hpp : `_k4` CellData emission path
#   Build verification below greps the driver binary for AgingLawPsi
#   (confirms F01+F02 linked) and checks for the endpoint re-eval
#   comment marker (confirms K01 compiled in).
#
# Plan references:
#   - debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_debug_plan.md
#     Step 8 (12 s confirmation), §20.5 closure gate
#   - tpv102_200m_p1_12.0s_400rank_v91.sbatch       (pre-fix baseline)
#   - tpv102_200m_p1_4.0s_400r_v92_f01f02_dev.sbatch (4 s dev counterpart)
#   - feedback_sbatch_modules.md
#   - feedback_frontera_approval.md — EXPLICIT APPROVAL REQUIRED.
#
# WALL-CLOCK budget: ~4.5 hr (v91 4 hr + ~25% endpoint re-eval overhead).
# 5 hr cap leaves margin.
# SU COST          : ~5 hr × 8 nodes × 32 SU/node-hr ≈ 1280 SU.
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_12.0s_400r_v92_postfix.sbatch

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

# Full rebuild so the v9.2.0 driver + flux + IO changes are picked up.
rm -f seas_tpv102_driver drivers/tpv102_driver.o dynamic/*.o io/*.o fault/*.o
make seas_tpv102_driver -j"${JOBS}" SEAS_EXTRA_CPPFLAGS="" \
   || { echo "BUILD FAILED"; exit 1; }

# Build verification 1 — F01+F02: AgingLawPsi::Rate must be linked.
if ! nm -C seas_tpv102_driver 2>/dev/null | grep -q "AgingLawPsi"; then
   echo "ERROR: seas_tpv102_driver does NOT contain AgingLawPsi symbol."
   echo "       F01+F02 patch in drivers/tpv102_driver.cpp is not in the"
   echo "       build.  Check git log + rebuild order."
   exit 2
fi
echo "BUILD OK (1/2) — AgingLawPsi symbol present (F01+F02 linked)."

# Build verification 2 — R-V92-K01: endpoint re-eval marker in source.
# The endpoint re-eval is a call site, not a new symbol, so we grep the
# driver source for the R-V92-K01 marker that names the fix.  This
# confirms the checked-out tree contains the endpoint re-eval block.
if ! grep -q "R-V92-K01" drivers/tpv102_driver.cpp 2>/dev/null; then
   echo "ERROR: R-V92-K01 endpoint re-evaluation marker missing from"
   echo "       drivers/tpv102_driver.cpp.  The post-step re-evaluation"
   echo "       of dof_data at Q(t+dt) is the second prerequisite of"
   echo "       this confirmation run."
   exit 2
fi
echo "BUILD OK (2/2) — R-V92-K01 endpoint re-eval marker present."

RESULT_DIR="tpv102/results_200m_p1_12.0s_400r_v92_postfix_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"
mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase v9.2.0 Step 8: 400 ranks / 8 nodes / 200 m / P1 / tfinal=12.0 s / cfl=0.5"
echo "Purpose  : F01+F02 + R-V92-K01 endpoint re-eval 12 s confirmation"
echo "Baseline : v91 job 7668434 — |σ_n−120|_peak ≈ 218 MPa, |slip_dip|_peak ≈ 1.83 m (VTU)"
echo "Closure  : |σ_n−120 MPa|_peak ≤ 1 MPa AND |slip_dip|_peak ≤ 10^-3 m at every station"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 12.0 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv \
      --debug-qnorm

# -------------------------------------------------------------------------
# Post-run classification.  Apply the §20.5 closure gate, then the tiered
# loose/partial/no-close decision matrix from the header.
# -------------------------------------------------------------------------
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_v92_postfix_12s_400r_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_v92_postfix_12s_400r_${SLURM_JOB_ID}.err"

if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}')." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData\|ERROR: NaN" "${JOB_ERR}" 2>/dev/null; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData|ERROR: NaN" "${JOB_ERR}" | head -1)
   echo "FAIL/NAN: run aborted — coupled-RK4 stability or endpoint-reeval issue?" > "${RESULT_FILE}"
   echo "      ${ABORT_LINE}" >> "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

SIGMA_N0_PA=120000000

# Tier thresholds (absolute physics, applied to GLOBAL station peaks).
CLOSES_SPEC_SIGMA_N_DEV_PA=1000000       # 1   MPa
CLOSES_SPEC_SLIP_DIP_M=0.001             # 1 mm
CLOSES_LOOSE_SIGMA_N_DEV_PA=45000000     # 45  MPa (Poisson bound)
CLOSES_LOOSE_SLIP_DIP_M=0.005            # 5 mm
PARTIAL_SIGMA_N_DEV_PA=80000000          # 80  MPa (v91 @ t=4 s VTU peak)
PARTIAL_SLIP_DIP_M=0.05                  # 50 mm

# Station .dat columns per TPV102StationWriter::WriteStep
#   1=t  2=slip1(dip)  3=slip2(strike)  4=V1(dip)  5=V2(strike)
#   6=tau1  7=tau2  8=sigma_n  9=log10_theta
GLOBAL_DIP_PEAK=0
GLOBAL_SN_DEV_PEAK_PA=0
HYPO_DIP_PEAK=0
HYPO_SN_DEV_PEAK_PA=0
PER_STATION=""
SLIP_COUNT=0
SLIPPED_STATIONS=""
SILENT_STATIONS=""
SPEC_PASS_STATIONS=""
SPEC_FAIL_STATIONS=""

for sta in flt_0_3 flt_0_7.5 flt_0_12 flt_9_7.5 flt_12_3 flt_12_12 \
           flt_n9_7.5 flt_n12_3 flt_n12_12; do
   STA_FILE="${RESULT_DIR}/tpv102_station_${sta}.dat"
   if [ ! -f "${STA_FILE}" ]; then
      STA_FILE="${RESULT_DIR}/${sta}.dat"
   fi
   if [ ! -f "${STA_FILE}" ]; then
      PER_STATION="${PER_STATION}  ${sta}: missing"$'\n'
      SILENT_STATIONS="${SILENT_STATIONS} ${sta}(missing)"
      continue
   fi
   DIP_PK=$(awk '!/^#/ && NF>=2 {v=($2<0?-$2:$2); if (v>m) m=v} END{print m+0}' "${STA_FILE}")
   SN_DEV_PK=$(awk -v sn0="${SIGMA_N0_PA}" '!/^#/ && NF>=8 {d=$8-sn0; if (d<0) d=-d; if (d>m) m=d} END{print m+0}' "${STA_FILE}")
   V2=$(awk '!/^#/ && NF>=5 {v=$5} END {print (v=="" ? "0" : v)}' "${STA_FILE}")
   DIP_PK_FMT=$(awk -v x=${DIP_PK} 'BEGIN{printf "%.4e", x}')
   SN_DEV_PK_MPA=$(awk -v x=${SN_DEV_PK} 'BEGIN{printf "%.3f", x/1e6}')
   PER_STATION="${PER_STATION}  ${sta}: |slip_dip|_peak=${DIP_PK_FMT} m, |σ_n−120|_peak=${SN_DEV_PK_MPA} MPa"$'\n'

   awk "BEGIN {exit !( ${DIP_PK} > ${GLOBAL_DIP_PEAK} )}" && GLOBAL_DIP_PEAK="${DIP_PK}"
   awk "BEGIN {exit !( ${SN_DEV_PK} > ${GLOBAL_SN_DEV_PEAK_PA} )}" && GLOBAL_SN_DEV_PEAK_PA="${SN_DEV_PK}"
   if [ "${sta}" = "flt_0_7.5" ]; then
      HYPO_DIP_PEAK="${DIP_PK}"
      HYPO_SN_DEV_PEAK_PA="${SN_DEV_PK}"
   fi

   if awk "BEGIN {exit !( ($V2)*($V2) > 1e-12 )}"; then
      SLIP_COUNT=$((SLIP_COUNT + 1))
      SLIPPED_STATIONS="${SLIPPED_STATIONS} ${sta}"
   else
      SILENT_STATIONS="${SILENT_STATIONS} ${sta}"
   fi

   # Per-station spec-tier pass check (for §20.5 closure gate).
   if awk "BEGIN {exit !( ${SN_DEV_PK} <= ${CLOSES_SPEC_SIGMA_N_DEV_PA} && ${DIP_PK} <= ${CLOSES_SPEC_SLIP_DIP_M} )}"; then
      SPEC_PASS_STATIONS="${SPEC_PASS_STATIONS} ${sta}"
   else
      SPEC_FAIL_STATIONS="${SPEC_FAIL_STATIONS} ${sta}"
   fi
done

GLOBAL_SN_DEV_PEAK_MPA=$(awk -v x=${GLOBAL_SN_DEV_PEAK_PA} 'BEGIN{printf "%.3f", x/1e6}')
HYPO_SN_DEV_PEAK_MPA=$(awk -v x=${HYPO_SN_DEV_PEAK_PA} 'BEGIN{printf "%.3f", x/1e6}')

V_MAX_PEAK=$(grep "V_max" "${JOB_LOG}" \
             | awk '{for(i=1;i<=NF;i++) if ($i ~ /^[0-9.e+-]+$/) v=$i; print v}' \
             | sort -g | tail -1)
FINAL_WATCH=$(grep "\[qnorm:watch\]" "${JOB_LOG}" | tail -1)

# Count spec-pass stations (need "9" for plan-closure gate).
SPEC_PASS_COUNT=$(echo "${SPEC_PASS_STATIONS}" | awk '{print NF}')

{
   echo "TPV102 v9.2.0 Step 8 — F01+F02 + R-V92-K01 confirmation (t=12 s)"
   echo "================================================================"
   echo "  v91 pre-fix baseline (VTU per-triangle max, job 7668434):"
   echo "    |σ_n − 120 MPa|_peak ≈ 218 MPa"
   echo "    |slip_dip|_peak      ≈ 1.83 m"
   echo ""
   echo "  Post-fix station peaks (this run, 0 ≤ t ≤ 12 s):"
   echo "    peak |σ_n − 120 MPa| (global 9 stations) = ${GLOBAL_SN_DEV_PEAK_MPA} MPa"
   echo "    peak |slip_dip|      (global 9 stations) = ${GLOBAL_DIP_PEAK} m"
   echo "    peak |σ_n − 120 MPa| (hypo flt_0_7.5)    = ${HYPO_SN_DEV_PEAK_MPA} MPa"
   echo "    peak |slip_dip|      (hypo flt_0_7.5)    = ${HYPO_DIP_PEAK} m"
   echo "    peak V_max (global)                      = ${V_MAX_PEAK}"
   echo ""
   echo "  Stations that slipped (|V2| > 1e-6): ${SLIP_COUNT} / 9"
   echo "    Slipped:${SLIPPED_STATIONS}"
   echo "    Silent :${SILENT_STATIONS}"
   echo ""
   echo "  Per-station plan-closure gate (|σ_n−120 MPa| ≤ 1 MPa AND |slip_dip| ≤ 1e-3 m):"
   echo "    PASS (${SPEC_PASS_COUNT}):${SPEC_PASS_STATIONS}"
   echo "    FAIL:${SPEC_FAIL_STATIONS}"
   echo ""
   echo "Per-station peaks:"
   printf "%s" "${PER_STATION}"
   echo ""
   echo "Decision (§20.5 closure gate is CLOSES-spec, 9/9 per-station PASS):"
   if [ "${SPEC_PASS_COUNT}" = "9" ]; then
      echo "  CLOSES-spec — 9/9 stations meet |σ_n−120 MPa| ≤ 1 MPa AND"
      echo "                |slip_dip| ≤ 10^-3 m.  R-V92 CLOSED per plan §20.5."
      echo "  → write tpv102_debug_v9.2.0_fix.md and open"
      echo "    tpv102_debug_v9.2.0_check.md."
   elif awk "BEGIN {exit !( ${GLOBAL_SN_DEV_PEAK_PA} <= ${CLOSES_LOOSE_SIGMA_N_DEV_PA} && ${GLOBAL_DIP_PEAK} <= ${CLOSES_LOOSE_SLIP_DIP_M} )}"; then
      echo "  CLOSES-loose — global peaks within Poisson bound (45 MPa) AND"
      echo "                 |slip_dip| ≤ 5·10^-3 m (5× spec).  F01+F02 + K01"
      echo "                 eliminated the dominant anomaly; residual is"
      echo "                 near-physical.  Discuss with user whether to"
      echo "                 accept as closure or continue Step 3 bulk SYY."
   elif awk "BEGIN {exit !( ${GLOBAL_SN_DEV_PEAK_PA} <= ${PARTIAL_SIGMA_N_DEV_PA} && ${GLOBAL_DIP_PEAK} <= ${PARTIAL_SLIP_DIP_M} )}"; then
      echo "  PARTIAL — station peaks |σ_n−120 MPa| ≤ 80 MPa AND |slip_dip| ≤ 5·10^-2 m."
      echo "  F01+F02+K01 contributed substantial reduction but were not"
      echo "  the sole primary; continue Step 3 (§5.1 bulk SYY) for the"
      echo "  residual amplifier.  Consider also Step 7 [FAULT-INIT-V1] probe."
   else
      echo "  NO-CLOSE — post-fix peaks still exceed PARTIAL tier."
      echo "  F01+F02+K01 was not the primary mechanism.  Re-open with"
      echo "  H-V92-V (abs-BC per-channel phase) / H-V92-U (interior-fault"
      echo "  multi-step) / H-V92-W (init-V1).  Compare station trajectories"
      echo "  against v91 job 7668434 to localise which hypothesis fits best."
   fi
   echo ""
   if [ -n "${FINAL_WATCH}" ]; then
      echo "  Final qnorm watch: ${FINAL_WATCH}"
   fi
   echo ""
   echo "See tpv102_debug_v9.2.0_debug_plan.md §20.5 for closure gate"
   echo "and §20.3 for the post-closure action path."
} > "${RESULT_FILE}"

cat "${RESULT_FILE}"
