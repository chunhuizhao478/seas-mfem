#!/bin/bash
#SBATCH -J tpv102_v92_f01f02_dev
#SBATCH -o tpv102_v92_f01f02_dev_%j.out
#SBATCH -e tpv102_v92_f01f02_dev_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV102 v9.2.0 plan Step 5 dev-check — coupled RK4 on (Q, psi) post-fix run.
# (tpv102_debug_v9.2.0_debug_plan.md Step 5 F01+F02 verification.)
#
# Identical to v91 12 s baseline (job 7668434) EXCEPT:
#   1. Source tree includes the F01+F02 fix in drivers/tpv102_driver.cpp —
#      classical coupled RK4 on psi (replaces four operator-split
#      UpdateStateAnalytic calls with stage-wise dpsi/dt + RK4 linear
#      combination).  See plan Step 5 detail for the exact scope.
#   2. tfinal 12.0 s -> 4.0 s   (fits the 2 hr/8-node budget).
#      v91 data at t ≈ 4 s:
#        peak |sigma_n - 120 MPa| ≈ 80 MPa   (1.8× Poisson bound)
#        peak |slip_dip|         ≈ 0.26 m   (260× SCEC spec)
#      extracted from tpv102/r_v92_outlier_detection/outlier_report.json
#      (cycles 39-42, same mesh / order / CFL / 400 ranks as this run).
#
# Purpose — first-look Frontera signal on whether F01+F02 alone
# materially reduces the σ_n / slip_dip anomaly.  Decision matrix
# (REVIEW R-V92-H04+H05: all thresholds are ABSOLUTE physics on the
# GLOBAL peak across the 9 SCEC stations.  The v91 reference values
# 80 MPa / 0.26 m quoted in the header come from the VTU per-triangle
# max over ~38 k fault faces — NOT from station .dat peaks — and are
# UPPER BOUNDS on what any station can observe.  Classifying against
# VTU-reference "reduction factors" would be 4-orders-of-magnitude
# biased by the 9-station vs ~38 k-triangle sampling ratio, so the
# classifier below uses absolute thresholds tied to SCEC spec + the
# Poisson bound, not v91-relative reductions):
#
#   | post-fix station peaks at t=4 s                     | Verdict                                 |
#   |-----------------------------------------------------|-----------------------------------------|
#   | global |σ_n−120 MPa|_peak ≤ 45 MPa AND              | F01+F02 CLOSES R-V92 at t=4 s;          |
#   | global |slip_dip|_peak     ≤ 5e-3 m                 | schedule Step 8 (12 s) confirmation.    |
#   | global |σ_n−120 MPa|_peak ≤ 80 MPa AND              | F01+F02 PARTIAL — meaningful reduction; |
#   | global |slip_dip|_peak     ≤ 5e-2 m                 | psi-split contributed but was not       |
#   |                                                     | the sole primary; continue Step 3.      |
#   | otherwise                                           | F01+F02 DID NOT CLOSE — psi-split was   |
#   |                                                     | not primary.  Step 1 dt-halving (on     |
#   |                                                     | the PRE-fix HEAD per R-V92-H03) remains |
#   |                                                     | the cleanest next probe; failing that,  |
#   |                                                     | Step 3 §5.1 bulk SYY.                   |
#   | NaN / abort                                         | Coupled-RK4 stability issue; audit the  |
#   |                                                     | psi rate call or reduce CFL.            |
#
# Note on Step 1 interaction: the existing dt_half/dt_quarter sbatches
# rebuild from branch HEAD, so if F01+F02 is pushed before they are
# submitted, they measure dt scaling on the POST-FIX driver.  That is
# still informative (characterises coupled-RK4 convergence in dt), but
# it is NOT the pre-fix dt scaling the plan Step 1 originally designed
# to measure.  The v91 job-7668434 data (pre-fix cfl=0.5) remains the
# pre-fix reference; this dev sbatch establishes the post-fix cfl=0.5
# counterpart.  The dt_half/dt_quarter runs then answer "does the
# post-fix driver still show dt-dependence?" — a clean residual probe.
#
# PREREQUISITE:
#   F01+F02 source change committed + pushed to origin on the branch
#   Frontera pulls from.  Verification step below greps the driver .o
#   file for the new AgingLawPsi symbol to catch an unpatched build.
#
# Plan references:
#   - debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_debug_plan.md
#     Step 5 (F01+F02 coupled RK4)
#   - tpv102_200m_p1_12.0s_400rank_v91.sbatch    (pre-fix baseline)
#   - tpv102/r_v92_outlier_detection/outlier_report.json (t=4s v91 values)
#   - feedback_sbatch_modules.md
#   - feedback_frontera_approval.md — EXPLICIT APPROVAL REQUIRED.
#
# WALL-CLOCK budget: ~1.3-1.6 hr (matches v92_diag 4 s throughput at
#                    400 ranks on 8 nodes; 2 hr cap leaves margin).
# SU COST          : ~512 SU (8 nodes × 2 hr × 32 SU/node-hr).
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_4.0s_400r_v92_f01f02_dev.sbatch

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

# Full rebuild of the driver and touched TUs (dynamic + io + fault)
# so the pushed F01+F02 patch to drivers/tpv102_driver.cpp is picked up.
rm -f seas_tpv102_driver drivers/tpv102_driver.o dynamic/*.o io/*.o fault/*.o
make seas_tpv102_driver -j"${JOBS}" SEAS_EXTRA_CPPFLAGS="" \
   || { echo "BUILD FAILED"; exit 1; }

# Build verification: the F01+F02 fix uses AgingLawPsi::Rate inside the
# driver's RK4 loop.  A binary lacking that symbol indicates the push/
# checkout cycle did not land the patch — bail early.
if ! nm -C seas_tpv102_driver 2>/dev/null | grep -q "AgingLawPsi"; then
   echo "ERROR: seas_tpv102_driver does NOT contain AgingLawPsi symbol"
   echo "       F01+F02 patch in drivers/tpv102_driver.cpp is not in"
   echo "       the build.  Check git log + rebuild order."
   exit 2
fi
echo "BUILD OK — AgingLawPsi symbol present in driver binary."

RESULT_DIR="tpv102/results_200m_p1_4.0s_400r_v92_f01f02_dev_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"
mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase v9.2.0 Step 5 dev: 400 ranks / 8 nodes / 200 m / P1 / tfinal=4.0 s / cfl=0.5"
echo "Purpose  : F01+F02 coupled-RK4-on-psi first-look verification"
echo "Baseline : v91 job 7668434 at t=4s — |σ_n−120|_peak ≈ 80 MPa, |slip_dip|_peak ≈ 0.26 m"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 4.0 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv \
      --debug-qnorm

# -------------------------------------------------------------------------
# Post-run classification.  Compare post-fix t=4 s peaks against v91
# pre-fix t=4 s (80 MPa σ_n dev, 0.26 m slip_dip) and classify into the
# decision matrix in the header.
# -------------------------------------------------------------------------
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_v92_f01f02_dev_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_v92_f01f02_dev_${SLURM_JOB_ID}.err"

if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}')." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData\|ERROR: NaN" "${JOB_ERR}" 2>/dev/null; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData|ERROR: NaN" "${JOB_ERR}" | head -1)
   echo "FAIL/NAN: run aborted — coupled-RK4 stability issue?  audit Rate() + CFL." > "${RESULT_FILE}"
   echo "      ${ABORT_LINE}" >> "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

# v91 pre-fix t=4 s VTU-peak (per-triangle max over ~38 k fault faces)
# reference values from outlier_report.json cycles 39-42.  These are
# UPPER BOUNDS on station peaks; recorded here for narrative context
# only — the classifier thresholds below are absolute physics
# (Poisson bound + SCEC spec), not v91-relative.  Per REVIEW R-V92-H05.
V91_VTU_PEAK_SIGMA_N_DEV_MPA=80.0
V91_VTU_PEAK_SLIP_DIP_M=0.26
SIGMA_N0_PA=120000000

# Absolute physics thresholds (apply to GLOBAL peaks across stations).
# CLOSES = below Poisson bound + within 5× SCEC spec (1 mm) for slip_dip.
# PARTIAL = still within v91-VTU envelope (80 MPa) AND station-level
# slip_dip ≤ 5e-2 m (50× spec) — a ~5× improvement over v91 VTU peak
# translated to station context.
CLOSES_SIGMA_N_DEV_PA=45000000
CLOSES_SLIP_DIP_M=0.005
PARTIAL_SIGMA_N_DEV_PA=80000000
PARTIAL_SLIP_DIP_M=0.05

# Station .dat columns per TPV102StationWriter::WriteStep
# (tpv102_setup.hpp:305-314):
#   1=t  2=slip1(dip)  3=slip2(strike)  4=V1(dip)  5=V2(strike)
#   6=tau1  7=tau2  8=sigma_n  9=log10_theta
GLOBAL_DIP_PEAK=0
GLOBAL_SN_DEV_PEAK_PA=0
HYPO_DIP_PEAK=0
HYPO_SN_DEV_PEAK_PA=0
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
   # Max |slip_dip| over t
   DIP_PK=$(awk '!/^#/ && NF>=2 {v=($2<0?-$2:$2); if (v>m) m=v} END{print m+0}' "${STA_FILE}")
   # Max |sigma_n - sigma_n0| over t
   SN_DEV_PK=$(awk -v sn0="${SIGMA_N0_PA}" '!/^#/ && NF>=8 {d=$8-sn0; if (d<0) d=-d; if (d>m) m=d} END{print m+0}' "${STA_FILE}")
   PER_STATION="${PER_STATION}  ${sta}: |slip_dip|_peak=${DIP_PK} m, |σ_n−120|_peak=$(awk -v x=${SN_DEV_PK} 'BEGIN{printf \"%.2f\", x/1e6}') MPa"$'\n'

   awk "BEGIN {exit !( ${DIP_PK} > ${GLOBAL_DIP_PEAK} )}" && GLOBAL_DIP_PEAK="${DIP_PK}"
   awk "BEGIN {exit !( ${SN_DEV_PK} > ${GLOBAL_SN_DEV_PEAK_PA} )}" && GLOBAL_SN_DEV_PEAK_PA="${SN_DEV_PK}"
   if [ "${sta}" = "flt_0_7.5" ]; then
      HYPO_DIP_PEAK="${DIP_PK}"
      HYPO_SN_DEV_PEAK_PA="${SN_DEV_PK}"
   fi
done

GLOBAL_SN_DEV_PEAK_MPA=$(awk -v x=${GLOBAL_SN_DEV_PEAK_PA} 'BEGIN{printf "%.2f", x/1e6}')
HYPO_SN_DEV_PEAK_MPA=$(awk -v x=${HYPO_SN_DEV_PEAK_PA} 'BEGIN{printf "%.2f", x/1e6}')

V_MAX_PEAK=$(grep "V_max" "${JOB_LOG}" \
             | awk '{for(i=1;i<=NF;i++) if ($i ~ /^[0-9.e+-]+$/) v=$i; print v}' \
             | sort -g | tail -1)
FINAL_WATCH=$(grep "\[qnorm:watch\]" "${JOB_LOG}" | tail -1)

{
   echo "TPV102 v9.2.0 Step 5 F01+F02 dev-check result (t=4 s)"
   echo "================================================================"
   echo "  v91 pre-fix context (narrative only — VTU per-triangle max,"
   echo "  NOT directly comparable to station peaks — see R-V92-H05):"
   echo "    v91 VTU peak |σ_n − 120 MPa| ≈ ${V91_VTU_PEAK_SIGMA_N_DEV_MPA} MPa"
   echo "    v91 VTU peak |slip_dip|      ≈ ${V91_VTU_PEAK_SLIP_DIP_M} m"
   echo ""
   echo "  Post-fix station peaks (this run, at t=4 s):"
   echo "    peak |σ_n − 120 MPa| (global 9 stations) = ${GLOBAL_SN_DEV_PEAK_MPA} MPa"
   echo "    peak |slip_dip|      (global 9 stations) = ${GLOBAL_DIP_PEAK} m"
   echo "    peak |σ_n − 120 MPa| (hypo flt_0_7.5)    = ${HYPO_SN_DEV_PEAK_MPA} MPa"
   echo "    peak |slip_dip|      (hypo flt_0_7.5)    = ${HYPO_DIP_PEAK} m"
   echo "    peak V_max (global)                      = ${V_MAX_PEAK}"
   echo ""
   echo "  Absolute-physics thresholds (Poisson bound + SCEC spec):"
   echo "    CLOSES  : |σ_n−120|_peak ≤ $(awk -v x=${CLOSES_SIGMA_N_DEV_PA} 'BEGIN{printf \"%.0f\", x/1e6}') MPa AND |slip_dip|_peak ≤ ${CLOSES_SLIP_DIP_M} m"
   echo "    PARTIAL : |σ_n−120|_peak ≤ $(awk -v x=${PARTIAL_SIGMA_N_DEV_PA} 'BEGIN{printf \"%.0f\", x/1e6}') MPa AND |slip_dip|_peak ≤ ${PARTIAL_SLIP_DIP_M} m"
   echo ""
   echo "Per-station peaks:"
   printf "%s" "${PER_STATION}"
   echo ""
   echo "Decision (matrix in header — uses GLOBAL peaks consistently per R-V92-H04):"
   if awk "BEGIN {exit !( ${GLOBAL_SN_DEV_PEAK_PA} <= ${CLOSES_SIGMA_N_DEV_PA} && ${GLOBAL_DIP_PEAK} <= ${CLOSES_SLIP_DIP_M} )}"; then
      echo "  F01+F02 CLOSES R-V92 at t=4 s (below Poisson bound AND slip_dip within 5× spec)"
      echo "  → schedule Step 8 (12 s / 200 m / 400 ranks) confirmation run."
   elif awk "BEGIN {exit !( ${GLOBAL_SN_DEV_PEAK_PA} <= ${PARTIAL_SIGMA_N_DEV_PA} && ${GLOBAL_DIP_PEAK} <= ${PARTIAL_SLIP_DIP_M} )}"; then
      echo "  F01+F02 PARTIAL — station peaks within |σ_n−120| ≤ $(awk -v x=${PARTIAL_SIGMA_N_DEV_PA} 'BEGIN{printf \"%.0f\", x/1e6}') MPa AND"
      echo "  |slip_dip| ≤ ${PARTIAL_SLIP_DIP_M} m.  psi-split contributed but was not the sole"
      echo "  primary; continue Step 3 (§5.1 bulk SYY) for the residual amplifier."
   else
      echo "  F01+F02 DID NOT CLOSE — station peaks exceed PARTIAL thresholds."
      echo "  psi-split was not primary.  Next: per R-V92-H03, wait for PRE-fix"
      echo "  Step 1 dt-halving results; if those also don't close, go to Step 3"
      echo "  §5.1 bulk SYY probe."
   fi
   echo ""
   if [ -n "${FINAL_WATCH}" ]; then
      echo "  Final qnorm watch: ${FINAL_WATCH}"
   fi
   echo ""
   echo "See tpv102_debug_v9.2.0_debug_plan.md Step 5 for the full fix scope."
} > "${RESULT_FILE}"

cat "${RESULT_FILE}"
