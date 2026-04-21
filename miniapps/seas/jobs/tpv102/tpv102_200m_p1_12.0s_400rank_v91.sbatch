#!/bin/bash
#SBATCH -J tpv102_v91_12s_400r
#SBATCH -o tpv102_v91_12s_400r_%j.out
#SBATCH -e tpv102_v91_12s_400r_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 04:00:00
#SBATCH -A EAR20006

# TPV102 v9.1.0 Phase 9.1 — extended 12 s rerun after Pelties-9
# (v9.0.0) + VTU-CellData / dip-orthonormality (v9.1.0) fixes.
#
# Decides R-V91-B: does the hypocentre dip-channel transient
# (flt_0_7.5: |slip_dip| ~ 5e-4 m, |shear_dip| ~ 0.3 MPa at t ~ 1.5 s)
# decay into the SeisSol reference envelope (|slip_dip| <= 2e-3 m,
# |v_dip| RMS within 2x reference) by t ~ 5 s?  If yes, H-V91-B3
# (friction-iteration transient sensitivity) is the mechanism and
# R-V91-B closes with no further code change.  If no, open v9.1.1
# targeting H-V91-B1/B2 with production data.
#
# Secondary: R-V91-A VTU-speckle visual confirmation.  The v9.1.0 C0
# (R-001) CellData rewrite must render speckle-free fault-surface
# VTU from t=0 onward.  Any residual element-scale colour variation
# would be genuine DG-jump (H-V91-A1), not speckle.
#
# Configuration: match job 7667881 exactly — 200 m mesh, P1, 400 ranks
# on 8 nodes — with tfinal extended from 2.0 s to 12.0 s.  No DIAG
# flags; keep output volume bounded (~3.6 GB VTU + ~50 MB stations).
# V_max halt at 400 m/s (driver default filename suffix vmax_halt).
#
# Plan references:
#   - debug_document/tpv102_debug_document/tpv102_debug_v9.1.0_debug_plan.md
#     §6 Frontera 12 s re-run specification
#   - v9.1.0_fix.md — C0/C1/C2 changes that must be in-tree on the build
#   - feedback_sbatch_modules.md — module list + LD_LIBRARY_PATH template
#   - feedback_frontera_approval.md — submission approved 2026-04-20
#
# WALL-CLOCK budget: ~4 hr (~6x the 2.0 s job 7667881 wallclock).
# SU COST: ~1024 SU (8 nodes x 4 hr x 32 SU/node-hr, budget).
#
# Usage:
#   sbatch jobs/tpv102/tpv102_200m_p1_12.0s_400rank_v91.sbatch
#
# Submission gated on tree state: v9.1.0 C0 + C1 + C2 must be committed
# and pushed to origin so the Frontera checkout picks them up.  Confirm
# via the branch/commit echo in the build log below.

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
# Force full rebuild so v9.1.0 C0 source changes in paraview_output.hpp
# and fault_basis.hpp are unconditionally picked up.
rm -f seas_tpv102_driver drivers/tpv102_driver.o dynamic/*.o io/*.o fault/*.o
make seas_tpv102_driver -j"${JOBS}" SEAS_EXTRA_CPPFLAGS="" \
   || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_200m_p1_12.0s_400r_v91_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase 9.1: 400 ranks / 8 nodes / 200 m / P1 / tfinal=12.0 s"
echo "Purpose  : R-V91-B dip-transient decay + R-V91-A speckle visual"

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
# Post-run classification for R-V91-B (dip channel decay).
# -------------------------------------------------------------------------
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_v91_12s_400r_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_v91_12s_400r_${SLURM_JOB_ID}.err"
if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}'). Hint: SLURM_SUBMIT_DIR='${SLURM_SUBMIT_DIR:-<unset>}'." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

# Early bailout for MFEM abort.
if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData" "${JOB_ERR}" 2>/dev/null; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData" "${JOB_ERR}" | head -1)
   echo "FAIL: MFEM abort during Phase 9.1 run" > "${RESULT_FILE}"
   echo "      ${ABORT_LINE}" >> "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

# Station-by-station slip check on the 9 standard stations (same panel
# as job 7667881) — for v9.1.0 we ALSO compute per-station dip peak and
# tail magnitudes so the decision table can classify R-V91-B directly.
#
# Station .dat columns (per tpv102_setup.hpp::TPV102StationWriter):
#   1=t, 2=slip1(dip), 3=slip2(strike), 4=V1(dip), 5=V2(strike),
#   6=tau1(dip), 7=tau2(strike), 8=sigma_n, 9=psi
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

# R-V91-B metrics at flt_0_7.5 (hypocentre):
#   slip_dip_peak      = max_t |slip1|
#   slip_dip_tail      = |slip1| at the LAST recorded t (>= 5 s)
#   v_dip_rms_tail     = RMS of V1 over t >= 5 s
HYPO_FILE="${RESULT_DIR}/flt_0_7.5.dat"
DIP_PEAK=""
DIP_TAIL=""
V_DIP_RMS_TAIL=""
if [ -f "${HYPO_FILE}" ]; then
   DIP_PEAK=$(awk '!/^#/ && NF>=2 {v=($2<0?-$2:$2); if (v>m) m=v} END{print m+0}' "${HYPO_FILE}")
   DIP_TAIL=$(awk '!/^#/ && NF>=2 {v=$2} END{print (v==""?0:(v<0?-v:v))}' "${HYPO_FILE}")
   V_DIP_RMS_TAIL=$(awk '!/^#/ && NF>=4 && $1>=5.0 { s+=$4*$4; n+=1 } END { if (n>0) print sqrt(s/n); else print "NA" }' "${HYPO_FILE}")
fi

# V_max timeline
V_MAX_INIT=$(grep "V_max" "${JOB_LOG}" | head -1)
V_MAX_FINAL=$(grep "V_max" "${JOB_LOG}" | tail -1)
V_MAX_PEAK=$(grep "V_max" "${JOB_LOG}" \
             | awk '{for(i=1;i<=NF;i++) if ($i ~ /^[0-9.e+-]+$/) v=$i; print v}' \
             | sort -g | tail -1)
FINAL_WATCH=$(grep "\[qnorm:watch\]" "${JOB_LOG}" | tail -1)

{
   echo "Phase 9.1 (v9.1.0, 400-rank / tfinal=12.0 s) score: ${SLIP_COUNT} / 9 stations slipped"
   echo "  Slipped:${SLIPPED_STATIONS}"
   echo "  Silent :${SILENT_STATIONS}"
   if [ -n "${V_MAX_INIT}" ]; then
      echo "  First  V_max line: ${V_MAX_INIT}"
   fi
   if [ -n "${V_MAX_PEAK}" ]; then
      echo "  Peak   V_max value: ${V_MAX_PEAK}"
   fi
   if [ -n "${V_MAX_FINAL}" ]; then
      echo "  Final  V_max line: ${V_MAX_FINAL}"
   fi
   if [ -n "${FINAL_WATCH}" ]; then
      echo "  Final qnorm watch: ${FINAL_WATCH}"
   fi
   echo "  Hypo dip peak  (|slip_dip| max over all t): ${DIP_PEAK}"
   echo "  Hypo dip tail  (|slip_dip| at last sample): ${DIP_TAIL}"
   echo "  Hypo V_dip RMS (t >= 5 s)                 : ${V_DIP_RMS_TAIL}"

   # R-V91-B pass bar per plan §6.3:
   #   pass-1: |slip_dip| <= 2e-3 m for t >= 5 s  (tail within 2x reference)
   # Pre-v9.1.0 (2 s): peak = O(5e-4 m).  Transient only if tail << peak.
   if [ -n "${DIP_TAIL}" ] && awk "BEGIN {exit !( ${DIP_TAIL} <= 2e-3 )}"; then
      echo "  R-V91-B-pass-1: |slip_dip| tail <= 2e-3 m  ==> PASS"
      RVB1=PASS
   else
      echo "  R-V91-B-pass-1: |slip_dip| tail > 2e-3 m   ==> FAIL (open v9.1.1)"
      RVB1=FAIL
   fi

   # Decision echo — reproduce plan §6.3 A0-like table.
   if [ "${SLIP_COUNT}" -eq 0 ]; then
      echo "  Decision: 0/9 slipped — rupture didn't reach stations; investigate build/mesh"
   elif [ "${SLIP_COUNT}" -ge 6 ] && [ "${RVB1}" = "PASS" ]; then
      echo "  Decision: ${SLIP_COUNT}/9 slipped AND R-V91-B-pass-1 PASS  ==> close R-V91 with tpv102_debug_v9.1.0_check.md"
   elif [ "${SLIP_COUNT}" -ge 6 ] && [ "${RVB1}" = "FAIL" ]; then
      echo "  Decision: ${SLIP_COUNT}/9 slipped but R-V91-B-pass-1 FAIL ==> open v9.1.1 targeting H-V91-B1/B2 with production data"
   else
      echo "  Decision: ${SLIP_COUNT}/9 slipped — partial rupture; investigate station-by-station"
   fi
} > "${RESULT_FILE}"

cat "${RESULT_FILE}"
