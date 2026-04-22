#!/bin/bash
#SBATCH -J tpv102_v92_4s_diag
#SBATCH -o tpv102_v92_4s_diag_%j.out
#SBATCH -e tpv102_v92_4s_diag_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 01:40:00
#SBATCH -A EAR20006

# TPV102 v9.2.0 plan §5.1 — RANK-1 diagnostic run with
# SEAS_DIAG_FAULT_SIGMA enabled.  Instruments every Mult() call with
# MPI-reduced `max|Q_self[SYY]|`, `max|Q_self[VY]|`, `max|Q_self[SYZ]|`,
# `max|Q_self[VZ]|` on fault QPs (global-frame, before rotation).
#
# Plan hook — decision tree §7:
#   max|SYY_g|  ~ 0          → H-V92-K (RK4 station-output bug)
#   max|SYY_g|  ~ 12 MPa     → H-V92-G (bulk mode-II amplification)
#   max|SYY_g|  >> 50 MPa    → something else (investigate with §6.3)
#   max|SYY_g|  rank-dep     → H-V92-M compounds
#
# PREREQUISITE — source patch application (NOT applied by this script):
#   The SEAS_DIAG_FAULT_SIGMA flag requires source edits in
#     dynamic/wave_operator.inl   (Mult + ComputeFaceFluxRHS + ComputeSharedFaceFluxRHS)
#     dynamic/wave_operator.hpp   (member accumulators)
#   The patch file
#     debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_diag_fault_sigma.patch
#   is additive only: all edits are inside `#ifdef SEAS_DIAG_FAULT_SIGMA`
#   guards, so the default build (without the flag) is unchanged.
#
#   User must approve + apply on the local laptop:
#     cd /Users/chunhuizhao/projects/seas-mfem
#     patch -p0 < miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_diag_fault_sigma.patch
#     cd miniapps/seas && make test             # sanity check (no DIAG flag)
#   Then commit + push to origin so Frontera checkout sees the patched
#   source on the workpath.
#
# After application: this sbatch compiles with -DSEAS_DIAG_FAULT_SIGMA
# so the #ifdef blocks become live.  Without application the sbatch
# still compiles and runs (flag becomes a no-op), but the [FAULT-DIAG-SYY]
# lines never appear in the log — so the run is diagnostic-dead.
#
# Plan references:
#   - debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_debug_plan.md
#     §5.1 (C0 diagnostic patch), §7 (decision matrix)
#   - debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_diag_fault_sigma.patch
#   - feedback_sbatch_modules.md
#   - feedback_frontera_approval.md — EXPLICIT APPROVAL REQUIRED
#
# WALL-CLOCK budget : ~1.4 hr (4 s tfinal at 400 ranks, same throughput
#                    as v91 job 7668434's first 4 s).
# SU COST           : ~358 SU (8 nodes × 1.4 hr × 32 SU/node-hr).

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

# Full rebuild with the diagnostic flag.  Without application of the
# v9.2.0_diag_fault_sigma.patch, the `#ifdef` blocks are empty and the
# flag is a no-op — use `tpv102_200m_p1_4.0s_400rank_v92_diag_check.sh`
# below to detect this case.
rm -f seas_tpv102_driver drivers/tpv102_driver.o dynamic/*.o io/*.o fault/*.o
make seas_tpv102_driver -j"${JOBS}" \
   SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_SIGMA" \
   || { echo "BUILD FAILED"; exit 1; }

# Sanity: grep the built driver binary for the diag signature string.
if ! strings seas_tpv102_driver 2>/dev/null | grep -q "FAULT-DIAG-SYY"; then
   echo "ERROR: binary does NOT contain [FAULT-DIAG-SYY] — patch was not applied?"
   echo "       Re-apply debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_diag_fault_sigma.patch"
   echo "       and re-push before re-submitting this sbatch."
   exit 2
fi

RESULT_DIR="tpv102/results_200m_p1_4.0s_400r_v92_diag_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"
mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase v9.2.0 §5.1 diagnostic: 400 ranks / 8 nodes / 200 m / P1 / tfinal=4.0 s"
echo "Purpose  : H-V92-G vs H-V92-K discrimination via max|Q_self[SYY]| trace"

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
      --debug-qnorm

# -------------------------------------------------------------------------
# Post-run classification for §7 decision matrix.
# Extracts peak max|SYY_g| across all Mult calls, peak max|VY_g|, and
# decides H-V92-G vs H-V92-K vs compound.
# -------------------------------------------------------------------------
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_v92_4s_diag_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_v92_4s_diag_${SLURM_JOB_ID}.err"

if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData" "${JOB_ERR}" 2>/dev/null; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData" "${JOB_ERR}" | head -1)
   echo "FAIL: MFEM abort during §5.1 diagnostic run" > "${RESULT_FILE}"
   echo "      ${ABORT_LINE}" >> "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

# Parse max|SYY_g| (column 3 in the `[FAULT-DIAG-SYY] max|SYY|=...` lines).
# stderr → *.err; also duplicated on stdout if we used 2>&1.  Search both.
DIAG_LOG="${JOB_LOG}"
if [ ! -s "${DIAG_LOG}" ] && [ -s "${JOB_ERR}" ]; then DIAG_LOG="${JOB_ERR}"; fi

N_DIAG=$(grep -c "\[FAULT-DIAG-SYY\]" "${DIAG_LOG}" 2>/dev/null || echo 0)
PEAK_SYY=$(grep "\[FAULT-DIAG-SYY\]" "${DIAG_LOG}" 2>/dev/null \
           | awk '{ for (i=1;i<=NF;i++) if ($i ~ /^max\|SYY\|=/) { v=$i; sub(/max\|SYY\|=/,"",v); if (v+0 > m) m=v+0 } } END { print m+0 }')
PEAK_VY=$(grep "\[FAULT-DIAG-SYY\]" "${DIAG_LOG}" 2>/dev/null \
          | awk '{ for (i=1;i<=NF;i++) if ($i ~ /^max\|VY\|=/) { v=$i; sub(/max\|VY\|=/,"",v); if (v+0 > m) m=v+0 } } END { print m+0 }')

{
   echo "§5.1 diagnostic run (SEAS_DIAG_FAULT_SIGMA) summary"
   echo "  Mult-call diagnostic lines : ${N_DIAG}"
   echo "  Peak max|Q_self[SYY]|      : ${PEAK_SYY} Pa"
   echo "  Peak max|Q_self[VY]|       : ${PEAK_VY} m/s"
   echo
   echo "Decision tree (plan §7):"
   if awk "BEGIN { exit !( ${PEAK_SYY} < 1e6 ) }"; then
      echo "  PEAK_SYY < 1 MPa      → bulk σ_yy is clean"
      echo "  → If station σ_n_corr still drifts, H-V92-K (RK4 stage-averaging) PRIMARY"
      echo "  → Else investigate H-V92-M (shared-fault bit-disagreement) with §4.5"
   elif awk "BEGIN { exit !( ${PEAK_SYY} < 2e7 ) }"; then
      echo "  1 MPa ≤ PEAK_SYY < 20 MPa → mode-II bulk coupling present but bounded"
      echo "  → Investigate H-V92-C (corner) via §6.3 ParaView slice"
   else
      echo "  PEAK_SYY ≥ 20 MPa     → bulk σ_yy drives the observed drift"
      echo "  → H-V92-G CONFIRMED as primary; identify amplifier (a/b/c/d in §3 H-V92-G)"
   fi
   echo
   echo "Next step: run §4.5 rankA sbatch (tpv102_200m_p1_4.0s_56rank_v92_rankA.sbatch)"
   echo "           + off-Frontera harness to co-classify H-V92-M alongside H-V92-G."
} > "${RESULT_FILE}"

cat "${RESULT_FILE}"
