#!/bin/bash
#SBATCH -J tpv102_v92_4s_56r
#SBATCH -o tpv102_v92_4s_56r_%j.out
#SBATCH -e tpv102_v92_4s_56r_%j.err
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 56
#SBATCH -t 01:30:00
#SBATCH -A EAR20006

# TPV102 v9.2.0 plan §4.5 — rank-consistency diagnostic "run A"
# (1-node / 56-rank partition).  Companion job for the existing
# v91 400-rank dataset (job 7668434) which serves as "run B".
# Together they exercise the §4.5 1-rank-vs-4-rank classification
# at production mesh resolution (200 m, P1, identical build).
#
# Plan hook — decision matrix §7 / §13.4 option B:
#   σ_n agree  (|Δ σ_n| ≤ 0.1 MPa) → H-V92-M ELIMINATED; bulk/corner/driver bug
#   σ_n differ (|Δ σ_n| > 0.1 MPa) → H-V92-M CONFIRMED; audit shared-fault ctor
#
# Tfinal = 4.0 s matches the onset window (σ_n drift starts ~3 s at
# hypocentre per job 7668434 flt_0_7.5.dat) and is the minimum
# duration at which the disagreement (if any) becomes measurable.
#
# Rank choice = 56 (1 node).  Reasons:
#   - 1 rank on 200 m mesh is infeasible (~533 h wallclock for 4 s).
#   - 56 ranks still exercises MPI shared-fault code — within a node
#     ParMETIS cuts the domain and shared faces exist across ranks.
#   - ParMETIS cut-line topology at 56 ranks differs from 400 ranks,
#     so a bit-disagreement on `elem1_on_plus` will manifest as
#     different σ_n(t) at stations.
#
# Output dir names include the Slurm job id so run A and run B are
# easy to pair in the post-processing.
#
# Plan references:
#   - debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_debug_plan.md
#     §4.5 (rank-1 diagnostic), §7 (decision matrix),
#     §13.4 (option B), §13.5 (implications)
#   - feedback_sbatch_modules.md
#   - feedback_frontera_approval.md — EXPLICIT APPROVAL REQUIRED BEFORE
#     `sbatch` (this file is generated, not yet submitted)
#
# WALL-CLOCK budget : ~1.5 hr (~70 min expected wallclock).
# SU COST           : ~48 SU (1 node × 1.5 hr × 32 SU/node-hr).
#
# Usage (AFTER user approval):
#   sbatch jobs/tpv102/tpv102_200m_p1_4.0s_56rank_v92_rankA.sbatch

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
# Force full rebuild so the build for this rank-consistency job is
# provably identical to the build used by run B (the 400-rank v91 job
# 7668434 at commit 99f355a).  This is critical: ANY build difference
# silently invalidates the comparison.
rm -f seas_tpv102_driver drivers/tpv102_driver.o dynamic/*.o io/*.o fault/*.o
make seas_tpv102_driver -j"${JOBS}" SEAS_EXTRA_CPPFLAGS="" \
   || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_200m_p1_4.0s_56r_v92_rankA_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"
mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase v9.2.0 §4.5 rank A: 56 ranks / 1 node / 200 m / P1 / tfinal=4.0 s"
echo "Purpose  : H-V92-M rank-consistency classification (A vs existing v91 400-rank run B)"

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
# Post-run: one-sided σ_n sanity report at t=3.5s (plan §4.5 sample time).
# The cross-rank comparison with run B (the v91 job 7668434 dataset) is
# done offline by `tests/scripts/test_sigma_n_rank_consistency.py`.
# -------------------------------------------------------------------------
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_v92_4s_56r_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv102_v92_4s_56r_${SLURM_JOB_ID}.err"

if grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData" "${JOB_ERR}" 2>/dev/null; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData" "${JOB_ERR}" | head -1)
   echo "FAIL: MFEM abort during §4.5 rank-A run" > "${RESULT_FILE}"
   echo "      ${ABORT_LINE}" >> "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

HYPO_FILE="${RESULT_DIR}/flt_0_7.5.dat"
# Columns (tpv102_setup.hpp::TPV102StationWriter):
#   1=t 2=slip1(dip) 3=slip2(strike) 4=V1(dip) 5=V2(strike)
#   6=tau1(dip) 7=tau2(strike) 8=sigma_n 9=log10_theta
SIGMA_N_AT_35=""
SIGMA_N_AT_20=""
SIGMA_N_INIT=""
if [ -f "${HYPO_FILE}" ]; then
   SIGMA_N_INIT=$(awk '!/^#/ && NF>=8 {print $8; exit}' "${HYPO_FILE}")
   SIGMA_N_AT_20=$(awk '!/^#/ && NF>=8 && $1>=2.0 {print $8; exit}' "${HYPO_FILE}")
   SIGMA_N_AT_35=$(awk '!/^#/ && NF>=8 && $1>=3.5 {print $8; exit}' "${HYPO_FILE}")
fi

{
   echo "§4.5 rank A (56 ranks, 1 node, tfinal=4.0s) hypocentre σ_n snapshots"
   echo "  σ_n(t=0.0)  (initial) : ${SIGMA_N_INIT}"
   echo "  σ_n(t=2.0s)           : ${SIGMA_N_AT_20}"
   echo "  σ_n(t=3.5s)           : ${SIGMA_N_AT_35}"
   echo
   echo "NEXT STEP — off-Frontera post-processing:"
   echo "  python3 tests/scripts/test_sigma_n_rank_consistency.py \\"
   echo "      --skip-run \\"
   echo "      --output-root tpv102/r_v92_r4_rank_consistency \\"
   echo "      --run-a ${RESULT_DIR} \\"
   echo "      --run-b tpv102/results_200m_p1_12.0s_400r_v91_job7668434 \\"
   echo "      --tol-mpa 0.1 --t-check 2.0 3.5"
   echo
   echo "CLASSIFICATION (plan §4.5 decision matrix):"
   echo "  |Δ σ_n| ≤ 0.1 MPa at all stations → H-V92-M ELIMINATED"
   echo "  |Δ σ_n| > 0.1 MPa at any station  → H-V92-M CONFIRMED"
} > "${RESULT_FILE}"

cat "${RESULT_FILE}"
