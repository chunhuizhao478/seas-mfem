#!/bin/bash
#SBATCH -J g0_factor
#SBATCH -o g0_factor_%j.out
#SBATCH -e g0_factor_%j.err
#SBATCH -p gg
#SBATCH -A EAR20006

# GPU plan G0: MUMPS factor size for one mesh and order.
#
#   sbatch -N nodes -n ranks -t hh:mm:ss --export=ALL,CASE=500m_p1 \
#          jobs/vista/g0_factor_size.sbatch
#   (submit from miniapps/seas; jobs/vista/g0_submit.sh sets nodes, ranks,
#   wall time and the chain)
#
# CASE = prod_p1, prod_p2, 500m_p1, 250m_p1, 500m_p2 or 250m_p2.
# One assembly and one factorization with --matrix-stats
# --exit-after-factorize (as jobs/bp5/bp5_matrix_stats_dev_2hr.sbatch did on
# Frontera), v63 flags.  --matrix-stats raises MUMPS's print level to 2, so
# the log has MUMPS's analysis estimate (written before it factorizes, so a
# job stopped by its wall time or memory still reports it) and the measured
# factor entries and memory.  prod_p1 and prod_p2 run on the timing job's
# layout (1 node, 144 ranks) and stand in for its factor statistics.
#
# g0_submit.sh chains 500m_p1 -> 250m_p1 -> 500m_p2 -> 250m_p2 with
# --dependency=afterok and --kill-on-invalid-dep=yes: once one case fails or
# runs out of time, the larger ones are cancelled and their factor size stays
# unknown until G3a's cuDSS analysis.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"
g0_check_binary

CASE="${CASE:?set CASE}"
case "${CASE}" in
   prod_p1) MESH_FILE="${G0_MESH_PROD}"; ORDER=1 ;;
   prod_p2) MESH_FILE="${G0_MESH_PROD}"; ORDER=2 ;;
   500m_p1) MESH_FILE="${G0_MESH_500M}"; ORDER=1 ;;
   500m_p2) MESH_FILE="${G0_MESH_500M}"; ORDER=2 ;;
   250m_p1) MESH_FILE="${G0_MESH_250M}"; ORDER=1 ;;
   250m_p2) MESH_FILE="${G0_MESH_250M}"; ORDER=2 ;;
   *) echo "ERROR: unknown CASE '${CASE}'" >&2; exit 1 ;;
esac

cd "${SEAS_DIR}"
OUT="bp5/results_g0_vista/factor_size/${CASE}_n${SLURM_JOB_NUM_NODES}_r${SLURM_NTASKS}_job${SLURM_JOB_ID}"
PFX="g0_factor_${CASE}"
mkdir -p "${OUT}"
g0_provenance "${OUT}"

echo "=== ${CASE}: ${SLURM_JOB_NUM_NODES} nodes, ${SLURM_NTASKS} ranks, start $(date -u +%T) ==="
ibrun ./seas_bp5_full --mesh "${MESH_FILE}" "${G0_V63_FLAGS[@]}" --qd-operators loops \
      --order "${ORDER}" --matrix-stats --exit-after-factorize \
      --output-dir "${OUT}" --output-prefix "${PFX}" > "${OUT}/run.log" 2>&1
rc=$?
echo "${rc}" > "${OUT}/exit_code"
echo "=== ${CASE}: exit ${rc} at $(date -u +%T) ==="

# MUMPS's analysis estimate and measured factor, as in
# jobs/bp5/bp5_matrix_stats_dev_2hr.sbatch.
{
   echo "case ${CASE}, ${SLURM_JOB_NUM_NODES} nodes, ${SLURM_NTASKS} ranks, exit ${rc}"
   grep -E " (ANALYSIS|FACTORIZATION) STEP " "${OUT}/run.log" || true
   grep -E "ESTIMATED|estimated" "${OUT}/run.log" | head -n 40 || true
   grep -E "(REAL|INTEGER) space for factors" "${OUT}/run.log" || true
   grep -E "Maximum estimated front size" "${OUT}/run.log" || true
   grep -E "Time for factorization|Elapsed time in factorization" "${OUT}/run.log" || true
   grep -E "INFOG\(|RINFOG\(" "${OUT}/run.log" | head -n 40 || true
   grep -E "MUMPS.*(retry|workspace|INFO)" "${OUT}/run.log" | head -n 20 || true
} | tee "${OUT}/mumps_summary.txt"
exit "${rc}"
