#!/bin/bash
#SBATCH -J g0_mcount
#SBATCH -o g0_mcount_%j.out
#SBATCH -e g0_mcount_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -A EAR20006

# GPU plan G0: dimensions and nonzero counts of K, no factorization.
#
#   sbatch [-N nodes] --export=ALL,MESH=prod|500m|250m,ORDER=P[,RANKS="1 2 4 8"] \
#          jobs/vista/g0_matrix_counts.sbatch
#   (submit from miniapps/seas; jobs/vista/g0_submit.sh does this)
#
# For each rank count in RANKS (the GPU counts G6 plans), one
# --exit-after-assembly run with the v63 flags: K is assembled and its global,
# diag, offd, per-rank and per-row nonzero counts go to a JSON file, then the
# driver exits before any solver is set up.  The table at the end says
# whether each rank's matrix fits 32-bit indices.
#
# The 250 m mesh at p=2 (about 4.3e9 nonzeros by the upper bound) runs on 4
# and 8 ranks only: on 1 or 2 ranks a rank's local matrix would not fit the
# 32-bit indices MFEM assembles with.  Its 4 ranks need about 50 GB each
# while assembling (estimate, unmeasured); g0_submit.sh gives it 4 nodes.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"
g0_check_binary

MESH="${MESH:?set MESH=prod, 500m or 250m}"
ORDER="${ORDER:?set ORDER=1 or 2}"
RANKS="${RANKS:-1 2 4 8}"
case "${MESH}" in
   prod) MESH_FILE="${G0_MESH_PROD}" ;;
   500m) MESH_FILE="${G0_MESH_500M}" ;;
   250m) MESH_FILE="${G0_MESH_250M}" ;;
   *) echo "ERROR: MESH must be prod, 500m or 250m" >&2; exit 1 ;;
esac

cd "${SEAS_DIR}"
OUT="bp5/results_g0_vista/matrix_counts/${MESH}_p${ORDER}"
mkdir -p "${OUT}"
g0_provenance "${OUT}"

status=0
for r in ${RANKS}; do
   if [ "${MESH}" = 250m ] && [ "${ORDER}" -ge 2 ] && [ "${r}" -lt 4 ]; then
      echo "skip ${MESH} p=${ORDER} on ${r} ranks (local nnz above int32)"
      continue
   fi
   if [ "${r}" -gt "${SLURM_NTASKS}" ]; then
      echo "ERROR: ${r} ranks requested but the job has ${SLURM_NTASKS} tasks" >&2
      status=1; continue
   fi
   d="${OUT}/r${r}"
   mkdir -p "${d}"
   echo "=== ${MESH} p=${ORDER} ranks=${r}: $(date -u +%T) ==="
   ibrun -n "${r}" ./seas_bp5_full --mesh "${MESH_FILE}" "${G0_V63_FLAGS[@]}" \
         --order "${ORDER}" --exit-after-assembly \
         --output-dir "${d}" --output-prefix "g0_${MESH}_p${ORDER}_r${r}" \
         > "${d}/run.log" 2>&1
   rc=$?
   echo "${rc}" > "${d}/exit_code"
   echo "=== exit ${rc} at $(date -u +%T) ==="
   grep -E 'exit-after-assembly|\[matrix-stats\]|local nnz/rank' "${d}/run.log" | head -n 20
   [ "${rc}" -eq 0 ] || status=1
done

python3 jobs/vista/g0_matrix_table.py "${OUT}"/r*/*_matrix_stats.json \
   | tee "${OUT}/table.txt"
exit "${status}"
