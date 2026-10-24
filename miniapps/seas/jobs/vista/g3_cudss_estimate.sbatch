#!/bin/bash
#SBATCH -J g3_cudss_estimate
#SBATCH -o g3_cudss_estimate_%j.out
#SBATCH -e g3_cudss_estimate_%j.err
#SBATCH -p gh
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# GPU plan G3a: cuDSS's analysis estimates (LU nonzeros, permanent and peak
# device and host memory per rank) for one mesh and order, with no
# factorization, on N GPUs (one rank per gh node), in the CUDA clone after
# g3_gpu_gate.sbatch built seas_test_cudss_vs_mumps:
#
#   sbatch --export=ALL,RES=250m,ORDER=2 jobs/vista/g3_cudss_estimate.sbatch
#   (from miniapps/seas of the CUDA clone; RES = 1000m, 500m or 250m;
#   -N N -n N for another rank count)
#
# 250 m at p=2 has 4.1e9 nonzeros, more than a 32-bit index counts on one
# rank (G0), so it needs at least 4 ranks.  A measurement, not a gate.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g1_cuda_env.sh"
g1_check_checkout
cd "${SEAS_DIR}"
[ -x ./seas_test_cudss_vs_mumps ] || {
   echo "ERROR: seas_test_cudss_vs_mumps missing; run g3_gpu_gate.sbatch first." >&2
   exit 1
}

RES="${RES:?set RES=1000m, 500m or 250m}"
ORDER="${ORDER:?set ORDER=1 or 2}"
case "${RES}" in
   1000m) MESH="${G0_MESH_PROD}" ;;
   500m)  MESH="${G0_MESH_500M}" ;;
   250m)  MESH="${G0_MESH_250M}" ;;
   *) echo "ERROR: RES must be 1000m, 500m or 250m" >&2; exit 1 ;;
esac
NP="${SLURM_NTASKS:-4}"
OUT="bp5/goldens/g3_cudss_estimate/job${SLURM_JOB_ID}"
LOG="${OUT}/${RES}_p${ORDER}_r${NP}.log"
mkdir -p "${OUT}"
g0_provenance "${OUT}"
echo "=== cuDSS estimate ${RES} p=${ORDER} on ${NP} GPU(s): start $(date -u +%T) ==="
ibrun -n "${NP}" ./seas_test_cudss_vs_mumps --mesh "${MESH}" --order "${ORDER}" \
      --device cuda --estimate-only > "${LOG}" 2>&1
rc=$?
grep -E "^G3a|cuDSS analysis|^  rank|^  [0-9]|^  \(" "${LOG}"
echo "ESTIMATE ${RES}_p${ORDER}_r${NP}: exit ${rc}"
exit "${rc}"
