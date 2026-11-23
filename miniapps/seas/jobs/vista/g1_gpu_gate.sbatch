#!/bin/bash
#SBATCH -J g1_gpu_gate
#SBATCH -o g1_gpu_gate_%j.out
#SBATCH -e g1_gpu_gate_%j.err
#SBATCH -p gh
#SBATCH -N 1
#SBATCH -t 04:00:00
#SBATCH -A EAR20006

# GPU plan G1, first GPU run of the unported code: seas_bp5_full from the
# CUDA build with --device cuda and nothing ported, so PETSc's RK vector
# updates run on the GPU and everything else (assembly, MUMPS, friction,
# output) on the host.  A stale host or device copy shows up as an error of
# order one, not as rounding.
#
#   sbatch -N N -n N jobs/vista/g1_gpu_gate.sbatch   (N = 1, 2, 4, 8 GPUs;
#   one GPU per gh node; from miniapps/seas of the CUDA clone)
#
# The 50-step p=1 gate against G0's golden at N ranks with
# TOL = max(1e-8, 10 floor_N) (floor.json key G1_first_horizon_run), judged
# by check_gate.sh.  The measured difference is the CPU-to-GPU floor the
# later GPU gates add.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g1_cuda_env.sh"
g1_check_checkout
g0_check_binary

cd "${SEAS_DIR}"
NP="${SLURM_NTASKS:-1}"
TAG="p1_r${NP}"
GOLD="${G1_GOLDEN_ROOT}/${TAG}"
GPFX="g0_${TAG}"
OUT="bp5/goldens/g1_gpu_gate/job${SLURM_JOB_ID}/${TAG}"
PFX="g1gpu_${TAG}"
[ -f "${GOLD}/floor.json" ] && [ -d "${GOLD}/run_a" ] || {
   echo "ERROR: no G0 golden in ${GOLD}" >&2
   exit 1
}
TOL="$(python3 -c "import json; print(json.load(open('${GOLD}/floor.json'))['tol']['G1_first_horizon_run'])")"
mkdir -p "${OUT}"
g0_provenance "${OUT}"
echo "=== ${TAG} on ${NP} GPU(s): TOL ${TOL}, start $(date -u +%T) ==="
ibrun -n "${NP}" ./seas_bp5_full --mesh "${G0_MESH_PROD}" "${G0_V63_FLAGS[@]}" --qd-operators loops \
      --order 1 --device cuda \
      --max-steps 50 --write-every-step --checkpoint-interval 0 \
      --output-dir "${OUT}" --output-prefix "${PFX}" \
      --ref-dir "${GOLD}/run_a" --ref-prefix "${GPFX}" \
      --regression-tolerance "${TOL}" > "${OUT}/run.log" 2>&1
rc=$?
echo "${rc}" > "${OUT}/exit_code"
echo "=== ${TAG}: exit ${rc} at $(date -u +%T) ==="
tail -n 25 "${OUT}/run.log"
bash jobs/horizon/check_gate.sh "${OUT}" "${PFX}" "${GOLD}/run_a" "${GPFX}" \
     "${TOL}" "${OUT}/run.log" "${rc}"
