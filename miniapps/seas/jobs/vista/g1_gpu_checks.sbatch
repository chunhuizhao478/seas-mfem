#!/bin/bash
#SBATCH -J g1_gpu_checks
#SBATCH -o g1_gpu_checks_%j.out
#SBATCH -e g1_gpu_checks_%j.err
#SBATCH -p gh-dev
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -A EAR20006

# GPU plan G1 acceptance on Vista gh nodes (one GH200 per node), before
# any seas run on the GPU:
#   1. ex1p -d cuda and ex2p -d cuda on one GPU (ex2p without -elast);
#   2. a device-buffer MPI_Allreduce on 2 GPUs (CUDA-aware MPI);
#   3. ex1p -d cuda -cudss on 1 GPU and on 2 GPUs (cuDSS distributed mode),
#      relative residual ||AX-B||/||B|| below 1e-10.
#
#   sbatch jobs/vista/g1_gpu_checks.sbatch   (from miniapps/seas of the CUDA clone)

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g1_cuda_env.sh"
g1_check_checkout
EX="${SEAS_MFEM_ROOT}/examples"
fail=0
check() {  # check NAME COMMAND...: run, keep the log, record pass/fail
   local name="$1"; shift
   echo "=== ${name}: $(date -u +%T) ==="
   "$@" > "g1_check_${name}_${SLURM_JOB_ID}.log" 2>&1
   local rc=$?
   tail -n 5 "g1_check_${name}_${SLURM_JOB_ID}.log"
   echo "=== ${name}: exit ${rc} ==="
   [ "${rc}" = 0 ] || fail=1
   return "${rc}"
}
residual_ok() {  # residual_ok LOG: the cuDSS residual line exists and is < 1e-10
   python3 - "$1" <<'PY'
import re, sys
m = re.findall(r"cuDSS relative residual \|\|AX-B\|\|/\|\|B\|\| = (\S+)", open(sys.argv[1]).read())
ok = bool(m) and float(m[-1]) < 1e-10
print(f"cuDSS residual: {m[-1] if m else 'missing'} -> {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
PY
}

cd "${EX}"
check ex1p_cuda  ibrun -n 1 -o 0 ./ex1p -d cuda -no-vis
check ex2p_cuda  ibrun -n 1 -o 0 ./ex2p -d cuda -no-vis
cd "${SEAS_DIR}"
check mpi_allreduce_2gpu ibrun -n 2 -o 0 ./jobs/vista/g1_mpi_cuda_allreduce
cd "${EX}"
check ex1p_cudss_1gpu ibrun -n 1 -o 0 ./ex1p -d cuda -cudss -no-vis
residual_ok "g1_check_ex1p_cudss_1gpu_${SLURM_JOB_ID}.log" || fail=1
check ex1p_cudss_2gpu ibrun -n 2 -o 0 ./ex1p -d cuda -cudss -no-vis
residual_ok "g1_check_ex1p_cudss_2gpu_${SLURM_JOB_ID}.log" || fail=1
echo "G1 GPU checks overall: $([ "${fail}" = 0 ] && echo PASS || echo FAIL)"
exit "${fail}"
