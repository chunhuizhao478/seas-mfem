#!/bin/bash
#SBATCH -J h_cuda_build
#SBATCH -o h_cuda_build_%j.out
#SBATCH -e h_cuda_build_%j.err
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -A EAR26024

# GPU plan on Horizon: the CUDA build on a GB200 node (sm_100), as
# jobs/vista/g1_cuda_build.sbatch on Vista, with the settings of
# jobs/horizon/h_env.sh.  Run in the Horizon clone
# ($SCRATCH/seas-mfem-gpu-cuda).  First, on a login node (internet):
#   source jobs/horizon/h_env.sh && MODE=cuda STAGE=fetch bash ../../build_horizon.sh
# (from miniapps/seas), then
#   sbatch jobs/horizon/h_cuda_build.sbatch
#
# Builds PETSc 3.24.6 with CUDA (hypre on the GPU, MUMPS on the CPU),
# unpacks cuDSS, builds MFEM with CUDA and cuDSS and seas_bp5_full, then
# ex1p, ex2p and the device-buffer MPI test for h_gpu_checks.sbatch.  A
# build error stops the job; it is reported, not worked around here.

set -euo pipefail
export LC_ALL=C LANG=C
ROOT="$(cd "${SLURM_SUBMIT_DIR}/../.." && pwd)"
cd "${ROOT}"
if [ -n "$(git status --porcelain --untracked-files=no)" ] ||
   [ "$(git rev-parse HEAD)" != "$(git rev-parse '@{u}' 2>/dev/null)" ]; then
   echo "ERROR: ${ROOT} is not clean at its upstream commit; build only pushed commits." >&2
   exit 1
fi
echo "Checkout: ${ROOT} at $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD))"
source "${SLURM_SUBMIT_DIR}/jobs/horizon/h_env.sh" || exit 1

MODE=cuda STAGE="deps mfem seas" JOBS=144 bash build_horizon.sh

# Acceptance programs (same modules as the build).
# shellcheck disable=SC2086
module load ${MODULES}
make -C examples -j16 ex1p ex2p
nvcc -ccbin mpicxx -arch="sm_${GPU_ARCH}" -o miniapps/seas/jobs/vista/g1_mpi_cuda_allreduce \
     miniapps/seas/jobs/vista/g1_mpi_cuda_allreduce.cu
ls -la examples/ex1p examples/ex2p miniapps/seas/seas_bp5_full \
       miniapps/seas/jobs/vista/g1_mpi_cuda_allreduce
