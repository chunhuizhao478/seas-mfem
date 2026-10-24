#!/bin/bash
#SBATCH -J g1_cuda_build
#SBATCH -o g1_cuda_build_%j.out
#SBATCH -e g1_cuda_build_%j.err
#SBATCH -p gh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -A EAR20006

# GPU plan G1: the CUDA build on a Vista Grace-Hopper node (GH200, sm_90).
#
# Run in a clone of its own (MFEM is configured in place, so the CPU
# checkout $WORK/seas-mfem-gpu cannot hold it), e.g. $WORK/seas-mfem-gpu-cuda
# on the same pushed commit.  First, on a login node (internet):
#   MODE=cuda GPU_ARCH=90 STAGE=fetch bash build_horizon.sh
# then, from that clone's miniapps/seas:
#   sbatch jobs/vista/g1_cuda_build.sbatch
#
# Builds PETSc 3.24.6 with CUDA (hypre on the GPU, MUMPS on the CPU),
# unpacks cuDSS, builds MFEM with CUDA and cuDSS and seas_bp5_full (through
# nvcc), then ex1p and ex2p and the device-buffer MPI test for
# g1_gpu_checks.sbatch.  A seas file nvcc cannot compile stops the job; it
# is reported, not worked around here.

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

MODE=cuda GPU_ARCH=90 STAGE="deps mfem seas" JOBS=72 bash build_horizon.sh

# Acceptance programs (same modules as the build).
module load gcc/13.2.0 "${CUDA_MODULE:-cuda/13.0}" openmpi/5.0.11 nvpl/26.9 cmake/3.31.5
make -C examples -j8 ex1p ex2p
nvcc -ccbin mpicxx -arch=sm_90 -o miniapps/seas/jobs/vista/g1_mpi_cuda_allreduce \
     miniapps/seas/jobs/vista/g1_mpi_cuda_allreduce.cu
ls -la examples/ex1p examples/ex2p miniapps/seas/seas_bp5_full \
       miniapps/seas/jobs/vista/g1_mpi_cuda_allreduce
