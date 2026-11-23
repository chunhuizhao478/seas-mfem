#!/bin/bash
#SBATCH -J d1_cuda_build
#SBATCH -o d1_cuda_build_%j.out
#SBATCH -e d1_cuda_build_%j.err
#SBATCH -p gh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -A EAR20006

# GPU dynamic-rupture plan D1 (plan rev2 §4.2, Build): the CUDA build of
# seas_nitsche_dyn_driver on a Vista Grace-Hopper node (GH200, sm_90), with
# build_horizon.sh MODE=cuda, in a clone of its own (MFEM is configured in
# place): $WORK/seas-mfem-gpu-dyn-cuda on gpu-dyn-nitsche.
#
# It links the quasi-dynamic thread's CUDA PETSc (hypre, MUMPS) and cuDSS in
# DEPS_ROOT=$WORK/seas-gpu-deps without writing there (the build provenance
# goes to $OUT), and builds the parallel HDF5 the driver needs in
# HDF5_ROOT=$WORK/seas-gpu-dyn-deps (hdf5-1.14.6-cuda; skipped when present).
# Login node first (internet), from the clone's root:
#   MODE=cuda STAGE=hdf5-fetch DEPS_ROOT=$WORK/seas-gpu-deps \
#      HDF5_ROOT=$WORK/seas-gpu-dyn-deps PROVENANCE_FILE=/dev/null bash build_horizon.sh
#   cd miniapps/seas && sbatch jobs/vista/d1_cuda_build.sbatch
# Then: the clone's meshes are D0a's (copied from $WORK/seas-mfem-gpu-dyn-nit
# when absent, sha256 checked), and BUILD_COMMIT_CUDA records the commit for
# d1_gh_check.sbatch.  A file nvcc cannot compile stops the job; it is
# reported, not worked around here.
# Cost: 1 gh node x <= 4 h (1 SU per node-hour, plan rev2 §1).
set -euo pipefail
export LC_ALL=C LANG=C
ROOT="$(cd "${SLURM_SUBMIT_DIR}/../.." && pwd)"
NIT="${NIT:-${WORK}/seas-mfem-gpu-dyn-nit}"
# shellcheck source=d1_env.sh
source "${ROOT}/miniapps/seas/jobs/vista/d1_env.sh"
d1_guard "${ROOT}"
[ "$(git -C "${ROOT}" rev-parse HEAD)" = "$(git -C "${ROOT}" rev-parse '@{u}' 2>/dev/null)" ] \
   || { echo "ERROR: ${ROOT} HEAD is not its upstream; pull first"; exit 2; }
OUT="${SCRATCH}/gpu_dyn/d1_cuda_build_${SLURM_JOB_ID}"
mkdir -p "${OUT}"
echo "job ${SLURM_JOB_ID} on $(hostname), $(date -u +%FT%TZ); checkout ${ROOT} at $(git -C "${ROOT}" rev-parse HEAD)"

cd "${ROOT}"
MODE=cuda GPU_ARCH=90 STAGE="hdf5 mfem seas" WITH_HDF5=1 \
   DEPS_ROOT="${WORK}/seas-gpu-deps" HDF5_ROOT="${WORK}/seas-gpu-dyn-deps" \
   PROVENANCE_FILE="${OUT}/build_provenance_cuda.txt" \
   SEAS_TARGETS="seas_nitsche_dyn_driver" JOBS=72 bash build_horizon.sh
git rev-parse HEAD > miniapps/seas/BUILD_COMMIT_CUDA
grep -E '^MFEM_USE_(CUDA|HDF5|PETSC|CUDSS) ' config/config.mk

cd "${ROOT}/miniapps/seas"
for m in tpv205/mesh/tpv2053d_1000m.msh tpv104/mesh/tpv104_1000m.msh; do
   [ -f "${m}" ] || cp -v "${NIT}/miniapps/seas/${m}" "${m}"
done
d1_check_meshes "${ROOT}/miniapps/seas"
echo "D1 CUDA build: done; commit $(cat BUILD_COMMIT_CUDA); provenance ${OUT}/build_provenance_cuda.txt"
