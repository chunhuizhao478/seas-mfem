#!/bin/bash
#SBATCH -J tandem_petsc
#SBATCH -o tandem_petsc_%j.out
#SBATCH -e tandem_petsc_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 72
#SBATCH -t 03:00:00
#SBATCH -A EAR20006

# Tandem vs MFEM-SEAS speed test: PETSc 3.22.3 for Tandem only (upstream
# Tandem b75f6669 does not compile against our 3.24.6: job 1113226).  Same
# script and configure as our own PETSc (build_horizon.sh MODE=cpu
# STAGE=deps: gcc 13, OpenMPI 5, nvpl BLAS/LAPACK, -mcpu=neoverse-v2,
# MUMPS, ScaLAPACK, METIS, ParMETIS, PT-Scotch, hypre at the versions
# PETSc 3.22.3 pins), installed in ${TANDEM_ROOT}/petsc-3.22.3-cpu.  Our
# MFEM build and ${WORK}/seas-gpu-deps are not touched.
#
# First, on a login node (internet), from the checkout root:
#   PETSC_VERSION=3.22.3 DEPS_ROOT=$WORK/seas-mfem-tandem STAGE=fetch bash build_horizon.sh
# Then, from miniapps/seas:
#   sbatch jobs/vista_qd_speed/tandem/tandem_petsc.sbatch

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista_qd_speed/tandem/tandem_env.sh"
seas_guard

cd "${SEAS_MFEM_ROOT}" || exit 1
MODE=cpu STAGE=deps JOBS=64 MODULES=none \
   PETSC_VERSION="${TANDEM_PETSC_VERSION}" DEPS_ROOT="${TANDEM_ROOT}" \
   bash build_horizon.sh \
   || { echo "ERROR: PETSc ${TANDEM_PETSC_VERSION} build failed" >&2; exit 1; }
[ "${TANDEM_PETSC}" = "${TANDEM_ROOT}/petsc-${TANDEM_PETSC_VERSION}-cpu" ] \
   && [ -f "${TANDEM_PETSC}/lib/petsc/conf/petscvariables" ] \
   || { echo "ERROR: no PETSc in ${TANDEM_PETSC}" >&2; exit 1; }
echo "=== PETSc for Tandem: ${TANDEM_PETSC} ==="
tail -n 20 "${TANDEM_ROOT}/build_provenance_cpu.txt"
