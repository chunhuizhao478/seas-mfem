#!/bin/bash
#SBATCH -J bp5_fs_cg
#SBATCH -o bp5_fs_cg_%j.out
#SBATCH -e bp5_fs_cg_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# BP5 v59 test 1: First-step dump with CG+AMG solver instead of MUMPS.
# If displacement matches MUMPS run → solver is not the issue, K is different.
# If displacement differs → MUMPS interface may be the issue.

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${LD_LIBRARY_PATH}"

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }

cd miniapps/seas
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="bp5/results_v59_fs_cg_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_v59_fs_cg_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo "TEST: CG+AMG solver (instead of MUMPS)"

export SEAS_DEBUG_FIRST_STEP_DUMP=1
export SEAS_DEBUG_FIRST_STEP_TARGET_RANK=96
echo "SEAS_DEBUG_FIRST_STEP_DUMP=${SEAS_DEBUG_FIRST_STEP_DUMP}"
echo "SEAS_DEBUG_FIRST_STEP_TARGET_RANK=${SEAS_DEBUG_FIRST_STEP_TARGET_RANK}"

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver cg \
      --tfinal 1.0 \
      --ref-dir bp5/benchmark_data \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --checkpoint-interval 0 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg

echo "First-step dump files:"
ls -lh "${RESULT_DIR}"/first_step_*.csv 2>/dev/null || echo "  (none found)"
