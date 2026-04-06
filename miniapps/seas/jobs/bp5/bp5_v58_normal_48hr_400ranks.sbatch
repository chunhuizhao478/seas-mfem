#!/bin/bash
#SBATCH -J bp5_v58
#SBATCH -o bp5_v58_%j.out
#SBATCH -e bp5_v58_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5 v58: Full production run (normal queue, 48hr).
# Sign chain aligned with Tandem, skeleton Dirichlet orientation fix,
# output negation removed.

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module list

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${LD_LIBRARY_PATH}"

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

export JOBS=8

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || {
    echo "ERROR: MFEM was not built with PETSc enabled"
    exit 1
}
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || {
    echo "ERROR: MFEM was not built with MUMPS enabled"
    exit 1
}

cd miniapps/seas
make clean
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="bp5/results_v58_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_v58_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Output prefix: ${OUTPUT_PREFIX}"

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --checkpoint-interval 0 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg
