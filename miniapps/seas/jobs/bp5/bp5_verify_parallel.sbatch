#!/bin/bash
#SBATCH -J bp5_vp4
#SBATCH -o bp5_verify_parallel_%j.out
#SBATCH -e bp5_verify_parallel_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# BP5 verification: parallel (400 ranks) test.
# Runs --verify diagnostics + short simulation (~50 steps).
# Compare output against bp5_verify_serial to detect shared-face bugs.
#
# Key comparison:
#   grep "\[VERIFY\]" bp5_verify_serial_*.out > serial_verify.txt
#   grep "\[VERIFY\]" bp5_verify_parallel_*.out > parallel_verify.txt
#   diff serial_verify.txt parallel_verify.txt
# If ||b_dir|| or ||b_total|| differ → shared-face Dirichlet bug confirmed.

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
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
make clean
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="bp5/results_verify_parallel_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_verify_parallel_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Verification: parallel (400 ranks) — 30 min dev queue"

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --ref-dir bp5/benchmark_data \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --checkpoint-interval 0 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --max-steps 50 \
      --verify
