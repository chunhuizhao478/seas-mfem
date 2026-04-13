#!/bin/bash
#SBATCH -J bp5_p1p
#SBATCH -o bp5_phase1_verify_parallel_%j.out
#SBATCH -e bp5_phase1_verify_parallel_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# Phase 1 verification: parallel (400 ranks) regression check.
# Runs 50-step simulation on refactor/phase1-mpi branch,
# compares against Phase 0 golden (golden_parallel_400r_50step).
#
# Expected: byte-for-byte match (regression-tolerance 0).
# If near-zero fields differ by FP noise, relax to 1e-12.
#
# Compare serial vs parallel:
#   grep "\[VERIFY\]" bp5_phase1_verify_serial_*.out > p1_serial_verify.txt
#   grep "\[VERIFY\]" bp5_phase1_verify_parallel_*.out > p1_parallel_verify.txt
#   diff p1_serial_verify.txt p1_parallel_verify.txt

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

# Checkout Phase 1 branch
git fetch origin
git checkout refactor/phase1-mpi
git pull origin refactor/phase1-mpi

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
make clean
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="bp5/results_phase1_verify_parallel_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_phase1_verify_parallel_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Output dir: ${RESULT_DIR}"
echo "Phase 1 parallel verification: 400 ranks, 50 steps, regression against Phase 0 golden"

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --ref-dir bp5/benchmark_data/golden_parallel_400r_50step \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --checkpoint-interval 0 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --max-steps 50 \
      --regression-tolerance 0 \
      --verify
