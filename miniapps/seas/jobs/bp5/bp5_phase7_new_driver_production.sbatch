#!/bin/bash
#SBATCH -J bp5_p7p
#SBATCH -o bp5_phase7_production_%j.out
#SBATCH -e bp5_phase7_production_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# Phase 7: Full production run using NEW TOML-based driver with PETSc TS.
# Tests all refactored code paths: ConstitutiveModel, BoundaryConfig,
# DomainConfig, SEASConfigParser, and PETSc TS time integration.
# Compares against seas-mfem parallel golden (golden_parallel_400r_p1_1000m).
#
# Before submitting, on the login node:
#   cd /scratch2/10024/zhaochun/seas-project/seas-mfem
#   git fetch origin && git pull
#   git submodule update --init miniapps/seas/extern/toml11
#   cd miniapps/seas && sbatch jobs/bp5/bp5_phase7_new_driver_production.sbatch
#
# Post-run comparison:
#   python3 scripts/regression_check.py \
#       --sim-dir $RESULT_DIR --sim-prefix $OUTPUT_PREFIX \
#       --ref-dir bp5/benchmark_data/golden_parallel_400r_p1_1000m \
#       --tolerance 1e-2

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas

# Verify toml11 submodule is present
test -f extern/toml11/include/toml.hpp || test -f extern/toml11/toml.hpp || \
   { echo "ERROR: toml11 not initialized. Run: git submodule update --init miniapps/seas/extern/toml11"; exit 1; }

make clean
make seas_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="bp5/results_phase7_production_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_phase7_prod_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase 7: New TOML driver + PETSc TS, full production (400 ranks, 48hr)"

ibrun ./seas_driver config/bp5_production_new_driver.toml \
      --override output.output_dir="${RESULT_DIR}" \
      --override output.output_prefix="${OUTPUT_PREFIX}" \
      --write-every-step
