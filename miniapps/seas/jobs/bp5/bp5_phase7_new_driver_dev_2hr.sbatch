#!/bin/bash
#SBATCH -J bp5_p7v
#SBATCH -o bp5_phase7_verify_%j.out
#SBATCH -e bp5_phase7_verify_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# Phase 7: Verification run using NEW TOML-based driver.
# Tests refactored code paths: ConstitutiveModel, BoundaryConfig,
# DomainConfig, SEASConfigParser.

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
make clean
make seas_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="bp5/results_phase7_verify_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase 7 dev test: 400 ranks, 2hr, 50 steps"

ibrun ./seas_driver config/bp5_verification_50step.toml \
      --override output.output_dir="${RESULT_DIR}" \
      --override output.output_prefix="bp5_phase7_verify_job${SLURM_JOB_ID}" \
      --write-every-step
