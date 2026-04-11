#!/bin/bash
#SBATCH -J bp5_v63d
#SBATCH -o bp5_v63_diag_sliprate_%j.out
#SBATCH -e bp5_v63_diag_sliprate_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# BP5 v63 diagnostic: trace WHY SolveSlipRatePsi returns V=0 at VW-core DOFs.
#
# Added stderr diagnostics at two levels:
#   1. [FRIC-FALLBACK] in SolveSlipRatePsi: logs when primary bracket fails
#      at VW-core DOFs (a<0.01, psi/a>100) and whether the fallback succeeds.
#   2. [RECOMPUTE-V0] in RecomputeSlipRate: logs when V_abs=0 despite
#      tau_abs>1kPa at VW-core DOFs.
#
# The issue reproduces within ~2 days of simulation time (~5000 steps).
# A 2hr dev run is sufficient.
#
# Check stderr for [FRIC-FALLBACK] and [RECOMPUTE-V0] messages.

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

RESULT_DIR="bp5/results_v63_diag_sliprate_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_v63_diag_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Diagnostic run: tracing V=0 at VW-core DOFs in RecomputeSlipRate"

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
      --verify
