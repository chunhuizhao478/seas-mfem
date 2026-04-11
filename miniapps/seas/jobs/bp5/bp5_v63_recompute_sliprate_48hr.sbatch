#!/bin/bash
#SBATCH -J bp5_v63_rV
#SBATCH -o bp5_v63_recompute_sliprate_%j.out
#SBATCH -e bp5_v63_recompute_sliprate_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5 v63: Test RecomputeSlipRate fix for stale slip-rate output.
#
# Root cause: station output used cached GetSlipRate() from an intermediate
# RK stage.  At locked fault DOFs (a=0.004, VW core), V ~ exp(-psi/a) is
# exponentially sensitive to psi — stale cache produced V=0 at dp+10
# stations during interseismic.
#
# Fix: RecomputeSlipRate() re-solves the friction equation from the accepted
# state at every output time, matching Tandem's RateAndState::state().
#
# What to check:
#   1. dp+10 stations: log10(V_strike) should stay ~-14 to -36 during
#      interseismic, NOT drop to -300
#   2. dp+00 stations: should be unchanged from v62
#   3. Recurrence interval and coseismic slip: should match v62

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

RESULT_DIR="bp5/results_v63_recompute_sliprate_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_v63_recompute_sliprate_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Fix: RecomputeSlipRate — re-solve friction at output time (commit e2732d7)"

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
      --checkpoint-interval 5000 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --paraview-dt 1.1e8 \
      --verify
