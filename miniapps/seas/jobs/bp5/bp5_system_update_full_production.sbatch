#!/bin/bash
#SBATCH -J bp5_suf
#SBATCH -o bp5_system_update_full_production_%j.out
#SBATCH -e bp5_system_update_full_production_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# system_update branch: full production run (400 ranks, ~1800 yr).
# Final verification after all 7 refactoring phases (0-6).
# Compares against Phase 0 full production golden (golden_tandem_p1_1000m).
#
# Before submitting, on the login node:
#   cd /scratch2/10024/zhaochun/seas-project/seas-mfem
#   git fetch origin && git checkout system_update && git pull
#   cd miniapps/seas && sbatch jobs/bp5/bp5_system_update_full_production.sbatch
#
# What to check:
#   1. Recurrence interval ~240 yr (same as v62/Phase 0 golden)
#   2. Coseismic slip ~3-5 m per event
#   3. Station time series match golden within recompilation noise
#   4. No blowup, no dt→0, no bracket failures
#
# Post-run comparison:
#   python3 scripts/regression_check.py \
#       --sim-dir $RESULT_DIR --sim-prefix $OUTPUT_PREFIX \
#       --ref-dir bp5/benchmark_data/golden_tandem_p1_1000m \
#       --tolerance 1e-2

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

RESULT_DIR="bp5/results_system_update_full_production_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_system_update_full_production_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "system_update full production: 400 ranks, 48hr, ~1800 yr earthquake cycles"
echo "Phases 0-6 complete. Final verification run."

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --ref-dir bp5/benchmark_data/golden_tandem_p1_1000m \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --checkpoint-interval 5000 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --regression-tolerance 1e-2 \
      --verify
