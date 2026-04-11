#!/bin/bash
#SBATCH -J bp5_v63f
#SBATCH -o bp5_v63_ftz_fix_%j.out
#SBATCH -e bp5_v63_ftz_fix_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5 v63 FTZ fix: Replace nextafter(0,1) fallback bracket with -300.
#
# Root cause: Intel's FTZ mode flushes nextafter(0.0,1.0) to 0.0, making
# Va_min = log10(0) = -inf.  Brent solver with [-inf, Vb] converges to
# V=0 instead of the true root at ~1e-35.
#
# Fix: Va_min = -300.0 (hard-coded, within normal double range).
#
# What to check:
#   1. dp+10 stations: V_strike should stay ~-30 to -40 during interseismic
#   2. No -300 sentinel values in interseismic V output
#   3. dp+00, dp+22 stations: unchanged from v62

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

RESULT_DIR="bp5/results_v63_ftz_fix_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_v63_ftz_fix_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Fix: Va_min=-300 (replaces -inf from FTZ-flushed nextafter)"

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
