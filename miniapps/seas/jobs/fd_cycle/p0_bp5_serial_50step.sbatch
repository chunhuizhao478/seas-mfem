#!/bin/bash
#SBATCH -J fdc_p0_bp5_vs1
#SBATCH -o fdc_p0_bp5_serial_50step_%j.out
#SBATCH -e fdc_p0_bp5_serial_50step_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# Fully dynamic BP5 cycle, phase P0 (plan rev 3): re-check the QD gold
# configuration on the branch point.  A guarded copy of
# jobs/bp5/bp5_verify_serial.sbatch (the 50-step golden, Frontera job 7643780):
# same modules, same seas_bp5_full flags, 1 rank, 50 steps; then a check
# against bp5/benchmark_data/golden_serial_1r_50step.
#
# The branch carries the station sigma_n fix ef70517 (gpu-qd-bp5), which the
# golden predates: station log10 V and tau now use the solver's elastic
# sigma_n, so they differ from the golden by design (up to about 2.4e-2 in
# log10 V at -16dp+10 in 50 steps); time, slip, state and the global file do
# not change.  So the pass/fail check is tests/cycle/compare_runs.py with
# relative L2 <= 1e-12 (the measure of scripts/regression_check.py) on
# station columns time, slip strike, slip dip and state and on every column
# of the global file.  The shared regression check then runs informational
# (all columns, always exits 0) to print the V and tau differences.
# Submit only on Chunhui Zhao's typed go for this job.
#
# Differences from the original, none of which touch the numerics: the
# GitHub-first guard (CLAUDE.md), the build through
# jobs/fd_cycle/build_frontera.sh (cycle/cycle.mk), output under $SCRATCH,
# and the regression check run here.

export LC_ALL=C
export LANG=C
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

ROOT="${FDC_ROOT:-$WORK/seas-mfem-fdcycle}"

# ---- GitHub-first guard (CLAUDE.md) ----
cd "$ROOT" || { echo "ERROR: no checkout at $ROOT"; exit 2; }
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "commit $(git rev-parse HEAD)"
echo "branch $(git rev-parse --abbrev-ref HEAD)"

# ---- build from this commit (modules and guard repeated inside) ----
(cd miniapps/seas && make clean > /dev/null)
bash miniapps/seas/jobs/fd_cycle/build_frontera.sh || { echo "BUILD FAILED"; exit 1; }

module load intel/19.1.1 impi/19.0.9 hypre/2.31.0 mumps/5.3 parmetis petsc/3.15 fftw3/3.3.8
export LD_LIBRARY_PATH="$ROOT/extern/hdf5/install/lib:$ROOT/extern/h5z-zfp/install/lib:${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"

cd "$ROOT/miniapps/seas"
RESULT_DIR="$SCRATCH/fd_cycle/p0_bp5_serial_50step_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="fdc_p0_bp5_serial_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"
echo "Output dir: ${RESULT_DIR}"

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
rc=$?
echo "seas_bp5_full exit ${rc}"
[ "${rc}" = 0 ] || exit "${rc}"

python3 tests/cycle/compare_runs.py "${RESULT_DIR}" "${OUTPUT_PREFIX}" \
        bp5/benchmark_data/golden_serial_1r_50step bp5_verify_serial_job7643780 \
        --tol 1e-12 --l2 --station-columns 0,1,2,7
rc=$?
echo "compare_runs exit ${rc} (0 = time, slip, state and global within 1e-12)"
make -f cycle/cycle.mk regression-check SIM_DIR="${RESULT_DIR}" SIM_PREFIX="${OUTPUT_PREFIX}" \
     REF_DIR=bp5/benchmark_data/golden_serial_1r_50step
exit "${rc}"
