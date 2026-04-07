#!/bin/bash
#SBATCH -J bp5_tss
#SBATCH -o bp5_v59_tssolve_%j.out
#SBATCH -e bp5_v59_tssolve_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# BP5 v59 TSSolve: Replace manual TSStep loop with Tandem-style single
# TSSolve call.  All per-step monitoring (I/O, earthquake detection,
# console output) now runs inside a PETSc TSMonitor callback.
#
# Root cause fix: the manual TSStep loop called TSSetTimeStep before
# every step, which interfered with PETSc's FSAL optimisation and
# adaptive dt management, producing earthquake cycles 5x slower than
# Tandem.  TSSolve lets PETSc manage everything internally.
#
# Expected: first earthquake at ~225 yr (matching Tandem), not ~1191 yr.

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

RESULT_DIR="bp5/results_v59_tssolve_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_v59_tssolve_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo "Fix: TSSolve replaces TSStep loop (Tandem-style time stepping)"

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --tfinal 1.0 \
      --max-steps 1 \
      --ref-dir bp5/benchmark_data \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --checkpoint-interval 0 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --paraview-every 1
