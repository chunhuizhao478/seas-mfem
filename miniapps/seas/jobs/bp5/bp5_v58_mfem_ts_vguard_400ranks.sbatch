#!/bin/bash
#SBATCH -J bp5_v58_mfem_vg
#SBATCH -o bp5_v58_mfem_vg_%j.out
#SBATCH -e bp5_v58_mfem_vg_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5 v58: 1000 m, IP, p=1, MFEM time stepper with V-guard, MUMPS
#
# Purpose:
# - use MFEM's own RK45 time stepper (not PETSc TS)
# - V-guard enabled (default factor 100) to prevent dt overshoot
# - compare with PETSc TS + dt_max run
#
# Note: --tandem-time-stepping is NOT used (it disables V-guard)
#       --petsc-ts is NOT used (uses MFEM stepper instead)

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module list

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${LD_LIBRARY_PATH}"

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

export JOBS=8

grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || {
    echo "ERROR: MFEM was not built with MUMPS enabled"
    exit 1
}

cd miniapps/seas
make clean
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }
mkdir -p bp5/results_v58_mfem_ts_vguard_400ranks

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v58_mfem_ts_vguard_400ranks \
      --output-prefix bp5_v58_mfem_vg \
      --checkpoint-interval 0
