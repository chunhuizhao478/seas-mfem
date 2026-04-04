#!/bin/bash
#SBATCH -J bp5_tndtq
#SBATCH -o bp5_tndtq_%j.out
#SBATCH -e bp5_tndtq_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# BP5 v58: PETSc TS with [MFEM-TQ] diagnostic at tip faces.
# Outputs per-QP traction data matching Tandem [TND-TQ] format
# for side-by-side comparison at fault tip faces (|x| > 49 km).

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

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || {
    echo "ERROR: MFEM was not built with PETSc enabled"
    exit 1
}
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || {
    echo "ERROR: MFEM was not built with MUMPS enabled"
    exit 1
}

cd miniapps/seas
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }
mkdir -p bp5/results_v58_tndtq

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v58_tndtq \
      --output-prefix bp5_v58_tndtq \
      --checkpoint-interval 0 \
      --diag-tip-monitor \
      --diag-tnd-tq \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg
