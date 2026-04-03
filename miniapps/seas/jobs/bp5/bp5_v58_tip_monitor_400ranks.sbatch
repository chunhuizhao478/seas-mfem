#!/bin/bash
#SBATCH -J bp5_tipmon
#SBATCH -o bp5_tipmon_%j.out
#SBATCH -e bp5_tipmon_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# BP5 v58: PETSc TS production with tip DOF monitor
# Tracks psi, V, traction, sigma_n_eff at fault tip DOFs
# to catch the feedback loop that causes blowup at ~0.5yr.

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load parmetis
module load petsc/3.15
module list

# CG solver is built into MFEM — no external module needed.
# MUMPS lib path kept in LD_LIBRARY_PATH in case the PETSc build links it.
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_PETSC_LIB}:${LD_LIBRARY_PATH}"

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

export JOBS=8

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || {
    echo "ERROR: MFEM was not built with PETSc enabled"
    exit 1
}

cd miniapps/seas
make clean
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }
mkdir -p bp5/results_v58_tipmon

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver cg \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v58_tipmon \
      --output-prefix bp5_v58_tipmon \
      --checkpoint-interval 0 \
      --diag-tip-monitor \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg
