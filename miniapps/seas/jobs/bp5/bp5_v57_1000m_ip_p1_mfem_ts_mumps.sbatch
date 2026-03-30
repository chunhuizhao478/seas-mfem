#!/bin/bash
#SBATCH -J bp5_v57_mts_p1_mumps
#SBATCH -o bp5_v57_mts_p1_mumps_%j.out
#SBATCH -e bp5_v57_mts_p1_mumps_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# BP5 v57: 1000 m, IP, p=1, MFEM time stepper with MUMPS BLR
#
# Purpose:
# - use the existing MFEM time stepper, not PETSc TS
# - keep the BP5 elasticity solve on MUMPS BLR
#
# Assumption:
# - MFEM was already built before job submission with MFEM_USE_MUMPS=YES
# - only the `seas_bp5_full` miniapp may need an incremental rebuild

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
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }
mkdir -p bp5/results_v57_1000m_ip_p1_mfem_ts_mumps

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v57_1000m_ip_p1_mfem_ts_mumps \
      --output-prefix bp5_v57_mfem_ts_p1_mumps \
      --checkpoint-interval 1000 \
      --write-every-step
