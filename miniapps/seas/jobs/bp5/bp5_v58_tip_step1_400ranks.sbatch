#!/bin/bash
#SBATCH -J bp5_tip_step1
#SBATCH -o bp5_tip_step1_%j.out
#SBATCH -e bp5_tip_step1_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# One-step tip reproducer: isolate first-solve asymmetry at fault tip
# Runs: init → apply first nonzero slip → one K*u=b solve → dump → exit
# Should complete in minutes, not hours.

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

cd miniapps/seas
make clean
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }
mkdir -p bp5/results_v58_tip_step1

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v58_tip_step1 \
      --output-prefix bp5_v58_tip_step1 \
      --checkpoint-interval 0 \
      --diag-tip-step1
