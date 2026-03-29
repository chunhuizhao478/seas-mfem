#!/bin/bash
#SBATCH -J bp5_v55_base_p1           # Job name
#SBATCH -o bp5_v55_base_p1_%j.out
#SBATCH -e bp5_v55_base_p1_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v55 baseline: NO D8/D4 changes — reproduces v53/v54 behavior
#
# Purpose: confirm the pre-change code still produces the same nucleation-dying
# trajectory as v53 (V peaks at 0.0136 then dies to 0.0080 by t=38s).
# This is the control run for comparison with D8-only and D8+D4 runs.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v55_1000m_ip_p1_baseline

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v55_1000m_ip_p1_baseline \
      --output-prefix bp5_v55_base_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
