#!/bin/bash
#SBATCH -J bp5_v41a_ip_p1      # Job name
#SBATCH -o bp5_v41a_ip_p1_%j.out
#SBATCH -e bp5_v41a_ip_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v41a: Baseline — IP p=1, 1000m mesh, full 1800 yr
#
# Same code as v40 (all IP traction sign bugs fixed).
# Full 1800yr run to capture 7-8 earthquake cycles for benchmark comparison.
# v40 only reached ~402yr before walltime.
#
# Expected: recurrence ~250yr, 7 events, 20% strike slip deficit (p=1 effect)

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

mkdir -p bp5/results_v41a_1000m_ip_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v41a_1000m_ip_p1_full \
      --output-prefix bp5_v41a_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
