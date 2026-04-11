#!/bin/bash
#SBATCH -J bp5_v33_p1_base     # Job name
#SBATCH -o bp5_v33_p1_base_%j.out
#SBATCH -e bp5_v33_p1_base_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes (same as v32)
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v33: Baseline — p=1 on 1000m mesh, BR2, full 1800 yr
#
# Purpose: Regression check. Confirm v33 code at p=1 reproduces v32 result
#   (~435 yr recurrence). This proves the general-order changes are no-op
#   at order=1 and any improvement at p=2 is due to higher-order accuracy,
#   not a code change that affects p=1.
#
# Expected: identical to v32_1000m_br2_full result.

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

mkdir -p bp5/results_v33_1000m_br2_p1_baseline

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v33_1000m_br2_p1_baseline \
      --output-prefix bp5_v33_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
