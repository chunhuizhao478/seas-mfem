#!/bin/bash
#SBATCH -J bp5_v34_p2_full     # Job name
#SBATCH -o bp5_v34_p2_full_%j.out
#SBATCH -e bp5_v34_p2_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes
#SBATCH -n 800                       # 50 tasks/node x 16 nodes
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v34: Full run — p=2, BR2, 1000m mesh, 1800 yr
#
# Purpose: Test v34 BC fix at p=2. v33 p=2 showed VS zone slower than p=1.
#   With the corrected BC (u=Vp*t on non-fault Y=0), p=2 may now converge
#   toward Tandem's reference.
#
# Fix: Non-fault Y=0 interior faces now get u_D = (Vp*t, 0, 0) instead
#   of u_D = (0, 0, 0). This matches Tandem's bp5.lua boundary().

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

mkdir -p bp5/results_v34_1000m_br2_p2_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v34_1000m_br2_p2_full \
      --output-prefix bp5_v34_p2 \
      --checkpoint-interval 1000 \
      --write-every-step
