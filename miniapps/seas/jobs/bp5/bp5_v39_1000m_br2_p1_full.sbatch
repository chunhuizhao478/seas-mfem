#!/bin/bash
#SBATCH -J bp5_v39_br2_p1     # Job name
#SBATCH -o bp5_v39_br2_p1_%j.out
#SBATCH -e bp5_v39_br2_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v39: BR2 method — p=1, 1000m mesh, 1800 yr
#
# Same code as v39 IP run but with BR2, to measure the effect of the v38b
# shared-face Dirichlet loading fix on BR2 recurrence.
#
# v38b fix: shared Y=0 non-fault faces on processor boundaries were missing
# Dirichlet loading since v34. This weakened tectonic driving near those faces.
#
# The IP traction sign fixes (v38c, v39) only affect the IP code path;
# BR2 traction is invariant to element ordering by construction.
#
# Previous: v37 (BR2) = 295 yr, Tandem (IP) = 240 yr

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

mkdir -p bp5/results_v39_1000m_br2_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v39_1000m_br2_p1_full \
      --output-prefix bp5_v39_br2_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
