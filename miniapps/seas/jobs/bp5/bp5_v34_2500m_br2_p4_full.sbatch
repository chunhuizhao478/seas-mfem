#!/bin/bash
#SBATCH -J bp5_v34_p4_full     # Job name
#SBATCH -o bp5_v34_p4_full_%j.out
#SBATCH -e bp5_v34_p4_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v34: Full run — p=4, BR2, 2500m mesh, 1800 yr
#
# Purpose: Match Tandem's SCEC "uphoff.2" submission (dx=2.5km, P4).
#   Test whether P4 on a moderate mesh resolves the recurrence interval
#   and dip stress mismatch.
#
# DOFs: P4 on 13503 tets → 35 DOFs/tet × 3 components = 1.42M DOFs
#   (smaller than our 1km P2 run at 1.90M DOFs)
#
# Key physics: BR2 traction correction scales as O(h^p).
#   At P4 with h=2.5km: the DG residual is greatly reduced vs P1/P2.

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

mkdir -p bp5/results_v34_2500m_br2_p4_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v34_2500m_br2_p4_full \
      --output-prefix bp5_v34_p4 \
      --checkpoint-interval 1000 \
      --write-every-step
