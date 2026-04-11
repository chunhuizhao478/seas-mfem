#!/bin/bash
#SBATCH -J bp5_v37_p1_full     # Job name
#SBATCH -o bp5_v37_p1_full_%j.out
#SBATCH -e bp5_v37_p1_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v37: Full run — p=1, BR2, 1000m mesh, 1800 yr
#
# Fix: Interior Dirichlet face_int2 sign error. The BR2 lifting of the
#   prescribed jump u_D = Vp*t into element 2 used the wrong sign (-=
#   instead of +=), causing:
#   - v34 (no cross-element): elem2 pushed in wrong direction → 277 yr
#   - v36 (cross-element + 0.5): eval1+eval2 cancelled → ~0 penalty → 363+ yr
#   - v37 (corrected sign): both elements driven correctly → expected < 240 yr
#
# Previous: v34 = 277 yr, v36 = 363+ yr, Tandem = 240 yr

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

mkdir -p bp5/results_v37_1000m_br2_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v37_1000m_br2_p1_full \
      --output-prefix bp5_v37_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
