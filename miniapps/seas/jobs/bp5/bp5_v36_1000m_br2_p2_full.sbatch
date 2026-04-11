#!/bin/bash
#SBATCH -J bp5_v36_p2_full     # Job name
#SBATCH -o bp5_v36_p2_full_%j.out
#SBATCH -e bp5_v36_p2_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes
#SBATCH -n 800                       # 50 tasks/node x 16 nodes
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v36: Full run — p=2, BR2, 1000m mesh, 1800 yr
#
# Purpose: KEY TEST for combined v36 fixes at p=2.
#   v34 p=2 showed ~396 yr recurrence (significantly delayed vs p=1 ~271 yr),
#   caused by centroid-only stress evaluation and avg_shapes BR2 correction
#   (both exact for p=1 but wrong for p>=2), PLUS the 2x interior Dirichlet
#   penalty overshoot affecting all polynomial orders.
#
# v36 fixes:
#   1. Interior Dirichlet BR2 skeleton: cross-element lifting + 0.5 factor
#   2. ComputeTraction: per-quad-point {sigma.n} (gradient varies for p=2)
#   3. BR2 correction uses per-point shapes instead of avg_shapes
#
# Expected: Significantly shorter recurrence than v34's 396 yr. With all
#   three fixes, p=2 should converge TOWARD p=1, both moving closer to
#   Tandem's ~197 yr reference.

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

mkdir -p bp5/results_v36_1000m_br2_p2_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v36_1000m_br2_p2_full \
      --output-prefix bp5_v36_p2 \
      --checkpoint-interval 1000 \
      --write-every-step
