#!/bin/bash
#SBATCH -J bp5_v36_p4_full     # Job name
#SBATCH -o bp5_v36_p4_full_%j.out
#SBATCH -e bp5_v36_p4_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes
#SBATCH -n 800                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v36: Full run — p=4, BR2, 2500m mesh, 1800 yr
#
# Purpose: CRITICAL TEST — v34 p=4 produced NO earthquakes in 1800 yr
#   (fault locked, 17-20% of expected interseismic slip). This was caused
#   by the combination of centroid-only traction evaluation (large error at
#   p=4), avg_shapes BR2 correction (Minv amplifies high-degree modes),
#   and 2x interior Dirichlet penalty overshoot.
#
# Matches Tandem's SCEC "uphoff.2" submission (dx=2.5km, P4).
# DOFs: P4 on 13503 tets -> 35 DOFs/tet x 3 components = 1.42M DOFs
#
# v36 fixes:
#   1. Interior Dirichlet BR2 skeleton: cross-element lifting + 0.5 factor
#   2. ComputeTraction: per-quad-point {sigma.n} (large effect at p=4)
#   3. BR2 correction uses per-point shapes (large effect — Minv amplifies
#      high-degree modes that avg_shapes completely misrepresents)
#
# Expected: Earthquakes should nucleate. With all fixes, p=4 at 2500m
#   should produce recurrence intervals approaching Tandem's ~197 yr.

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

mkdir -p bp5/results_v36_2500m_br2_p4_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v36_2500m_br2_p4_full \
      --output-prefix bp5_v36_p4 \
      --checkpoint-interval 1000 \
      --write-every-step
