#!/bin/bash
#SBATCH -J bp5_v35_p1_full     # Job name
#SBATCH -o bp5_v35_p1_full_%j.out
#SBATCH -e bp5_v35_p1_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v35: Full run — p=1, BR2, 1000m mesh, 1800 yr
#
# Purpose: Regression check after v35 per-quadrature-point traction fix.
#   v34 p=1 baseline: ~271 yr recurrence. For p=1 (linear elements), the
#   centroid evaluation is exact (constant gradient), so v35 should give
#   similar results. The interior Dirichlet skeleton fix (cross-element
#   lifting + 0.5 factor correction) may cause minor changes.
#
# v35 fixes:
#   1. ComputeTraction: per-quad-point {sigma.n} + BR2 correction (exact for p=1)
#   2. BR2 correction uses per-point shapes instead of avg_shapes (exact for p=1)
#   3. Interior Dirichlet skeleton: cross-element lifting + 0.5 factor fix
#
# Expected: Similar recurrence to v34 p=1 (~271 yr). Any difference comes
#   from fix #3 (interior Dirichlet scaling correction).

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

mkdir -p bp5/results_v35_1000m_br2_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v35_1000m_br2_p1_full \
      --output-prefix bp5_v35_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
