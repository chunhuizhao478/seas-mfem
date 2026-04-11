#!/bin/bash
#SBATCH -J bp5_v36_p1_full     # Job name
#SBATCH -o bp5_v36_p1_full_%j.out
#SBATCH -e bp5_v36_p1_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v36: Full run — p=1, BR2, 1000m mesh, 1800 yr
#
# Purpose: PRIMARY TEST for v36 interior Dirichlet BR2 skeleton fix.
#   v34 p=1 showed ~271 yr recurrence (vs Tandem ~197 yr) and growing
#   antisymmetric tau_dip. The tau_dip growth was traced to the interior
#   Dirichlet faces introduced in v34 — the BR2 skeleton loading had:
#     (a) Missing cross-element lifts (elem1 saw only its own lift)
#     (b) Missing 0.5 factor (2x penalty overshoot vs bilinear form)
#
# v36 fixes (re-applied from v35 with correct analysis):
#   1. Interior Dirichlet BR2 skeleton: cross-element lifting + 0.5 factor
#      (CRITICAL — this is the only fix that affects p=1)
#   2. ComputeTraction: per-quad-point {sigma.n} + BR2 correction (exact for p=1)
#   3. BR2 correction uses per-point shapes instead of avg_shapes (exact for p=1)
#
# Expected: tau_dip growth eliminated. Recurrence should shift toward
#   Tandem's ~197 yr reference (v34 was ~271 yr with 2x penalty overshoot).

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

mkdir -p bp5/results_v36_1000m_br2_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v36_1000m_br2_p1_full \
      --output-prefix bp5_v36_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
