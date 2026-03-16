#!/bin/bash
#SBATCH -J bp5_v35_p4_full     # Job name
#SBATCH -o bp5_v35_p4_full_%j.out
#SBATCH -e bp5_v35_p4_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes
#SBATCH -n 800                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v35: Full run — p=4, BR2, 2500m mesh, 1800 yr
#
# Purpose: CRITICAL TEST — v34 p=4 produced NO earthquakes in 1800 yr
#   (fault locked, 17-20% of expected interseismic slip). This was the
#   primary evidence of formulation bugs that grow with polynomial order.
#
# Matches Tandem's SCEC "uphoff.2" submission (dx=2.5km, P4).
# DOFs: P4 on 13503 tets → 35 DOFs/tet × 3 components = 1.42M DOFs
#
# v35 fixes:
#   1. ComputeTraction: per-quad-point {sigma.n} (large effect at p=4)
#   2. BR2 correction uses per-point shapes (large effect — Minv amplifies
#      high-degree modes that avg_shapes completely misrepresents)
#   3. Interior Dirichlet skeleton: cross-element lifting + 0.5 factor fix
#
# Expected: Earthquakes should nucleate. If the per-quadrature-point fix
#   is correct, p=4 at 2500m should produce recurrence intervals in the
#   ballpark of Tandem's ~197 yr reference. The DG residual at P4 with
#   h=2.5km is small, so the solution should be accurate.

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

mkdir -p bp5/results_v35_2500m_br2_p4_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v35_2500m_br2_p4_full \
      --output-prefix bp5_v35_p4 \
      --checkpoint-interval 1000 \
      --write-every-step
