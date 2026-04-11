#!/bin/bash
#SBATCH -J bp5_v36_p6_full     # Job name
#SBATCH -o bp5_v36_p6_full_%j.out
#SBATCH -e bp5_v36_p6_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes
#SBATCH -n 800                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v36: Full run — p=6, BR2, 4000m mesh, 1800 yr
#
# Purpose: CRITICAL TEST — v34 p=6 produced NO earthquakes in 1291 yr
#   (fault locked). This configuration matches Tandem's SCEC "uphoff"
#   submission (dx=4km, P6), which produces excellent benchmark results.
#
# DOFs: P6 on 6724 tets -> 84 DOFs/tet x 3 components = 1.69M DOFs
#
# v36 fixes:
#   1. Interior Dirichlet BR2 skeleton: cross-element lifting + 0.5 factor
#   2. ComputeTraction: per-quad-point {sigma.n} (very large effect at p=6)
#   3. BR2 correction uses per-point shapes (very large effect — Minv
#      eigenvalues grow as O(p^3/h^3), massively amplifying the error
#      from avg_shapes at p=6)
#
# Expected: Earthquakes should nucleate. At P6 with h=4km, the DG residual
#   is negligible (spectral convergence). This should match Tandem's ~197 yr
#   recurrence closely — any remaining mismatch indicates a different bug.

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

mkdir -p bp5/results_v36_4000m_br2_p6_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_4000m.msh \
      --mesh-scale 1000 \
      --order 6 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v36_4000m_br2_p6_full \
      --output-prefix bp5_v36_p6 \
      --checkpoint-interval 1000 \
      --write-every-step
