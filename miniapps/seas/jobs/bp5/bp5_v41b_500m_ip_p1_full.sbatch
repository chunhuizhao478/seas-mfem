#!/bin/bash
#SBATCH -J bp5_v41b_500m_ip      # Job name
#SBATCH -o bp5_v41b_500m_ip_%j.out
#SBATCH -e bp5_v41b_500m_ip_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 32                        # 32 nodes (234k tets, ~4x 1000m)
#SBATCH -n 1600                      # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v41b: h-refinement — IP p=1, 500m mesh, full 1800 yr
#
# Tests whether halving mesh size reduces the 20% strike slip deficit.
# Expected: O(h) error halved → deficit reduced from ~20% to ~10%.
#
# Mesh: bp5_tandem_500m.msh (233,783 tets, ~4x the 1000m mesh)
# DOFs: ~2.8M (3 components × 4 DOFs/tet × 234k tets)
# Memory: ~4x baseline, requires 32 nodes.
#
# Key metrics to compare with v41a (1000m):
#   1. Strike slip at z=22km: should be ~90% of Tandem (vs 80% at 1000m)
#   2. tau_dip at center: should be ~2x smaller (vs ~5x larger than Tandem at 1000m)
#   3. Recurrence interval: should be closer to Tandem's 240yr

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

mkdir -p bp5/results_v41b_500m_ip_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_500m.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v41b_500m_ip_p1_full \
      --output-prefix bp5_v41b_500m_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
