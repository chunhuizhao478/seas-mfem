#!/bin/bash
#SBATCH -J bp5_v48f_1k_exact     # Job name
#SBATCH -o bp5_v48f_1k_exact_%j.out
#SBATCH -e bp5_v48f_1k_exact_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 4                         # 4 nodes (more memory for exact MUMPS)
#SBATCH -n 8                         # 8 ranks (spread across nodes for memory)
#SBATCH -t 04:00:00                  # 4 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v48f: IP p=2, 1000m mesh, 8 ranks, EXACT MUMPS
#
# The 1000m mesh at p=2 has ~1.9M vector DOFs. With 8 ranks, each rank
# handles ~240K DOFs. Exact MUMPS needs ~O(n^1.5) memory for 3D problems.
# Using 4 nodes (8 ranks spread across them) for maximum memory per rank.
#
# This tests whether the ×3 penalty works at production resolution
# with an accurate solver.
#
# Previous 1000m results:
#   v47b (400 rank, BLR): IMMEDIATE CRASH
#   v47l was on 2500m not 1000m (serial 1000m too large)
#
# If V grows → ×3 penalty validated at 1000m, BLR is the blocker
# If V decays → 1000m resolution may have additional issues
# If OOM → need more nodes or use BLR with tight tolerance

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

mkdir -p bp5/results_v48f_1000m_mumps_exact_8rank

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v48f_1000m_mumps_exact_8rank \
      --output-prefix bp5_v48f_1k_exact \
      --write-every-step
