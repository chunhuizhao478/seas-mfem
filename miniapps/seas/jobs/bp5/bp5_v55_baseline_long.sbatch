#!/bin/bash
#SBATCH -J bp5_v55_base_long         # Job name
#SBATCH -o bp5_v55_base_long_%j.out
#SBATCH -e bp5_v55_base_long_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v55 baseline LONG run: pre-D5 code, run to t>200s
#
# Purpose: determine if the pre-fix code (v53/v54 era) also recovers
# and nucleates when run long enough. Previous v53 run only reached
# t=38s at station output before wall-time kill.
#
# Uses git worktree to build from the baseline commit (d27c9cf)
# without modifying the main working tree.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

SEAS_DIR=/scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas
BASELINE_COMMIT=d27c9cf

# Create a temporary worktree for the baseline commit
WORKTREE_DIR=/scratch2/10024/zhaochun/seas-project/seas-mfem-baseline-worktree
cd /scratch2/10024/zhaochun/seas-project/seas-mfem

# Clean up any previous worktree
git worktree remove --force "$WORKTREE_DIR" 2>/dev/null || true
git worktree add "$WORKTREE_DIR" "$BASELINE_COMMIT"

cd "$WORKTREE_DIR/miniapps/seas"

make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

# Use mesh and benchmark data from the main tree
mkdir -p bp5/results_v55_baseline_long

ibrun ./seas_bp5_full \
      --mesh "$SEAS_DIR/bp5/mesh/reference/bp5_tandem.msh" \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir "$SEAS_DIR/bp5/benchmark_data" \
      --output-dir bp5/results_v55_baseline_long \
      --output-prefix bp5_v55_base_long \
      --checkpoint-interval 1000 \
      --write-every-step

# Clean up worktree after run
cd /scratch2/10024/zhaochun/seas-project/seas-mfem
git worktree remove --force "$WORKTREE_DIR" 2>/dev/null || true
