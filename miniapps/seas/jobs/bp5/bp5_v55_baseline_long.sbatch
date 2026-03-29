#!/bin/bash
#SBATCH -J bp5_v55_base_long         # Job name
#SBATCH -o bp5_v55_base_long_%j.out
#SBATCH -e bp5_v55_base_long_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v55 baseline LONG run: pre-D5 code to t>200s
#
# Strategy: checkout baseline seas files, build, run, then restore.
# The MFEM library itself doesn't change between commits.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

REPO_DIR=/scratch2/10024/zhaochun/seas-project/seas-mfem
SEAS_DIR=$REPO_DIR/miniapps/seas
BASELINE_COMMIT=d27c9cf
CURRENT_BRANCH=$(cd $REPO_DIR && git rev-parse HEAD)

cd $SEAS_DIR

# Save current seas files
git stash push -m "baseline-long-run" -- .

# Checkout baseline seas files only
cd $REPO_DIR
git checkout $BASELINE_COMMIT -- miniapps/seas/

cd $SEAS_DIR
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; cd $REPO_DIR; git checkout $CURRENT_BRANCH -- miniapps/seas/; git stash pop; exit 1; }

mkdir -p bp5/results_v55_baseline_long

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v55_baseline_long \
      --output-prefix bp5_v55_base_long \
      --checkpoint-interval 1000 \
      --write-every-step

# Restore current code
cd $REPO_DIR
git checkout $CURRENT_BRANCH -- miniapps/seas/
git stash pop 2>/dev/null || true
