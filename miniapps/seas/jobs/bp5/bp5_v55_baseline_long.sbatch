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
# Strategy: checkout baseline commit into a temp dir, copy only the
# miniapps/seas source files, and build against the existing MFEM library.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

SEAS_DIR=/scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas
REPO_DIR=/scratch2/10024/zhaochun/seas-project/seas-mfem
BASELINE_COMMIT=d27c9cf
BUILD_DIR=/tmp/seas_baseline_$$

# Extract only miniapps/seas source from the baseline commit
mkdir -p "$BUILD_DIR"
cd "$REPO_DIR"
git archive "$BASELINE_COMMIT" -- miniapps/seas | tar -x -C "$BUILD_DIR"

cd "$BUILD_DIR/miniapps/seas"

# Symlink to the existing MFEM library (two levels up from miniapps/seas)
ln -sf "$REPO_DIR/libmfem.a" ../../libmfem.a
ln -sf "$REPO_DIR/config" ../../config
ln -sf "$REPO_DIR/fem" ../../fem
ln -sf "$REPO_DIR/general" ../../general
ln -sf "$REPO_DIR/linalg" ../../linalg
ln -sf "$REPO_DIR/mesh" ../../mesh

# Build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

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

# Copy results back to main tree for easy access
cp -r bp5/results_v55_baseline_long "$SEAS_DIR/bp5/" 2>/dev/null || true

# Cleanup
rm -rf "$BUILD_DIR"
