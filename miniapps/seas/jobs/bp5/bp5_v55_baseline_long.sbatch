#!/bin/bash
#SBATCH -J bp5_v55_base_long         # Job name
#SBATCH -o bp5_v55_base_long_%j.out
#SBATCH -e bp5_v55_base_long_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v55 baseline LONG: pre-D5 code, 48hr run
#
# Copies the existing MFEM build (libmfem.a + headers + config) from
# the main tree, then checks out only miniapps/seas at the baseline
# commit and builds the seas binary against the existing library.

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
BUILD_DIR=$SCRATCH/seas_baseline_build_$$

mkdir -p "$BUILD_DIR"

# Copy the MFEM library and all headers from the existing build
# (the MFEM library itself is identical across our commits —
# we only changed miniapps/seas/)
cp "$REPO_DIR/libmfem.a" "$BUILD_DIR/"
cp "$REPO_DIR/mfem.hpp" "$BUILD_DIR/"
cp "$REPO_DIR/mfem-performance.hpp" "$BUILD_DIR/" 2>/dev/null || true
cp -a "$REPO_DIR/config" "$BUILD_DIR/"
for dir in general linalg mesh fem; do
    cp -a "$REPO_DIR/$dir" "$BUILD_DIR/"
done

# Extract miniapps/seas from the baseline commit
cd "$REPO_DIR"
git archive "$BASELINE_COMMIT" -- miniapps/seas | tar -x -C "$BUILD_DIR"

# Build only the seas binary
cd "$BUILD_DIR/miniapps/seas"
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

# Copy results back
cp -r bp5/results_v55_baseline_long "$SEAS_DIR/bp5/" 2>/dev/null || true
echo "Results in: $BUILD_DIR/miniapps/seas/bp5/results_v55_baseline_long/"
