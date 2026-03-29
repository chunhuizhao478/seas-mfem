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
# Clones the repo at the baseline commit into $SCRATCH/tmp, builds
# everything from scratch, and runs. Does NOT touch the main working tree.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

REPO_DIR=/scratch2/10024/zhaochun/seas-project/seas-mfem
BASELINE_COMMIT=d27c9cf
BUILD_DIR=$SCRATCH/seas_baseline_build_$$

# Clone at baseline commit
git clone --single-branch --branch feature/elasticity "$REPO_DIR" "$BUILD_DIR"
cd "$BUILD_DIR"
git checkout "$BASELINE_COMMIT"

# Build MFEM library
cp config/defaults.mk config/config.mk
# Use the same config as the main build
cp "$REPO_DIR/config/config.mk" config/config.mk 2>/dev/null || true
make lib -j8 || { echo "MFEM LIB BUILD FAILED"; exit 1; }

# Build seas
cd miniapps/seas
make seas_bp5_full || { echo "SEAS BUILD FAILED"; exit 1; }

# Use mesh/data from main tree
SEAS_DIR=$REPO_DIR/miniapps/seas
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
echo "Also copied to: $SEAS_DIR/bp5/results_v55_baseline_long/"
