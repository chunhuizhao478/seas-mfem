#!/bin/bash
# Build MFEM + SEAS miniapp locally on macOS with conda (mfem-dev env)
#
# Prerequisites:
#   conda activate mfem-dev
#   (env should have: openmpi, hypre, mumps-mpi, parmetis, metis, lapack, scalapack)
#
# Usage:
#   conda activate mfem-dev
#   bash build_local.sh          # full build (MFEM lib + seas miniapp)
#   bash build_local.sh seas     # rebuild seas miniapp only (after MFEM is built)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Verify conda env has what we need
if [ -z "$CONDA_PREFIX" ]; then
    echo "ERROR: No conda environment active. Run: conda activate mfem-dev"
    exit 1
fi

for lib in libHYPRE libdmumps libparmetis libmetis; do
    if ! ls "$CONDA_PREFIX/lib/${lib}"* &>/dev/null; then
        echo "ERROR: ${lib} not found in $CONDA_PREFIX/lib/"
        echo "Install with: conda install -c conda-forge hypre mumps-mpi parmetis metis"
        exit 1
    fi
done

echo "=== Environment ==="
echo "  CONDA_PREFIX: $CONDA_PREFIX"
echo "  CXX: $(which mpicxx)"
echo ""

# Skip MFEM lib build if "seas" argument is given and lib already exists
if [ "$1" = "seas" ] && [ -f libmfem.a ]; then
    echo "=== Skipping MFEM lib build (libmfem.a exists) ==="
else
    echo "=== Configuring MFEM ==="
    make config \
      MFEM_USE_MPI=YES \
      MFEM_USE_METIS=YES \
      MFEM_USE_METIS_5=YES \
      MFEM_USE_LAPACK=YES \
      MFEM_USE_MUMPS=YES \
      HYPRE_OPT="-I${CONDA_PREFIX}/include" \
      HYPRE_LIB="-L${CONDA_PREFIX}/lib -lHYPRE" \
      METIS_OPT="-I${CONDA_PREFIX}/include -DMETIS_EXPORT=" \
      METIS_LIB="-L${CONDA_PREFIX}/lib -lparmetis -lmetis" \
      MUMPS_OPT="-I${CONDA_PREFIX}/include" \
      MUMPS_LIB="-L${CONDA_PREFIX}/lib -ldmumps -lmumps_common -lpord -lscalapack" \
      LAPACK_OPT="" \
      LAPACK_LIB="-framework Accelerate"

    echo ""
    echo "=== Building MFEM library ==="
    make -j8
fi

echo ""
echo "=== Building SEAS BP1 miniapp ==="
cd miniapps/seas
make seas_bp1_bdrload

echo ""
echo "=== Build complete ==="
echo "Binary: $(pwd)/seas_bp1_bdrload"
echo ""
echo "Quick test:"
echo "  mpirun -np 2 ./seas_bp1_bdrload --mesh bp1/mesh/bp1_ss_100m.msh --tfinal 1e5 --checkpoint-interval 0 --max-steps 20"
