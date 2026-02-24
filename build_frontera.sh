#!/bin/bash
# Build MFEM + SEAS miniapp on Frontera (TACC)
#
# Usage:
#   bash build_frontera.sh
#
# Prerequisites:
#   module load hypre/2.31.0
#   module load mumps/5.3
#   module load parmetis

set -e

# Load modules if not already loaded
module load hypre/2.31.0 2>/dev/null || true
module load mumps/5.3 2>/dev/null || true
module load parmetis 2>/dev/null || true

# Verify modules
echo "=== Checking environment ==="
for var in TACC_HYPRE_INC TACC_HYPRE_LIB TACC_MUMPS_INC TACC_MUMPS_LIB TACC_PARMETIS_INC TACC_PARMETIS_LIB MKLROOT; do
    val=$(eval echo \$$var)
    if [ -z "$val" ]; then
        echo "ERROR: $var is not set. Check module loads."
        exit 1
    fi
    echo "  $var = $val"
done

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo ""
echo "=== Configuring MFEM ==="
make config \
  MFEM_USE_MPI=YES \
  MFEM_USE_METIS=YES \
  MFEM_USE_METIS_5=YES \
  MFEM_USE_LAPACK=YES \
  MFEM_USE_MUMPS=YES \
  HYPRE_OPT="-I${TACC_HYPRE_INC}" \
  HYPRE_LIB="-L${TACC_HYPRE_LIB} -lHYPRE" \
  METIS_OPT="-I${TACC_PARMETIS_INC} -DMETIS_EXPORT=" \
  METIS_LIB="-L${TACC_PARMETIS_LIB} -lparmetis -lmetis" \
  MUMPS_OPT="-I${TACC_MUMPS_INC}" \
  MUMPS_LIB="-L${TACC_MUMPS_LIB} -ldmumps -lmumps_common -lpord -lesmumps -lptscotch -lscotch -lptscotcherr -lscotcherr -L${TACC_MKL_LIB} -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lmpifort -lifcore" \
  LAPACK_OPT="-I${MKLROOT}/include" \
  LAPACK_LIB="-L${TACC_MKL_LIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread"

echo ""
echo "=== Building MFEM library ==="
make -j8

echo ""
echo "=== Building SEAS BP1 miniapp ==="
cd miniapps/seas
make seas_bp1_full -j8

echo ""
echo "=== Build complete ==="
echo "Binary: $(pwd)/seas_bp1_full"
