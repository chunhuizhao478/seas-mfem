#!/bin/bash
# setup_mfem.sh - Configure and build MFEM with MUMPS support for SEAS miniapp
#
# Usage:
#   ./setup_mfem.sh                     # Auto-detect environment
#   ./setup_mfem.sh /path/to/mfem       # Use specific MFEM source directory
#
# This script:
#   1. Detects MPI Fortran library name (OpenMPI vs MPICH)
#   2. Detects SuiteSparse KLU library
#   3. Configures MFEM with MUMPS + SuiteSparse + Hypre + METIS
#   4. Builds MFEM
#   5. Copies config files and library to seas-mfem
#
# Prerequisites (install via conda, spack, or module load):
#   - MPI (openmpi or mpich)
#   - hypre
#   - metis (v5)
#   - mumps (with -ldmumps -lmumps_common -lpord)
#   - scalapack
#   - suitesparse (umfpack, klu)
#   - lapack/blas

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MFEM_SRC="${1:-/Users/chunhuizhao/projects/mfem}"

if [ ! -f "$MFEM_SRC/makefile" ] && [ ! -f "$MFEM_SRC/Makefile" ]; then
    echo "ERROR: MFEM source not found at $MFEM_SRC"
    echo "Usage: $0 /path/to/mfem"
    exit 1
fi

echo "=== MFEM Setup for SEAS Miniapp ==="
echo "  MFEM source: $MFEM_SRC"
echo "  SEAS-MFEM:   $SCRIPT_DIR"

# --- Detect library prefix (conda env, spack, system, module) ---
# Try to find Hypre header to determine the library prefix
LIB_PREFIX=""
for candidate in \
    "$CONDA_PREFIX" \
    "$MFEM_DEV_PREFIX" \
    "/usr/local" \
    "/usr"; do
    if [ -n "$candidate" ] && [ -f "$candidate/include/HYPRE.h" ]; then
        LIB_PREFIX="$candidate"
        break
    fi
done

if [ -z "$LIB_PREFIX" ]; then
    echo "ERROR: Cannot find Hypre headers. Set CONDA_PREFIX or MFEM_DEV_PREFIX."
    exit 1
fi
echo "  Library prefix: $LIB_PREFIX"

# --- Detect MPI Fortran library ---
MPI_FORT_LIB=""
if [ -f "$LIB_PREFIX/lib/libmpi_mpifh.dylib" ] || [ -f "$LIB_PREFIX/lib/libmpi_mpifh.so" ]; then
    MPI_FORT_LIB="-lmpi_mpifh"  # OpenMPI
    echo "  MPI Fortran: OpenMPI (libmpi_mpifh)"
elif [ -f "$LIB_PREFIX/lib/libmpifort.dylib" ] || [ -f "$LIB_PREFIX/lib/libmpifort.so" ]; then
    MPI_FORT_LIB="-lmpifort"    # MPICH
    echo "  MPI Fortran: MPICH (libmpifort)"
elif [ -f "$LIB_PREFIX/lib/libmpi_f90.dylib" ] || [ -f "$LIB_PREFIX/lib/libmpi_f90.so" ]; then
    MPI_FORT_LIB="-lmpi_f90"    # Older OpenMPI
    echo "  MPI Fortran: OpenMPI legacy (libmpi_f90)"
else
    echo "WARNING: Cannot detect MPI Fortran library. Trying without it."
    MPI_FORT_LIB=""
fi

# --- Detect SuiteSparse include path ---
SUITESPARSE_OPT="-I$LIB_PREFIX/include"
if [ -f "$LIB_PREFIX/include/suitesparse/umfpack.h" ]; then
    SUITESPARSE_OPT="-I$LIB_PREFIX/include/suitesparse"
    echo "  SuiteSparse: $LIB_PREFIX/include/suitesparse"
elif [ -f "$LIB_PREFIX/include/umfpack.h" ]; then
    echo "  SuiteSparse: $LIB_PREFIX/include"
else
    echo "WARNING: Cannot find umfpack.h. SuiteSparse may not link."
fi

# --- Build SuiteSparse library list (include KLU if available) ---
SUITESPARSE_LIBS="-lumfpack"
if [ -f "$LIB_PREFIX/lib/libklu.dylib" ] || [ -f "$LIB_PREFIX/lib/libklu.so" ] || [ -f "$LIB_PREFIX/lib/libklu.a" ]; then
    SUITESPARSE_LIBS="$SUITESPARSE_LIBS -lklu"
    echo "  KLU: found"
fi
SUITESPARSE_LIBS="$SUITESPARSE_LIBS -lcholmod -lcolamd -lamd -lcamd -lccolamd -lsuitesparseconfig"

# --- Detect MUMPS ---
MUMPS_AVAILABLE="NO"
if [ -f "$LIB_PREFIX/include/dmumps_c.h" ]; then
    MUMPS_AVAILABLE="YES"
    MUMPS_LIBS="-ldmumps -lmumps_common -lpord -lscalapack -llapack -lblas $MPI_FORT_LIB"
    echo "  MUMPS: found"
else
    echo "  MUMPS: not found (will use HypreILU fallback)"
fi

# --- Configure MFEM ---
echo ""
echo "=== Configuring MFEM ==="
cd "$MFEM_SRC"

CONFIG_ARGS=(
    "MFEM_USE_MPI=YES"
    "MFEM_USE_METIS=YES"
    "MFEM_USE_METIS_5=YES"
    "MFEM_USE_LAPACK=YES"
    "MFEM_USE_SUITESPARSE=YES"
    "HYPRE_OPT=-I$LIB_PREFIX/include"
    "HYPRE_LIB=-L$LIB_PREFIX/lib -lHYPRE"
    "METIS_OPT=-I$LIB_PREFIX/include"
    "METIS_LIB=-L$LIB_PREFIX/lib -lmetis"
    "LAPACK_OPT=-I$LIB_PREFIX/include"
    "LAPACK_LIB=-L$LIB_PREFIX/lib -llapack -lblas"
    "SUITESPARSE_OPT=$SUITESPARSE_OPT"
    "SUITESPARSE_LIB=-L$LIB_PREFIX/lib $SUITESPARSE_LIBS"
)

if [ "$MUMPS_AVAILABLE" = "YES" ]; then
    CONFIG_ARGS+=(
        "MFEM_USE_MUMPS=YES"
        "MUMPS_OPT=-I$LIB_PREFIX/include"
        "MUMPS_LIB=-L$LIB_PREFIX/lib $MUMPS_LIBS"
    )
fi

echo "  Running: make config ${CONFIG_ARGS[*]}"
make config "${CONFIG_ARGS[@]}"

# --- Build MFEM ---
echo ""
echo "=== Building MFEM ==="
NPROC=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
make -j"$NPROC"

# --- Copy to seas-mfem ---
echo ""
echo "=== Copying to seas-mfem ==="
cp "$MFEM_SRC/config/config.mk" "$SCRIPT_DIR/config/config.mk"
cp "$MFEM_SRC/config/_config.hpp" "$SCRIPT_DIR/config/_config.hpp"
cp "$MFEM_SRC/libmfem.a" "$SCRIPT_DIR/libmfem.a"
echo "  Copied config.mk, _config.hpp, libmfem.a"

# --- Verify ---
echo ""
echo "=== Verification ==="
if grep -q "MFEM_USE_MUMPS.*YES" "$SCRIPT_DIR/config/config.mk"; then
    echo "  MUMPS:        enabled"
else
    echo "  MUMPS:        disabled (using HypreILU fallback)"
fi
if grep -q "MFEM_USE_SUITESPARSE.*YES" "$SCRIPT_DIR/config/config.mk"; then
    echo "  SuiteSparse:  enabled"
fi
if grep -q "MFEM_USE_MPI.*YES" "$SCRIPT_DIR/config/config.mk"; then
    echo "  MPI:          enabled"
fi

echo ""
echo "=== Done ==="
echo "Now build the SEAS miniapp:"
echo "  cd miniapps/seas && make seas_bp1_full -j$NPROC"
