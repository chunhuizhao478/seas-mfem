#!/bin/bash
# setup_mfem.sh - Configure and build MFEM (in-tree) with the SEAS-flavour
# patches for the seas miniapp.
#
# Usage:
#   ./setup_mfem.sh                          # Build IN PLACE (default;
#                                            # uses the in-tree patched MFEM
#                                            # source under this directory).
#   ./setup_mfem.sh /path/to/mfem            # Build an external MFEM tree
#                                            # and copy libmfem.a + configs
#                                            # over (legacy "vanilla MFEM"
#                                            # path; SetHDFCompression and
#                                            # other seas-specific patches
#                                            # are absent — VTKHDF compression
#                                            # paths will fail to link).
#
# Why in-tree?  The repository root (`/Users/chunhuizhao/projects/
# seas-mfem-paraview/`) is itself a patched MFEM source tree — `fem/`,
# `mesh/`, `linalg/`, etc. contain seas-specific modifications (the most
# critical of which adds `ParaViewHDFDataCollection::SetHDFCompression`
# to support H5Z-ZFP volume PV output, plan §Phase 2d.3 / §Phase 6.2).
# Building IN PLACE preserves those patches; building upstream MFEM and
# copying the resulting libmfem.a silently overwrites them and breaks
# `seas_bp5_full` linking under `MFEM_USE_HDF5=YES`.
#
# This script:
#   1. Detects MPI Fortran library name (OpenMPI vs MPICH)
#   2. Detects SuiteSparse KLU library
#   3. Detects optional HDF5 + PETSc and toggles MFEM_USE_HDF5 / MFEM_USE_PETSC
#   4. Configures MFEM with MUMPS + SuiteSparse + Hypre + METIS (+ HDF5/PETSc
#      when present)
#   5. Builds MFEM in place (in-tree mode) OR builds external tree + copies
#      libmfem.a/_config.hpp/config.mk over (legacy mode)
#
# Prerequisites (install via conda, spack, or module load):
#   - MPI (openmpi or mpich)
#   - hypre
#   - metis (v5)
#   - mumps (with -ldmumps -lmumps_common -lpord)
#   - scalapack
#   - suitesparse (umfpack, klu)
#   - lapack/blas
#   - hdf5 (optional, for VTKHDF / ParaView 5.11+ output paths)
#   - petsc (optional, for --petsc-ts BP5 driver path + V2 restart)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Default: in-tree build (preserves the seas patches in fem/, mesh/, etc.).
# Pass an alternate path as $1 to fall back to the legacy "build external
# MFEM tree + copy libmfem.a over" workflow (caveat: external = unpatched).
MFEM_SRC="${1:-$SCRIPT_DIR}"

if [ ! -f "$MFEM_SRC/makefile" ] && [ ! -f "$MFEM_SRC/Makefile" ]; then
    echo "ERROR: MFEM source not found at $MFEM_SRC"
    echo "Usage: $0 [/path/to/mfem]   (default: in-tree build in $SCRIPT_DIR)"
    exit 1
fi

# Distinguish in-tree vs external (legacy) mode.  In-tree mode skips the
# libmfem.a copy step (it's already in the right place after `make`); the
# legacy mode copies the build artefacts into SCRIPT_DIR so existing
# downstream consumers keep working.
if [ "$(cd "$MFEM_SRC" && pwd)" = "$SCRIPT_DIR" ]; then
    BUILD_MODE="in-tree"
else
    BUILD_MODE="external"
fi

echo "=== MFEM Setup for SEAS Miniapp ==="
echo "  MFEM source: $MFEM_SRC"
echo "  SEAS-MFEM:   $SCRIPT_DIR"
echo "  Build mode:  $BUILD_MODE"

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

# --- Detect HDF5 ---
# Enables MFEM_USE_HDF5 (VTKHDF output, ParaViewHDFDataCollection).  When
# combined with MFEM_USE_MPI, MFEM auto-detects parallel HDF5 by probing
# `H5_HAVE_PARALLEL` in the linked hdf5 library.  The seas Phase 6
# `DefaultVolumeOutputMode()` switch keys off `MFEM_PARALLEL_HDF5`,
# which MFEM sets automatically when both MPI and a parallel HDF5 are
# present.  No extra config knob is required here.
HDF5_AVAILABLE="NO"
if [ -f "$LIB_PREFIX/include/hdf5.h" ]; then
    HDF5_AVAILABLE="YES"
    HDF5_LIBS="-lhdf5_hl -lhdf5"
    echo "  HDF5: found"
else
    echo "  HDF5: not found (VTKHDF output paths will be VTU-only)"
fi

# --- Detect PETSc ---
# Enables MFEM_USE_PETSC: required by the BP5 `--petsc-ts` driver path
# (Tandem-style PETSc TS time stepper) and the V2 checkpoint restart
# extension (`io/petsc_ts_checkpoint.hpp`, plan
# petsc_ts_restart_plan_2026-05-16.md).  PETSc 3.15+ is required.
PETSC_AVAILABLE="NO"
if [ -f "$LIB_PREFIX/include/petsc.h" ] && [ -f "$LIB_PREFIX/include/petscconf.h" ]; then
    PETSC_AVAILABLE="YES"
    PETSC_LIBS="-lpetsc"
    # The conda-forge petsc bundles its own PETSC_DIR/PETSC_ARCH-aware
    # layout under $LIB_PREFIX, so a single -I/-L pair is sufficient.
    # If you ever switch to a source build with a separate $PETSC_ARCH
    # directory, also append `-I$PETSC_DIR/$PETSC_ARCH/include` and
    # `-L$PETSC_DIR/$PETSC_ARCH/lib` here.
    echo "  PETSc: found (header at $LIB_PREFIX/include/petsc.h)"
else
    echo "  PETSc: not found (--petsc-ts driver path + V2 restart will "
    echo "          be compiled out; build the V1-only paths still work)"
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

if [ "$HDF5_AVAILABLE" = "YES" ]; then
    CONFIG_ARGS+=(
        "MFEM_USE_HDF5=YES"
        "HDF5_OPT=-I$LIB_PREFIX/include"
        "HDF5_LIB=-L$LIB_PREFIX/lib $HDF5_LIBS"
    )
fi

if [ "$PETSC_AVAILABLE" = "YES" ]; then
    CONFIG_ARGS+=(
        "MFEM_USE_PETSC=YES"
        "PETSC_OPT=-I$LIB_PREFIX/include"
        "PETSC_LIB=-L$LIB_PREFIX/lib $PETSC_LIBS"
    )
fi

echo "  Running: make config ${CONFIG_ARGS[*]}"
make config "${CONFIG_ARGS[@]}"

# --- Build MFEM ---
echo ""
echo "=== Building MFEM ==="
NPROC=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
make -j"$NPROC"

# --- Copy artefacts (legacy external-build mode only) ---
if [ "$BUILD_MODE" = "external" ]; then
    echo ""
    echo "=== Copying external build artefacts into $SCRIPT_DIR ==="
    cp "$MFEM_SRC/config/config.mk" "$SCRIPT_DIR/config/config.mk"
    cp "$MFEM_SRC/config/_config.hpp" "$SCRIPT_DIR/config/_config.hpp"
    cp "$MFEM_SRC/libmfem.a" "$SCRIPT_DIR/libmfem.a"
    echo "  Copied config.mk, _config.hpp, libmfem.a"
    echo "  WARNING (external mode): the in-tree seas patches under"
    echo "           fem/, mesh/, linalg/ are NOT in the copied"
    echo "           libmfem.a.  HDF5 compression / VTKHDF / other"
    echo "           seas-specific MFEM APIs will be missing at link time."
    echo "           Re-run without arguments for an in-tree build that"
    echo "           preserves the patches."
else
    echo ""
    echo "=== In-tree build — no copy needed ==="
    echo "  config.mk, _config.hpp, libmfem.a already live in $SCRIPT_DIR"
fi

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
if grep -q "MFEM_USE_HDF5.*YES" "$SCRIPT_DIR/config/config.mk"; then
    echo "  HDF5:         enabled"
    # `MFEM_PARALLEL_HDF5` is preprocessor-defined inside
    # `mesh/vtkhdf.hpp` when `MFEM_USE_MPI` is YES AND the HDF5 header
    # exposes `H5_HAVE_PARALLEL` — NOT a config.mk knob.  Probe the
    # HDF5 header directly so the verification reflects what the
    # compile-time macro will actually see.
    if grep -q "define H5_HAVE_PARALLEL" "$LIB_PREFIX/include/H5pubconf.h" 2>/dev/null \
       && grep -q "MFEM_USE_MPI.*YES" "$SCRIPT_DIR/config/config.mk"; then
        echo "                (parallel HDF5 auto-detected at compile time — VTKHDF default ON for ParMesh)"
    else
        echo "                (serial HDF5 only — VTKHDF default OFF for ParMesh)"
    fi
else
    echo "  HDF5:         disabled (VTKHDF paths compile out)"
fi
if grep -q "MFEM_USE_PETSC.*YES" "$SCRIPT_DIR/config/config.mk"; then
    echo "  PETSc:        enabled"
else
    echo "  PETSc:        disabled (--petsc-ts path + V2 restart compile out)"
fi

echo ""
echo "=== Done ==="
echo "Now rebuild the SEAS miniapp against the freshly-built MFEM:"
echo "  cd miniapps/seas && make clean && make seas_bp5_full \\"
echo "      seas_test_bp5_petsc_ts_restart -j$NPROC"
echo ""
echo "Then verify restart:"
echo "  ./seas_test_bp5_petsc_ts_restart    # V2 restart unit tests"
