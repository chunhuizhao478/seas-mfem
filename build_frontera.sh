#!/bin/bash
# Build MFEM + SEAS miniapp on Frontera (TACC), with PETSc enabled for
# `seas_bp5_full --petsc-ts`.
#
# Usage:
#   bash build_frontera.sh
#
# Optional environment overrides:
#   PETSC_MODULE=petsc/3.23   # Frontera PETSc module to load
#   JOBS=8                    # Parallel build jobs

set -euo pipefail

PETSC_MODULE="${PETSC_MODULE:-petsc/3.23}"
JOBS="${JOBS:-8}"

# Intel classic compiler wrappers on Frontera can abort during config probes when
# the login shell carries an empty or invalid locale. Match the working batch
# scripts in this repo and force the portable C locale for the build.
export LC_ALL=C
export LANG=C

if ! type module >/dev/null 2>&1; then
    echo "ERROR: environment modules are not available in this shell."
    exit 1
fi

# Load modules if not already loaded.
module load intel/19.1.1 2>/dev/null || true
module load impi/19.0.9 2>/dev/null || true
module load hypre/2.31.0 2>/dev/null || true
module load mumps/5.3 2>/dev/null || true
module load parmetis 2>/dev/null || true
module load "${PETSC_MODULE}" 2>/dev/null || true

resolve_petsc_dir() {
    local candidates=(
        "${PETSC_DIR:-}"
        "${TACC_PETSC_DIR:-}"
        "${PETSC_ROOT:-}"
        "${PETSC_HOME:-}"
    )
    local d=""
    for d in "${candidates[@]}"; do
        if [ -n "${d}" ] && [ -f "${d}/lib/petsc/conf/petscvariables" ]; then
            printf '%s\n' "${d}"
            return 0
        fi
    done

    if [ -n "${TACC_PETSC_LIB:-}" ]; then
        d="$(cd "${TACC_PETSC_LIB}/.." 2>/dev/null && pwd || true)"
        if [ -n "${d}" ] && [ -f "${d}/lib/petsc/conf/petscvariables" ]; then
            printf '%s\n' "${d}"
            return 0
        fi
    fi

    d="$(module show "${PETSC_MODULE}" 2>&1 \
        | awk '
            /TACC_PETSC_DIR/ {
                gsub(/"/, "", $0);
                print $NF;
                exit;
            }
        ')"
    if [ -n "${d}" ] && [ -f "${d}/lib/petsc/conf/petscvariables" ]; then
        printf '%s\n' "${d}"
        return 0
    fi

    return 1
}

PETSC_DIR_RESOLVED="$(resolve_petsc_dir || true)"

# Verify modules
echo "=== Checking environment ==="
for var in \
    TACC_HYPRE_INC TACC_HYPRE_LIB \
    TACC_MUMPS_INC TACC_MUMPS_LIB \
    TACC_PARMETIS_INC TACC_PARMETIS_LIB \
    MKLROOT TACC_MKL_LIB; do
    val="$(eval echo \$$var)"
    if [ -z "${val}" ]; then
        echo "ERROR: ${var} is not set. Check module loads."
        exit 1
    fi
    echo "  ${var} = ${val}"
done

if [ -z "${PETSC_DIR_RESOLVED}" ]; then
    echo "ERROR: could not resolve PETSC_DIR after loading ${PETSC_MODULE}."
    echo "Set PETSC_DIR explicitly or choose a different PETSC_MODULE."
    exit 1
fi

echo "  PETSC_MODULE = ${PETSC_MODULE}"
echo "  PETSC_DIR    = ${PETSC_DIR_RESOLVED}"
if [ -n "${TACC_PETSC_LIB:-}" ]; then
    echo "  TACC_PETSC_LIB = ${TACC_PETSC_LIB}"
fi
if [ -n "${TACC_PETSC_INC:-}" ]; then
    echo "  TACC_PETSC_INC = ${TACC_PETSC_INC}"
fi

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
  MFEM_USE_PETSC=YES \
  PETSC_DIR="${PETSC_DIR_RESOLVED}" \
  HYPRE_OPT="-I${TACC_HYPRE_INC}" \
  HYPRE_LIB="-L${TACC_HYPRE_LIB} -Wl,-rpath,${TACC_HYPRE_LIB} -lHYPRE" \
  METIS_OPT="-I${TACC_PARMETIS_INC} -DMETIS_EXPORT=" \
  METIS_LIB="-L${TACC_PARMETIS_LIB} -lparmetis -lmetis" \
  MUMPS_OPT="-I${TACC_MUMPS_INC}" \
  MUMPS_LIB="-L${TACC_MUMPS_LIB} -ldmumps -lmumps_common -lpord -lesmumps -lptscotch -lscotch -lptscotcherr -lscotcherr -L${TACC_MKL_LIB} -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lmpifort -lifcore" \
  LAPACK_OPT="-I${MKLROOT}/include" \
  LAPACK_LIB="-L${TACC_MKL_LIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread"

echo ""
echo "=== Building MFEM library ==="
make -j"${JOBS}"

echo ""
echo "=== Building SEAS miniapps ==="
cd miniapps/seas
make seas_bp1_full seas_bp5_full seas_test_parallel_elasticity seas_test_bp5_parallel_smoke -j"${JOBS}"

echo ""
echo "=== Build complete ==="
echo "Binaries:"
echo "  $(pwd)/seas_bp1_full"
echo "  $(pwd)/seas_bp5_full"
echo "  $(pwd)/seas_test_parallel_elasticity"
echo "  $(pwd)/seas_test_bp5_parallel_smoke"
echo ""
echo "PETSc-enabled run example:"
echo "  ibrun ./seas_bp5_full --petsc-ts"
