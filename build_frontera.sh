#!/bin/bash
# Build MFEM + SEAS miniapp on Frontera (TACC), with PETSc enabled for
# `seas_bp5_full --petsc-ts`.
#
# Usage:
#   bash build_frontera.sh
#
# Optional environment overrides:
#   PETSC_MODULE=petsc/3.23   # Frontera PETSc module to load
#   USE_MUMPS=0               # Disable MFEM's direct MUMPS integration
#   JOBS=8                    # Parallel build jobs

set -euo pipefail

PETSC_MODULE="${PETSC_MODULE:-petsc/3.15}"
USE_MUMPS="${USE_MUMPS:-YES}"
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
USE_MUMPS_RESOLVED="YES"
MUMPS_OPT_RESOLVED="-I${TACC_MUMPS_INC:-}"
PETSC_OPT_RESOLVED=""

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

# Frontera's PETSc modules do not always expose include flags in the exact form
# MFEM's make-based PETSc autodetection expects, especially for older PETSc
# trees such as 3.15. Provide the PETSc include flags explicitly so headers like
# petscversion.h and petscconf.h are always found during MFEM compilation.
if [ -n "${TACC_PETSC_INC:-}" ] && [ -f "${TACC_PETSC_INC}/petscversion.h" ]; then
    PETSC_OPT_RESOLVED="-I${TACC_PETSC_INC}"
elif [ -f "${PETSC_DIR_RESOLVED}/include/petscversion.h" ]; then
    PETSC_OPT_RESOLVED="-I${PETSC_DIR_RESOLVED}/include"
else
    echo "ERROR: could not locate petscversion.h under TACC_PETSC_INC or PETSC_DIR."
    exit 1
fi

# MUMPS is always enabled. Both the petsc and mumps/5.3 modules on Frontera
# are under the same PETSc 3.15 tree, so there is no ABI conflict.
# If MUMPS includes contain PETSc headers (petscconf.h), use -isystem to
# deprioritize them so the PETSc module's own headers win.
if [ "${USE_MUMPS}" = "0" ] || [ "${USE_MUMPS}" = "NO" ] || [ "${USE_MUMPS}" = "no" ]; then
    USE_MUMPS_RESOLVED="NO"
else
    USE_MUMPS_RESOLVED="YES"
fi

if [ "${USE_MUMPS_RESOLVED}" = "YES" ] &&
   [ -n "${TACC_MUMPS_INC:-}" ] && [ -f "${TACC_MUMPS_INC}/petscconf.h" ]; then
    MUMPS_OPT_RESOLVED="-isystem ${TACC_MUMPS_INC}"
    echo "  NOTE: TACC_MUMPS_INC contains PETSc headers; using '-isystem' to avoid header conflicts."
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Initialize toml11 submodule (header-only, v3.8.1 for GCC 8.3 compat)
echo ""
echo "=== Initializing toml11 submodule ==="
git submodule update --init miniapps/seas/extern/toml11 2>/dev/null || true
if [ -d miniapps/seas/extern/toml11 ]; then
    (cd miniapps/seas/extern/toml11 && git checkout v3.8.1 2>/dev/null || true)
    if [ -f miniapps/seas/extern/toml11/toml.hpp ]; then
        echo "  toml11 v3.8.1 ready (header-only, no build needed)"
    else
        echo "  WARNING: toml11 header not found. TOML driver will not build."
    fi
else
    echo "  WARNING: toml11 submodule not available. TOML driver will not build."
fi

echo ""
echo "=== Configuring MFEM ==="
CONFIG_ARGS=(
  MFEM_USE_MPI=YES
  MFEM_USE_METIS=YES
  MFEM_USE_METIS_5=YES
  MFEM_USE_LAPACK=YES
  MFEM_USE_MUMPS="${USE_MUMPS_RESOLVED}"
  MFEM_USE_PETSC=YES
  PETSC_DIR="${PETSC_DIR_RESOLVED}"
  PETSC_OPT="${PETSC_OPT_RESOLVED}"
  HYPRE_OPT="-I${TACC_HYPRE_INC}"
  HYPRE_LIB="-L${TACC_HYPRE_LIB} -Wl,-rpath,${TACC_HYPRE_LIB} -lHYPRE"
  METIS_OPT="-I${TACC_PARMETIS_INC} -DMETIS_EXPORT="
  METIS_LIB="-L${TACC_PARMETIS_LIB} -lparmetis -lmetis"
  LAPACK_OPT="-I${MKLROOT}/include"
  LAPACK_LIB="-L${TACC_MKL_LIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread"
)

if [ "${USE_MUMPS_RESOLVED}" = "YES" ]; then
    CONFIG_ARGS+=(
      MUMPS_OPT="${MUMPS_OPT_RESOLVED}"
      MUMPS_LIB="-L${TACC_MUMPS_LIB} -ldmumps -lmumps_common -lpord -lesmumps -lptscotch -lscotch -lptscotcherr -lscotcherr -L${TACC_MKL_LIB} -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lmpifort -lifcore"
    )
fi

make config "${CONFIG_ARGS[@]}"

echo ""
echo "=== Building MFEM library ==="
make -j"${JOBS}"

echo ""
echo "=== Building SEAS miniapps ==="
cd miniapps/seas
make seas_bp1_full seas_bp5_full seas_test_parallel_elasticity seas_test_bp5_parallel_smoke -j"${JOBS}"

# Build TOML driver if toml11 is available
if [ -f extern/toml11/toml.hpp ] || [ -f extern/toml11/include/toml.hpp ]; then
    echo ""
    echo "=== Building TOML driver ==="
    make seas_driver -j"${JOBS}"
    DRIVER_BUILT=1
else
    DRIVER_BUILT=0
fi

echo ""
echo "=== Build complete ==="
echo "Binaries:"
echo "  $(pwd)/seas_bp1_full"
echo "  $(pwd)/seas_bp5_full"
echo "  $(pwd)/seas_test_parallel_elasticity"
echo "  $(pwd)/seas_test_bp5_parallel_smoke"
if [ "${DRIVER_BUILT}" = "1" ]; then
    echo "  $(pwd)/seas_driver          (TOML-based, new code paths)"
fi
echo ""
echo "Run examples:"
echo "  ibrun ./seas_bp5_full --petsc-ts"
if [ "${DRIVER_BUILT}" = "1" ]; then
    echo "  ibrun ./seas_driver config/bp5_example.toml"
fi
