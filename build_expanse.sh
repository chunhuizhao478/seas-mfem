#!/bin/bash
# Build MFEM + SEAS miniapp on Expanse (SDSC, UCSD).
#
# TEMPORARY first-iteration port of build_frontera.sh — expect to iterate.
#
# ---------------------------------------------------------------------------
# Expanse vs Frontera, and the choices made here
# ---------------------------------------------------------------------------
#   * OS / CPU:  Rocky Linux 8 (glibc 2.28), AMD EPYC 7742 "Rome"
#     (2 x 64 cores / node, 256 GB).  We use GCC, not Intel: SDSC built
#     most of the Expanse stack with gcc/10.2.0 (-march=znver2) and the
#     Intel 19 modules there are NOT the toolchain the math stack was
#     built against (unlike Frontera).
#   * Module system: Lmod over a Spack tree.  Two CPU stacks exist
#     (cpu/0.15.4 and cpu/0.17.3b); we use the newer cpu/0.17.3b:
#         module load cpu/0.17.3b gcc/10.2.0 openmpi/4.1.x cmake/3.21.4
#     (verified working combination for large C++ builds on Expanse).
#   * hypre / metis: Expanse has no TACC-style hypre/parmetis modules on
#     the gcc stack, so BOTH are built from source below:
#         hypre v2.31.0  (same version as the Frontera production stack)
#         metis 5.1.0    (MFEM's required METIS; serial metis suffices)
#     Static libs + explicit -L/-rpath — this also structurally avoids
#     the Frontera "undefined symbol: HYPRE_Initialize" class of
#     LD_LIBRARY_PATH shadowing failures.
#   * MUMPS / PETSc: default ON (the QD solvers need MUMPS), but
#     Expanse has NO mumps/petsc modules on the gcc stack — so PETSc
#     ${PETSC_VERSION} is built from source with
#         --download-mumps --download-scalapack
#         --download-metis --download-parmetis
#     which delivers MUMPS + its whole dependency chain in one
#     validated configure.  This mirrors the Frontera architecture
#     (TACC's mumps/5.3 module lives inside the same PETSc 3.15 tree);
#     the default v3.15.5 pins MUMPS at 5.3.x = Frontera parity.
#     PETSc is built STATIC (--with-shared-libraries=0) so the
#     Frontera "libHYPRE shadowing" class of LD_LIBRARY_PATH failures
#     cannot recur.  When PETSc is on, MFEM's METIS comes from the
#     PETSc prefix too (one metis copy, not two).
#     NOTE: `./configure --download-*` fetches tarballs — Expanse
#     compute nodes have outbound internet, but if a fetch fails run
#     the script once on a login node up through the PETSc configure
#     (or pre-stage with --with-packages-download-dir).
#   * HDF5: built from source (1.14.6), same reason as Frontera — MFEM
#     mesh/vtkhdf.cpp needs HDF5 >= 1.14 (H5S_BLOCK etc.) and Expanse's
#     hdf5 modules are 1.10.x.
#   * gmsh: Expanse's glibc 2.28 >= 2.23, so the prebuilt gmsh.info
#     Linux64 tarball WORKS here (it did not on CentOS-7 Frontera).
#     Default = binary download (~1 min); GMSH_BINARY=NO falls back to
#     the Frontera-style source build.
#   * LAPACK: from the `openblas` module when it resolves; otherwise
#     MFEM_USE_LAPACK=NO with a warning (MFEM has internal fallbacks).
#   * MPI launcher: srun --mpi=pmi2 (or OpenMPI's mpirun).  There is NO
#     ibrun on Expanse — sbatch templates ported from Frontera must
#     replace `ibrun` with `srun --mpi=pmi2 -n $SLURM_NTASKS`.
#   * Filesystems: build in $HOME (100 GB, backed up) or
#     /expanse/lustre/scratch/$USER/temp_project (purged after 90 days).
#     Jobs touching Lustre need `#SBATCH --constraint="lustre"`.
#
# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
#   SDSC asks that compiles run in an interactive session, not on the
#   login nodes (gcc has been OOM-killed on login/GPU nodes: "Killed
#   signal terminated program cc1plus").  Get a shared-partition shell:
#
#     srun --partition=shared --pty --account=<abc123> --nodes=1 \
#          --ntasks-per-node=1 --cpus-per-task=16 --mem=64G \
#          -t 04:00:00 --wait=0 --export=ALL /bin/bash
#
#   then:
#     bash build_expanse.sh
#     QUICK=1 bash build_expanse.sh        # (default) reuse built deps
#     FORCE_REBUILD=1 bash build_expanse.sh
#
# Optional environment overrides:
#   ACCOUNT=abc123              # only used in the printed srun examples
#   CPU_MODULE=cpu/0.17.3b      # Expanse Spack stack to use
#   GCC_MODULE=gcc/10.2.0
#   OPENMPI_MODULES="openmpi/4.1.5 openmpi/4.1.3 openmpi/4.1.1 openmpi"
#                               # tried in order; first that loads wins
#   CMAKE_MODULES="cmake/3.21.4 cmake"
#   OPENBLAS_MODULES="openblas/0.3.18 openblas/0.3.17 openblas"
#   USE_LAPACK=YES              # NO skips the openblas/LAPACK probe
#   USE_MUMPS=YES               # MUMPS direct solver for the QD drivers;
#                               #   delivered by the PETSc source build
#                               #   (requires USE_PETSC=YES)
#   USE_PETSC=YES               # build PETSc from source (also provides
#                               #   MUMPS/ScaLAPACK/METIS/ParMETIS and
#                               #   enables seas_bp5_full --petsc-ts)
#   PETSC_VERSION=v3.15.5       # pinned for Frontera parity (petsc/3.15
#                               #   + mumps 5.3.x); if it fights the
#                               #   gcc-10 toolchain try v3.21.6
#   PETSC_EXTRA_ARGS="..."      # appended verbatim to the PETSc
#                               #   configure (word-split), e.g.
#                               #   "--download-fblaslapack" to swap BLAS
#                               # NOTE: BLAS comes from PETSc's
#                               #   --download-openblas (LP64).  The
#                               #   Expanse openblas MODULE is ILP64 and
#                               #   was rejected by configure ("Cannot
#                               #   use scalapack with 64 bit BLAS/
#                               #   Lapack indices") — OPENBLAS_MODULES
#                               #   is only consulted when USE_PETSC=NO.
#   USE_HDF5=YES                # NO disables HDF5 + H5Z-ZFP integration
#   HDF5_VERSION=1.14.6
#   HYPRE_VERSION=v2.31.0       # hypre tag to build from source
#   METIS_VERSION=5.1.0         # metis tarball version (mfem/tpls mirror)
#   ZFP_VERSION=1.0.1
#   H5Z_ZFP_VERSION=v1.1.1
#   USE_GMSH=YES                # NO skips gmsh entirely
#   GMSH_BINARY=YES             # NO = build gmsh from source instead of
#                               #   downloading the Linux64 binary tarball
#   GMSH_VERSION=4.13.1
#   USE_CALIPER=NO              # default OFF on Expanse iteration 1
#                               #   (perfgraph not needed for first
#                               #   validation; YES builds Caliper from
#                               #   source like build_frontera.sh and is
#                               #   non-fatal on failure)
#   CALIPER_VERSION=v2.11.0
#   MARCH_FLAG=-march=znver2    # AMD Rome tuning; set empty to disable
#   JOBS=8                      # raise to match --cpus-per-task
#   QUICK=1 / FORCE_REBUILD=1   # same caching semantics as Frontera

set -euo pipefail

ACCOUNT="${ACCOUNT:-<abc123>}"
CPU_MODULE="${CPU_MODULE:-cpu/0.17.3b}"
GCC_MODULE="${GCC_MODULE:-gcc/10.2.0}"
OPENMPI_MODULES="${OPENMPI_MODULES:-openmpi/4.1.5 openmpi/4.1.3 openmpi/4.1.1 openmpi}"
CMAKE_MODULES="${CMAKE_MODULES:-cmake/3.21.4 cmake}"
OPENBLAS_MODULES="${OPENBLAS_MODULES:-openblas/0.3.18 openblas/0.3.17 openblas}"
USE_LAPACK="${USE_LAPACK:-YES}"
USE_MUMPS="${USE_MUMPS:-YES}"
USE_PETSC="${USE_PETSC:-YES}"
PETSC_VERSION="${PETSC_VERSION:-v3.15.5}"
PETSC_EXTRA_ARGS="${PETSC_EXTRA_ARGS:-}"
PYTHON_MODULES="${PYTHON_MODULES:-python/3.8.12 python}"
USE_HDF5="${USE_HDF5:-YES}"
HDF5_VERSION="${HDF5_VERSION:-1.14.6}"
HYPRE_VERSION="${HYPRE_VERSION:-v2.31.0}"
METIS_VERSION="${METIS_VERSION:-5.1.0}"
ZFP_VERSION="${ZFP_VERSION:-1.0.1}"
H5Z_ZFP_VERSION="${H5Z_ZFP_VERSION:-v1.1.1}"
USE_GMSH="${USE_GMSH:-YES}"
GMSH_BINARY="${GMSH_BINARY:-YES}"
GMSH_VERSION="${GMSH_VERSION:-4.13.1}"
USE_CALIPER="${USE_CALIPER:-NO}"
CALIPER_DIR="${CALIPER_DIR:-}"
CALIPER_VERSION="${CALIPER_VERSION:-v2.11.0}"
ADIAK_DIR="${ADIAK_DIR:-}"
MARCH_FLAG="${MARCH_FLAG:--march=znver2}"
JOBS="${JOBS:-8}"
QUICK="${QUICK:-1}"
FORCE_REBUILD="${FORCE_REBUILD:-0}"

if [ "${FORCE_REBUILD}" = "1" ] || [ "${FORCE_REBUILD}" = "YES" ] || [ "${FORCE_REBUILD}" = "yes" ]; then
    QUICK=0
fi

is_yes() {
    case "${1:-}" in
        YES|yes|1) return 0 ;;
        *)         return 1 ;;
    esac
}

# MUMPS is delivered by the PETSc source build (--download-mumps); a
# standalone MUMPS would need its own ScaLAPACK/orderings stack.
if is_yes "${USE_MUMPS}" && ! is_yes "${USE_PETSC}"; then
    echo "ERROR: USE_MUMPS=YES requires USE_PETSC=YES on Expanse."
    echo "  MUMPS is built via the PETSc configure (--download-mumps"
    echo "  --download-scalapack); there is no standalone MUMPS path here."
    exit 1
fi

# Match build_frontera.sh: force the portable C locale for config probes.
export LC_ALL=C
export LANG=C

if ! type module >/dev/null 2>&1; then
    # Expanse provides the Lmod init here when the shell missed it.
    if [ -f /etc/profile.d/modules.sh ]; then
        # shellcheck disable=SC1091
        source /etc/profile.d/modules.sh
    fi
fi
if ! type module >/dev/null 2>&1; then
    echo "ERROR: environment modules (Lmod) are not available in this shell."
    exit 1
fi

# ---------------------------------------------------------------------------
# Modules.  Start from a clean slate; cpu/* selects the Spack tree, and
# the spack-generated modules autoload their own dependents.
# ---------------------------------------------------------------------------
echo "=== Loading Expanse modules ==="
module purge 2>/dev/null || true
module load slurm 2>/dev/null || true
module load "${CPU_MODULE}"
module load "${GCC_MODULE}"

# Load the first module of a candidate list that actually resolves.
load_first_of() {
    local label="$1"; shift
    local m
    for m in "$@"; do
        if module load "${m}" 2>/dev/null && module is-loaded "${m}" 2>/dev/null; then
            echo "  ${label}: ${m}"
            LOADED_MODULE="${m}"
            return 0
        fi
    done
    LOADED_MODULE=""
    return 1
}

LOADED_MODULE=""
# shellcheck disable=SC2086
if ! load_first_of "MPI" ${OPENMPI_MODULES}; then
    echo "ERROR: none of '${OPENMPI_MODULES}' loaded under ${CPU_MODULE} + ${GCC_MODULE}."
    echo "  Run 'module spider openmpi' on Expanse and set OPENMPI_MODULES."
    exit 1
fi
MPI_MODULE_RESOLVED="${LOADED_MODULE}"

# shellcheck disable=SC2086
if ! load_first_of "cmake" ${CMAKE_MODULES}; then
    echo "ERROR: no cmake module resolved (tried: ${CMAKE_MODULES})."
    exit 1
fi

# python3 for the PETSc configure (Rocky 8's system /usr/bin/python3
# also works — the module is a nice-to-have, so failure is non-fatal).
if is_yes "${USE_PETSC}"; then
    # shellcheck disable=SC2086
    load_first_of "python" ${PYTHON_MODULES} || true
    if ! command -v python3 >/dev/null 2>&1; then
        echo "ERROR: no python3 available — PETSc configure needs it."
        exit 1
    fi
fi

echo "  cpu stack: ${CPU_MODULE}"
echo "  compiler : ${GCC_MODULE} ($(gcc -dumpversion 2>/dev/null || echo '?'))"
REQUIRED_TOOLS=(mpicc mpicxx cmake)
if is_yes "${USE_MUMPS}"; then
    # MUMPS is Fortran — the PETSc build needs the Fortran MPI wrapper.
    REQUIRED_TOOLS+=(mpif90)
fi
for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "ERROR: ${tool} not on PATH after module loads."
        exit 1
    fi
done
echo "  mpicxx   : $(command -v mpicxx)"
echo "  cmake    : $(cmake --version | head -1)"

# Resolve a module's install prefix: SDSC spack modules export
# <NAME>HOME (e.g. OPENBLASHOME); fall back to parsing `module show`.
prefix_from_module_show() {
    module show "$1" 2>&1 \
        | sed -n 's/.*prepend_path("LD_LIBRARY_PATH", *"\([^"]*\)".*/\1/p' \
        | head -1 | sed -e 's:/lib64$::' -e 's:/lib$::'
}

# ---------------------------------------------------------------------------
# BLAS/LAPACK.
#
# VERIFIED ON THE MACHINE (first build attempt, 2026-06): Expanse's
# openblas module under cpu/0.17.3b + gcc/10.2.0 is built ILP64
# (64-bit integer indices) — PETSc's configure rejected it with
#   "Cannot use scalapack with 64 bit BLAS/Lapack indices".
# An ILP64 BLAS is also silently incompatible with MFEM's LAPACK
# interface (32-bit ints), so it must not be linked ANYWHERE in this
# build.  When PETSc is on (the default) we therefore skip the module
# probe entirely and let PETSc --download-openblas build a standard
# LP64 copy that ScaLAPACK, MUMPS, and MFEM's LAPACK all share.
# The module probe below only runs in a USE_PETSC=NO configuration —
# and may hand MFEM an ILP64 library; verify before trusting results.
# ---------------------------------------------------------------------------
USE_LAPACK_RESOLVED="NO"
LAPACK_LIBDIR_RESOLVED=""
if is_yes "${USE_PETSC}"; then
    echo "  LAPACK   : deferred to PETSc --download-openblas (LP64);"
    echo "             Expanse's openblas module is ILP64 and rejected"
    echo "             by ScaLAPACK/MFEM (see comment in this script)."
elif is_yes "${USE_LAPACK}"; then
    # shellcheck disable=SC2086
    if load_first_of "openblas" ${OPENBLAS_MODULES}; then
        OPENBLAS_MODULE_RESOLVED="${LOADED_MODULE}"
        OPENBLAS_PREFIX="${OPENBLASHOME:-$(prefix_from_module_show "${OPENBLAS_MODULE_RESOLVED}")}"
        for d in "${OPENBLAS_PREFIX}/lib64" "${OPENBLAS_PREFIX}/lib"; do
            if [ -f "${d}/libopenblas.so" ] || [ -f "${d}/libopenblas.a" ]; then
                LAPACK_LIBDIR_RESOLVED="${d}"
                break
            fi
        done
        if [ -n "${LAPACK_LIBDIR_RESOLVED}" ]; then
            USE_LAPACK_RESOLVED="YES"
            echo "  LAPACK   : openblas at ${LAPACK_LIBDIR_RESOLVED}"
        fi
    fi
    if [ "${USE_LAPACK_RESOLVED}" = "NO" ]; then
        echo "  WARNING: no usable openblas module found — building with"
        echo "           MFEM_USE_LAPACK=NO (MFEM internal dense kernels)."
        echo "           Set OPENBLAS_MODULES after 'module spider openblas'."
    fi
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}"

HYPRE_PREFIX="${HYPRE_PREFIX:-${SCRIPT_DIR}/extern/hypre/install}"
METIS_PREFIX="${METIS_PREFIX:-${SCRIPT_DIR}/extern/metis/install}"
PETSC_PREFIX="${PETSC_PREFIX:-${SCRIPT_DIR}/extern/petsc/install}"
HDF5_PREFIX="${HDF5_PREFIX:-${SCRIPT_DIR}/extern/hdf5/install}"
ZFP_PREFIX="${ZFP_PREFIX:-${SCRIPT_DIR}/extern/zfp/install}"
H5Z_ZFP_PREFIX="${H5Z_ZFP_PREFIX:-${SCRIPT_DIR}/extern/h5z-zfp/install}"
GMSH_PREFIX="${GMSH_PREFIX:-${SCRIPT_DIR}/extern/gmsh}"
CALIPER_PREFIX="${CALIPER_PREFIX:-${SCRIPT_DIR}/extern/caliper/install}"

download() {
    local url="$1" dest="$2"
    if command -v wget >/dev/null 2>&1; then
        wget -q -O "${dest}" "${url}"
    elif command -v curl >/dev/null 2>&1; then
        curl -sL -o "${dest}" "${url}"
    else
        echo "ERROR: neither wget nor curl available."
        exit 1
    fi
}

# ---------------------------------------------------------------------------
# hypre ${HYPRE_VERSION} from source (REQUIRED by MFEM_USE_MPI).
# Static lib; hypre's internal blas/lapack (no external dependency).
# ---------------------------------------------------------------------------
build_hypre() {
    local hypre_lib="${HYPRE_PREFIX}/lib/libHYPRE.a"
    if [ "${QUICK}" = "1" ] || [ "${QUICK}" = "YES" ]; then
        if [ -f "${hypre_lib}" ] && [ -f "${HYPRE_PREFIX}/include/HYPRE.h" ]; then
            echo ""
            echo "=== Reusing existing hypre at ${HYPRE_PREFIX} ==="
            echo "    (set FORCE_REBUILD=1 to clear and rebuild)"
            return 0
        fi
        echo ""
        echo "=== Cache miss: ${hypre_lib} not found — building hypre ==="
    fi

    echo ""
    echo "=== Building hypre ${HYPRE_VERSION} from source ==="
    mkdir -p "${SCRIPT_DIR}/extern/hypre"
    local hypre_src="${SCRIPT_DIR}/extern/hypre/src"
    if [ ! -d "${hypre_src}/.git" ]; then
        rm -rf "${hypre_src}"
        git clone --quiet --depth 1 --branch "${HYPRE_VERSION}" \
            https://github.com/hypre-space/hypre.git "${hypre_src}"
    fi
    (
        cd "${hypre_src}"
        git fetch --tags --quiet || true
        git checkout --quiet "${HYPRE_VERSION}" 2>/dev/null \
            || git checkout --quiet "tags/${HYPRE_VERSION}" 2>/dev/null || true
        cd src
        make distclean >/dev/null 2>&1 || true
        CC="$(command -v mpicc)" CXX="$(command -v mpicxx)" \
        ./configure \
            --prefix="${HYPRE_PREFIX}" \
            --disable-fortran \
            --without-superlu
        make -j"${JOBS}"
        make install
    )

    if [ ! -f "${hypre_lib}" ] || [ ! -f "${HYPRE_PREFIX}/include/HYPRE.h" ]; then
        echo "ERROR: hypre build finished but ${hypre_lib} or HYPRE.h is missing."
        exit 1
    fi
    echo "  hypre ${HYPRE_VERSION} installed at ${HYPRE_PREFIX}"
}

# ---------------------------------------------------------------------------
# metis 5.1.0 from source (REQUIRED by MFEM_USE_METIS / parallel meshes).
# Tarball from MFEM's tpls mirror (the original UMN site is unreliable).
# Serial metis suffices — MFEM's ParMesh partitioner calls METIS, not
# ParMETIS.
# ---------------------------------------------------------------------------
build_metis() {
    local metis_lib="${METIS_PREFIX}/lib/libmetis.a"
    if [ "${QUICK}" = "1" ] || [ "${QUICK}" = "YES" ]; then
        if [ -f "${metis_lib}" ] && [ -f "${METIS_PREFIX}/include/metis.h" ]; then
            echo ""
            echo "=== Reusing existing metis at ${METIS_PREFIX} ==="
            echo "    (set FORCE_REBUILD=1 to clear and rebuild)"
            return 0
        fi
        echo ""
        echo "=== Cache miss: ${metis_lib} not found — building metis ==="
    fi

    echo ""
    echo "=== Building metis ${METIS_VERSION} from source ==="
    mkdir -p "${SCRIPT_DIR}/extern/metis"
    local tarball="${SCRIPT_DIR}/extern/metis/metis-${METIS_VERSION}.tar.gz"
    local url="https://github.com/mfem/tpls/raw/gh-pages/metis-${METIS_VERSION}.tar.gz"
    if [ ! -f "${tarball}" ]; then
        echo "  Downloading ${url}"
        download "${url}" "${tarball}"
    fi
    local metis_src="${SCRIPT_DIR}/extern/metis/src"
    rm -rf "${metis_src}"
    mkdir -p "${metis_src}"
    tar -xzf "${tarball}" -C "${metis_src}" --strip-components=1
    (
        cd "${metis_src}"
        make config prefix="${METIS_PREFIX}" cc=gcc
        make -j"${JOBS}"
        make install
    )

    if [ ! -f "${metis_lib}" ] || [ ! -f "${METIS_PREFIX}/include/metis.h" ]; then
        echo "ERROR: metis build finished but ${metis_lib} or metis.h is missing."
        exit 1
    fi
    echo "  metis ${METIS_VERSION} installed at ${METIS_PREFIX}"
}

# ---------------------------------------------------------------------------
# PETSc ${PETSC_VERSION} from source — the MUMPS delivery vehicle.
# --download-mumps --download-scalapack --download-metis
# --download-parmetis builds the whole direct-solver chain with the same
# mpicc/mpicxx/mpif90 toolchain in one configure.  Static
# (--with-shared-libraries=0): every dependent lib is linked into the
# final binaries, so no LD_LIBRARY_PATH ordering can break them.
# BLAS/LAPACK: the openblas module when resolved, else
# --download-openblas (and MFEM's LAPACK then reuses that copy).
# NO --download-hypre: MFEM links our own extern hypre 2.31.0; a second
# hypre inside PETSc is exactly the Frontera shadowing-bug setup.
# ---------------------------------------------------------------------------
build_petsc() {
    if [ "${USE_PETSC_RESOLVED}" != "YES" ]; then
        echo ""
        echo "=== Skipping PETSc build (USE_PETSC=${USE_PETSC}) ==="
        return 0
    fi

    local petsc_conf="${PETSC_PREFIX}/lib/petsc/conf/petscvariables"
    local mumps_lib="${PETSC_PREFIX}/lib/libdmumps.a"
    if [ "${QUICK}" = "1" ] || [ "${QUICK}" = "YES" ]; then
        local cache_ok=1
        [ -f "${petsc_conf}" ] || cache_ok=0
        [ -f "${PETSC_PREFIX}/include/petscversion.h" ] || cache_ok=0
        if is_yes "${USE_MUMPS}"; then
            [ -f "${mumps_lib}" ] || cache_ok=0
        fi
        if [ "${cache_ok}" = "1" ]; then
            echo ""
            echo "=== Reusing existing PETSc at ${PETSC_PREFIX} ==="
            echo "    (set FORCE_REBUILD=1 to clear and rebuild)"
            return 0
        fi
        echo ""
        echo "=== Cache miss: PETSc install incomplete — building ==="
    fi

    echo ""
    echo "=== Building PETSc ${PETSC_VERSION} from source (+MUMPS/ScaLAPACK/ParMETIS) ==="
    mkdir -p "${SCRIPT_DIR}/extern/petsc"
    local petsc_src="${SCRIPT_DIR}/extern/petsc/src"
    if [ ! -d "${petsc_src}/.git" ]; then
        rm -rf "${petsc_src}"
        git clone --quiet --depth 1 --branch "${PETSC_VERSION}" \
            https://github.com/petsc/petsc.git "${petsc_src}"
    fi

    local mumps_args=()
    if is_yes "${USE_MUMPS}"; then
        mumps_args=( "--download-mumps" "--download-scalapack" )
    fi
    # Free-form extra configure args for on-machine iteration, e.g.
    #   PETSC_EXTRA_ARGS="--download-fblaslapack" (swap the BLAS)
    local extra_args=()
    if [ -n "${PETSC_EXTRA_ARGS}" ]; then
        # shellcheck disable=SC2206
        extra_args=( ${PETSC_EXTRA_ARGS} )
    fi

    local arch="arch-expanse-opt"
    (
        cd "${petsc_src}"
        git fetch --tags --quiet || true
        git checkout --quiet "${PETSC_VERSION}" 2>/dev/null \
            || git checkout --quiet "tags/${PETSC_VERSION}" 2>/dev/null || true
        rm -rf "${arch}"
        # The spack openmpi module exports F77/F90/FC/CC/CXX; PETSc's
        # configure warns about them and we pass --with-cc/cxx/fc
        # explicitly anyway — drop them to keep the probe deterministic.
        unset CC CXX FC F77 F90 || true
        # --download-openblas: standard LP64 OpenBLAS shared by
        # ScaLAPACK/MUMPS/MFEM.  Do NOT substitute Expanse's openblas
        # module — it is ILP64 and configure aborts with "Cannot use
        # scalapack with 64 bit BLAS/Lapack indices".
        ./configure PETSC_ARCH="${arch}" \
            --prefix="${PETSC_PREFIX}" \
            --with-debugging=0 \
            --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
            COPTFLAGS="-O2 ${MARCH_FLAG}" \
            CXXOPTFLAGS="-O2 ${MARCH_FLAG}" \
            FOPTFLAGS="-O2 ${MARCH_FLAG}" \
            --with-shared-libraries=0 \
            --with-fortran-bindings=0 \
            --with-x=0 \
            --download-openblas \
            --download-metis \
            --download-parmetis \
            "${mumps_args[@]}" \
            "${extra_args[@]}"
        make PETSC_DIR="${petsc_src}" PETSC_ARCH="${arch}" all -j"${JOBS}"
        make PETSC_DIR="${petsc_src}" PETSC_ARCH="${arch}" install
    )

    if [ ! -f "${petsc_conf}" ] || [ ! -f "${PETSC_PREFIX}/include/petscversion.h" ]; then
        echo "ERROR: PETSc build finished but ${petsc_conf} or petscversion.h is missing."
        exit 1
    fi
    if is_yes "${USE_MUMPS}" && [ ! -f "${mumps_lib}" ]; then
        echo "ERROR: PETSc installed but ${mumps_lib} is missing — check the"
        echo "       --download-mumps step in ${petsc_src}/${arch}/lib/petsc/conf/configure.log"
        exit 1
    fi
    echo "  PETSc ${PETSC_VERSION} installed at ${PETSC_PREFIX}"
}

# ---------------------------------------------------------------------------
# HDF5 ${HDF5_VERSION} from source with --enable-parallel (same rationale
# as build_frontera.sh: MFEM mesh/vtkhdf.cpp needs HDF5 >= 1.14 and the
# Expanse hdf5 modules are 1.10.x).
# ---------------------------------------------------------------------------
build_hdf5() {
    if [ "${USE_HDF5_RESOLVED}" != "YES" ]; then
        echo ""
        echo "=== Skipping HDF5 build (USE_HDF5=${USE_HDF5}) ==="
        return 0
    fi

    local hdf5_lib="${HDF5_PREFIX}/lib/libhdf5.so"
    if [ "${QUICK}" = "1" ] || [ "${QUICK}" = "YES" ]; then
        if [ -f "${hdf5_lib}" ] && [ -f "${HDF5_PREFIX}/include/hdf5.h" ]; then
            echo ""
            echo "=== Reusing existing HDF5 at ${HDF5_PREFIX} ==="
            echo "    (set FORCE_REBUILD=1 to clear and rebuild)"
            return 0
        fi
        echo ""
        echo "=== Cache miss: ${hdf5_lib} not found — building HDF5 ==="
    fi

    echo ""
    echo "=== Building HDF5 ${HDF5_VERSION} ==="
    mkdir -p "${SCRIPT_DIR}/extern"
    local hdf5_src="${SCRIPT_DIR}/extern/hdf5/src"
    if [ ! -d "${hdf5_src}/.git" ]; then
        rm -rf "${hdf5_src}"
        mkdir -p "${SCRIPT_DIR}/extern/hdf5"
        git clone --quiet --depth 1 --branch "hdf5_${HDF5_VERSION}" \
            https://github.com/HDFGroup/hdf5.git "${hdf5_src}"
    fi
    (
        cd "${hdf5_src}"
        git fetch --tags --quiet || true
        git checkout --quiet "hdf5_${HDF5_VERSION}" 2>/dev/null \
            || git checkout --quiet "tags/hdf5_${HDF5_VERSION}" 2>/dev/null || true
        rm -rf build_${HDF5_VERSION}
        mkdir -p build_${HDF5_VERSION}
        cd build_${HDF5_VERSION}
        CC="$(command -v mpicc)" CXX="$(command -v mpicxx)" \
        ../configure \
            --prefix="${HDF5_PREFIX}" \
            --enable-parallel \
            --enable-shared \
            --enable-hl \
            --disable-fortran \
            --disable-cxx \
            --disable-tests \
            --disable-tools \
            --without-szlib
        make -j"${JOBS}"
        make install
    )

    if [ ! -f "${hdf5_lib}" ] || [ ! -f "${HDF5_PREFIX}/include/hdf5.h" ]; then
        echo "ERROR: HDF5 build finished but ${hdf5_lib} or hdf5.h is missing."
        exit 1
    fi
    echo "  HDF5 ${HDF5_VERSION} installed at ${HDF5_PREFIX}"
}

# ---------------------------------------------------------------------------
# ZFP + H5Z-ZFP plugin (lossy floating-point compression for VTKHDF
# output).  Identical recipe to build_frontera.sh; the plugin is loaded
# at runtime via HDF5_PLUGIN_PATH.
# ---------------------------------------------------------------------------
build_zfp_and_h5z_zfp() {
    if [ "${USE_HDF5_RESOLVED}" != "YES" ]; then
        echo ""
        echo "=== Skipping zfp / h5z-zfp build (USE_HDF5=${USE_HDF5}) ==="
        return 0
    fi

    local zfp_lib_so="${ZFP_PREFIX}/lib64/libzfp.so"
    local zfp_lib_so_alt="${ZFP_PREFIX}/lib/libzfp.so"
    local plugin_so="${H5Z_ZFP_PREFIX}/plugin/libh5zzfp.so"
    if [ "${QUICK}" = "1" ] || [ "${QUICK}" = "YES" ]; then
        if { [ -f "${zfp_lib_so}" ] || [ -f "${zfp_lib_so_alt}" ]; } && [ -f "${plugin_so}" ]; then
            echo ""
            echo "=== Reusing existing ZFP + H5Z-ZFP at ${H5Z_ZFP_PREFIX} ==="
            echo "    (set FORCE_REBUILD=1 to clear and rebuild — necessary"
            echo "     when HDF5_VERSION has changed)"
            return 0
        fi
        echo ""
        echo "=== Cache miss: ZFP or H5Z-ZFP install incomplete — building ==="
    fi

    echo ""
    echo "=== Building ZFP ${ZFP_VERSION} ==="
    mkdir -p "${SCRIPT_DIR}/extern"
    if [ ! -d "${SCRIPT_DIR}/extern/zfp/.git" ]; then
        git clone --quiet --depth 1 --branch "${ZFP_VERSION}" \
            https://github.com/LLNL/zfp.git "${SCRIPT_DIR}/extern/zfp"
    fi
    (
        cd "${SCRIPT_DIR}/extern/zfp"
        git fetch --tags --quiet || true
        git checkout --quiet "${ZFP_VERSION}" 2>/dev/null \
            || git checkout --quiet "tags/${ZFP_VERSION}" 2>/dev/null || true
        mkdir -p build
        cd build
        # ZFP_BIT_STREAM_WORD_SIZE=8 is REQUIRED by the H5Z-ZFP plugin's
        # can_apply callback (same as the Frontera build — see
        # build_frontera.sh for the full explanation).
        cmake .. \
            -DCMAKE_INSTALL_PREFIX="${ZFP_PREFIX}" \
            -DCMAKE_C_COMPILER="$(command -v mpicc)" \
            -DBUILD_SHARED_LIBS=ON \
            -DBUILD_TESTING=OFF \
            -DBUILD_EXAMPLES=OFF \
            -DZFP_WITH_OPENMP=OFF \
            -DZFP_BIT_STREAM_WORD_SIZE=8
        cmake --build . --target install -j "${JOBS}"
    )

    echo ""
    echo "=== Building H5Z-ZFP ${H5Z_ZFP_VERSION} ==="
    if [ ! -d "${SCRIPT_DIR}/extern/H5Z-ZFP/.git" ]; then
        git clone --quiet --depth 1 --branch "${H5Z_ZFP_VERSION}" \
            https://github.com/LLNL/H5Z-ZFP.git "${SCRIPT_DIR}/extern/H5Z-ZFP"
    fi
    (
        cd "${SCRIPT_DIR}/extern/H5Z-ZFP"
        git fetch --tags --quiet || true
        git checkout --quiet "${H5Z_ZFP_VERSION}" 2>/dev/null \
            || git checkout --quiet "tags/${H5Z_ZFP_VERSION}" 2>/dev/null || true
        # FC= (empty) DISABLES H5Z-ZFP's optional Fortran interface.
        # H5Z-ZFP/src/Makefile builds H5Zzfp_props_f.o whenever $(FC) is
        # non-empty ("ifneq ($(FC),)").  Our HDF5 is built --disable-fortran
        # (no hdf5.mod), and the Expanse spack openmpi module EXPORTS FC into
        # the environment — so without this the Fortran object compiles and
        # dies on "USE HDF5: Cannot open module file 'hdf5.mod'".  It must be
        # a command-line assignment, NOT `unset FC`: GNU Make's built-in
        # default is FC=f77, so unsetting the env var would let the gate fire
        # again with f77.  MFEM consumes only the C plugin (libh5zzfp.so) via
        # HDF5_PLUGIN_PATH, so the Fortran interface is not needed here.
        make FC= \
             HDF5_HOME="${HDF5_PREFIX}" \
             ZFP_HOME="${ZFP_PREFIX}" \
             PREFIX="${H5Z_ZFP_PREFIX}" \
             clean || true
        # SERIAL build (no -j): H5Z-ZFP's Makefile updates libh5zzfp.a via
        # two separate archive-member rules
        #   libh5zzfp.a(H5Zzfp_lib.o):   ... ; $(AR) cr libh5zzfp.a $<
        #   libh5zzfp.a(H5Zzfp_props.o): ... ; $(AR) cr libh5zzfp.a $<
        # with NO .NOTPARALLEL.  Under -j both `ar cr` run concurrently on
        # the same archive and race -> "ar: libh5zzfp.a: File format not
        # recognized" + "unable to copy file ...: No such file or directory".
        # Only 3 source files, so serial is instant.
        make CC="$(command -v mpicc)" FC= \
             HDF5_HOME="${HDF5_PREFIX}" \
             ZFP_HOME="${ZFP_PREFIX}" \
             PREFIX="${H5Z_ZFP_PREFIX}" \
             install
    )

    if [ ! -f "${plugin_so}" ]; then
        echo "ERROR: H5Z-ZFP build finished but ${plugin_so} is missing."
        exit 1
    fi
    echo "  H5Z-ZFP plugin installed at ${H5Z_ZFP_PREFIX}/plugin"
}

# ---------------------------------------------------------------------------
# gmsh: binary tarball by default (Rocky 8 glibc 2.28 satisfies gmsh's
# 2.23 requirement — the reason Frontera had to source-build does not
# apply here).  GMSH_BINARY=NO uses the Frontera-style source build.
# ---------------------------------------------------------------------------
build_gmsh() {
    if ! is_yes "${USE_GMSH}"; then
        echo ""
        echo "=== Skipping gmsh (USE_GMSH=${USE_GMSH}) ==="
        return 0
    fi

    local gmsh_bin="${GMSH_PREFIX}/bin/gmsh"
    if [ "${QUICK}" = "1" ] || [ "${QUICK}" = "YES" ]; then
        if [ -x "${gmsh_bin}" ] && "${gmsh_bin}" -version >/dev/null 2>&1; then
            echo ""
            echo "=== Reusing existing gmsh at ${gmsh_bin} ==="
            echo "    (set FORCE_REBUILD=1 to rebuild)"
            return 0
        fi
        echo ""
        echo "=== Cache miss or stale install: fetching gmsh ==="
    fi

    mkdir -p "${SCRIPT_DIR}/extern"
    if is_yes "${GMSH_BINARY}"; then
        echo ""
        echo "=== Installing gmsh ${GMSH_VERSION} (prebuilt Linux64 binary) ==="
        local tarball="${SCRIPT_DIR}/extern/gmsh-${GMSH_VERSION}-Linux64.tgz"
        local url="https://gmsh.info/bin/Linux/gmsh-${GMSH_VERSION}-Linux64.tgz"
        if [ ! -f "${tarball}" ]; then
            echo "  Downloading ${url}"
            download "${url}" "${tarball}"
        fi
        rm -rf "${GMSH_PREFIX}"
        mkdir -p "${GMSH_PREFIX}"
        tar -xzf "${tarball}" -C "${GMSH_PREFIX}" --strip-components=1
        if [ ! -x "${gmsh_bin}" ] || ! "${gmsh_bin}" -version >/dev/null 2>&1; then
            echo "ERROR: prebuilt gmsh failed to run (check 'ldd ${gmsh_bin}')."
            echo "  Re-run with GMSH_BINARY=NO to build gmsh from source."
            exit 1
        fi
        echo "  gmsh $("${gmsh_bin}" -version 2>&1) installed at ${gmsh_bin}"
        return 0
    fi

    echo ""
    echo "=== Building gmsh ${GMSH_VERSION} from source ==="
    local src_tarball="${SCRIPT_DIR}/extern/gmsh-${GMSH_VERSION}-source.tgz"
    local src_url="https://gmsh.info/src/gmsh-${GMSH_VERSION}-source.tgz"
    if [ ! -f "${src_tarball}" ]; then
        echo "  Downloading ${src_url}"
        download "${src_url}" "${src_tarball}"
    fi
    local gmsh_src="${SCRIPT_DIR}/extern/gmsh-src"
    rm -rf "${gmsh_src}" "${GMSH_PREFIX}"
    mkdir -p "${gmsh_src}"
    tar -xzf "${src_tarball}" -C "${gmsh_src}" --strip-components=1
    (
        cd "${gmsh_src}"
        mkdir -p build_${GMSH_VERSION}
        cd build_${GMSH_VERSION}
        # CMAKE_POLICY_VERSION_MINIMUM=3.5: same modern-CMake workaround
        # as the Frontera build (gmsh still declares 2.8.12).
        cmake .. \
            -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
            -DCMAKE_INSTALL_PREFIX="${GMSH_PREFIX}" \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_C_COMPILER=gcc \
            -DCMAKE_CXX_COMPILER=g++ \
            -DENABLE_FLTK=OFF \
            -DENABLE_OCC=OFF \
            -DENABLE_MED=OFF \
            -DENABLE_PETSC=OFF \
            -DENABLE_SLEPC=OFF \
            -DENABLE_MPI=OFF \
            -DENABLE_WRAP_PYTHON=OFF \
            -DENABLE_BUILD_LIB=OFF \
            -DENABLE_BUILD_DYNAMIC=OFF
        cmake --build . --target install -j "${JOBS}"
    )
    if [ ! -x "${gmsh_bin}" ]; then
        echo "ERROR: gmsh build finished but ${gmsh_bin} missing."
        exit 1
    fi
    echo "  gmsh ${GMSH_VERSION} built at ${gmsh_bin}"
}

# ---------------------------------------------------------------------------
# Caliper (optional, default OFF on Expanse iteration 1).  Same
# from-source recipe as build_frontera.sh; MUST use the same
# gcc + openmpi toolchain as MFEM (libcaliper is linked in).
# Non-fatal: failure WARNs and the build continues without Caliper.
# ---------------------------------------------------------------------------
build_caliper() {
    if ! is_yes "${USE_CALIPER}"; then
        echo ""
        echo "=== Skipping Caliper build (USE_CALIPER=${USE_CALIPER}) ==="
        return 0
    fi
    if [ -n "${CALIPER_DIR}" ] && [ -f "${CALIPER_DIR}/include/caliper/cali.h" ]; then
        echo ""
        echo "=== Using provided Caliper at ${CALIPER_DIR} (not building from source) ==="
        return 0
    fi

    local cali_so="${CALIPER_PREFIX}/lib64/libcaliper.so"
    local cali_so_alt="${CALIPER_PREFIX}/lib/libcaliper.so"
    local cali_hdr="${CALIPER_PREFIX}/include/caliper/cali.h"
    if [ "${QUICK}" = "1" ] || [ "${QUICK}" = "YES" ]; then
        if { [ -f "${cali_so}" ] || [ -f "${cali_so_alt}" ]; } && [ -f "${cali_hdr}" ]; then
            echo ""
            echo "=== Reusing existing Caliper at ${CALIPER_PREFIX} ==="
            CALIPER_DIR="${CALIPER_PREFIX}"
            return 0
        fi
        echo ""
        echo "=== Cache miss: Caliper install incomplete — building ==="
    fi

    echo ""
    echo "=== Building Caliper ${CALIPER_VERSION} from source (gcc + openmpi) ==="
    mkdir -p "${SCRIPT_DIR}/extern/caliper"
    local cali_src="${SCRIPT_DIR}/extern/caliper/src"
    if [ ! -d "${cali_src}/.git" ]; then
        rm -rf "${cali_src}"
        if ! git clone --quiet --depth 1 --branch "${CALIPER_VERSION}" \
                https://github.com/LLNL/Caliper.git "${cali_src}"; then
            echo "WARNING: Caliper clone (${CALIPER_VERSION}) failed — continuing WITHOUT Caliper."
            return 1
        fi
    fi
    if ! (
        cd "${cali_src}"
        git fetch --tags --quiet || true
        git checkout --quiet "${CALIPER_VERSION}" 2>/dev/null \
            || git checkout --quiet "tags/${CALIPER_VERSION}" 2>/dev/null || true
        rm -rf "build_${CALIPER_VERSION}"
        mkdir -p "build_${CALIPER_VERSION}"
        cd "build_${CALIPER_VERSION}"
        cmake .. \
            -DCMAKE_INSTALL_PREFIX="${CALIPER_PREFIX}" \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_C_COMPILER="$(command -v mpicc)" \
            -DCMAKE_CXX_COMPILER="$(command -v mpicxx)" \
            -DBUILD_SHARED_LIBS=ON \
            -DWITH_MPI=ON \
            -DWITH_TESTS=OFF \
            -DWITH_FORTRAN=OFF \
            -DWITH_DOCS=OFF \
            -DBUILD_TESTING=OFF
        make -j"${JOBS}"
        make install
    ); then
        echo "WARNING: Caliper build failed — continuing WITHOUT Caliper."
        return 1
    fi

    if [ ! -f "${cali_hdr}" ] || { [ ! -f "${cali_so}" ] && [ ! -f "${cali_so_alt}" ]; }; then
        echo "WARNING: Caliper build finished but cali.h / libcaliper missing — continuing WITHOUT Caliper."
        return 1
    fi
    CALIPER_DIR="${CALIPER_PREFIX}"
    echo "  Caliper ${CALIPER_VERSION} installed at ${CALIPER_PREFIX}"
}

# ---------------------------------------------------------------------------
# Run the dependency builds.
# ---------------------------------------------------------------------------
if is_yes "${USE_HDF5}"; then
    USE_HDF5_RESOLVED="YES"
else
    USE_HDF5_RESOLVED="NO"
fi
if is_yes "${USE_PETSC}"; then
    USE_PETSC_RESOLVED="YES"
else
    USE_PETSC_RESOLVED="NO"
fi
if is_yes "${USE_MUMPS}"; then
    USE_MUMPS_RESOLVED="YES"
else
    USE_MUMPS_RESOLVED="NO"
fi

build_hypre
build_petsc

# Reuse PETSc's --download-openblas copy (LP64) for MFEM's LAPACK and
# the MUMPS link line.  Runs here (not inside build_petsc) so QUICK
# cache-hit re-runs resolve it too.
if [ "${USE_PETSC_RESOLVED}" = "YES" ] && [ -z "${LAPACK_LIBDIR_RESOLVED}" ] && \
   { [ -f "${PETSC_PREFIX}/lib/libopenblas.a" ] || [ -f "${PETSC_PREFIX}/lib/libopenblas.so" ]; }; then
    LAPACK_LIBDIR_RESOLVED="${PETSC_PREFIX}/lib"
    USE_LAPACK_RESOLVED="YES"
    echo "  LAPACK   : PETSc-downloaded openblas at ${LAPACK_LIBDIR_RESOLVED}"
fi
if [ "${USE_MUMPS_RESOLVED}" = "YES" ] && [ -z "${LAPACK_LIBDIR_RESOLVED}" ]; then
    echo "ERROR: MUMPS is enabled but no BLAS/LAPACK resolved (no openblas"
    echo "       module and no PETSc-downloaded openblas) — the MUMPS link"
    echo "       line would be incomplete.  Check the PETSc build."
    exit 1
fi

# METIS for MFEM: when PETSc is on, its --download-metis/--download-parmetis
# install into ${PETSC_PREFIX} — use that single copy (linking a second
# extern metis alongside MUMPS' would risk mixed-version symbols).
if [ "${USE_PETSC_RESOLVED}" = "YES" ]; then
    METIS_INC_RESOLVED="${PETSC_PREFIX}/include"
    METIS_LIBDIR_RESOLVED="${PETSC_PREFIX}/lib"
    METIS_LIB_RESOLVED="-L${METIS_LIBDIR_RESOLVED} -lparmetis -lmetis"
    if [ ! -f "${METIS_INC_RESOLVED}/metis.h" ]; then
        echo "ERROR: ${METIS_INC_RESOLVED}/metis.h missing — PETSc was expected"
        echo "       to install METIS (--download-metis).  Check configure.log."
        exit 1
    fi
    echo "  METIS    : from PETSc prefix (${METIS_LIBDIR_RESOLVED})"
else
    build_metis
    METIS_INC_RESOLVED="${METIS_PREFIX}/include"
    METIS_LIBDIR_RESOLVED="${METIS_PREFIX}/lib"
    METIS_LIB_RESOLVED="-L${METIS_LIBDIR_RESOLVED} -lmetis"
fi

build_hdf5
build_zfp_and_h5z_zfp
build_gmsh

if [ "${USE_HDF5_RESOLVED}" = "YES" ]; then
    if grep -q '^#define *H5_HAVE_PARALLEL' "${HDF5_PREFIX}/include/H5pubconf.h" 2>/dev/null; then
        echo "  HDF5         = parallel (H5_HAVE_PARALLEL defined)"
    else
        echo "  WARNING: HDF5 build does NOT define H5_HAVE_PARALLEL — Phase 6"
        echo "           ParMesh runs will fail at the static_assert in"
        echo "           paraview_output.hpp.  Rebuild HDF5 with --enable-parallel."
    fi
fi

# Initialize toml11 submodule (header-only, v3.8.1)
echo ""
echo "=== Initializing toml11 submodule ==="
git submodule update --init miniapps/seas/extern/toml11 2>/dev/null || true
if [ -d miniapps/seas/extern/toml11 ]; then
    (cd miniapps/seas/extern/toml11 && git checkout v3.8.1 2>/dev/null || true)
    if [ -f miniapps/seas/extern/toml11/toml.hpp ]; then
        echo "  toml11 v3.8.1 ready (header-only, no build needed)"
    else
        echo "  WARNING: toml11 header not found. TOML drivers will not build."
    fi
else
    echo "  WARNING: toml11 submodule not available. TOML drivers will not build."
fi

# Caliper resolution: explicit env > spack > built-from-source.
USE_CALIPER_RESOLVED="NO"
if is_yes "${USE_CALIPER}"; then
    if [ -z "${CALIPER_DIR}" ] && command -v spack >/dev/null 2>&1; then
        CALIPER_DIR="$(spack location -i caliper 2>/dev/null || true)"
    fi
    build_caliper || true
    if [ -n "${CALIPER_DIR}" ] && [ -f "${CALIPER_DIR}/include/caliper/cali.h" ]; then
        USE_CALIPER_RESOLVED="YES"
        echo "  CALIPER_DIR  = ${CALIPER_DIR}"
        if [ -n "${ADIAK_DIR}" ]; then echo "  ADIAK_DIR    = ${ADIAK_DIR}"; fi
    else
        echo "  WARNING: USE_CALIPER=YES but no usable Caliper — continuing"
        echo "           WITHOUT Caliper (perfgraph annotations stay no-ops)."
    fi
fi

# ---------------------------------------------------------------------------
# Configure + build MFEM
# ---------------------------------------------------------------------------
echo ""
echo "=== Configuring MFEM ==="
CONFIG_ARGS=(
  MFEM_USE_MPI=YES
  MFEM_USE_METIS=YES
  MFEM_USE_METIS_5=YES
  MFEM_USE_LAPACK="${USE_LAPACK_RESOLVED}"
  MFEM_USE_MUMPS="${USE_MUMPS_RESOLVED}"
  MFEM_USE_PETSC="${USE_PETSC_RESOLVED}"
  HYPRE_OPT="-I${HYPRE_PREFIX}/include"
  HYPRE_LIB="-L${HYPRE_PREFIX}/lib -lHYPRE"
  METIS_OPT="-I${METIS_INC_RESOLVED}"
  METIS_LIB="${METIS_LIB_RESOLVED}"
  OPTIM_FLAGS="-O3 ${MARCH_FLAG} -std=c++17"
)

if [ "${USE_LAPACK_RESOLVED}" = "YES" ]; then
    CONFIG_ARGS+=(
      LAPACK_OPT=""
      LAPACK_LIB="-L${LAPACK_LIBDIR_RESOLVED} -Wl,-rpath,${LAPACK_LIBDIR_RESOLVED} -lopenblas"
    )
fi

if [ "${USE_PETSC_RESOLVED}" = "YES" ]; then
    # Prefix install: MFEM reads ${PETSC_DIR}/lib/petsc/conf/petscvariables
    # directly (no PETSC_ARCH), same as the Frontera TACC module layout.
    CONFIG_ARGS+=(
      PETSC_DIR="${PETSC_PREFIX}"
      PETSC_OPT="-I${PETSC_PREFIX}/include"
    )
fi

if [ "${USE_MUMPS_RESOLVED}" = "YES" ]; then
    # Mirrors MFEM's own MUMPS_LIB recipe (config/defaults.mk:357-364):
    # dmumps + mumps_common + pord + ScaLAPACK + LAPACK + MPI-Fortran.
    # parmetis/metis are MUMPS' ordering libs (PETSc built it with them);
    # -lmpi_mpifh is OpenMPI 4.x's Fortran MPI library and -lgfortran the
    # GNU Fortran runtime (MUMPS is Fortran, the final link is mpicxx).
    #
    # -Wl,--allow-multiple-definition: netlib ScaLAPACK (PETSc
    # --download-scalapack) BUNDLES a private copy of the LAPACK 3.7+ aux
    # routine dcombssq (in pdtreecomb.o), which collides with the same symbol
    # in OpenBLAS's LAPACK (dcombssq.o) at static-link time:
    #   "multiple definition of `dcombssq_'".  Both are the identical
    # reference-LAPACK routine, so taking the first is harmless; this flag
    # (a global ld mode toggle) tells the linker to keep it and drop the
    # duplicate.  Only needed when MUMPS is on (the sole config that links
    # BOTH -lscalapack and -lopenblas).  Placed at the FRONT of MUMPS_LIB so
    # it is present on every MFEM miniapp/test link line.
    MUMPS_LAPACK_PART=""
    if [ -n "${LAPACK_LIBDIR_RESOLVED}" ]; then
        MUMPS_LAPACK_PART="-L${LAPACK_LIBDIR_RESOLVED} -Wl,-rpath,${LAPACK_LIBDIR_RESOLVED} -lopenblas"
    fi
    CONFIG_ARGS+=(
      MUMPS_OPT="-I${PETSC_PREFIX}/include"
      MUMPS_LIB="-Wl,--allow-multiple-definition -L${PETSC_PREFIX}/lib -ldmumps -lmumps_common -lpord -lscalapack -lparmetis -lmetis ${MUMPS_LAPACK_PART} -lgfortran -lmpi_mpifh -lpthread"
    )
fi

if [ "${USE_HDF5_RESOLVED}" = "YES" ]; then
    CONFIG_ARGS+=(
      MFEM_USE_HDF5=YES
      MFEM_USE_H5Z_ZFP=YES
      HDF5_OPT="-I${HDF5_PREFIX}/include"
      HDF5_LIB="-L${HDF5_PREFIX}/lib -Wl,-rpath,${HDF5_PREFIX}/lib -lhdf5_hl -lhdf5 -lz"
    )
fi

if [ "${USE_CALIPER_RESOLVED}" = "YES" ]; then
    CONFIG_ARGS+=( MFEM_USE_CALIPER=YES CALIPER_DIR="${CALIPER_DIR}" )
    if [ -n "${ADIAK_DIR}" ]; then
        CONFIG_ARGS+=( ADIAK_DIR="${ADIAK_DIR}" )
    fi
fi

make config "${CONFIG_ARGS[@]}"

echo ""
echo "=== Building MFEM library ==="
make -j"${JOBS}"

# ---------------------------------------------------------------------------
# Build the SEAS miniapps (same target filtering as build_frontera.sh:
# skip the targets that instantiate SEASQuasiDynamicOperator with
# AntiplaneDomainOperator — see miniapps/seas/debug_document/
# bp5_debug_document/antiplane_template_instantiation_2026-05-11.md).
# ---------------------------------------------------------------------------
echo ""
echo "=== Building SEAS miniapps ==="
cd miniapps/seas
if [ "${USE_CALIPER_RESOLVED}" = "YES" ]; then
    echo "  Caliper build: 'make clean' in miniapps/seas (no-header-deps stale-object trap)"
    make clean
fi
EXPLICIT_SKIPS=(
    seas_test_br2_consistency
    seas_test_parallel_fault
    seas_test_serial_parallel_consistency
    seas_test_scaling
    seas_test_checkpoint
    seas_test_io
    seas_test_quasi_dynamic
)
MINIAPPS_LIST="$(make -s -f Makefile -f - print-miniapps <<'PRINT_MAKEFILE'
.PHONY: print-miniapps
print-miniapps:
	@echo $(MINIAPPS)
PRINT_MAKEFILE
)"
is_explicit_skip() {
    local target="$1"
    for s in "${EXPLICIT_SKIPS[@]}"; do
        [ "${s}" = "${target}" ] && return 0
    done
    return 1
}
FILTERED_TARGETS=""
SKIPPED_TARGETS=""
for t in ${MINIAPPS_LIST}; do
    case "${t}" in
        *bp1*|*bp2*|*antiplane*|*pseas*)
            SKIPPED_TARGETS="${SKIPPED_TARGETS} ${t}"
            continue
            ;;
    esac
    if is_explicit_skip "${t}"; then
        SKIPPED_TARGETS="${SKIPPED_TARGETS} ${t}"
    else
        FILTERED_TARGETS="${FILTERED_TARGETS} ${t}"
    fi
done
if [ -n "${SKIPPED_TARGETS}" ]; then
    echo "  Skipping antiplane-template-instantiating targets:${SKIPPED_TARGETS}"
fi
# seas_spatial_dyn_driver is NOT in MINIAPPS (own target, like the tpv
# drivers) and is the primary production driver — build it explicitly.
make ${FILTERED_TARGETS} \
     seas_tpv102_driver seas_tpv104_driver seas_tpv205_driver \
     seas_spatial_dyn_driver \
     -j"${JOBS}"

if [ -f extern/toml11/toml.hpp ] || [ -f extern/toml11/include/toml.hpp ]; then
    echo ""
    echo "=== Building TOML driver ==="
    make seas_driver -j"${JOBS}"
    DRIVER_BUILT=1
else
    DRIVER_BUILT=0
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "=== Build complete (Expanse) ==="
echo "Modules used: ${CPU_MODULE} ${GCC_MODULE} ${MPI_MODULE_RESOLVED}"
echo ""
echo "Key binaries (full list under $(pwd)):"
echo "  ./seas_spatial_dyn_driver                  # SAFS dynamic-rupture driver"
echo "  ./seas_tpv102_driver"
echo "  ./seas_tpv104_driver"
echo "  ./seas_tpv205_driver"
if [ "${USE_MUMPS_RESOLVED}" = "YES" ]; then
    echo "  ./seas_bp5_full                            # MUMPS direct solver enabled"
    if [ "${USE_PETSC_RESOLVED}" = "YES" ]; then
        echo "                                             #   (--petsc-ts available)"
    fi
else
    echo "  ./seas_bp5_full                            # NO MUMPS: QD solve falls"
    echo "                                             #   back to CG + HypreILU"
fi
echo "  ./seas_test_parallel_elasticity            # (and ~50 other unit tests)"
if [ "${DRIVER_BUILT}" = "1" ]; then
    echo "  $(pwd)/seas_driver          (TOML-based)"
fi
if [ "${USE_HDF5_RESOLVED}" = "YES" ]; then
    echo ""
    echo "H5Z-ZFP plugin (set at runtime in your sbatch):"
    echo "  export HDF5_PLUGIN_PATH=${H5Z_ZFP_PREFIX}/plugin"
fi
if is_yes "${USE_GMSH}"; then
    echo ""
    echo "gmsh: ${GMSH_PREFIX}/bin/gmsh  (emit v2.2 meshes: gmsh ... -format msh22)"
fi
echo ""
echo "Running on Expanse (NO ibrun here — use srun or mpirun):"
echo "  srun --mpi=pmi2 -n \$SLURM_NTASKS ./seas_spatial_dyn_driver --config <cfg.toml>"
echo ""
echo "sbatch header sketch:"
echo "  #SBATCH --account=${ACCOUNT}"
echo "  #SBATCH --partition=compute          # or shared / debug"
echo "  #SBATCH --nodes=N --ntasks-per-node=128 --mem=0"
echo "  #SBATCH --constraint=\"lustre\"        # REQUIRED if the job touches"
echo "                                       # /expanse/lustre/* paths"
echo "  module purge && module load slurm ${CPU_MODULE} ${GCC_MODULE} ${MPI_MODULE_RESOLVED}"
echo ""
echo "Output dirs belong on /expanse/lustre/scratch/\$USER/temp_project"
echo "(purged 90 days after creation; \$HOME is 100 GB and slow)."
