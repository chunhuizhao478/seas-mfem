#!/bin/bash
# Build MFEM + SEAS miniapp on Frontera (TACC), with PETSc enabled for
# `seas_bp5_full --petsc-ts` and HDF5 + H5Z-ZFP enabled for the Phase 6
# ParaView output (VTKHDF + lossy floating-point compression).
#
# Usage:
#   bash build_frontera.sh
#   QUICK=1 bash build_frontera.sh          # skip h5z-zfp / zfp rebuild
#
# Optional environment overrides:
#   PETSC_MODULE=petsc/3.23     # Frontera PETSc module to load
#   HDF5_USE_MODULE=NO          # By default, build HDF5 1.14.x from
#                                 source under extern/hdf5/install rather
#                                 than load Frontera's phdf5 module.
#                                 Reason: MFEM mesh/vtkhdf.cpp requires
#                                 HDF5 >= 1.14 (uses H5S_BLOCK and assumes
#                                 hsize_t == unsigned long long).  The
#                                 phdf5/1.14.x module on Frontera requires
#                                 intel/23.1.0 + impi/21.9.0, incompatible
#                                 with the intel/19.1.1 + impi/19.0.9
#                                 hypre/mumps/petsc stack we use.  Older
#                                 modules (1.10.x / 1.12.x) DO load with
#                                 intel/19 but lack the symbols MFEM uses.
#   HDF5_VERSION=1.14.6         # HDF5 release tag to build from source
#                                 (matches local conda-env's 1.14.3 API).
#   PHDF5_MODULE=phdf5/1.12.2   # Only used when HDF5_USE_MODULE=YES (the
#                                 fallback path that loads a TACC phdf5
#                                 module instead of building one).
#   USE_MUMPS=0                 # Disable MFEM's direct MUMPS integration
#   USE_HDF5=NO                 # Disable HDF5 + H5Z-ZFP integration
#   QUICK=1   (default)         # Reuse already-installed deps when their
#                               # key files are present (hdf5: libhdf5.so
#                               # + hdf5.h; zfp/h5z-zfp: libh5zzfp.{so,
#                               # dylib}).  Run with FORCE_REBUILD=1 (or
#                               # QUICK=0) to clear the cache and build
#                               # everything from scratch.
#   FORCE_REBUILD=1             # Force fresh rebuild of HDF5, ZFP, and
#                               # H5Z-ZFP even if existing installs look
#                               # complete.  Useful after changing
#                               # HDF5_VERSION or after switching from
#                               # HDF5_USE_MODULE=YES to NO (the cached
#                               # H5Z-ZFP plugin is linked against the
#                               # old HDF5).
#   ZFP_VERSION=1.0.1           # ZFP base library tag to check out
#   H5Z_ZFP_VERSION=v1.1.1      # H5Z-ZFP plugin tag to check out
#   ZFP_PREFIX=...              # Override the ZFP install prefix
#                               # (default: $SCRIPT_DIR/extern/zfp/install)
#   H5Z_ZFP_PREFIX=...          # Override the H5Z-ZFP install prefix
#                               # (default: $SCRIPT_DIR/extern/h5z-zfp/install).
#                               # NOTE: the Phase 6 sbatch defaults look
#                               # for the plugin at $WORK/h5z-zfp/install/
#                               # plugin — if you change this prefix you
#                               # also need to update HDF5_PLUGIN_PATH in
#                               # miniapps/seas/jobs/{bp5,tpv*}/*phase6*.sbatch
#                               # (or set HDF5_PLUGIN_PATH at runtime).
#   USE_GMSH=YES                # Build gmsh from source into extern/gmsh.
#                               # Frontera doesn't ship a gmsh module and
#                               # the 200 m TPV .msh files (~125 MB each)
#                               # are too large to commit, so Phase 6
#                               # sbatch regenerate .msh from tracked .geo
#                               # at submit time using this gmsh.  Source
#                               # build (rather than the prebuilt Linux64
#                               # tarball) is required because gmsh's
#                               # binary releases need glibc 2.23+ and
#                               # Frontera ships CentOS 7 / glibc 2.17.
#                               # Build is ~5 min, install ~50 MB (no GUI,
#                               # no OpenCASCADE, no MPI, no Python).
#                               # USE_GMSH=NO skips the build.
#   GMSH_VERSION=4.13.1         # gmsh release tag to build from source.
#   JOBS=8                      # Parallel build jobs

set -euo pipefail

PETSC_MODULE="${PETSC_MODULE:-petsc/3.15}"
PHDF5_MODULE="${PHDF5_MODULE:-phdf5/1.12.2}"
USE_MUMPS="${USE_MUMPS:-YES}"
USE_HDF5="${USE_HDF5:-YES}"
HDF5_USE_MODULE="${HDF5_USE_MODULE:-NO}"
HDF5_VERSION="${HDF5_VERSION:-1.14.6}"
QUICK="${QUICK:-1}"
FORCE_REBUILD="${FORCE_REBUILD:-0}"
USE_GMSH="${USE_GMSH:-YES}"
GMSH_VERSION="${GMSH_VERSION:-4.13.1}"
JOBS="${JOBS:-8}"

# FORCE_REBUILD=1 wins over QUICK=1: clear caches and rebuild from source.
if [ "${FORCE_REBUILD}" = "1" ] || [ "${FORCE_REBUILD}" = "YES" ] || [ "${FORCE_REBUILD}" = "yes" ]; then
    QUICK=0
fi
ZFP_VERSION="${ZFP_VERSION:-1.0.1}"
H5Z_ZFP_VERSION="${H5Z_ZFP_VERSION:-v1.1.1}"

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
module load fftw3/3.3.8 2>/dev/null || true   # PETSc links against libfftw3_mpi
# Phase 6 ParaView output uses VTKHDF + H5Z-ZFP.  We build HDF5 1.14.x
# from source by default (see build_hdf5() below); only load Frontera's
# phdf5 module when HDF5_USE_MODULE=YES.  Older Frontera modules
# (1.10/1.12) load fine against intel/19 but lack symbols MFEM's
# mesh/vtkhdf.cpp needs (H5S_BLOCK, the unsigned-long-long hsize_t typedef
# both arrived in HDF5 1.14).
if { [ "${USE_HDF5}" = "YES" ] || [ "${USE_HDF5}" = "1" ] || [ "${USE_HDF5}" = "yes" ]; } &&
   { [ "${HDF5_USE_MODULE}" = "YES" ] || [ "${HDF5_USE_MODULE}" = "1" ] || [ "${HDF5_USE_MODULE}" = "yes" ]; }; then
    if ! module load "${PHDF5_MODULE}" 2>&1; then
        echo "ERROR: failed to load ${PHDF5_MODULE}."
        echo "  Run 'module spider ${PHDF5_MODULE}' to see compatible prerequisites."
        echo "  Pick a phdf5 version compatible with intel/19.1.1 + impi/19.0.9"
        echo "  (e.g. phdf5/1.12.2, phdf5/1.12.0, or phdf5/1.10.4) and set"
        echo "  PHDF5_MODULE accordingly, or re-run with USE_HDF5=NO to skip"
        echo "  Phase 6 VTKHDF output entirely."
        exit 1
    fi
    echo "  NOTE: HDF5_USE_MODULE=YES — using ${PHDF5_MODULE}, which is"
    echo "        1.10.x/1.12.x on intel-19.  MFEM mesh/vtkhdf.cpp build"
    echo "        will fail unless you patch the H5S_BLOCK / hsize_t uses."
fi

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
REQUIRED_VARS=(
    TACC_HYPRE_INC TACC_HYPRE_LIB
    TACC_MUMPS_INC TACC_MUMPS_LIB
    TACC_PARMETIS_INC TACC_PARMETIS_LIB
    TACC_FFTW3_LIB
    MKLROOT TACC_MKL_LIB
)
if { [ "${USE_HDF5}" = "YES" ] || [ "${USE_HDF5}" = "1" ] || [ "${USE_HDF5}" = "yes" ]; } &&
   { [ "${HDF5_USE_MODULE}" = "YES" ] || [ "${HDF5_USE_MODULE}" = "1" ] || [ "${HDF5_USE_MODULE}" = "yes" ]; }; then
    # Only required when sourcing HDF5 from a module; the from-source
    # path resolves its own HDF5_DIR / INC / LIB below.
    REQUIRED_VARS+=(TACC_HDF5_INC TACC_HDF5_LIB)
fi
for var in "${REQUIRED_VARS[@]}"; do
    # Indirect expansion with a default keeps the check working under
    # `set -u` — `eval echo \$$var` would itself trip the unbound-var
    # error before we could print the intended diagnostic.
    val="${!var:-}"
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

# Install prefixes for the hdf5 / zfp / h5z-zfp builds.  Default to a
# tree inside the MFEM checkout so the build is self-contained.
# Override via env vars if you want them installed elsewhere (e.g.
# $WORK to match the path the Phase 6 sbatch use).
HDF5_PREFIX="${HDF5_PREFIX:-${SCRIPT_DIR}/extern/hdf5/install}"
ZFP_PREFIX="${ZFP_PREFIX:-${SCRIPT_DIR}/extern/zfp/install}"
H5Z_ZFP_PREFIX="${H5Z_ZFP_PREFIX:-${SCRIPT_DIR}/extern/h5z-zfp/install}"
GMSH_PREFIX="${GMSH_PREFIX:-${SCRIPT_DIR}/extern/gmsh}"

# Build HDF5 ${HDF5_VERSION} from source, configured with --enable-parallel
# against the loaded intel/19 + impi/19 stack.  HDF5 1.14.x is required by
# MFEM mesh/vtkhdf.cpp (H5S_BLOCK + unsigned-long-long hsize_t typedef);
# Frontera's only intel-19-compatible phdf5 modules cap at 1.12.2.
build_hdf5() {
    if [ "${USE_HDF5_RESOLVED}" != "YES" ]; then
        echo ""
        echo "=== Skipping HDF5 build (USE_HDF5=${USE_HDF5}) ==="
        return 0
    fi
    if [ "${HDF5_USE_MODULE}" = "YES" ] || [ "${HDF5_USE_MODULE}" = "1" ] || [ "${HDF5_USE_MODULE}" = "yes" ]; then
        echo ""
        echo "=== Skipping HDF5 build (HDF5_USE_MODULE=YES — using ${PHDF5_MODULE}) ==="
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
        # Out-of-tree build to keep the source clean across runs.
        rm -rf build_${HDF5_VERSION}
        mkdir -p build_${HDF5_VERSION}
        cd build_${HDF5_VERSION}
        CC="$(which mpicc)" CXX="$(which mpicxx)" \
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

# Resolve HDF5_DIR / INC / LIB.  Two paths:
#   (a) HDF5_USE_MODULE=YES: use the loaded phdf5 module (TACC_HDF5_*).
#   (b) default: build HDF5 from source under ${HDF5_PREFIX}.
if [ "${USE_HDF5}" = "YES" ] || [ "${USE_HDF5}" = "1" ] || [ "${USE_HDF5}" = "yes" ]; then
    USE_HDF5_RESOLVED="YES"

    build_hdf5

    if [ "${HDF5_USE_MODULE}" = "YES" ] || [ "${HDF5_USE_MODULE}" = "1" ] || [ "${HDF5_USE_MODULE}" = "yes" ]; then
        HDF5_DIR_RESOLVED="$(cd "${TACC_HDF5_LIB}/.." 2>/dev/null && pwd || true)"
        HDF5_INC_RESOLVED="${TACC_HDF5_INC}"
        HDF5_LIB_RESOLVED="${TACC_HDF5_LIB}"
    else
        HDF5_DIR_RESOLVED="${HDF5_PREFIX}"
        HDF5_INC_RESOLVED="${HDF5_PREFIX}/include"
        HDF5_LIB_RESOLVED="${HDF5_PREFIX}/lib"
    fi

    if [ -z "${HDF5_DIR_RESOLVED}" ] || [ ! -f "${HDF5_INC_RESOLVED}/hdf5.h" ]; then
        echo "ERROR: could not locate hdf5.h under HDF5_INC=${HDF5_INC_RESOLVED:-<unset>}."
        exit 1
    fi
    echo "  HDF5_DIR     = ${HDF5_DIR_RESOLVED}"
    # Detect parallel-HDF5 support so the user sees the same flag MFEM
    # is going to check (`H5_HAVE_PARALLEL` in `<hdf5.h>` ⇒
    # `MFEM_PARALLEL_HDF5`, see mesh/vtkhdf.hpp:25).  Required for
    # `ParaViewOutput<ParMesh>::DefaultVolumeOutputMode() == Hdf5`.
    if grep -q '^#define *H5_HAVE_PARALLEL' "${HDF5_INC_RESOLVED}/H5pubconf.h" 2>/dev/null; then
        echo "  HDF5         = parallel (H5_HAVE_PARALLEL defined)"
    else
        echo "  WARNING: HDF5 build does NOT define H5_HAVE_PARALLEL — Phase 6"
        echo "           ParMesh runs will fail at the static_assert in"
        echo "           paraview_output.hpp:126.  Rebuild HDF5 with --enable-parallel."
    fi
else
    USE_HDF5_RESOLVED="NO"
fi

# Build gmsh ${GMSH_VERSION} from source into ${GMSH_PREFIX}.  Frontera
# does not ship a gmsh module, and the Phase 6 TPV meshes (~125 MB
# each) are too large to commit to git; the .geo sources are tracked and
# the Phase 6 sbatch use ${GMSH_PREFIX}/bin/gmsh to regenerate the
# matching .msh at submit time.
#
# Why source build instead of binary releases:
#   * gmsh.info Linux64 tarball requires glibc 2.23+ (Frontera = 2.17).
#   * The PyPI manylinux1 wheel claims glibc 2.5 baseline but the actual
#     libgmsh.so.X has glibc 2.23 symbols anyway (manylinux1 violation).
# Both fail on CentOS 7 with `GLIBC_2.23 not found`.  Building from
# source against the loaded intel/19 toolchain gives a binary that uses
# CentOS 7's glibc 2.17, which actually works.
#
# Minimal CMake config: no FLTK (no GUI), no OpenCASCADE (the SCEC .geo
# files only use Gmsh's built-in geometry kernel), no MED, no HDF5
# (Gmsh ships .msh files, we don't need HDF5 mesh output), no PETSc/
# SLEPc, no MPI, no Python bindings.  This keeps the build under ~5 min
# and the install under ~50 MB.
build_gmsh() {
    if [ "${USE_GMSH}" != "YES" ] && [ "${USE_GMSH}" != "1" ] && [ "${USE_GMSH}" != "yes" ]; then
        echo ""
        echo "=== Skipping gmsh build (USE_GMSH=${USE_GMSH}) ==="
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
        echo "=== Cache miss or stale install: rebuilding gmsh ==="
    fi

    echo ""
    echo "=== Building gmsh ${GMSH_VERSION} from source (CentOS 7 / glibc 2.17) ==="
    mkdir -p "${SCRIPT_DIR}/extern"

    local tarball="${SCRIPT_DIR}/extern/gmsh-${GMSH_VERSION}-source.tgz"
    local url="https://gmsh.info/src/gmsh-${GMSH_VERSION}-source.tgz"
    if [ ! -f "${tarball}" ]; then
        echo "  Downloading source tarball from ${url}"
        if command -v wget >/dev/null 2>&1; then
            wget -q -O "${tarball}" "${url}"
        elif command -v curl >/dev/null 2>&1; then
            curl -sL -o "${tarball}" "${url}"
        else
            echo "ERROR: neither wget nor curl available."
            exit 1
        fi
    fi

    local gmsh_src="${SCRIPT_DIR}/extern/gmsh-src"
    rm -rf "${gmsh_src}" "${GMSH_PREFIX}"
    mkdir -p "${gmsh_src}"
    tar -xzf "${tarball}" -C "${gmsh_src}" --strip-components=1

    # Out-of-tree build.  Force the intel/19 + impi/19 MPI wrappers as
    # the C/C++ compilers — they share libstdc++ with the rest of the
    # MFEM build so there's no ABI mismatch when downstream tools mix
    # binaries.  All optional deps OFF to keep the build minimal.
    (
        cd "${gmsh_src}"
        mkdir -p build_${GMSH_VERSION}
        cd build_${GMSH_VERSION}
        cmake .. \
            -DCMAKE_INSTALL_PREFIX="${GMSH_PREFIX}" \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_C_COMPILER="$(which mpicc)" \
            -DCMAKE_CXX_COMPILER="$(which mpicxx)" \
            -DENABLE_FLTK=OFF \
            -DENABLE_OCC=OFF \
            -DENABLE_MED=OFF \
            -DENABLE_HDF5=OFF \
            -DENABLE_PETSC=OFF \
            -DENABLE_SLEPC=OFF \
            -DENABLE_MPI=OFF \
            -DENABLE_PYTHON=OFF \
            -DENABLE_BUILD_LIB=OFF \
            -DENABLE_BUILD_DYNAMIC=OFF
        cmake --build . --target install -j "${JOBS}"
    )

    if [ ! -x "${gmsh_bin}" ]; then
        echo "ERROR: gmsh build finished but ${gmsh_bin} missing."
        exit 1
    fi
    if ! "${gmsh_bin}" -version >/dev/null 2>&1; then
        echo "WARNING: ${gmsh_bin} -version failed.  Check ldd output."
    fi
    echo "  gmsh ${GMSH_VERSION} built at ${gmsh_bin}"
}

build_gmsh

# Build the H5Z-ZFP plugin (+ its dependency, ZFP) for Phase 6 lossy
# floating-point compression in VTKHDF output.  The plugin lands at
# ${H5Z_ZFP_PREFIX}/plugin/libh5zzfp.so and is discovered at runtime
# via HDF5_PLUGIN_PATH (the Phase 6 sbatch do this conditionally).
build_zfp_and_h5z_zfp() {
    if [ "${USE_HDF5_RESOLVED}" != "YES" ]; then
        echo ""
        echo "=== Skipping zfp / h5z-zfp build (USE_HDF5=${USE_HDF5}) ==="
        return 0
    fi

    local zfp_lib_so="${ZFP_PREFIX}/lib64/libzfp.so"
    local zfp_lib_so_alt="${ZFP_PREFIX}/lib/libzfp.so"
    local plugin_so="${H5Z_ZFP_PREFIX}/plugin/libh5zzfp.so"
    local plugin_dylib="${H5Z_ZFP_PREFIX}/plugin/libh5zzfp.dylib"
    if [ "${QUICK}" = "1" ] || [ "${QUICK}" = "YES" ]; then
        # Reuse only when BOTH ZFP and the H5Z-ZFP plugin look complete —
        # otherwise a half-finished extern/ tree would silently succeed.
        if { [ -f "${zfp_lib_so}" ] || [ -f "${zfp_lib_so_alt}" ]; } &&
           { [ -f "${plugin_so}" ] || [ -f "${plugin_dylib}" ]; }; then
            echo ""
            echo "=== Reusing existing ZFP + H5Z-ZFP at ${H5Z_ZFP_PREFIX} ==="
            echo "    (set FORCE_REBUILD=1 to clear and rebuild — necessary"
            echo "     when HDF5_VERSION or HDF5_USE_MODULE has changed)"
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
        # Fetch the tag if the shallow clone landed on a different ref
        # (e.g. when the directory was pre-populated).
        git fetch --tags --quiet || true
        git checkout --quiet "${ZFP_VERSION}" 2>/dev/null \
            || git checkout --quiet "tags/${ZFP_VERSION}" 2>/dev/null || true
        mkdir -p build
        cd build
        cmake .. \
            -DCMAKE_INSTALL_PREFIX="${ZFP_PREFIX}" \
            -DCMAKE_C_COMPILER="$(which mpicc)" \
            -DBUILD_SHARED_LIBS=ON \
            -DBUILD_TESTING=OFF \
            -DBUILD_EXAMPLES=OFF \
            -DZFP_WITH_OPENMP=OFF
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
        # The H5Z-ZFP Makefile reads HDF5_HOME / ZFP_HOME and writes
        # the plugin to ${PREFIX}/plugin.  Use the Makefile build path
        # (not the cmake one) because it is the canonical install
        # layout the runtime probe in paraview_output.hpp expects.
        # Force a clean build so an existing plugin that was linked
        # against a different HDF5 (e.g. an older Frontera module from a
        # prior run) gets rebuilt against ${HDF5_DIR_RESOLVED}.
        make HDF5_HOME="${HDF5_DIR_RESOLVED}" \
             ZFP_HOME="${ZFP_PREFIX}" \
             PREFIX="${H5Z_ZFP_PREFIX}" \
             clean || true
        make CC="$(which mpicc)" \
             HDF5_HOME="${HDF5_DIR_RESOLVED}" \
             ZFP_HOME="${ZFP_PREFIX}" \
             PREFIX="${H5Z_ZFP_PREFIX}" \
             install -j "${JOBS}"
    )

    if [ ! -f "${plugin_so}" ] && [ ! -f "${plugin_dylib}" ]; then
        echo "ERROR: H5Z-ZFP build finished but ${plugin_so} is missing."
        exit 1
    fi
    echo "  H5Z-ZFP plugin installed at ${H5Z_ZFP_PREFIX}/plugin"
}

build_zfp_and_h5z_zfp

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

# Phase 6: VTKHDF + H5Z-ZFP.  MFEM_USE_H5Z_ZFP is a pure preprocessor
# flag (no link line — the plugin is loaded at runtime via
# HDF5_PLUGIN_PATH); the only link flags are HDF5_OPT / HDF5_LIB
# (the parallel-HDF5 module's libhdf5{_hl}).
if [ "${USE_HDF5_RESOLVED}" = "YES" ]; then
    CONFIG_ARGS+=(
      MFEM_USE_HDF5=YES
      MFEM_USE_H5Z_ZFP=YES
      HDF5_OPT="-I${HDF5_INC_RESOLVED}"
      HDF5_LIB="-L${HDF5_LIB_RESOLVED} -Wl,-rpath,${HDF5_LIB_RESOLVED} -lhdf5_hl -lhdf5 -lz"
    )
fi

make config "${CONFIG_ARGS[@]}"

echo ""
echo "=== Building MFEM library ==="
make -j"${JOBS}"

echo ""
echo "=== Building SEAS miniapps ==="
cd miniapps/seas
# Build target list:
#   - Everything `make all` would build (MINIAPPS = SEQ_MINIAPPS +
#     PAR_MINIAPPS + SEQ_LONG + PAR_LONG, see Makefile:570-700) ...
#   - ... minus all targets that instantiate SEASQuasiDynamicOperator
#     with AntiplaneDomainOperator.  That template calls
#     ComputeTractionDiagnostics / IsFirstStepDebugEnabled on the domain
#     operator (added April 2026 in commits 8cb7c00 / b631c63 for BP5
#     face tracing) — methods only ElasticityDomainOperator exposes.
#     Deferred fix; see
#     miniapps/seas/debug_document/bp5_debug_document/
#       antiplane_template_instantiation_2026-05-11.md
#     The pattern (*bp1*|*bp2*|*antiplane*|*pseas*) catches most; the
#     EXPLICIT_SKIPS list catches the seven the pattern misses.
#   - Plus the tpv102/104/205 dynamic-rupture drivers, which the seas
#     Makefile keeps outside MINIAPPS.
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
    echo "  (see miniapps/seas/debug_document/bp5_debug_document/"
    echo "   antiplane_template_instantiation_2026-05-11.md)"
fi
make ${FILTERED_TARGETS} \
     seas_tpv102_driver seas_tpv104_driver seas_tpv205_driver \
     -j"${JOBS}"

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
echo "Key binaries (full list under $(pwd)):"
echo "  ./seas_bp5_full                            # BP5 verification driver"
echo "  ./seas_tpv102_driver"
echo "  ./seas_tpv104_driver"
echo "  ./seas_tpv205_driver"
echo "  ./seas_test_parallel_elasticity"
echo "  ./seas_test_bp5_parallel_smoke"
echo "  ./seas_test_bp5_integration"
echo "  ./seas_test_bp5_fault_operator"
echo "  ./seas_test_friction_law                   # (and ~50 other unit/parallel tests)"
if [ "${DRIVER_BUILT}" = "1" ]; then
    echo "  $(pwd)/seas_driver          (TOML-based, new code paths)"
fi
if [ "${USE_HDF5_RESOLVED}" = "YES" ]; then
    echo ""
    echo "H5Z-ZFP plugin:"
    echo "  ${H5Z_ZFP_PREFIX}/plugin/libh5zzfp.so"
    echo "  (HDF5_PLUGIN_PATH=${H5Z_ZFP_PREFIX}/plugin at runtime — the"
    echo "   Phase 6 sbatch in miniapps/seas/jobs/{bp5,tpv*}/*phase6*"
    echo "   already point at \${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin,"
    echo "   so no manual override is needed if you submit from this checkout.)"
fi
echo ""
echo "Run examples:"
echo "  ibrun ./seas_bp5_full --petsc-ts"
if [ "${DRIVER_BUILT}" = "1" ]; then
    echo "  ibrun ./seas_driver config/bp5_example.toml"
fi
