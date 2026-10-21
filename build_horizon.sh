#!/bin/bash
# Build the dependency stack, MFEM and the SEAS BP5 driver for the GPU plan
# (GPU BP5 on Horizon, phases G0 and G1).
#
#   MODE=cpu (default)  CPU-only build for Vista's Grace-Grace (gg) nodes:
#                       G0's CPU baseline and goldens.  No CUDA, no cuDSS.
#   MODE=cuda           CUDA build (G1): PETSc with CUDA (hypre on the GPU;
#                       MUMPS stays on the CPU), cuDSS, and MFEM with
#                       MFEM_USE_CUDA and MFEM_USE_CUDSS.  Same versions as
#                       MODE=cpu, so a CPU golden and a GPU run differ only
#                       by the GPU work, the CUDA libraries and nvcc's host
#                       compilation.  GPU_ARCH picks the GPU (90: Vista gh,
#                       GH200; 100: Horizon gb, GB200).  MFEM is configured
#                       in place, so a MODE=cuda build needs its own clone
#                       (one checkout holds one mode; the script checks).
#
# Every library is built from source at a pinned version (Vista's PETSc and
# MUMPS modules cannot be used: petsc/3.24 does not load with gcc/13.2.0 +
# openmpi/5.0.11, and mumps/5.8.0 is the sequential build), so Vista and
# Horizon get the same versions.  Only the compiler, MPI, NVPL (BLAS/LAPACK)
# and cmake come from modules (MODULES below).
#
# Stages (STAGE=..., several separated by spaces; default "deps mfem seas"):
#   fetch  On a LOGIN node, which has internet (compute nodes may not):
#          download the PETSc release and every package tarball its configure
#          asks for into ${DEPS_ROOT}/downloads (MODE=cuda: also the cuDSS
#          wheel from PyPI).
#   deps   PETSc with MUMPS, ScaLAPACK, METIS, ParMETIS, PT-Scotch and hypre,
#          built from those tarballs, installed in ${PETSC_PREFIX}
#          (MODE=cuda: also unpacks cuDSS into ${CUDSS_DIR}).
#   mfem   Configure and build MFEM (this checkout, in place) against it.
#   seas   Build seas_bp5_full and the unit tests G0 runs
#          (seas_test_elasticity_operator, seas_test_g0_matrix_only).
#
# Vista login nodes give one core per user, so run deps, mfem and seas in a
# compute job (jobs/vista/g0_build.sbatch).  Each stage writes what it used
# to ${DEPS_ROOT}/build_provenance_${MODE}.txt.
#
# Environment overrides:
#   MODE=cpu|cuda               (default cpu)
#   STAGE="deps mfem seas"      stages to run, in order
#   DEPS_ROOT=$WORK/seas-gpu-deps
#   PETSC_VERSION=3.24.6        PETSc release (its configure pins the package
#                               versions: MUMPS, hypre, METIS, ...)
#   MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
#                               set MODULES=none to use the current shell
#   BLASLAPACK_LIB="..."        override the NVPL BLAS/LAPACK link line
#                               (-L/-l flags only; rpaths are added for MFEM)
#   CPU_FLAGS="-mcpu=neoverse-v2"   (Grace on both Vista and Horizon)
#   JOBS=<n>                    parallel make (default: nproc)
#   MODE=cuda only:
#   CUDA_MODULE=cuda/13.0       CUDA toolkit module added to MODULES
#   GPU_ARCH=90                 CUDA compute capability (90 GH200, 100 GB200)
#   CUDSS_VERSION=0.7.1.6       cuDSS release, from the PyPI wheel
#                               nvidia-cudss-cu<CUDA major> (aarch64)
#   CUDSS_DIR=<dir>             cuDSS install to use (default: the wheel
#                               unpacked in ${DEPS_ROOT}/cudss-${CUDSS_VERSION})

set -euo pipefail

MODE="${MODE:-cpu}"
STAGE="${STAGE:-deps mfem seas}"
PETSC_VERSION="${PETSC_VERSION:-3.24.6}"
DEPS_ROOT="${DEPS_ROOT:-${WORK:-$HOME}/seas-gpu-deps}"
CPU_MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
# cuda/13.x is the CUDA line whose openmpi/5.0.11 build on Vista is
# CUDA-aware (with cuda/12.x Lmod loads the CPU-only OpenMPI).
CUDA_MODULE="${CUDA_MODULE:-cuda/13.0}"
GPU_ARCH="${GPU_ARCH:-90}"
CUDSS_VERSION="${CUDSS_VERSION:-0.7.1.6}"
if [ "${MODE}" = "cuda" ]; then
   # Load the CUDA module before openmpi so Lmod picks its CUDA-aware build.
   MODULES="${MODULES:-gcc/13.2.0 ${CUDA_MODULE} openmpi/5.0.11 nvpl/26.9 cmake/3.31.5}"
else
   MODULES="${MODULES:-${CPU_MODULES}}"
fi
CPU_FLAGS="${CPU_FLAGS:--mcpu=neoverse-v2}"
JOBS="${JOBS:-$(nproc)}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOWNLOADS="${DEPS_ROOT}/downloads"
PETSC_SRC="${DEPS_ROOT}/src/petsc-${PETSC_VERSION}"
PETSC_PREFIX="${DEPS_ROOT}/petsc-${PETSC_VERSION}-${MODE}"
# Build directory inside the shared PETSc source tree, one per mode (cpu keeps
# the name PETSc's configure picks by default).
if [ "${MODE}" = "cuda" ]; then
   PETSC_BUILD_ARCH=arch-linux-cuda-opt
else
   PETSC_BUILD_ARCH=arch-linux-c-opt
fi
PROVENANCE="${DEPS_ROOT}/build_provenance_${MODE}.txt"
PETSC_TARBALL="petsc-${PETSC_VERSION}.tar.gz"
PETSC_URLS=(
   "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/${PETSC_TARBALL}"
   "https://gitlab.com/petsc/petsc/-/archive/v${PETSC_VERSION}/petsc-v${PETSC_VERSION}.tar.gz"
)

if [ "${MODE}" != "cpu" ] && [ "${MODE}" != "cuda" ]; then
   echo "ERROR: MODE must be cpu or cuda (got '${MODE}')." >&2
   exit 1
fi
CUDSS_DIR="${CUDSS_DIR:-${DEPS_ROOT}/cudss-${CUDSS_VERSION}}"
CUDSS_CUDA=""   # CUDA major version of the cuDSS wheel, from nvcc (below)
# One checkout holds one MFEM configuration (MFEM is configured in place).
if [ -f "${SCRIPT_DIR}/config/config.mk" ]; then
   _cfg_cuda="$(grep -m1 '^MFEM_USE_CUDA ' "${SCRIPT_DIR}/config/config.mk" | awk '{print $3}')"
   if { [ "${MODE}" = "cuda" ] && [ "${_cfg_cuda}" = "NO" ]; } ||
      { [ "${MODE}" = "cpu" ] && [ "${_cfg_cuda}" = "YES" ]; }; then
      echo "ERROR: this checkout is configured with MFEM_USE_CUDA=${_cfg_cuda};" \
           "use a separate clone for MODE=${MODE}." >&2
      exit 1
   fi
   unset _cfg_cuda
fi

log() { echo "[build_horizon ${MODE}] $*"; }

# --------------------------------------------------------------------------
# Toolchain
# --------------------------------------------------------------------------
if [ "${MODULES}" != "none" ]; then
   if ! type module >/dev/null 2>&1; then
      echo "ERROR: environment modules are not available in this shell." >&2
      exit 1
   fi
   # Lmod swaps the default nvidia compiler for gcc and reloads MPI/NVPL.
   # shellcheck disable=SC2086
   module load ${MODULES}
fi

if [ "${MODE}" = "cuda" ]; then
   command -v nvcc >/dev/null || { echo "ERROR: MODE=cuda needs nvcc (CUDA_MODULE)." >&2; exit 1; }
   CUDSS_CUDA="$(nvcc --version | sed -n 's/.*release \([0-9]*\)\..*/\1/p')"
   [ -n "${CUDSS_CUDA}" ] || { echo "ERROR: cannot read the CUDA version from nvcc." >&2; exit 1; }
   # hypre, cuDSS and the device-buffer test pass GPU buffers to MPI.
   # (No grep -q: it stops reading at the first match, and under pipefail a
   # SIGPIPE'd ompi_info would fail the test.)
   ompi_info --parsable --all 2>/dev/null | grep 'built_with_cuda_support:value:true' >/dev/null || {
      echo "ERROR: the loaded MPI is not CUDA-aware (ompi_info); check MODULES." >&2
      exit 1
   }
fi
CUDSS_WHEEL="nvidia_cudss_cu${CUDSS_CUDA}-${CUDSS_VERSION}-py3-none-manylinux_2_17_aarch64.whl"

BLASLAPACK_LIB="${BLASLAPACK_LIB:-}"
if [ -z "${BLASLAPACK_LIB}" ]; then
   if [ -z "${TACC_NVPL_LIB:-}" ]; then
      echo "ERROR: TACC_NVPL_LIB is not set (nvpl module not loaded) and" \
           "BLASLAPACK_LIB is empty." >&2
      exit 1
   fi
   # Sequential NVPL: one MPI rank per core, no BLAS threads (threaded BLAS
   # under MPI oversubscribes cores and makes results run-dependent).
   BLASLAPACK_LIB="-L${TACC_NVPL_LIB} -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core"
fi
# The same libraries for MFEM's link line, with an rpath for every -L dir.
MFEM_LAPACK_LIB="${BLASLAPACK_LIB}"
for flag in ${BLASLAPACK_LIB}; do
   case "${flag}" in
      -L*) MFEM_LAPACK_LIB="-Wl,-rpath,${flag#-L} ${MFEM_LAPACK_LIB}" ;;
   esac
done

record_toolchain() {
   {
      echo "# build_horizon.sh MODE=${MODE} stage '$1' on $(hostname) at $(date -u +%FT%TZ)"
      echo "source: $(git -C "${SCRIPT_DIR}" rev-parse HEAD 2>/dev/null || echo unknown)" \
           "($(git -C "${SCRIPT_DIR}" rev-parse --abbrev-ref HEAD 2>/dev/null || echo ?))," \
           "dirty files: $(git -C "${SCRIPT_DIR}" status --porcelain 2>/dev/null | wc -l)"
      echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
      echo "mpicc: $(mpicc --version 2>&1 | head -1)"
      echo "mpif90: $(mpif90 --version 2>&1 | head -1)"
      echo "mpi: $(ompi_info 2>/dev/null | grep -m1 'Open MPI:' | xargs || mpichversion 2>/dev/null | head -1)"
      echo "blas/lapack: ${BLASLAPACK_LIB}"
      echo "cpu flags: ${CPU_FLAGS}"
      if [ "${MODE}" = "cuda" ]; then
         echo "nvcc: $(nvcc --version 2>&1 | grep -m1 release)"
         echo "gpu arch: sm_${GPU_ARCH}"
         echo "cudss: ${CUDSS_DIR} (CUDSS_VERSION ${CUDSS_VERSION})"
      fi
      echo "petsc: ${PETSC_VERSION} at ${PETSC_PREFIX}"
   } >> "${PROVENANCE}"
}

mkdir -p "${DEPS_ROOT}" "${DOWNLOADS}" "${DEPS_ROOT}/src"

# --------------------------------------------------------------------------
# PETSc configure arguments (shared by fetch and deps)
# --------------------------------------------------------------------------
PETSC_CONFIGURE_ARGS=(
   --prefix="${PETSC_PREFIX}"
   --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90
   --with-debugging=0
   "COPTFLAGS=-O3 ${CPU_FLAGS}"
   "CXXOPTFLAGS=-O3 ${CPU_FLAGS}"
   "FOPTFLAGS=-O3 ${CPU_FLAGS}"
   --with-shared-libraries=1
   --with-64-bit-indices=0
   --with-x=0
   "--with-blaslapack-lib=${BLASLAPACK_LIB}"
   --download-scalapack
   --download-mumps
   --download-metis
   --download-parmetis
   --download-ptscotch
   --download-hypre
   --with-packages-download-dir="${DOWNLOADS}"
)
if [ "${MODE}" = "cuda" ]; then
   # hypre (downloaded by PETSc) is then built with CUDA too; MUMPS and the
   # rest stay on the CPU.
   PETSC_CONFIGURE_ARGS+=(--with-cuda=1 --with-cuda-arch="${GPU_ARCH}")
else
   PETSC_CONFIGURE_ARGS+=(--with-cuda=0)
fi

unpack_petsc() {
   if [ -f "${PETSC_SRC}/configure" ]; then return; fi
   [ -f "${DOWNLOADS}/${PETSC_TARBALL}" ] || {
      echo "ERROR: ${DOWNLOADS}/${PETSC_TARBALL} missing; run STAGE=fetch on a login node." >&2
      exit 1
   }
   log "unpacking ${PETSC_TARBALL}"
   rm -rf "${PETSC_SRC}.tmp"
   mkdir -p "${PETSC_SRC}.tmp"
   tar -xzf "${DOWNLOADS}/${PETSC_TARBALL}" -C "${PETSC_SRC}.tmp" --strip-components=1
   mv "${PETSC_SRC}.tmp" "${PETSC_SRC}"
}

# --------------------------------------------------------------------------
# Stage: fetch (login node)
# --------------------------------------------------------------------------
stage_fetch() {
   if [ ! -f "${DOWNLOADS}/${PETSC_TARBALL}" ]; then
      local ok=0
      for url in "${PETSC_URLS[@]}"; do
         log "downloading ${url}"
         if curl -fL --retry 3 -o "${DOWNLOADS}/${PETSC_TARBALL}.part" "${url}"; then
            mv "${DOWNLOADS}/${PETSC_TARBALL}.part" "${DOWNLOADS}/${PETSC_TARBALL}"
            ok=1
            break
         fi
      done
      [ "${ok}" = 1 ] || { echo "ERROR: could not download ${PETSC_TARBALL}" >&2; exit 1; }
   fi
   unpack_petsc

   # With --with-packages-download-dir, PETSc's configure lists the package
   # tarballs missing from that directory and stops.  Each line reads:
   # name ['url1', 'url2', ...].  When none is missing it goes on to build
   # the packages instead, so it is asked with an empty directory (it then
   # lists them all), a prefix inside it and a PETSC_ARCH of its own, which
   # leaves the build directories, the install prefix and DOWNLOADS alone.
   log "asking PETSc's configure which package tarballs it needs"
   local listing="${DEPS_ROOT}/petsc_download_list.txt"
   local probe probe_arch=arch-fetch-probe
   probe="$(mktemp -d "${DEPS_ROOT}/fetch-probe.XXXXXX")"
   local probe_args=() a
   for a in "${PETSC_CONFIGURE_ARGS[@]}"; do
      case "${a}" in
         --prefix=*) probe_args+=(--prefix="${probe}/prefix") ;;
         --with-packages-download-dir=*)
            probe_args+=(--with-packages-download-dir="${probe}") ;;
         *) probe_args+=("${a}") ;;
      esac
   done
   ( cd "${PETSC_SRC}" && ./configure PETSC_ARCH="${probe_arch}" "${probe_args[@]}" ) \
      > "${listing}" 2>&1 || true
   rm -rf "${probe}" "${PETSC_SRC:?}/${probe_arch}"
   python3 - "${listing}" "${DOWNLOADS}" <<'PY'
import ast, os, re, subprocess, sys
listing, dest = sys.argv[1], sys.argv[2]
text = open(listing).read()
pkgs = re.findall(r"^\s*([A-Za-z0-9_\-]+)\s+(\[.*\])\s*$", text, re.M)
if not pkgs:
    sys.exit("ERROR: no package list in PETSc's configure output; see " + listing)
failed = []
for name, urls in pkgs:
    urls = [u for u in ast.literal_eval(urls) if u.startswith("http")]
    have = [f for f in os.listdir(dest)
            if any(os.path.basename(u.split("?")[0]) == f for u in urls)]
    if have:
        print(f"  {name}: have {have[0]}")
        continue
    for u in urls:
        out = os.path.join(dest, os.path.basename(u.split("?")[0]))
        print(f"  {name}: {u}")
        if subprocess.call(["curl", "-fL", "--retry", "3", "-o", out + ".part", u]) == 0:
            os.replace(out + ".part", out)
            break
        if os.path.exists(out + ".part"):
            os.remove(out + ".part")
    else:
        failed.append(name)
if failed:
    sys.exit("ERROR: could not download: " + ", ".join(failed))
PY
   if [ "${MODE}" = "cuda" ] && [ ! -f "${DOWNLOADS}/${CUDSS_WHEEL}" ]; then
      log "downloading cuDSS ${CUDSS_VERSION} (PyPI wheel nvidia-cudss-cu${CUDSS_CUDA}, aarch64)"
      python3 -m pip download --no-deps --only-binary=:all: \
         --platform manylinux_2_17_aarch64 --python-version 3 \
         -d "${DOWNLOADS}" "nvidia-cudss-cu${CUDSS_CUDA}==${CUDSS_VERSION}"
      [ -f "${DOWNLOADS}/${CUDSS_WHEEL}" ] || {
         echo "ERROR: ${CUDSS_WHEEL} not in ${DOWNLOADS} after pip download" >&2
         exit 1
      }
   fi
   log "downloads in ${DOWNLOADS}:"
   ls -la "${DOWNLOADS}"
   record_toolchain fetch
}

# cuDSS from its wheel: nvidia/cuNN/{include,lib} -> ${CUDSS_DIR}/{include,lib},
# plus the unversioned .so names the link line and MFEM's cuDSS settings use.
unpack_cudss() {
   if [ -f "${CUDSS_DIR}/include/cudss.h" ]; then return; fi
   [ -f "${DOWNLOADS}/${CUDSS_WHEEL}" ] || {
      echo "ERROR: ${DOWNLOADS}/${CUDSS_WHEEL} missing; run STAGE=fetch on a login node," \
           "or set CUDSS_DIR to an existing cuDSS install." >&2
      exit 1
   }
   log "unpacking ${CUDSS_WHEEL} into ${CUDSS_DIR}"
   python3 - "${DOWNLOADS}/${CUDSS_WHEEL}" "${CUDSS_DIR}" "${CUDSS_CUDA}" <<'PY'
import os, sys, zipfile
whl, dest, cuda = sys.argv[1], sys.argv[2], sys.argv[3]
with zipfile.ZipFile(whl) as z:
    for name in z.namelist():
        for sub in ("include", "lib"):
            pre = f"nvidia/cu{cuda}/{sub}/"
            if name.startswith(pre) and not name.endswith("/"):
                out = os.path.join(dest, sub, name[len(pre):])
                os.makedirs(os.path.dirname(out), exist_ok=True)
                with open(out, "wb") as f:
                    f.write(z.read(name))
libdir = os.path.join(dest, "lib")
for so in os.listdir(libdir):
    if so.endswith(".so.0"):
        link = os.path.join(libdir, so[:-2])
        if not os.path.exists(link):
            os.symlink(so, link)
PY
   [ -f "${CUDSS_DIR}/include/cudss.h" ] && [ -e "${CUDSS_DIR}/lib/libcudss.so" ] || {
      echo "ERROR: cuDSS unpack into ${CUDSS_DIR} incomplete" >&2
      exit 1
   }
}

# --------------------------------------------------------------------------
# Stage: deps (compute node)
# --------------------------------------------------------------------------
stage_deps() {
   unpack_petsc
   if [ "${MODE}" = "cuda" ]; then unpack_cudss; fi
   if [ -f "${PETSC_PREFIX}/lib/petsc/conf/petscvariables" ]; then
      log "PETSc already installed in ${PETSC_PREFIX}; skipping (remove it to rebuild)"
   else
      log "configuring PETSc ${PETSC_VERSION}"
      ( cd "${PETSC_SRC}" && ./configure PETSC_ARCH="${PETSC_BUILD_ARCH}" "${PETSC_CONFIGURE_ARGS[@]}" )
      ( cd "${PETSC_SRC}" && make PETSC_DIR="${PETSC_SRC}" PETSC_ARCH="${PETSC_BUILD_ARCH}" MAKE_NP="${JOBS}" all )
      ( cd "${PETSC_SRC}" && make PETSC_DIR="${PETSC_SRC}" PETSC_ARCH="${PETSC_BUILD_ARCH}" install )
   fi
   record_toolchain deps
   {
      echo "petsc configure options: $(grep -m1 '^CONFIGURE_OPTIONS' "${PETSC_PREFIX}/lib/petsc/conf/petscvariables" | cut -d= -f2-)"
      for pkg in mumps hypre metis parmetis ptscotch scalapack; do
         f="${PETSC_PREFIX}/lib/petsc/conf/pkg.conf.${pkg}"
         [ -f "${f}" ] && echo "${pkg}: $(grep -m1 -i 'version' "${f}" 2>/dev/null || echo '(no version line)')"
      done
      grep -h -m1 'define PETSC_PKG_MUMPS_VERSION_\(MAJOR\|MINOR\|SUBMINOR\)' "${PETSC_PREFIX}/include/petscconf.h" 2>/dev/null || true
      grep -h -m1 'HYPRE_RELEASE_VERSION' "${PETSC_PREFIX}/include/HYPRE_config.h" 2>/dev/null || true
      grep -h -E 'define HYPRE_USING_(CUDA|GPU|UNIFIED_MEMORY|GPU_AWARE_MPI|DEVICE_MALLOC_ASYNC)' \
         "${PETSC_PREFIX}/include/HYPRE_config.h" 2>/dev/null || true
   } >> "${PROVENANCE}"
   if [ "${MODE}" = "cuda" ]; then
      # The plan builds hypre for the GPU without unified memory.
      grep -q 'define HYPRE_USING_CUDA' "${PETSC_PREFIX}/include/HYPRE_config.h" || {
         echo "ERROR: PETSc's hypre was not built with CUDA (HYPRE_config.h)." >&2
         exit 1
      }
      if grep -q 'define HYPRE_USING_UNIFIED_MEMORY' "${PETSC_PREFIX}/include/HYPRE_config.h"; then
         echo "ERROR: PETSc's hypre was built with unified memory; the plan needs it off." >&2
         exit 1
      fi
   fi
}

petsc_var() {
   grep -m1 "^$1 *=" "${PETSC_PREFIX}/lib/petsc/conf/petscvariables" | cut -d= -f2-
}

# --------------------------------------------------------------------------
# Stage: mfem (compute node)
# --------------------------------------------------------------------------
stage_mfem() {
   [ -f "${PETSC_PREFIX}/lib/petsc/conf/petscvariables" ] || {
      echo "ERROR: no PETSc in ${PETSC_PREFIX}; run STAGE=deps first." >&2
      exit 1
   }
   local P="${PETSC_PREFIX}"
   local ext_libs
   ext_libs="$(petsc_var PETSC_EXTERNAL_LIB_BASIC)"
   local mfem_cuda_args=()
   local mfem_cxx_args=(MFEM_CXX=mpicxx OPTIM_FLAGS="-O3 -std=c++17 ${CPU_FLAGS}")
   local xl="-Wl,"
   if [ "${MODE}" = "cuda" ]; then
      # nvcc compiles MFEM (and seas, which builds with MFEM's compiler) with
      # mpicxx as its host compiler; host-only flags go through -Xcompiler.
      mfem_cxx_args=(MFEM_CXX=nvcc MFEM_HOST_CXX=mpicxx
                     OPTIM_FLAGS="-O3 -std=c++17 -Xcompiler=${CPU_FLAGS}")
      [ -f "${CUDSS_DIR}/include/cudss.h" ] || {
         echo "ERROR: no cuDSS in ${CUDSS_DIR}; run STAGE=deps first." >&2
         exit 1
      }
      local cuda_dir
      cuda_dir="$(dirname "$(dirname "$(command -v nvcc)")")"
      mfem_cuda_args=(
         MFEM_USE_CUDA=YES CUDA_DIR="${cuda_dir}" CUDA_ARCH="sm_${GPU_ARCH}"
         MFEM_USE_CUDSS=YES CUDSS_DIR="${CUDSS_DIR}"
      )
      # nvcc links too and rejects -Wl,: linker flags go through -Xlinker=
      # (MFEM's CUDA_XLINKER), as MFEM already does for PETSc's own libs.
      xl="-Xlinker="
   fi
   ext_libs="${ext_libs//-Wl,/${xl}}"
   local lapack_lib="${MFEM_LAPACK_LIB//-Wl,/${xl}}"
   log "configuring MFEM in ${SCRIPT_DIR}"
   cd "${SCRIPT_DIR}"
   make config \
      MFEM_USE_MPI=YES \
      "${mfem_cxx_args[@]}" \
      MFEM_USE_METIS=YES MFEM_USE_METIS_5=YES \
      METIS_DIR="${P}" METIS_OPT="-I${P}/include" \
      METIS_LIB="${xl}-rpath,${P}/lib -L${P}/lib -lmetis" \
      HYPRE_DIR="${P}" HYPRE_OPT="-I${P}/include" \
      HYPRE_LIB="${xl}-rpath,${P}/lib -L${P}/lib -lHYPRE" \
      MFEM_USE_LAPACK=YES LAPACK_LIB="${lapack_lib}" \
      MFEM_USE_MUMPS=YES MUMPS_DIR="${P}" MUMPS_OPT="-I${P}/include" \
      MUMPS_LIB="${xl}-rpath,${P}/lib -L${P}/lib -ldmumps -lmumps_common -lpord ${ext_libs}" \
      MFEM_USE_PETSC=YES PETSC_DIR="${P}" PETSC_ARCH= \
      MFEM_USE_HDF5=NO \
      "${mfem_cuda_args[@]}"
   make -j"${JOBS}" lib
   record_toolchain mfem
   grep -E '^MFEM_USE_CUDSS ' config/config.mk >> "${PROVENANCE}" || true
   grep -E '^MFEM_USE_(MPI|PETSC|MUMPS|METIS|LAPACK|CUDA|HDF5) ' config/config.mk >> "${PROVENANCE}"
}

# --------------------------------------------------------------------------
# Stage: seas (compute node)
# --------------------------------------------------------------------------
stage_seas() {
   cd "${SCRIPT_DIR}/miniapps/seas"
   # Start from fresh objects for the binaries G0 uses.
   rm -f tests/verification/bp5_verification_full.o \
         tests/unit/test_elasticity_operator.o tests/unit/test_g0_matrix_only.o
   make -j"${JOBS}" seas_bp5_full seas_test_elasticity_operator \
        seas_test_g0_matrix_only
   record_toolchain seas
   ls -la seas_bp5_full seas_test_elasticity_operator \
          seas_test_g0_matrix_only >> "${PROVENANCE}"
}

for s in ${STAGE}; do
   case "${s}" in
      fetch) stage_fetch ;;
      deps)  stage_deps ;;
      mfem)  stage_mfem ;;
      seas)  stage_seas ;;
      *) echo "ERROR: unknown stage '${s}'" >&2; exit 1 ;;
   esac
done
log "done: ${STAGE}"
