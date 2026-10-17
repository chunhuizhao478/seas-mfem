#!/bin/bash
# jobs/fd_cycle/build_discovery.sh — Discovery (USC CARC) build of the fully
# dynamic cycle branch for the SCEC condo nodes (cluster condo, partition
# scec: Xeon E5-2640 v4 Broadwell), Chunhui Zhao, 2026-10-16.
#
# No Discovery module gives PETSc with MUMPS on gcc/13.3.0 + openmpi/5.0.5
# (petsc/3.22.2 has no MUMPS; mumps/5.6.2 needs hpcx-mpi), so this follows
# the Vista G0 recipe (build_horizon.sh MODE=cpu on gpu-qd-bp5): PETSc
# 3.24.6, the version of the Vista QD reference run and of the Mac, built
# from source with MUMPS, ScaLAPACK, METIS, ParMETIS, PT-Scotch and hypre.
# Compiler, MPI, BLAS/LAPACK (OpenBLAS, 32-bit integers), HDF5, cmake and
# python (PETSc's configure; the system python3 is 3.6) come from modules.
#
# Layout: the checkout and everything it builds stay in one folder, e.g.
#   /scratch1/$USER/fdcycle-discovery-20261016/seas-mfem   this checkout
#   /scratch1/$USER/fdcycle-discovery-20261016/deps        DEPS_ROOT (PETSc)
#
# Stages (STAGE=..., several separated by spaces; default "deps mfem seas"):
#   fetch  On the LOGIN node (internet): the PETSc tarball and every package
#          tarball its configure asks for, into ${DEPS_ROOT}/downloads.
#   deps   PETSc from those tarballs, into ${DEPS_ROOT}/petsc-<version>.
#   mfem   Configure and build MFEM (this checkout, in place) against it.
#   seas   cycle/cycle.mk cycle-all: seas_fd_cycle_driver, the cycle tests,
#          the API pins and the owned-paths check (rule M1).
# Run deps, mfem and seas on a condo node (jobs/fd_cycle/discovery_build
# .sbatch).  Each stage appends what it used to
# ${DEPS_ROOT}/build_provenance.txt.
#
# Every stage starts with the GitHub-first guard (CLAUDE.md): a clean tracked
# tree whose HEAD is on GitHub.
#
# Environment overrides:
#   STAGE="deps mfem seas"   DEPS_ROOT=<checkout>/../deps   PETSC_VERSION=3.24.6
#   MODULES="gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 hdf5/1.14.3 cmake/3.29.4 python/3.11.9"
#   BLASLAPACK_LIB="..."     (default: OpenBLAS from the openblas module)
#   OPENBLAS_LIBDIR=<dir>    (default: found from pkg-config or LD_LIBRARY_PATH)
#   CPU_FLAGS="-march=x86-64-v3" (AVX2, as the condo Nitsche build)
#   JOBS=<n> (default nproc)

set -euo pipefail
export LC_ALL=C LANG=C

STAGE="${STAGE:-deps mfem seas}"
PETSC_VERSION="${PETSC_VERSION:-3.24.6}"
MODULES="${MODULES:-gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 hdf5/1.14.3 cmake/3.29.4 python/3.11.9}"
CPU_FLAGS="${CPU_FLAGS:--march=x86-64-v3}"
JOBS="${JOBS:-$(nproc)}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
DEPS_ROOT="${DEPS_ROOT:-$(dirname "${ROOT}")/deps}"
DOWNLOADS="${DEPS_ROOT}/downloads"
PETSC_SRC="${DEPS_ROOT}/src/petsc-${PETSC_VERSION}"
PETSC_PREFIX="${DEPS_ROOT}/petsc-${PETSC_VERSION}"
PROVENANCE="${DEPS_ROOT}/build_provenance.txt"
PETSC_TARBALL="petsc-${PETSC_VERSION}.tar.gz"
PETSC_URLS=(
   "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/${PETSC_TARBALL}"
   "https://gitlab.com/petsc/petsc/-/archive/v${PETSC_VERSION}/petsc-v${PETSC_VERSION}.tar.gz"
)

log() { echo "[build_discovery] $*"; }
fail() { echo "ERROR: $*" >&2; exit 2; }

# ---- GitHub-first guard (CLAUDE.md) ----
cd "${ROOT}"
git fetch -q origin || fail "git fetch failed"
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { git status --short; fail "uncommitted changes in ${ROOT}"; }
[ -n "$(git branch -r --contains HEAD)" ] || fail "HEAD $(git rev-parse --short HEAD) is not on GitHub"
log "commit $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD)), stages: ${STAGE}"

# ---- toolchain: separate module statements (Lmod in a pipeline does nothing) ----
if ! type module >/dev/null 2>&1; then
   [ -f /etc/profile.d/z00_lmod.sh ] && source /etc/profile.d/z00_lmod.sh
   [ -f /usr/share/lmod/lmod/init/bash ] && source /usr/share/lmod/lmod/init/bash
fi
type module >/dev/null 2>&1 || fail "no 'module' command"
set +u   # Lmod's shell functions read unset variables
module purge
for m in ${MODULES}; do module load "${m}" || fail "module load ${m} failed"; done
set -u
modlist() { set +u; module -t list 2>&1 | tr '\n' ' '; set -u; }
log "modules: $(modlist)"
for tool in mpicc mpicxx mpif90 cmake python3 curl; do
   command -v "${tool}" >/dev/null || fail "${tool} not on PATH"
done
python3 -c 'import sys; sys.exit(sys.version_info < (3, 8))' || fail "python3 is $(python3 --version 2>&1); PETSc needs >= 3.8"
[ -n "${HDF5_ROOT:-}" ] && [ -f "${HDF5_ROOT}/include/hdf5.h" ] || fail "HDF5_ROOT not set by the hdf5 module"

BLASLAPACK_LIB="${BLASLAPACK_LIB:-}"
if [ -z "${BLASLAPACK_LIB}" ]; then
   # The openblas module sets no OPENBLAS_ROOT; its libraries are in lib64.
   OB_LIBDIR="${OPENBLAS_LIBDIR:-$(pkg-config --variable=libdir openblas 2>/dev/null || true)}"
   if [ -z "${OB_LIBDIR}" ] || [ ! -e "${OB_LIBDIR}/libopenblas.so" ]; then
      OB_LIBDIR=""
      IFS=: read -r -a dirs <<< "${LD_LIBRARY_PATH:-}"
      for d in "${dirs[@]}"; do
         case "${d}" in *openblas*) [ -e "${d}/libopenblas.so" ] && { OB_LIBDIR="${d}"; break; } ;; esac
      done
   fi
   [ -n "${OB_LIBDIR}" ] || fail "no libopenblas.so from the openblas module; set OPENBLAS_LIBDIR or BLASLAPACK_LIB"
   # PETSc and MUMPS need 32-bit BLAS integers.
   cfg="$(find "$(dirname "${OB_LIBDIR}")/include" -name openblas_config.h 2>/dev/null | head -1)"
   [ -n "${cfg}" ] || fail "no openblas_config.h next to ${OB_LIBDIR}"
   if grep -qE "^#define OPENBLAS_USE64BITINT" "${cfg}"; then
      fail "${OB_LIBDIR} is an ILP64 OpenBLAS (64-bit integers)"
   fi
   BLASLAPACK_LIB="-L${OB_LIBDIR} -lopenblas"
fi
MFEM_LAPACK_LIB="${BLASLAPACK_LIB}"
for flag in ${BLASLAPACK_LIB}; do
   case "${flag}" in -L*) MFEM_LAPACK_LIB="-Wl,-rpath,${flag#-L} ${MFEM_LAPACK_LIB}" ;; esac
done

mkdir -p "${DEPS_ROOT}" "${DOWNLOADS}" "${DEPS_ROOT}/src"
record() {
   {
      echo "# build_discovery.sh stage '$1' on $(hostname) at $(date -u +%FT%TZ)"
      echo "source: $(git -C "${ROOT}" rev-parse HEAD) ($(git -C "${ROOT}" rev-parse --abbrev-ref HEAD))," \
           "tracked changes: $(git -C "${ROOT}" status --porcelain --untracked-files=no | wc -l)"
      echo "modules: $(modlist)"
      echo "mpicc: $(mpicc --version 2>&1 | head -1)"
      echo "mpi: $(ompi_info 2>/dev/null | grep -m1 'Open MPI:' | xargs)"
      echo "blas/lapack: ${BLASLAPACK_LIB}; hdf5: ${HDF5_ROOT}; cpu flags: ${CPU_FLAGS}"
      echo "petsc: ${PETSC_VERSION} at ${PETSC_PREFIX}"
   } >> "${PROVENANCE}"
}

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
   --with-cuda=0
   "--with-blaslapack-lib=${BLASLAPACK_LIB}"
   --download-scalapack
   --download-mumps
   --download-metis
   --download-parmetis
   --download-ptscotch
   --download-hypre
   --with-packages-download-dir="${DOWNLOADS}"
)

unpack_petsc() {
   [ -f "${PETSC_SRC}/configure" ] && return
   [ -f "${DOWNLOADS}/${PETSC_TARBALL}" ] || fail "${DOWNLOADS}/${PETSC_TARBALL} missing; run STAGE=fetch on the login node"
   log "unpacking ${PETSC_TARBALL}"
   rm -rf "${PETSC_SRC}.tmp"
   mkdir -p "${PETSC_SRC}.tmp"
   tar -xzf "${DOWNLOADS}/${PETSC_TARBALL}" -C "${PETSC_SRC}.tmp" --strip-components=1
   mv "${PETSC_SRC}.tmp" "${PETSC_SRC}"
}

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
      [ "${ok}" = 1 ] || fail "could not download ${PETSC_TARBALL}"
   fi
   unpack_petsc
   # With --with-packages-download-dir, PETSc's configure lists the package
   # tarballs it needs and stops.  Each line reads: name ['url1', 'url2', ...]
   log "asking PETSc's configure which package tarballs it needs"
   local listing="${DEPS_ROOT}/petsc_download_list.txt"
   ( cd "${PETSC_SRC}" && ./configure "${PETSC_CONFIGURE_ARGS[@]}" ) > "${listing}" 2>&1 || true
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
   log "downloads in ${DOWNLOADS}:"
   ls -la "${DOWNLOADS}"
   record fetch
}

stage_deps() {
   unpack_petsc
   if [ -f "${PETSC_PREFIX}/lib/petsc/conf/petscvariables" ]; then
      log "PETSc already installed in ${PETSC_PREFIX}; skipping (remove it to rebuild)"
   else
      log "configuring PETSc ${PETSC_VERSION}"
      ( cd "${PETSC_SRC}" && ./configure "${PETSC_CONFIGURE_ARGS[@]}" )
      ( cd "${PETSC_SRC}" && make PETSC_DIR="${PETSC_SRC}" PETSC_ARCH=arch-linux-c-opt MAKE_NP="${JOBS}" all )
      ( cd "${PETSC_SRC}" && make PETSC_DIR="${PETSC_SRC}" PETSC_ARCH=arch-linux-c-opt install )
   fi
   grep -q "^#define PETSC_HAVE_MUMPS 1" "${PETSC_PREFIX}/include/petscconf.h" || fail "PETSc has no MUMPS"
   record deps
   {
      echo "petsc configure options: $(grep -m1 '^CONFIGURE_OPTIONS' "${PETSC_PREFIX}/lib/petsc/conf/petscvariables" | cut -d= -f2-)"
      grep -h 'define PETSC_PKG_MUMPS_VERSION_\(MAJOR\|MINOR\|SUBMINOR\)' "${PETSC_PREFIX}/include/petscconf.h" || true
      grep -h -m1 'HYPRE_RELEASE_VERSION' "${PETSC_PREFIX}/include/HYPRE_config.h" || true
   } >> "${PROVENANCE}"
}

petsc_var() {
   grep -m1 "^$1 *=" "${PETSC_PREFIX}/lib/petsc/conf/petscvariables" | cut -d= -f2-
}

stage_mfem() {
   [ -f "${PETSC_PREFIX}/lib/petsc/conf/petscvariables" ] || fail "no PETSc in ${PETSC_PREFIX}; run STAGE=deps first"
   local P="${PETSC_PREFIX}" ext_libs
   ext_libs="$(petsc_var PETSC_EXTERNAL_LIB_BASIC)"
   log "configuring MFEM in ${ROOT}"
   cd "${ROOT}"
   make config \
      MFEM_USE_MPI=YES MFEM_CXX=mpicxx \
      OPTIM_FLAGS="-O3 -std=c++17 ${CPU_FLAGS}" \
      MFEM_USE_METIS=YES MFEM_USE_METIS_5=YES \
      METIS_DIR="${P}" METIS_OPT="-I${P}/include" \
      METIS_LIB="-Wl,-rpath,${P}/lib -L${P}/lib -lmetis" \
      HYPRE_DIR="${P}" HYPRE_OPT="-I${P}/include" \
      HYPRE_LIB="-Wl,-rpath,${P}/lib -L${P}/lib -lHYPRE" \
      MFEM_USE_LAPACK=YES LAPACK_LIB="${MFEM_LAPACK_LIB}" \
      MFEM_USE_MUMPS=YES MUMPS_DIR="${P}" MUMPS_OPT="-I${P}/include" \
      MUMPS_LIB="-Wl,-rpath,${P}/lib -L${P}/lib -ldmumps -lmumps_common -lpord ${ext_libs}" \
      MFEM_USE_PETSC=YES PETSC_DIR="${P}" PETSC_ARCH= \
      MFEM_USE_HDF5=NO
   make -j"${JOBS}" lib
   [ -f libmfem.a ] || fail "libmfem.a missing after make"
   record mfem
   grep -E '^MFEM_USE_(MPI|PETSC|MUMPS|METIS|LAPACK|HDF5) ' config/config.mk >> "${PROVENANCE}"
}

stage_seas() {
   git -C "${ROOT}" submodule update --init miniapps/seas/extern/toml11 || fail "toml11 submodule"
   [ -f "${ROOT}/miniapps/seas/extern/toml11/toml.hpp" ] || fail "toml11/toml.hpp missing"
   cd "${ROOT}/miniapps/seas"
   # Fresh objects: the shared Makefile has no header dependencies.
   find . -name '*.o' -not -path './extern/*' -delete
   make -f cycle/cycle.mk -j"${JOBS}" cycle-all SEAS_FS_LIB= HDF5_DIR="${HDF5_ROOT}" \
      || fail "cycle-all failed"
   git -C "${ROOT}" rev-parse HEAD > BUILD_COMMIT_DISCOVERY
   record seas
   ls -la seas_fd_cycle_driver seas_test_cycle_p2 seas_test_cycle_p3 >> "${PROVENANCE}"
}

for s in ${STAGE}; do
   case "${s}" in
      fetch) stage_fetch ;;
      deps)  stage_deps ;;
      mfem)  stage_mfem ;;
      seas)  stage_seas ;;
      *) fail "unknown stage '${s}'" ;;
   esac
done
log "done: ${STAGE} at commit $(git -C "${ROOT}" rev-parse HEAD)"
