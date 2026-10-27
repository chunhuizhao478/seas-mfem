#!/bin/bash
# jobs/fd_cycle/build_vista.sh — Vista (TACC) build of the fully dynamic
# cycle branch for the Grace-Grace nodes (partition gg: 2 x 72 Neoverse V2
# cores), Chunhui Zhao, 2026-10-21.  build_discovery.sh, with the Vista
# toolchain of the QD reference runs (build_horizon.sh MODE=cpu on
# gpu-qd-bp5): gcc/13.2.0, openmpi/5.0.11, sequential NVPL BLAS/LAPACK and
# cmake from modules; PETSc 3.24.6 built from source with MUMPS, ScaLAPACK,
# METIS, ParMETIS, PT-Scotch and hypre (Vista's petsc/3.24 module does not
# load with gcc/13.2.0 + openmpi/5.0.11, and its mumps/5.8.0 is sequential).
# The cycle also links HDF5 (the snapshot files): PETSc builds it as well
# (--download-hdf5), so the whole stack comes from one configure and does
# not depend on which HDF5 modules load with openmpi/5.0.11.
#
# Layout: next to the checkout, e.g.
#   $WORK/seas-mfem-fdcycle        this checkout
#   $WORK/seas-mfem-fdcycle-deps   DEPS_ROOT (PETSc, its packages, HDF5)
#
# Stages (STAGE=..., several separated by spaces; default "deps mfem seas"):
#   fetch  On the LOGIN node (internet): the PETSc tarball and every package
#          tarball its configure asks for, into ${DEPS_ROOT}/downloads.
#   deps   PETSc from those tarballs, into ${DEPS_ROOT}/petsc-<version>.
#   mfem   Configure and build MFEM (this checkout, in place) against it.
#   seas   cycle/cycle.mk cycle-all: seas_fd_cycle_driver, the cycle tests,
#          the API pins and the owned-paths check (rule M1).
# Vista login nodes give one core per user: run deps, mfem and seas on a gg
# node (jobs/fd_cycle/vista_build.sbatch).  Each stage appends what it used
# to ${DEPS_ROOT}/build_provenance.txt.
#
# Every stage starts with the GitHub-first guard (CLAUDE.md): a clean tracked
# tree whose HEAD is in a remote-tracking branch.  It does not fetch, as
# compute nodes may have no internet: those refs hold only commits fetched
# from GitHub (clone or pull on the login node).
#
# Environment overrides:
#   STAGE="deps mfem seas"   DEPS_ROOT=<checkout>-deps   PETSC_VERSION=3.24.6
#   MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
#   BLASLAPACK_LIB="..."     (default: sequential NVPL from the nvpl module)
#   CPU_FLAGS="-mcpu=neoverse-v2" (Grace, as the QD build)
#   JOBS=<n> (default nproc)

set -euo pipefail
export LC_ALL=C LANG=C

STAGE="${STAGE:-deps mfem seas}"
PETSC_VERSION="${PETSC_VERSION:-3.24.6}"
MODULES="${MODULES:-gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5}"
CPU_FLAGS="${CPU_FLAGS:--mcpu=neoverse-v2}"
JOBS="${JOBS:-$(nproc)}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
DEPS_ROOT="${DEPS_ROOT:-${ROOT}-deps}"
DOWNLOADS="${DEPS_ROOT}/downloads"
PETSC_SRC="${DEPS_ROOT}/src/petsc-${PETSC_VERSION}"
PETSC_PREFIX="${DEPS_ROOT}/petsc-${PETSC_VERSION}"
PROVENANCE="${DEPS_ROOT}/build_provenance.txt"
PETSC_TARBALL="petsc-${PETSC_VERSION}.tar.gz"
PETSC_URLS=(
   "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/${PETSC_TARBALL}"
   "https://gitlab.com/petsc/petsc/-/archive/v${PETSC_VERSION}/petsc-v${PETSC_VERSION}.tar.gz"
)

log() { echo "[build_vista] $*"; }
fail() { echo "ERROR: $*" >&2; exit 2; }

# ---- GitHub-first guard (CLAUDE.md) ----
cd "${ROOT}"
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { git status --short; fail "uncommitted changes in ${ROOT}"; }
[ -n "$(git branch -r --contains HEAD)" ] || fail "HEAD $(git rev-parse --short HEAD) is not on GitHub"
log "commit $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD)), stages: ${STAGE}"

# ---- toolchain: separate module statements (Lmod in a pipeline does nothing) ----
# No purge, as build_horizon.sh: Lmod swaps the default nvidia compiler for
# gcc and reloads MPI and NVPL for it.
if ! type module >/dev/null 2>&1; then
   [ -f /etc/profile.d/z00_lmod.sh ] && source /etc/profile.d/z00_lmod.sh
   [ -f /usr/share/lmod/lmod/init/bash ] && source /usr/share/lmod/lmod/init/bash
fi
type module >/dev/null 2>&1 || fail "no 'module' command"
set +u   # Lmod's shell functions read unset variables
for m in ${MODULES}; do module load "${m}" || fail "module load ${m} failed"; done
set -u
modlist() { set +u; module -t list 2>&1 | tr '\n' ' '; set -u; }
log "modules: $(modlist)"
for tool in mpicc mpicxx mpif90 cmake python3 curl; do
   command -v "${tool}" >/dev/null || fail "${tool} not on PATH"
done
python3 -c 'import sys; sys.exit(sys.version_info < (3, 8))' || fail "python3 is $(python3 --version 2>&1); PETSc needs >= 3.8"

BLASLAPACK_LIB="${BLASLAPACK_LIB:-}"
if [ -z "${BLASLAPACK_LIB}" ]; then
   [ -n "${TACC_NVPL_LIB:-}" ] || fail "TACC_NVPL_LIB not set (nvpl module) and BLASLAPACK_LIB empty"
   # Sequential NVPL with 32-bit integers (PETSc and MUMPS need lp64): one
   # MPI rank per core, no BLAS threads.
   BLASLAPACK_LIB="-L${TACC_NVPL_LIB} -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core"
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
      echo "blas/lapack: ${BLASLAPACK_LIB}; hdf5: PETSc's, in ${PETSC_PREFIX}; cpu flags: ${CPU_FLAGS}"
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
   --download-hdf5
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
   # tarballs missing from that directory and stops.  Each line reads:
   # name ['url1', 'url2', ...].  When none is missing it goes on to build
   # them instead (on the login node), so, as build_horizon.sh, it is asked
   # with an empty directory (it then lists them all), a prefix inside it
   # and a PETSC_ARCH of its own, which leaves DOWNLOADS and the install alone.
   log "asking PETSc's configure which package tarballs it needs"
   local listing="${DEPS_ROOT}/petsc_download_list.txt"
   local probe probe_arch=arch-fetch-probe probe_args=() a
   probe="$(mktemp -d "${DEPS_ROOT}/fetch-probe.XXXXXX")"
   for a in "${PETSC_CONFIGURE_ARGS[@]}"; do
      case "${a}" in
         --prefix=*) probe_args+=(--prefix="${probe}/prefix") ;;
         --with-packages-download-dir=*) probe_args+=(--with-packages-download-dir="${probe}") ;;
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
      ( cd "${PETSC_SRC}" && ./configure PETSC_ARCH=arch-linux-c-opt "${PETSC_CONFIGURE_ARGS[@]}" )
      ( cd "${PETSC_SRC}" && make PETSC_DIR="${PETSC_SRC}" PETSC_ARCH=arch-linux-c-opt MAKE_NP="${JOBS}" all )
      ( cd "${PETSC_SRC}" && make PETSC_DIR="${PETSC_SRC}" PETSC_ARCH=arch-linux-c-opt install )
   fi
   grep -q "^#define PETSC_HAVE_MUMPS 1" "${PETSC_PREFIX}/include/petscconf.h" || fail "PETSc has no MUMPS"
   grep -q "^#define PETSC_HAVE_HDF5 1" "${PETSC_PREFIX}/include/petscconf.h" || fail "PETSc has no HDF5"
   [ -f "${PETSC_PREFIX}/include/hdf5.h" ] && [ -e "${PETSC_PREFIX}/lib/libhdf5.so" ] \
      || fail "no hdf5.h or libhdf5.so in ${PETSC_PREFIX}"
   record deps
   {
      echo "petsc configure options: $(grep -m1 '^CONFIGURE_OPTIONS' "${PETSC_PREFIX}/lib/petsc/conf/petscvariables" | cut -d= -f2-)"
      grep -h 'define PETSC_PKG_MUMPS_VERSION_\(MAJOR\|MINOR\|SUBMINOR\)' "${PETSC_PREFIX}/include/petscconf.h" || true
      grep -h -m1 'HYPRE_RELEASE_VERSION' "${PETSC_PREFIX}/include/HYPRE_config.h" || true
      grep -h -m1 'define H5_VERSION' "${PETSC_PREFIX}/include/H5pubconf.h" || true
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
   # The submodule comes with the clone (--recurse-submodules on the login
   # node); a compute node may not reach GitHub to fetch it.
   case "$(git -C "${ROOT}" submodule status miniapps/seas/extern/toml11)" in
      " "*) ;;
      *) fail "toml11 submodule missing or not at its recorded commit; on the login node: git submodule update --init miniapps/seas/extern/toml11" ;;
   esac
   [ -f "${ROOT}/miniapps/seas/extern/toml11/toml.hpp" ] || fail "toml11/toml.hpp missing"
   cd "${ROOT}/miniapps/seas"
   # Fresh objects: the shared Makefile has no header dependencies.
   find . -name '*.o' -not -path './extern/*' -delete
   make -f cycle/cycle.mk -j"${JOBS}" cycle-all SEAS_FS_LIB= HDF5_DIR="${PETSC_PREFIX}" \
      || fail "cycle-all failed"
   git -C "${ROOT}" rev-parse HEAD > BUILD_COMMIT_VISTA
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
