#!/bin/bash
#SBATCH -J tandem_build
#SBATCH -o tandem_build_%j.out
#SBATCH -e tandem_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 144
#SBATCH -t 01:00:00
#SBATCH -A EAR20006

# Tandem vs MFEM-SEAS speed test: build upstream Tandem (TANDEM_COMMIT in
# tandem_env.sh) for 3D at polynomial degree DEGREE (1, 2 or 3; Tandem fixes
# the degree at compile time) in ${TANDEM_ROOT}/build-d3-p<DEGREE>.
#
#   sbatch --export=ALL,DEGREE=1 jobs/vista_qd_speed/tandem/tandem_build.sbatch
#   (from miniapps/seas of the qd-speed checkout)
#
# Settings: Release, ARCH=noarch (Tandem has no aarch64 target; noarch
# compiles with -march=native and 8-byte alignment), HDF5 off (output is
# CSV probes), Tandem's default ageing law, mpif90 and nvpl LAPACK (named
# library by library, checked on the binary).  Needs tandem_petsc.sbatch
# first.  Before building it prints what it found (compilers, PETSc version
# and memalign, Eigen, Lua, ParMETIS, numpy) and stops with the first
# missing dependency named.  It never downloads or installs anything.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista_qd_speed/tandem/tandem_env.sh"
seas_guard
tandem_guard

DEGREE="${DEGREE:-}"
[[ "${DEGREE}" =~ ^[1-3]$ ]] || { echo "ERROR: DEGREE='${DEGREE}' must be 1, 2 or 3." >&2; exit 1; }
BUILD="${TANDEM_ROOT}/build-d3-p${DEGREE}"

echo "=== dependencies ==="
module -t list 2>&1 | tr '\n' ' '; echo
which mpicc mpicxx cmake python3
python3 -c 'import numpy; print("numpy", numpy.__version__)' \
   || { echo "ERROR: python3 has no numpy (yateto needs it)" >&2; exit 1; }
[ -f "${TANDEM_PETSC}/lib/pkgconfig/PETSc.pc" ] || [ -f "${TANDEM_PETSC}/lib/pkgconfig/petsc.pc" ] \
   || { echo "ERROR: no PETSc pkg-config file under ${TANDEM_PETSC}/lib/pkgconfig" >&2; exit 1; }
grep -E '#define PETSC_VERSION_(MAJOR|MINOR|SUBMINOR) |#define PETSC_MEMALIGN ' \
   "${TANDEM_PETSC}/include/petscversion.h" "${TANDEM_PETSC}/include/petscconf.h"
ls "${TANDEM_PETSC}"/lib/libparmetis* "${TANDEM_PETSC}"/lib/libmetis* \
   || { echo "ERROR: ParMETIS/METIS not in ${TANDEM_PETSC}/lib" >&2; exit 1; }
echo "lua modules:"; module spider lua 2>&1 | grep -iE "lua/|liblua" | head -n 5

echo "=== configure ${BUILD} ==="
# Always a fresh build folder: CMake keeps the compilers and libraries of an
# earlier configure in its cache.
case "${BUILD}" in "${TANDEM_ROOT}"/build-d3-p[1-3]) rm -rf "${BUILD}" ;; esac
mkdir -p "${BUILD}"
cd "${BUILD}" || exit 1
# nvpl named library by library: CMake's BLA_VENDOR=NVPL search finds none
# on Vista (job 1114529).  The same four libraries our PETSc links
# (build_horizon.sh BLASLAPACK_LIB).
[ -n "${TACC_NVPL_LIB:-}" ] || { echo "ERROR: TACC_NVPL_LIB unset (nvpl module)" >&2; exit 1; }
NVPL_LIBS=""
for l in nvpl_lapack_lp64_seq nvpl_lapack_core nvpl_blas_lp64_seq nvpl_blas_core; do
   [ -f "${TACC_NVPL_LIB}/lib${l}.so" ] || { echo "ERROR: no ${TACC_NVPL_LIB}/lib${l}.so" >&2; exit 1; }
   NVPL_LIBS="${NVPL_LIBS:+${NVPL_LIBS};}${TACC_NVPL_LIB}/lib${l}.so"
done
echo "nvpl: ${NVPL_LIBS}"
cmake "${TANDEM_SRC}" \
   -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
   -DCMAKE_BUILD_TYPE=Release \
   -DDOMAIN_DIMENSION=3 -DPOLYNOMIAL_DEGREE="${DEGREE}" \
   -DARCH=noarch -DENABLE_HDF5=OFF \
   -DPARMETIS_DIR="${TANDEM_PETSC}" -DMETIS_DIR="${TANDEM_PETSC}" \
   -DPython3_EXECUTABLE="$(command -v python3)" \
   -DCMAKE_Fortran_COMPILER=mpif90 \
   -DBLAS_LIBRARIES="${NVPL_LIBS}" -DLAPACK_LIBRARIES="${NVPL_LIBS}" \
   > cmake.log 2>&1 \
   || { echo "ERROR: cmake failed; tail of ${BUILD}/cmake.log:" >&2; tail -n 40 cmake.log >&2; exit 1; }
grep -iE 'petsc|lua|eigen|parmetis|arch|fortran|lapack|blas' cmake.log | head -n 40
# Tandem's LAPACK must be the nvpl our PETSc and MFEM use (job 1113226
# found the system OpenBLAS and gfortran 11 instead).
grep -q 'Found LAPACK' cmake.log \
   || { echo "ERROR: cmake did not find LAPACK (nvpl); see ${BUILD}/cmake.log" >&2; exit 1; }
grep -E '^CMAKE_Fortran_COMPILER:' CMakeCache.txt

echo "=== build $(date -u +%T) ==="
make -j 72 tandem > make.log 2>&1 \
   || { echo "ERROR: make failed; tail of ${BUILD}/make.log:" >&2; tail -n 60 make.log >&2; exit 1; }
echo "=== done $(date -u +%T) ==="
BIN="$(tandem_bin "${DEGREE}")"
[ -x "${BIN}" ] || { echo "ERROR: ${BIN} missing after make" >&2; exit 1; }
if ldd "${BIN}" | grep -qi openblas || ! ldd "${BIN}" | grep -q nvpl_lapack; then
   echo "ERROR: ${BIN} does not link nvpl LAPACK only:" >&2
   ldd "${BIN}" | grep -iE 'blas|lapack' >&2
   exit 1
fi
ldd "${BIN}" | grep -E 'petsc|mpi|nvpl|lua|parmetis'
{
   echo "job ${SLURM_JOB_ID:-none} on $(hostname), $(date -u +%FT%TZ)"
   echo "seas-mfem $(git -C "${SEAS_MFEM_ROOT}" rev-parse HEAD)"
   echo "tandem ${TANDEM_COMMIT}, degree ${DEGREE}, 3D, Release, ARCH=noarch, HDF5 off"
   echo "petsc ${TANDEM_PETSC}"
   echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
   echo "binary sha256 $(sha256sum "${BIN}" | cut -c1-16)"
} > "${BUILD}/provenance.txt"
cat "${BUILD}/provenance.txt"
