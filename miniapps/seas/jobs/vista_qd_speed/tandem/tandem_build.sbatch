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
# CSV probes), Tandem's default ageing law.  Before building it prints what
# it found (compilers, PETSc version and memalign, Eigen, Lua, ParMETIS,
# numpy) and stops with the first missing dependency named.  It never
# downloads or installs anything.

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
mkdir -p "${BUILD}"
cd "${BUILD}" || exit 1
cmake "${TANDEM_SRC}" \
   -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
   -DCMAKE_BUILD_TYPE=Release \
   -DDOMAIN_DIMENSION=3 -DPOLYNOMIAL_DEGREE="${DEGREE}" \
   -DARCH=noarch -DENABLE_HDF5=OFF \
   -DPARMETIS_DIR="${TANDEM_PETSC}" -DMETIS_DIR="${TANDEM_PETSC}" \
   > cmake.log 2>&1 \
   || { echo "ERROR: cmake failed; tail of ${BUILD}/cmake.log:" >&2; tail -n 40 cmake.log >&2; exit 1; }
grep -iE 'petsc|lua|eigen|parmetis|arch' cmake.log | head -n 30

echo "=== build $(date -u +%T) ==="
make -j 72 tandem > make.log 2>&1 \
   || { echo "ERROR: make failed; tail of ${BUILD}/make.log:" >&2; tail -n 60 make.log >&2; exit 1; }
echo "=== done $(date -u +%T) ==="
BIN="$(tandem_bin "${DEGREE}")"
[ -x "${BIN}" ] || { echo "ERROR: ${BIN} missing after make" >&2; exit 1; }
{
   echo "job ${SLURM_JOB_ID:-none} on $(hostname), $(date -u +%FT%TZ)"
   echo "seas-mfem $(git -C "${SEAS_MFEM_ROOT}" rev-parse HEAD)"
   echo "tandem ${TANDEM_COMMIT}, degree ${DEGREE}, 3D, Release, ARCH=noarch, HDF5 off"
   echo "petsc ${TANDEM_PETSC}"
   echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
   echo "binary sha256 $(sha256sum "${BIN}" | cut -c1-16)"
} > "${BUILD}/provenance.txt"
cat "${BUILD}/provenance.txt"
