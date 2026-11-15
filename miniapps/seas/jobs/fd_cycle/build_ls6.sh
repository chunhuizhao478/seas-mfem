#!/bin/bash
# jobs/fd_cycle/build_ls6.sh — Lonestar6 (TACC) build of the cycle driver's
# dependencies for the AMD nodes (partition normal: 2 x 64 EPYC 7763 cores,
# Zen 3), Chunhui Zhao, 2026-11-15: "let's try on lonestar6".
#
# build_vista.sh does the work (same stages, same PETSc 3.24.6 from source
# with MUMPS, ScaLAPACK, METIS, ParMETIS, PT-Scotch, hypre and HDF5, same
# GitHub-first guard); this file only gives it the Lonestar6 toolchain
# through the overrides build_vista.sh already takes:
#   MODULES         gcc/13.2.0 impi/21.18 mkl/26.1 cmake/3.31.5 python/3.12.11
#                   (Lonestar6 has no OpenMPI and no NVPL; mvapich2 exists
#                   only for gcc/11.2.0.  Checked on login2, 2026-11-15:
#                   mpicxx -O3 -march=znver3 -std=c++17 compiles; mpif90
#                   compiles "use mpi" and include 'mpif.h'; dgesv_ links and
#                   runs with the MKL line below)
#   BLASLAPACK_LIB  sequential MKL with 32-bit integers and the GNU Fortran
#                   interface (mkl_gf_lp64, mkl_sequential, mkl_core): one
#                   MPI rank per core, no BLAS threads
#   CPU_FLAGS       -march=znver3
# Usage (as build_vista.sh; its header describes the stages and the layout):
#   login node (internet):  STAGE=fetch bash miniapps/seas/jobs/fd_cycle/build_ls6.sh
#   compute node:           STAGE="deps mfem" bash .../build_ls6.sh
#                           (jobs/2d_inplane/ls6_build.sbatch does this)
# Log lines read "[build_vista]", and the provenance file's stage headers
# name build_discovery.sh: both come from build_vista.sh, unchanged.  Its
# "seas" stage (cycle-all, BUILD_COMMIT_VISTA) is not for this branch:
# ls6_build.sbatch builds the driver and writes BUILD_COMMIT_LS6.
set -euo pipefail
export LC_ALL=C LANG=C
export MODULES="${MODULES:-gcc/13.2.0 impi/21.18 mkl/26.1 cmake/3.31.5 python/3.12.11}"
export CPU_FLAGS="${CPU_FLAGS:--march=znver3}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

fail() { echo "ERROR: $*" >&2; exit 2; }
if [ -z "${BLASLAPACK_LIB:-}" ]; then
   # MKLROOT comes from the mkl module: load the toolchain here as well
   # (build_vista.sh loads the same list again, which changes nothing).
   # Separate module statements: Lmod in a pipeline does nothing.
   if ! type module >/dev/null 2>&1; then
      [ -f /etc/profile.d/z00_lmod.sh ] && source /etc/profile.d/z00_lmod.sh
      [ -f /usr/share/lmod/lmod/init/bash ] && source /usr/share/lmod/lmod/init/bash
   fi
   type module >/dev/null 2>&1 || fail "no 'module' command"
   set +u   # Lmod's shell functions read unset variables
   for m in ${MODULES}; do module load "${m}" || fail "module load ${m} failed"; done
   set -u
   [ -n "${MKLROOT:-}" ] || fail "MKLROOT not set (mkl module) and BLASLAPACK_LIB empty"
   MKL_LIBDIR="${MKLROOT}/lib/intel64"
   [ -e "${MKL_LIBDIR}/libmkl_gf_lp64.so" ] || MKL_LIBDIR="${MKLROOT}/lib"
   for l in mkl_gf_lp64 mkl_sequential mkl_core; do
      [ -e "${MKL_LIBDIR}/lib${l}.so" ] || fail "no lib${l}.so in ${MKLROOT}/lib/intel64 or ${MKLROOT}/lib"
   done
   export BLASLAPACK_LIB="-L${MKL_LIBDIR} -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
fi
echo "[build_ls6] modules ${MODULES}; cpu flags ${CPU_FLAGS}; blas/lapack ${BLASLAPACK_LIB}"
echo "[build_ls6] mpi: $(mpicc -v 2>&1 | grep -m1 -i 'mpi' || true)"
exec bash "${HERE}/build_vista.sh" "$@"
