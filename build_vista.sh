#!/bin/bash
# Build MFEM + the Nitsche dynamic-rupture driver on TACC Vista (Grace CPU,
# aarch64, partition gg) for the TPV205 / TPV104 benchmarks (Chunhui Zhao,
# 2026-10-16).  Run it INSIDE a batch job (jobs/vista_nitsche/build_verify_vista.sbatch),
# never on a login node: Vista login nodes cap processes and threads.
#
# Toolchain (TACC modules, no Spack, nothing built from source):
#   gcc/14.2.0  cuda/12.6 (hierarchy only)  openmpi/5.0.5  hypre/2.31.0 (32-bit ints, no GPU)  metis/5.1.0.3  phdf5/1.14.6
# MFEM features: MPI + METIS 5 + HDF5 (fault.vtkhdf).  No MUMPS, PETSc, LAPACK,
# Caliper or H5Z-ZFP: the Nitsche driver uses none of them (explicit time
# stepping; the decks set paraview_fault_zfp_tol = 0).
#
# Usage (from the repository root):  JOBS=64 bash build_vista.sh
# Exit status 0 only when libmfem.a AND seas_nitsche_dyn_driver exist after a
# successful make; every module and path is printed.
set -o pipefail
JOBS="${JOBS:-32}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT" || exit 2
fail () { echo "ERROR: $*"; exit 2; }

# --- modules: separate statements, never a pipeline (Lmod in a subshell does nothing)
if ! type module >/dev/null 2>&1; then
   [ -f /etc/profile.d/z00_lmod.sh ] && source /etc/profile.d/z00_lmod.sh
   [ -f /usr/share/lmod/lmod/init/bash ] && source /usr/share/lmod/lmod/init/bash
fi
type module >/dev/null 2>&1 || fail "no 'module' command"
module purge
module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
# cuda/12.6: Lmod hierarchy only. TACC lists hypre/2.31.0 under gcc + cuda + openmpi, but its HYPRE_config.h
# defines no HYPRE_USING_CUDA/GPU (checked on Vista 2026-10-16). TACC's hypre, metis, phdf5 and gmsh are built with
# gcc 14 (/home1/apps/gcc14/...); gmsh under gcc/13.2.0 fails on CXXABI_1.3.15, so the whole stack uses gcc/14.2.0.
module load cuda/12.6
module load openmpi/5.0.5
module load hypre/2.31.0
module load metis/5.1.0.3
module load phdf5/1.14.6
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"

HDF5_DIR_V="${TACC_PHDF5_DIR:-${TACC_HDF5_DIR:-}}"
HDF5_INC_V="${TACC_PHDF5_INC:-${TACC_HDF5_INC:-$HDF5_DIR_V/include}}"
HDF5_LIB_V="${TACC_PHDF5_LIB:-${TACC_HDF5_LIB:-$HDF5_DIR_V/lib}}"
for v in TACC_HYPRE_INC TACC_HYPRE_LIB TACC_METIS_INC TACC_METIS_LIB; do
   [ -n "${!v:-}" ] || fail "$v is not set by the modules"
   echo "  $v = ${!v}"
done
echo "  HDF5 dir/inc/lib = $HDF5_DIR_V | $HDF5_INC_V | $HDF5_LIB_V"
[ -f "$TACC_HYPRE_INC/HYPRE.h" ]  || fail "HYPRE.h not in $TACC_HYPRE_INC"
[ -f "$TACC_METIS_INC/metis.h" ]  || fail "metis.h not in $TACC_METIS_INC"
[ -f "$HDF5_INC_V/hdf5.h" ]       || fail "hdf5.h not in $HDF5_INC_V"
grep -q "define IDXTYPEWIDTH 32" "$TACC_METIS_INC/metis.h" || fail "metis.h is not 32-bit idx_t (MFEM needs IDXTYPEWIDTH 32)"
grep -E "HYPRE_BIGINT|HYPRE_MIXEDINT" "$TACC_HYPRE_INC/HYPRE_config.h" | grep -q "^#define" \
   && fail "this hypre uses 64-bit (big/mixed) ints; load the 32-bit hypre/2.31.0 module"
grep -qE "^#define HYPRE_USING_(CUDA|GPU|HIP)" "$TACC_HYPRE_INC/HYPRE_config.h" \
   && fail "this hypre is a GPU build; MFEM would need MFEM_USE_CUDA (linalg/hypre.hpp:56)"
command -v mpicxx >/dev/null || fail "mpicxx not on PATH"
echo "  mpicxx: $(command -v mpicxx); $(mpicxx --version | head -1)"

# --- toml11 (the driver reads the decks through it)
git submodule update --init miniapps/seas/extern/toml11 || fail "toml11 submodule update failed"
[ -f miniapps/seas/extern/toml11/toml.hpp ] || fail "miniapps/seas/extern/toml11/toml.hpp missing"

# --- MFEM
echo "=== configuring MFEM (commit $(git rev-parse HEAD)) ==="
make config \
   MFEM_USE_MPI=YES MFEM_USE_METIS=YES MFEM_USE_METIS_5=YES MFEM_USE_HDF5=YES \
   MPICXX=mpicxx OPTIM_FLAGS="-O3 -std=c++17 -mcpu=neoverse-v2" \
   HYPRE_OPT="-I$TACC_HYPRE_INC" \
   HYPRE_LIB="-L$TACC_HYPRE_LIB -Wl,-rpath,$TACC_HYPRE_LIB -lHYPRE" \
   METIS_OPT="-I$TACC_METIS_INC" \
   METIS_LIB="-L$TACC_METIS_LIB -Wl,-rpath,$TACC_METIS_LIB -lmetis" \
   HDF5_OPT="-I$HDF5_INC_V" \
   HDF5_LIB="-L$HDF5_LIB_V -Wl,-rpath,$HDF5_LIB_V -lhdf5_hl -lhdf5 -lz" \
   || fail "make config failed"
grep -E "^MFEM_USE_(MPI|METIS_5|HDF5|MUMPS|PETSC) " config/config.mk
echo "=== building libmfem.a ==="
make -j"$JOBS" lib || fail "MFEM library build failed"
[ -f libmfem.a ] || fail "libmfem.a missing after make"

# --- the driver and the gate tests (fresh objects: the Makefile has no header deps)
cd miniapps/seas || fail "no miniapps/seas"
find . -name '*.o' -not -path './extern/*' -delete
TARGETS="seas_nitsche_dyn_driver seas_test_nitsche_gates seas_test_nitsche_p5 seas_test_nitsche_pf \
         seas_test_nitsche_hook_spy seas_test_nitsche_lsw seas_test_nitsche_srw"
echo "=== building $TARGETS ==="
make -j"$JOBS" MFEM_DIR=../.. HDF5_DIR="$HDF5_DIR_V" SEAS_FS_LIB= $TARGETS || fail "SEAS build failed"
for t in $TARGETS; do [ -x "$t" ] || fail "$t missing after make"; done
git rev-parse HEAD > BUILD_COMMIT_VISTA
echo "build OK: commit $(cat BUILD_COMMIT_VISTA), tracked changes: $(git status --porcelain --untracked-files=no | wc -l)"
