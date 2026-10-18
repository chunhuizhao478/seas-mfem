#!/bin/bash
# Build MFEM + the Nitsche dynamic-rupture driver on USC CARC Discovery for the
# SCEC condo nodes (cluster condo, partition scec: Xeon E5-2640 v4 Broadwell,
# x86-64) for the TPV205 / TPV104 benchmarks (Chunhui Zhao, 2026-10-16).  Run it
# INSIDE a batch job on a condo node (jobs/condo_nitsche/build_verify_condo.sbatch).
#
# Toolchain (CARC Spack modules, nothing built from source):
#   gcc/13.3.0  openmpi/5.0.5  hypre/2.31.0 (32-bit ints, no GPU)  metis/5.1.0  hdf5/1.14.3 (parallel)
# MFEM features: MPI + METIS 5 + HDF5 (fault.vtkhdf).  No MUMPS, PETSc, LAPACK,
# Caliper or H5Z-ZFP: the Nitsche driver uses none of them (explicit time
# stepping; the decks set paraview_fault_zfp_tol = 0).  Same targets as
# build_vista.sh.  -march=x86-64-v3 (AVX2) as in the 2026-10-14 Discovery build,
# whose binaries ran on these nodes (condo job 6459052).
#
# Usage (from the repository root):  JOBS=20 bash build_discovery.sh
# Exit status 0 only when libmfem.a AND every target exist after a successful
# make; every module and path is printed.
set -o pipefail
JOBS="${JOBS:-20}"
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
module load gcc/13.3.0
module load openmpi/5.0.5
module load hypre/2.31.0
module load metis/5.1.0
module load hdf5/1.14.3
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"

for v in HYPRE_ROOT METIS_ROOT HDF5_ROOT; do
   [ -n "${!v:-}" ] || fail "$v is not set by the modules"
   echo "  $v = ${!v}"
done
# The seas Makefile links -L$(HDF5_DIR)/lib, so every library must sit in <root>/lib.
[ -e "$HYPRE_ROOT/lib/libHYPRE.so" ]  || fail "libHYPRE.so not in $HYPRE_ROOT/lib"
[ -e "$METIS_ROOT/lib/libmetis.so" ]  || fail "libmetis.so not in $METIS_ROOT/lib"
[ -e "$HDF5_ROOT/lib/libhdf5.so" ]    || fail "libhdf5.so not in $HDF5_ROOT/lib"
[ -e "$HDF5_ROOT/lib/libhdf5_hl.so" ] || fail "libhdf5_hl.so not in $HDF5_ROOT/lib"
[ -f "$HYPRE_ROOT/include/HYPRE.h" ]  || fail "HYPRE.h not in $HYPRE_ROOT/include"
[ -f "$METIS_ROOT/include/metis.h" ]  || fail "metis.h not in $METIS_ROOT/include"
[ -f "$HDF5_ROOT/include/hdf5.h" ]    || fail "hdf5.h not in $HDF5_ROOT/include"
grep -q "define IDXTYPEWIDTH 32" "$METIS_ROOT/include/metis.h" || fail "metis.h is not 32-bit idx_t (MFEM needs IDXTYPEWIDTH 32)"
grep -E "HYPRE_BIGINT|HYPRE_MIXEDINT" "$HYPRE_ROOT/include/HYPRE_config.h" | grep -q "^#define" \
   && fail "this hypre uses 64-bit (big/mixed) ints"
grep -qE "^#define HYPRE_USING_(CUDA|GPU|HIP)" "$HYPRE_ROOT/include/HYPRE_config.h" \
   && fail "this hypre is a GPU build; MFEM would need MFEM_USE_CUDA (linalg/hypre.hpp:56)"
# MFEM's VTKHDF writer needs HDF5 >= 1.14 and parallel I/O; the decks use deflate level 1.
grep -q '^#define H5_VERSION "1\.14\.' "$HDF5_ROOT/include/H5pubconf.h" || fail "HDF5 is not 1.14.x"
grep -q "^#define H5_HAVE_PARALLEL 1" "$HDF5_ROOT/include/H5pubconf.h" || fail "HDF5 is not a parallel build"
grep -q "^#define H5_HAVE_FILTER_DEFLATE 1" "$HDF5_ROOT/include/H5pubconf.h" || fail "HDF5 has no deflate filter (decks set paraview_fault_deflate_level = 1)"
command -v mpicxx >/dev/null || fail "mpicxx not on PATH"
echo "  mpicxx: $(command -v mpicxx); $(mpicxx --version | head -1)"
echo "  cpu: $(grep -m1 'model name' /proc/cpuinfo)"

# --- toml11 (the driver reads the decks through it).  The login-node clone uses
# --recurse-submodules; this is then a no-op that needs no network.
git submodule update --init miniapps/seas/extern/toml11 || fail "toml11 submodule update failed"
[ -f miniapps/seas/extern/toml11/toml.hpp ] || fail "miniapps/seas/extern/toml11/toml.hpp missing"

# --- MFEM
echo "=== configuring MFEM (commit $(git rev-parse HEAD)) ==="
make config \
   MFEM_USE_MPI=YES MFEM_USE_METIS=YES MFEM_USE_METIS_5=YES MFEM_USE_HDF5=YES \
   MPICXX=mpicxx OPTIM_FLAGS="-O3 -std=c++17 -march=x86-64-v3" \
   HYPRE_OPT="-I$HYPRE_ROOT/include" \
   HYPRE_LIB="-L$HYPRE_ROOT/lib -Wl,-rpath,$HYPRE_ROOT/lib -lHYPRE" \
   METIS_OPT="-I$METIS_ROOT/include" \
   METIS_LIB="-L$METIS_ROOT/lib -Wl,-rpath,$METIS_ROOT/lib -lmetis" \
   HDF5_OPT="-I$HDF5_ROOT/include" \
   HDF5_LIB="-L$HDF5_ROOT/lib -Wl,-rpath,$HDF5_ROOT/lib -lhdf5_hl -lhdf5" \
   || fail "make config failed"
grep -E "^MFEM_USE_(MPI|METIS_5|HDF5|MUMPS|PETSC) " config/config.mk
echo "=== building libmfem.a ==="
make -j"$JOBS" lib || fail "MFEM library build failed"
[ -f libmfem.a ] || fail "libmfem.a missing after make"

# --- the driver and the gate tests (fresh objects: the Makefile has no header deps)
cd miniapps/seas || fail "no miniapps/seas"
find . -name '*.o' -not -path './extern/*' -delete
TARGETS="seas_nitsche_dyn_driver seas_test_nitsche_gates seas_test_nitsche_pf seas_test_nitsche_hook_spy \
         seas_test_nitsche_lsw seas_test_nitsche_srw seas_test_nitsche_lsw_inline"
echo "=== building $TARGETS ==="
make -j"$JOBS" MFEM_DIR=../.. HDF5_DIR="$HDF5_ROOT" SEAS_FS_LIB= $TARGETS || fail "SEAS build failed"
for t in $TARGETS; do [ -x "$t" ] || fail "$t missing after make"; done
git rev-parse HEAD > BUILD_COMMIT_DISCOVERY
echo "build OK: commit $(cat BUILD_COMMIT_DISCOVERY), tracked changes: $(git status --porcelain --untracked-files=no | wc -l)"
