#!/bin/bash
# jobs/fd_cycle/build_frontera.sh — guarded Frontera build of the fully
# dynamic cycle branch (plan rev 3, rules M2 and M10; CLAUDE.md at 010b56e).
#
# Usage, in a Frontera checkout of fully-dynamic-cycle (default
# $WORK/seas-mfem-fdcycle), on a login node or inside a job:
#   bash miniapps/seas/jobs/fd_cycle/build_frontera.sh [extra cycle.mk targets]
#
# 1. Refuses the Frontera checkouts that hold uncommitted changes
#    (seas-mfem-nitsche-damp, -kv, -pc, -fixsn).
# 2. The GitHub-first guard: fetch, clean tracked tree, HEAD on GitHub,
#    commit printed.
# 3. The module stack of the golden 50-step job (jobs/bp5/bp5_verify_serial
#    .sbatch) plus fftw3, as the Nitsche jobs load it.
# 4. MFEM: if the checkout is not configured yet (no config/config.mk with
#    PETSc and MUMPS, or no extern/hdf5 install), runs the root
#    build_frontera.sh once (QUICK=1); otherwise rebuilds libmfem
#    incrementally.  Either way the library is built from this commit.
# 5. seas_bp5_full, seas_fd_cycle_driver and seas_test_cycle_api_pin through
#    cycle/cycle.mk (never the shared Makefile's own build scripts), then the
#    API pin test and the owned-paths check.
set -uo pipefail
export LC_ALL=C LANG=C
JOBS="${JOBS:-8}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
case "$(basename "$ROOT")" in
   *-damp|*-kv|*-pc|*-fixsn)
      echo "ERROR: $ROOT holds uncommitted changes (version audit 2026-10-15); never build there"
      exit 2 ;;
esac

# ---- GitHub-first guard (CLAUDE.md) ----
cd "$ROOT"
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "commit $(git rev-parse HEAD)"
echo "branch $(git rev-parse --abbrev-ref HEAD)"

# ---- modules (bp5_verify_serial.sbatch + fftw3) ----
module load intel/19.1.1 impi/19.0.9 hypre/2.31.0 mumps/5.3 parmetis petsc/3.15 fftw3/3.3.8 \
   || { echo "ERROR: module load failed"; exit 2; }
export LD_LIBRARY_PATH="$ROOT/extern/hdf5/install/lib:$ROOT/extern/h5z-zfp/install/lib:${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"

# ---- MFEM library from this commit ----
configured=1
grep -q '^MFEM_USE_PETSC *= YES' config/config.mk 2>/dev/null || configured=0
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk 2>/dev/null || configured=0
[ -d extern/hdf5/install/lib ] || configured=0
if [ "$configured" = 0 ]; then
   echo "=== First build of this checkout: root build_frontera.sh (QUICK=1) ==="
   QUICK=1 JOBS="$JOBS" bash build_frontera.sh || { echo "ERROR: root build_frontera.sh failed"; exit 1; }
else
   echo "=== libmfem (incremental) ==="
   make -j"$JOBS" lib || { echo "ERROR: libmfem build failed"; exit 1; }
fi
grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

# ---- seas targets through cycle.mk ----
git submodule update --init miniapps/seas/extern/toml11 || { echo "ERROR: toml11 submodule"; exit 1; }
cd miniapps/seas
make -f cycle/cycle.mk -j"$JOBS" seas_bp5_full seas_fd_cycle_driver seas_test_cycle_api_pin "$@" \
   || { echo "BUILD FAILED"; exit 1; }
make -f cycle/cycle.mk test-cycle-api-pin cycle-owned-check || { echo "CHECK FAILED"; exit 1; }
echo "=== build complete at commit $(git rev-parse HEAD) ==="
