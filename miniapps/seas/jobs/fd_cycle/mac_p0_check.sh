#!/bin/bash
# jobs/fd_cycle/mac_p0_check.sh — phase P0 checks on Chunhui Zhao's Mac
# (plan rev 3).  No cluster.  Run from anywhere inside the worktree, in the
# conda env mumps-test:
#   bash miniapps/seas/jobs/fd_cycle/mac_p0_check.sh [step ...]
# Steps, in order, default all of them except lib:
#   lib     rebuild libmfem incrementally (the worktree must be configured)
#   build   seas_bp5_full, seas_fd_cycle_driver, seas_test_cycle_api_pin
#           through cycle/cycle.mk
#   pin     the API pin test (rule M5) and the owned-paths check (rule M1)
#   stub    the driver stub on the BP5 mesh at np = 1 and np = 4
#   bp5     seas_bp5_full with the golden's flags, 50 steps, 1 rank, then the
#           shared regression check against golden_serial_1r_50step.  The
#           Mac's MUMPS differs from the golden's (Frontera job 7643780), so
#           this is a smoke test: the check's numbers show how far the Mac is
#           from the golden; the gate is the Frontera job
#           jobs/fd_cycle/p0_bp5_serial_50step.sbatch.  When MFEM is built
#           without PETSc (the mumps-test env has none), the run drops
#           --petsc-ts and its options file and uses seas_bp5_full's native
#           Tandem-style RK45, and says so: a different stepper from the
#           golden's, so the check is only a rough sanity number then.
#
# On macOS the shared Makefile's std::filesystem probe (SEAS_FS_LIB) fails
# under GNU make >= 4.3, which keeps the backslash of '\#include' in its
# $(shell) call, and falls back to -lstdc++fs, which libc++ does not have.
# The build step passes the Makefile's documented override SEAS_FS_LIB= on
# Darwin (libc++ carries std::filesystem in its main library).
#
# Mac rules (team memory "Mac run lock"):
#   - one run at a time: the lock /Users/chunhuizhao/projects/.seas-mac-run.lock
#     is taken with mkdir (atomic), holds the PID and the thread name, and is
#     released on exit, also after a failure or a guard kill;
#   - right before every launch, refuse if any seas_*, *bp5*, prterun, mpirun
#     or mpiexec process is running;
#   - the memory guard kills only this script's own process tree, by PID,
#     never by name; OMP_NUM_THREADS=1.
# Logs and outputs go to $RESULTS (default <worktree>/../seas-mfem-fdcycle-results/p0).
set -uo pipefail
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
SEAS="$ROOT/miniapps/seas"
RESULTS="${RESULTS:-$(dirname "$ROOT")/seas-mfem-fdcycle-results/p0}"
LOCK=/Users/chunhuizhao/projects/.seas-mac-run.lock
LOCK_WAIT_MAX="${LOCK_WAIT_MAX:-7200}"     # s
MEM_LIMIT_GB="${MEM_LIMIT_GB:-16}"         # this script's tree only
JOBS="${JOBS:-6}"
STAMP="$(date +%Y%m%d-%H%M%S)"
mkdir -p "$RESULTS"
LOG="$RESULTS/mac_p0_check_$STAMP.log"
exec > >(tee -a "$LOG") 2>&1

case "$ROOT" in
   */seas-mfem-spatial-dyn-driver) echo "ERROR: never work in the main checkout $ROOT"; exit 2 ;;
esac

# ---- lock ----
have_lock=0
release_lock()
{
   if [ "$have_lock" = 1 ] && grep -q "pid $$ " "$LOCK/owner" 2>/dev/null; then
      rm -f "$LOCK/owner" && rmdir "$LOCK" && echo "lock released"
   fi
}
trap release_lock EXIT
waited=0
until mkdir "$LOCK" 2>/dev/null; do
   if [ "$waited" -ge "$LOCK_WAIT_MAX" ]; then
      echo "ERROR: lock still held after ${waited} s by: $(cat "$LOCK/owner" 2>/dev/null)"; exit 3
   fi
   [ $((waited % 300)) = 0 ] && echo "waiting for the lock, held by: $(cat "$LOCK/owner" 2>/dev/null)"
   sleep 30; waited=$((waited + 30))
done
have_lock=1
echo "pid $$ thread fully-dynamic-cycle (P0 mac check) since $(date)" > "$LOCK/owner"
echo "lock taken: $(cat "$LOCK/owner")"

# ---- refuse to launch next to another solver run ----
others_running()
{
   ps -axo pid=,comm= | awk '
      { n = split($2, p, "/"); b = p[n] }
      b ~ /^seas_/ || b ~ /bp5/ || b == "prterun" || b == "mpirun" || b == "mpiexec" { print }'
}
check_idle()
{
   local busy; busy="$(others_running)"
   if [ -n "$busy" ]; then
      echo "ERROR: refusing to launch; solver processes are running:"; echo "$busy"; exit 4
   fi
}

# ---- run one command with the memory guard on its own tree ----
tree_pids()
{
   ps -axo pid=,ppid= | awk -v root="$1" '
      { kids[$2] = kids[$2] " " $1 }
      END { q = root; out = root
            while (q != "") { split(q, a, " "); q = ""
               for (i in a) if (a[i] in kids) { out = out kids[a[i]]; q = q kids[a[i]] } }
            print out }'
}
guarded()
{
   local name="$1"; shift
   check_idle
   echo "=== $name: $* ==="
   local t0; t0=$(date +%s)
   "$@" &
   local pid=$! peak=0
   while kill -0 "$pid" 2>/dev/null; do
      local pids rss
      pids="$(tree_pids "$pid")"
      rss=$(ps -o rss= -p "$(echo $pids | tr ' ' ',')" 2>/dev/null | awk '{s += $1} END {print s + 0}')
      [ "$rss" -gt "$peak" ] && peak=$rss
      if [ "$rss" -gt $((MEM_LIMIT_GB * 1024 * 1024)) ]; then
         echo "GUARD: $name uses $((rss / 1024)) MB > ${MEM_LIMIT_GB} GB; killing its own tree: $pids"
         kill -TERM $pids 2>/dev/null; sleep 5; kill -KILL $pids 2>/dev/null
      fi
      sleep 2
   done
   wait "$pid"; local rc=$?
   echo "=== $name: exit $rc after $(( $(date +%s) - t0 )) s, peak RSS $((peak / 1024)) MB ==="
   return $rc
}

steps=("$@"); [ ${#steps[@]} = 0 ] && steps=(build pin stub bp5)
cd "$SEAS"
echo "worktree $ROOT"
echo "commit $(git rev-parse HEAD)  dirty: $(git status --porcelain --untracked-files=no | wc -l | tr -d ' ') files"
echo "conda env ${CONDA_DEFAULT_ENV:-none}; mpirun $(command -v mpirun)"
[ "${CONDA_DEFAULT_ENV:-}" = mumps-test ] || echo "WARNING: not in conda env mumps-test"

MESH=bp5/mesh/reference/bp5_tandem_exact.msh
for s in "${steps[@]}"; do
   case "$s" in
      lib)
         [ -f "$ROOT/config/config.mk" ] || { echo "ERROR: $ROOT is not configured"; exit 1; }
         guarded lib make -C "$ROOT" -j"$JOBS" lib || exit 1 ;;
      build)
         fs_override=(); [ "$(uname -s)" = Darwin ] && fs_override=(SEAS_FS_LIB=)
         guarded build make -f cycle/cycle.mk -j"$JOBS" ${fs_override[@]+"${fs_override[@]}"} seas_bp5_full \
            seas_fd_cycle_driver seas_test_cycle_api_pin || exit 1 ;;
      pin)
         guarded pin ./seas_test_cycle_api_pin || exit 1
         make -f cycle/cycle.mk cycle-owned-check || exit 1 ;;
      stub)
         for np in 1 4; do
            guarded "stub np=$np" mpirun -np "$np" ./seas_fd_cycle_driver --mesh "$MESH" \
               --mesh-scale 1000 --order 1 --fault-attr 3 || exit 1
         done ;;
      bp5)
         out="$RESULTS/bp5_50step_$STAMP"; mkdir -p "$out"
         ts=(--petsc-ts --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg)
         if ! grep -q '^MFEM_USE_PETSC *= YES' "$ROOT/config/config.mk"; then
            ts=()
            echo "NOTE: MFEM has no PETSc here: running seas_bp5_full's native Tandem-style"
            echo "      RK45 instead of the golden's PETSc TS rk 5dp (a different stepper)"
         fi
         guarded "bp5 50 steps" mpirun -np 1 ./seas_bp5_full --mesh "$MESH" --mesh-scale 1000 \
            --order 1 --dg-method IP --tandem-time-stepping --solver mumps \
            --ref-dir bp5/benchmark_data --output-dir "$out" --output-prefix fdc_p0_mac \
            --checkpoint-interval 0 ${ts[@]+"${ts[@]}"} \
            --max-steps 50 --verify > "$out/seas_bp5_full.log" 2>&1
         rc=$?; tail -20 "$out/seas_bp5_full.log"; [ "$rc" = 0 ] || exit "$rc"
         make -f cycle/cycle.mk regression-check SIM_DIR="$out" SIM_PREFIX=fdc_p0_mac \
            REF_DIR=bp5/benchmark_data/golden_serial_1r_50step
         echo "regression-check exit $? (informational on the Mac)" ;;
      *) echo "ERROR: unknown step $s"; exit 2 ;;
   esac
done
echo "mac_p0_check done; log $LOG"
