#!/bin/bash
# locked_exec.sh <cmd...>: run <cmd> under the Mac run lock with its output
# passed straight through; wall time and the peak RSS summed over its process
# tree go to $STATS (one line per run), never into the command's output.
set -o pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LOCK=/Users/chunhuizhao/projects/.seas-mac-run.lock
STATS=${STATS:-/Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/g3-bp3/run_stats.log}
have=0
release() { [ "$have" = 1 ] && grep -q "pid $$ " "$LOCK/owner" 2>/dev/null && rm -f "$LOCK/owner" && rmdir "$LOCK"; }
trap release EXIT
# Take the lock; if another solver process runs (e.g. a run that does not
# use the lock), release the lock, wait and try again: never launch next to it.
while :; do
   until mkdir "$LOCK" 2>/dev/null; do echo "waiting for lock: $(cat "$LOCK/owner" 2>/dev/null)" >&2; sleep 30; done
   have=1
   echo "pid $$ thread 2d-inplane (g3-bp3) since $(date)" > "$LOCK/owner"
   busy="$(ps -axo pid=,comm= | awk '{ n = split($2, p, "/"); b = p[n] }
      b ~ /^seas_/ || b ~ /bp5/ || b == "prterun" || b == "mpirun" || b == "mpiexec" { print }')"
   [ -z "$busy" ] && break
   echo "solver processes running, not launching; released the lock, retrying in 60 s:" >&2; echo "$busy" >&2
   echo "$(date '+%F %T') WAIT-BUSY $(echo $busy) :: $*" >> "$STATS"
   release; have=0; sleep 60
done
tree_pids() {
   ps -axo pid=,ppid= | awk -v root="$1" '{ kids[$2] = kids[$2] " " $1 }
      END { q = root; out = root
            while (q != "") { split(q, a, " "); q = ""
               for (i in a) if (a[i] in kids) { out = out kids[a[i]]; q = q kids[a[i]] } }
            print out }'
}
t0=$(date +%s)
"$@" &
pid=$! peak=0
while kill -0 $pid 2>/dev/null; do
   rss=$(ps -o rss= -p "$(tree_pids $pid | tr ' ' ',')" 2>/dev/null | awk '{s += $1} END {print s + 0}')
   [ "$rss" -gt "$peak" ] && peak=$rss
   sleep 1
done
wait $pid; rc=$?
echo "$(date '+%F %T') wall $(( $(date +%s) - t0 )) s peak_tree_rss $((peak / 1024)) MB rc $rc :: $*" >> "$STATS"
exit $rc
