#!/bin/bash
# One line per $1 seconds (default 5) until killed: host, time (s since
# 1970), then CommitLimit, Committed_AS and MemAvailable in GiB.  Vista's
# compute nodes refuse an allocation past CommitLimit (vm.overcommit_memory
# = 2), so Committed_AS, not RSS, says how close a run is to failing.  The
# Vista jobs start it on every node of their allocation in the background:
#   srun --overlap -N "$SLURM_NNODES" -n "$SLURM_NNODES" --ntasks-per-node=1 \
#      jobs/fd_cycle/vista_memlog.sh 5 > memlog.txt &
h=$(hostname -s)
while :; do
   awk -v h="$h" -v t="$(date +%s)" '
      /^(CommitLimit|Committed_AS|MemAvailable):/ { v[$1] = $2 / 1048576 }
      END { printf "%s %s limit %.1f committed %.1f available %.1f\n", h, t,
                   v["CommitLimit:"], v["Committed_AS:"], v["MemAvailable:"]; fflush() }' /proc/meminfo
   sleep "${1:-5}"
done
