#!/bin/bash
# perf fallback of run_map_vista.sh (speed-next, Chunhui Zhao 2026-10-20): ranks named in $PERF_RANKS
# run under `perf record` (sampling at 499 Hz, no call graph) into $PERF_DIR/perf.r<rank>.data; the
# other ranks run the program unchanged.  $PERF_EVENTS, when set (S1, s1_hooks.sh), is perf's -e list.
r=${OMPI_COMM_WORLD_RANK:-${PMIX_RANK:-0}}
for x in $PERF_RANKS; do
   [ "$x" = "$r" ] && exec perf record -q -F 499 ${PERF_EVENTS:+-e "$PERF_EVENTS"} -o "$PERF_DIR/perf.r$r.data" -- "$@"
done
exec "$@"
