# Hooks of the S1 jobs (slow-parts study, Chunhui Zhao 2026-10-25), sourced by s1_*.sbatch before
# run_m0_vista.sh, which calls them around each variant with its name and output folder:
#   before_variant  before the first variant: where every rank sits (rank_bindings.sh through that
#                   variant's launcher, into $BASE/bindings.txt; per node the ranks on each Grace CPU,
#                   cores 0-71 and 72-143).  Before perf: the ranks to record (pick_perf_ranks.py on
#                   the timers variant's rank balance) and the event group the node offers: cycles,
#                   instructions and backend memory stalls, else cycles and instructions, else perf's
#                   default
#   after_variant   after perf: perf report per recorded rank (symbols with each event's share and
#                   count), the top lines printed
S1_FIRST=1
before_variant() {
   local v=$1 out=$2 ev
   if [ $S1_FIRST = 1 ]; then
      S1_FIRST=0
      ibrun $LAUNCH jobs/vista_profiling/rank_bindings.sh > $BASE/bindings.txt 2>&1
      echo "--- rank placement through [${LAUNCH:-ibrun alone}] ($BASE/bindings.txt) ---"
      sort -n $BASE/bindings.txt | awk '
         NF == 4 { n[$3]++; if (!($3 in order)) { order[$3] = ++k; name[k] = $3 }
                   if ($4 ~ /^[0-9]+$/) { one[$3]++; if ($4 + 0 < 72) lo[$3]++; else hi[$3]++ }
                   else wide[$3] = $4
                   if ($2 < 8) first[$3] = first[$3] " " $2 ":" $4 }
         END { for (i = 1; i <= k; i++) { h = name[i]
                  printf "  %s: %d ranks, %d bound to one core (%d on cores 0-71, %d on 72-143)%s; local rank:cores%s\n",
                         h, n[h], one[h], lo[h], hi[h], (h in wide) ? ", the others on " wide[h] : "", first[h] } }'
   fi
   if [ "$v" = perf ]; then
      export PERF_DIR=$out/perf
      mkdir -p $PERF_DIR
      PERF_RANKS=$(python3 -I jobs/vista_profiling/pick_perf_ranks.py $BASE/timers/${PROB}_rank_balance.csv 2> $out/perf_pick.txt)
      cat $out/perf_pick.txt
      export PERF_RANKS
      PERF_EVENTS=
      for ev in '{cycles,instructions,stall_backend_mem}' '{cycles,instructions}'; do
         perf stat -e "$ev" -x, -o /dev/null -- true 2> /dev/null && { PERF_EVENTS=$ev; break; }
      done
      export PERF_EVENTS
      echo "perf: ranks [$PERF_RANKS], events [${PERF_EVENTS:-the default of perf}]; $(perf --version 2>&1)"
      perf list 2> /dev/null | grep -i -E 'stall_backend|l2d_cache_refill|ll_cache_miss|mem_access ' | head -8
   fi
}
after_variant() {
   local v=$1 out=$2 f
   [ "$v" = perf ] || return 0
   ls $out/perf/perf.r*.data > /dev/null 2>&1 || { echo "ERROR: no perf data in $out/perf"; rc_all=2; return 0; }
   for f in $out/perf/perf.r*.data; do
      perf report -i $f --stdio --no-children --sort symbol --group --show-total-period --percent-limit 0.5 \
         > ${f%.data}.txt 2> /dev/null
      echo "--- perf $(basename $f .data):$(grep -E '^# (Samples|Event count)' ${f%.data}.txt | tr -d '#' | tr '\n' ';')"
      grep -v '^#' ${f%.data}.txt | grep -v '^ *$' | head -25 | cut -c1-200
   done
}
