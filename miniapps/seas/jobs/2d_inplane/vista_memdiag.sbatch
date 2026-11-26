#!/bin/bash
#SBATCH -J memdiag
#SBATCH -o memdiag_%j.out
#SBATCH -e memdiag_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 00:50:00
#SBATCH -A EAR20006
# Memory diagnosis of the H5 +V run (Chunhui Zhao, 2026-11-25: "can we debug
# why there is a memory issue?").  Vista jobs 1167994 and 1170201 (deck
# ..._cvmv_100ev_vol_fd.toml, 144 ranks) both died about 3 and 30 min into the
# QD stretch after event 1: the step fell to 2e-5 s, node Committed_AS grew
# from 182 to 226 GB (about 2 GB per minute), and rank 0 was killed
# (signal 9, a UCX xpmem attach error just before).
# This job:
#   1. copies the run folder $1 (with its checkpoint at the event-1
#      hand-back) to $1_memdiag_<job id>, so the original stays untouched;
#   2. runs vista_2d_cycle.sbatch on the copy (a restart from that
#      checkpoint) with PETSC_OPTIONS=-ts_adapt_monitor (every step's
#      accept / reject and its dt reason) and a 40 min wall limit;
#   3. samples every rank's RSS and VmData on the node every 30 s into
#      rss_<job id>.txt (columns: time (s since 1970), pid, rank, RSS kB,
#      VmData kB), next to the node's memlog.
# Usage (from jobs/2d_inplane, on Chunhui Zhao's typed go):
#   CYCLE_DECK=config/fd_cycle/<deck> sbatch vista_memdiag.sbatch <run folder>
# Cost: 1 node x <= 50 min on gg, <= 0.3 SU (EAR20006).
SRC=${1:?the run folder to copy}
[ -d "$SRC" ] || { echo "ERROR: no run folder $SRC"; exit 2; }
DIAG=${SRC}_memdiag_$SLURM_JOB_ID
cp -a "$SRC" "$DIAG" || { echo "ERROR: copy to $DIAG failed"; exit 2; }
rm -f "$DIAG/DONE"
echo "copied $SRC -> $DIAG ($(du -sh "$DIAG" | cut -f1))"
export PETSC_OPTIONS="-ts_adapt_monitor"
export CYCLE_WALL_S=2400
# per-rank memory: OMPI_COMM_WORLD_RANK from each driver's environment
(
   while :; do
      t=$(date +%s)
      for p in $(pgrep -f seas_fd_cycle_driver); do
         r=$(tr '\0' '\n' < /proc/$p/environ 2>/dev/null | sed -n 's/^OMPI_COMM_WORLD_RANK=//p')
         [ -n "$r" ] || continue
         awk -v t=$t -v p=$p -v r=$r '/^VmRSS/{rss=$2} /^VmData/{d=$2} END{print t, p, r, rss, d}' /proc/$p/status
      done
      sleep 30
   done
) > "$DIAG/rss_$SLURM_JOB_ID.txt" 2>/dev/null &
RSS=$!
trap 'kill $RSS 2>/dev/null' EXIT
bash "$SLURM_SUBMIT_DIR/vista_2d_cycle.sbatch" "$DIAG"
rc=$?
echo "vista_2d_cycle.sbatch exit $rc; per-rank samples in $DIAG/rss_$SLURM_JOB_ID.txt"
exit $rc
