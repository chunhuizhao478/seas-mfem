#!/bin/zsh
# heavylock.sh -- a mutex around the stages that can exhaust this machine.
#
# Two agents are building ShakeOut-D meshes on one 36 GB box (PREFERRED and ALT).
# The build guide's rule: refiners 5-13 GB, censuses 6-9 GB, mmg 4.5-11.7 GB, but
# **tetgen volume fill can want 20-50 GB on its own and must never overlap**.
# The one crash in this campaign came when SWAP passed ~90 % and macOS killed both
# Python processes -- RSS alone gave no warning.
#
# So: take this lock around tetgen and mmg. Everything else (LEB refinement, gate
# census, verification) is well-behaved and may overlap freely.
#
#   heavylock.sh acquire <name> [timeout_s]   wait for the lock, then hold it
#   heavylock.sh release <name>               drop it
#   heavylock.sh status                       who holds it, and the memory picture
#
# The lock is a directory (mkdir is atomic on macOS; `flock` is not available).
# A stale lock whose PID is gone is reclaimed automatically -- a crashed job must
# not wedge the other agent forever.
LOCK=/tmp/safs_shakeoutD_heavy.lock
CMD=${1:-status}
NAME=${2:-anon}
TIMEOUT=${3:-14400}

mem_line() {
  local sw=$(sysctl -n vm.swapusage | sed 's/.*used = \([0-9.]*\)M.*/\1/')
  local st=$(sysctl -n vm.swapusage | sed 's/total = \([0-9.]*\)M.*/\1/')
  local pct=$(echo "$sw $st" | awk '{printf "%.0f", 100*$1/$2}')
  echo "swap ${sw}M/${st}M (${pct}%)"
}

case $CMD in
  acquire)
    local_start=$(date +%s)
    while true; do
      if mkdir "$LOCK" 2>/dev/null; then
        echo "$$" > "$LOCK/pid"; echo "$NAME" > "$LOCK/name"; date > "$LOCK/since"
        echo "[lock] ACQUIRED by $NAME (pid $$)  $(mem_line)"
        exit 0
      fi
      # reclaim if the holder is gone
      if [ -f "$LOCK/pid" ]; then
        holder=$(cat "$LOCK/pid" 2>/dev/null)
        if [ -n "$holder" ] && ! kill -0 "$holder" 2>/dev/null; then
          echo "[lock] stale lock from pid $holder ($(cat $LOCK/name 2>/dev/null)) -- reclaiming"
          rm -rf "$LOCK"; continue
        fi
      fi
      now=$(date +%s)
      if [ $((now - local_start)) -ge $TIMEOUT ]; then
        echo "[lock] TIMEOUT after ${TIMEOUT}s waiting for $(cat $LOCK/name 2>/dev/null)"; exit 1
      fi
      sleep 20
    done ;;
  release)
    if [ -d "$LOCK" ]; then rm -rf "$LOCK"; echo "[lock] released by $NAME"; else echo "[lock] not held"; fi ;;
  status)
    if [ -d "$LOCK" ]; then
      echo "[lock] HELD by $(cat $LOCK/name 2>/dev/null) pid $(cat $LOCK/pid 2>/dev/null) since $(cat $LOCK/since 2>/dev/null)"
    else
      echo "[lock] free"
    fi
    echo "  $(mem_line)"
    ps -Ao rss=,command= | grep -E "tetgen|mmg3d|leb_|census|fill_domain" | grep -v grep \
      | awk '{printf "  %.1f GB  %s\n", $1/1048576, substr($0, index($0,$2), 90)}' ;;
  *) echo "usage: $0 {acquire <name> [timeout]|release <name>|status}"; exit 2 ;;
esac
