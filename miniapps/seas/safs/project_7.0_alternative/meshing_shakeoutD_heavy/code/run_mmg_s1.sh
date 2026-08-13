#!/bin/zsh
# run_mmg_s1.sh -- launch the stage-1 mmg pass DETACHED, logging to build_tmp.
#
# Two reasons this is a script and not an inline command:
#
#  1. The pass is a multi-hour job. Run as a managed background task it was killed
#     at ~52 min with no output; `setsid` + `nohup` detaches it from the caller's
#     lifecycle so only the machine decides when it ends.
#  2. It must NOT overlap another mmg3d. Measured on this box: ALT's mmg alongside
#     PREFERRED's (13.0 GB, 8 h in) was killed while PREFERRED's survived untouched
#     and swap never moved -- so neither process looked close to the limit
#     individually. tetgen's fill DOES coexist fine (55 s alongside it); mmg does
#     not. The guard below refuses to start rather than repeat that.
#
# Switches match PREFERRED's stage 1 exactly (-hgrad 1.3 -hausd 30 -hmin 115
# -hmax 5000 -m 20000), fault frozen as RequiredTriangles, --no-freeze so the
# whole mesh can have its shape repaired -- the base CDT is slivers by construction.
set -u
cd "${0:a:h}/.."
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
LOG=build_tmp/mmg_s1.log

if [ "${1:-}" != "--force" ]; then
  if pgrep -x mmg3d_O3 >/dev/null 2>&1; then
    echo "REFUSING: an mmg3d is already running -- concurrent mmg gets this one killed."
    ps -Ao pid,etime,rss,command | grep mmg3d_O3 | grep -v grep \
      | awk '{printf "  pid %s  %s  %.2f GB\n", $1, $2, $3/1048576}'
    echo "Re-run with --force only if you mean it."
    exit 1
  fi
fi

echo "[start] $(date)  swap $(sysctl -n vm.swapusage | sed 's/.*used = //;s/ .*//')"
setsid nohup $PY -u code/mmg_refine_sizemap.py \
    build_tmp/base.msh build_tmp/s1.msh build_tmp/metric_s1.sol \
    --hgrad 1.3 --hausd 30 --hmin 115 --hmax 5000 --mem-mb 20000 --no-freeze \
    >> $LOG 2>&1 &
echo "[launched] pid $! -> $LOG"
