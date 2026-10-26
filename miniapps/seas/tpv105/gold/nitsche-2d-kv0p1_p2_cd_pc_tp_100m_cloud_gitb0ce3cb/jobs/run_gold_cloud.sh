#!/bin/bash
# 2D p=2 gold runs (TPV205-2D, TPV105-2D) in the cloud container at the pushed 2d-inplane head.
# Guard of CLAUDE.md: clean tracked tree, HEAD on GitHub, binary built here at HEAD.
S=/tmp/claude-0/-home-user-seas-mfem/ce558891-cad7-5003-a322-13bbf23d3889/scratchpad
ROOT=/home/user/seas-mfem
O=$S/gold/runs
cd $ROOT
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
C=$(git rev-parse HEAD)
echo "commit $C"
cd $ROOT/miniapps/seas
find . -name '*.o' -delete; rm -f seas_nitsche_dyn_driver
s=$(date +%s)
$S/mk.sh seas_nitsche_dyn_driver > $S/gold/build.log 2>&1 || { echo "ERROR: build failed (gold/build.log)"; exit 2; }
echo "build: $(( $(date +%s) - s )) s, binary sha256 $(sha256sum seas_nitsche_dyn_driver | cut -d' ' -f1)"
[ "$(git -C $ROOT rev-parse HEAD)" = "$C" ] && [ -z "$(git -C $ROOT status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: tree changed during the build"; exit 2; }
mkdir -p $O; cp seas_nitsche_dyn_driver $O/
export OMP_NUM_THREADS=1
run() {   # name deck args...
   local n=$1 d=$2; shift 2
   rm -rf $O/$n; mkdir -p $O/$n
   local s=$(date +%s)
   { echo "commit $C"; echo "deck $d"; echo "args $*"; echo "host cloud container $(nproc) cores, mpirun -np 4"; date -u; } > $O/$n.info
   mpirun --allow-run-as-root --oversubscribe --bind-to none -np 4 $O/seas_nitsche_dyn_driver --config $d \
      --output-dir $O/$n "$@" < /dev/null > $O/$n.log 2>&1
   local rc=$?
   echo "$n: exit $rc in $(( $(date +%s) - s )) s; $(grep '=== done' $O/$n.log | cut -c1-120)"
   date -u >> $O/$n.info
}
run tpv205_2d_p2 tpv205/configs/tpv205_2d_nitsche.toml --order 2 --log-every 2000
run tpv105_2d_p2 tpv105/configs/tpv105_2d_nitsche.toml --order 2 --log-every 2000
echo "all done $(date -u +%H:%M:%S)"
