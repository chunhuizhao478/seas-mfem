#!/bin/bash
# The 2.5D note's Level 2 test, the map-view runs on the Mac (Chunhui Zhao,
# 2026-10-30: "all run must be run locally, and save data on local disk"):
# the four elongated_2d_*.toml decks with seas_nitsche_dyn_driver, one after
# the other, each into RESULTS/elongated_2d_<variant>_<commit>/ with run.log
# and RUN_INFO.txt.
#
# Usage (from miniapps/seas of a clean checkout at a pushed commit):
#   make -f tests/2d_inplane/2d_inplane.mk seas_nitsche_dyn_driver   (mfem-dev)
#   git rev-parse HEAD > BUILD_COMMIT_MAC_NITSCHE
#   elongated2p5d/scripts/run_2d_mac.sh RESULTS [NP]
# RESULTS the results folder (e.g.
# /Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/level2), NP the
# ranks (default 4).  The driver must be built at a commit whose code equals
# HEAD's (BUILD_COMMIT_MAC_NITSCHE; only elongated2p5d/, config/fd_cycle/,
# mapview_bp5/scripts/ and jobs/2d_inplane/ may differ).  Runs under the
# Mac's one-run lock (/Users/chunhuizhao/projects/.seas-mac-run.lock,
# mkdir), OMP_NUM_THREADS=1.
# Each run: p = 2 on the committed 200 m mesh to 40 s, 68,000 steps of
# 5.88e-4 s (the cloud's 0.2 s smoke on 4 ranks: 340 steps in 14.5 s, so
# about 48 min a run, 3.2 h for the four); stations every 2nd step
# (1.2 ms, the 3-D run's 4th step is about 1 ms), the fault CSV dump every
# 1700 steps (1 s).  Then compare_2p5d_3d.py with the 3-D run (README.md).
set -u
export LC_ALL=C OMP_NUM_THREADS=1
RES=${1:?usage: run_2d_mac.sh RESULTS [NP]}
NP=${2:-4}
LOCK=${SEAS_MAC_LOCK:-/Users/chunhuizhao/projects/.seas-mac-run.lock}
MESH=elongated2p5d/mesh/elongated_2d_200m.msh
MESH_SHA=5d74930bb73a9145302190657c36213da9342e8bf814b7a9f9fd8c472991afbe

# GitHub-first guard (CLAUDE.md), and the binary matches this code
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
HEAD=$(git rev-parse HEAD)
BUILT=$(cat BUILD_COMMIT_MAC_NITSCHE 2>/dev/null)
git cat-file -e "$BUILT^{commit}" 2>/dev/null \
   || { echo "ERROR: no BUILD_COMMIT_MAC_NITSCHE ('$BUILT'); build, then git rev-parse HEAD > BUILD_COMMIT_MAC_NITSCHE"; exit 2; }
CODE=$(git diff --name-only "$BUILT" HEAD | grep -v -E '^miniapps/seas/(elongated2p5d|config/fd_cycle|mapview_bp5/scripts|jobs/2d_inplane)/')
[ -z "$CODE" ] || { echo "ERROR: code changed since the build at $BUILT:"; echo "$CODE"; exit 2; }
[ -x ./seas_nitsche_dyn_driver ] || { echo "ERROR: seas_nitsche_dyn_driver not built"; exit 2; }
[ ./seas_nitsche_dyn_driver -nt BUILD_COMMIT_MAC_NITSCHE ] \
   && { echo "ERROR: seas_nitsche_dyn_driver is newer than BUILD_COMMIT_MAC_NITSCHE; rewrite the stamp at the build's commit"; exit 2; }
echo "commit $HEAD, seas_nitsche_dyn_driver built at $BUILT"
SHA=$( (sha256sum "$MESH" 2>/dev/null || shasum -a 256 "$MESH") | cut -d' ' -f1)
[ "$SHA" = "$MESH_SHA" ] || { echo "ERROR: $MESH sha256 '$SHA', expected $MESH_SHA"; exit 2; }

mkdir "$LOCK" 2>/dev/null || { echo "ERROR: another run holds $LOCK"; exit 2; }
trap 'rmdir "$LOCK"' EXIT
short=$(git rev-parse --short=8 HEAD)
fails=0
for d in wa2020 lehner liang plane_strain; do
   case $d in
      plane_strain) deck=elongated2p5d/configs/elongated_2d_plane_strain.toml ;;
      *) deck=elongated2p5d/configs/elongated_2d_2p5d_$d.toml ;;
   esac
   out=$RES/elongated_2d_${d}_$short
   [ ! -e "$out" ] || { echo "ERROR: $out exists"; fails=$((fails + 1)); continue; }
   mkdir -p "$out" || exit 2
   cp "$deck" "$out/"
   { echo "commit $HEAD"; echo "built $BUILT"; echo "deck $deck"; echo "mesh $MESH $SHA"
     echo "ranks $NP"; echo "start $(date -u +%Y-%m-%dT%H:%M:%SZ)"; } > "$out/RUN_INFO.txt"
   s=$(date +%s)
   mpirun -np "$NP" ./seas_nitsche_dyn_driver --config "$deck" --output-dir "$out" \
      --station-every 2 --dump-every 1700 --log-every 2000 > "$out/run.log" 2>&1
   rc=$?
   echo "exit $rc after $(( $(date +%s) - s )) s" >> "$out/RUN_INFO.txt"
   echo "$d: exit $rc in $(( $(date +%s) - s )) s; $(grep -m1 -F '2.5D finite width' "$out/run.log" | sed 's/^ *//')"
   grep "=== done" "$out/run.log" || tail -10 "$out/run.log"
   [ $rc -eq 0 ] || fails=$((fails + 1))
done
echo "failing runs: $fails"
[ $fails -eq 0 ]
