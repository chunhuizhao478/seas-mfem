#!/bin/bash
# The bend convergence test on the Mac (Chunhui Zhao, 2026-10-30: "prepare
# the bend convergence test"; "all run must be run locally, and save data on
# local disk"): the nine bend_convergence/configs/*.toml decks with
# seas_fd_cycle_driver, one after the other, each into
# RESULTS/<tag>_<commit>/ with run.log and RUN_INFO.txt.  Then
# analyse_bends.py (README.md).
#
# Usage (from miniapps/seas of a clean checkout at a pushed commit):
#   make -f tests/2d_inplane/2d_inplane.mk seas_fd_cycle_driver   (conda env mumps-petsc)
#   git rev-parse HEAD > BUILD_COMMIT_MAC_FDCYCLE
#   bend_convergence/scripts/run_bends_mac.sh RESULTS [NP [TAG ...]]
# RESULTS the results folder (e.g.
# /Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/bend_convergence),
# NP the ranks (default 8), TAG a subset (default all nine, sharp05_h500
# first: the check against the earlier 500 m sharp run, abort near 3.69 s).
# The driver must be built at a commit whose code equals HEAD's
# (BUILD_COMMIT_MAC_FDCYCLE; only bend_convergence/, elongated2p5d/,
# config/fd_cycle/, mapview_bp5/scripts/ and jobs/2d_inplane/ may differ).
# Runs under the Mac's one-run lock
# (/Users/chunhuizhao/projects/.seas-mac-run.lock, mkdir), OMP_NUM_THREADS=1.
# A run that stops at sigma_n' < 0 exits non-zero: that is a result here,
# so the script goes on to the next deck.
set -u
export LC_ALL=C OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
RES=${1:?usage: run_bends_mac.sh RESULTS [NP [TAG ...]]}
NP=${2:-8}
shift $(( $# < 2 ? $# : 2 ))
TAGS=${*:-sharp05_h500 sharp05_h250 sharp05_h125 arc05L05_h500 arc05L05_h250 arc05L05_h125 arc05L10_h500 arc05L10_h250 arc05L10_h125}
LOCK=${SEAS_MAC_LOCK:-/Users/chunhuizhao/projects/.seas-mac-run.lock}

mesh_sha() {
   case $1 in
   sharp05_h500) echo be86e9ce5b9740d2dea6c66c8c63527cae59eff70c70e362eedcf18043695c26 ;;
   sharp05_h250) echo 7610004069c69ef7e89a2e565e40304a4e0e762e8fd3ae83482670dd085de41d ;;
   sharp05_h125) echo 13b3d0d02930d2e295b1d01791e172e52357759dd04379919313ffd083de1e8c ;;
   arc05L05_h500) echo 0f8a78206be07cef7ab74f8ca28d823eefc8262a5d3a381ed13b43c13870259d ;;
   arc05L05_h250) echo c953ea83c280076ff90d0a8181628935995dd0fb1b766d2b6bb37c827dfef795 ;;
   arc05L05_h125) echo 9ba4d920668bff7691bf76052a9704c1c00ae7017b01b9d338ad62a53163dddc ;;
   arc05L10_h500) echo e00f734b9159b23f888183274c71e63fd220899c8728aed1bd2b26c3eeabf536 ;;
   arc05L10_h250) echo 7e9b649e2048a51dd9cacbf5045a77d300cd2d438ad47abf6aedf46a1279f781 ;;
   arc05L10_h125) echo 42673988d0974546e06da6fc6a433a9ab6924084a2947ad080e1ec1ebe5e2255 ;;
   *) echo unknown ;;
   esac
}

# GitHub-first guard (CLAUDE.md), and the binary matches this code
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
HEAD=$(git rev-parse HEAD)
BUILT=$(cat BUILD_COMMIT_MAC_FDCYCLE 2>/dev/null)
git cat-file -e "$BUILT^{commit}" 2>/dev/null \
   || { echo "ERROR: no BUILD_COMMIT_MAC_FDCYCLE ('$BUILT'); build, then git rev-parse HEAD > BUILD_COMMIT_MAC_FDCYCLE"; exit 2; }
CODE=$(git diff --name-only "$BUILT" HEAD | grep -v -E '^miniapps/seas/(bend_convergence|elongated2p5d|config/fd_cycle|mapview_bp5/scripts|jobs/2d_inplane)/')
[ -z "$CODE" ] || { echo "ERROR: code changed since the build at $BUILT:"; echo "$CODE"; exit 2; }
[ -x ./seas_fd_cycle_driver ] || { echo "ERROR: seas_fd_cycle_driver not built"; exit 2; }
[ ./seas_fd_cycle_driver -nt BUILD_COMMIT_MAC_FDCYCLE ] \
   && { echo "ERROR: seas_fd_cycle_driver is newer than BUILD_COMMIT_MAC_FDCYCLE; rewrite the stamp at the build's commit"; exit 2; }
echo "commit $HEAD, seas_fd_cycle_driver built at $BUILT"
for t in $TAGS; do
   m=bend_convergence/mesh/$t.msh
   [ -f bend_convergence/configs/$t.toml ] || { echo "ERROR: no deck for $t"; exit 2; }
   SHA=$( (sha256sum "$m" 2>/dev/null || shasum -a 256 "$m") | cut -d' ' -f1)
   [ "$SHA" = "$(mesh_sha $t)" ] || { echo "ERROR: $m sha256 '$SHA', expected $(mesh_sha $t)"; exit 2; }
done

mkdir "$LOCK" 2>/dev/null || { echo "ERROR: another run holds $LOCK"; exit 2; }
trap 'rmdir "$LOCK"' EXIT
short=$(git rev-parse --short=8 HEAD)
for t in $TAGS; do
   deck=bend_convergence/configs/$t.toml
   out=$RES/${t}_$short
   [ ! -e "$out" ] || { echo "$t: $out exists, skipped"; continue; }
   mkdir -p "$out" || exit 2
   cp "$deck" "$out/"
   { echo "commit $HEAD"; echo "built $BUILT"; echo "deck $deck"
     echo "mesh bend_convergence/mesh/$t.msh $(mesh_sha $t)"
     echo "ranks $NP"; echo "start $(date -u +%Y-%m-%dT%H:%M:%SZ)"; } > "$out/RUN_INFO.txt"
   s=$(date +%s)
   mpirun -np "$NP" ./seas_fd_cycle_driver --deck "$deck" --output-dir "$out" > "$out/run.log" 2>&1
   rc=$?
   echo "exit $rc after $(( $(date +%s) - s )) s" >> "$out/RUN_INFO.txt"
   echo "$t: exit $rc in $(( $(date +%s) - s )) s; $(grep -m1 -E 'MFEM abort|Verification failed' "$out/run.log" | cut -c1-160)"
done
echo "all done; next: python3 bend_convergence/scripts/analyse_bends.py $RES"
