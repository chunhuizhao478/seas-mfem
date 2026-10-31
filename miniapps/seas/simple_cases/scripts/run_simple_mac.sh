#!/bin/bash
# The simple cases on the Mac (Chunhui Zhao, 2026-10-31: "go ahead with the
# simple cases, run them on local Mac"; 2026-10-30: "all run must be run
# locally, and save data on local disk"): the twelve
# simple_cases/configs/*.toml decks with seas_fd_cycle_driver, one after the
# other, each into RESULTS/<case>_<commit>/ with run.log and RUN_INFO.txt.
# Then analyse_simple.py (README.md).
#
# Usage (from miniapps/seas of a clean checkout at a pushed commit):
#   make -f tests/2d_inplane/2d_inplane.mk seas_fd_cycle_driver   (conda env mumps-petsc)
#   git rev-parse HEAD > BUILD_COMMIT_MAC_FDCYCLE
#   simple_cases/scripts/run_simple_mac.sh RESULTS [NP [CASE ...]]
# RESULTS the results folder (e.g.
# /Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/simple_cases),
# NP the ranks (default 4), CASE a subset (default all twelve, in the order
# below: D1f, D2f and Q1 first).  WALL_LIMIT (env, default 2700 s) is each
# run's --wall-limit: past it the driver stops at the next QD step with a
# checkpoint (exit 8), and the script goes on to the next case; such a run
# can be continued later with --restart from its checkpoint.
# The driver must be built at a commit whose code equals HEAD's
# (BUILD_COMMIT_MAC_FDCYCLE; only simple_cases/, bend_convergence/,
# elongated2p5d/, config/fd_cycle/, mapview_bp5/scripts/ and jobs/2d_inplane/
# may differ).  Runs under the Mac's one-run lock
# (/Users/chunhuizhao/projects/.seas-mac-run.lock, mkdir), OMP_NUM_THREADS=1.
# A run that stops at sigma_n' < 0 (the cases without a floor) exits
# non-zero: that is a result here, so the script goes on to the next case.
# The friction law's "[FRIC-GUARD] sigma_n_eff <= 0" lines (one per DOF and
# stage once a DOF is in tension) are kept up to the first 200 in run.log;
# the rest are counted at its end.
set -u
export LC_ALL=C OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
RES=${1:?usage: run_simple_mac.sh RESULTS [NP [CASE ...]]}
NP=${2:-4}
shift $(( $# < 2 ? $# : 2 ))
ALL="D1f_sharp015_floor D2f_sharp05_floor Q1_sharp015 Q0_straight Q2_sharp05 Q3_arc05L10 Q4_sharp015vs D0_straight D1_sharp015 D2_sharp05 D3_arc05L10 D4_sharp015vs"
CASES=${*:-$ALL}
WALL=${WALL_LIMIT:-2700}
LOCK=${SEAS_MAC_LOCK:-/Users/chunhuizhao/projects/.seas-mac-run.lock}

mesh_of() {
   case $1 in
   Q0_*|D0_*) echo straight_h500 ;;
   Q1_*|D1_*|D1f_*) echo sharp015_h500 ;;
   Q2_*|D2_*|D2f_*) echo sharp05_h500 ;;
   Q3_*|D3_*) echo arc05L10_h500 ;;
   Q4_*|D4_*) echo sharp015vs_h500 ;;
   *) echo unknown ;;
   esac
}
mesh_sha() {
   case $1 in
   straight_h500) echo cc57a4cae781e1bed98da25f3b4e951624f5135e4b0c6028d6a4573cb63c7f04 ;;
   sharp015_h500) echo fb5b5526d68d74106b5146be086e55089aad18e03edcc1d5a7ded6be3d39319d ;;
   sharp05_h500) echo 4f3ab900ad8a46be406c270c69610d588a8c4780e7a35d2dc48b7310e4b269fd ;;
   arc05L10_h500) echo 05956834b5a5d81e72192ac56fb83ec1b02079f0ffbbe45eab33d7c92337a529 ;;
   sharp015vs_h500) echo 83c995f896981a81060594d26e7081841aa8e8b5158aaca6fbd8d3cfe7ae9d04 ;;
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
CODE=$(git diff --name-only "$BUILT" HEAD | grep -v -E '^miniapps/seas/(simple_cases|bend_convergence|elongated2p5d|config/fd_cycle|mapview_bp5/scripts|jobs/2d_inplane)/')
[ -z "$CODE" ] || { echo "ERROR: code changed since the build at $BUILT:"; echo "$CODE"; exit 2; }
[ -x ./seas_fd_cycle_driver ] || { echo "ERROR: seas_fd_cycle_driver not built"; exit 2; }
[ ./seas_fd_cycle_driver -nt BUILD_COMMIT_MAC_FDCYCLE ] \
   && { echo "ERROR: seas_fd_cycle_driver is newer than BUILD_COMMIT_MAC_FDCYCLE; rewrite the stamp at the build's commit"; exit 2; }
python3 simple_cases/scripts/make_decks.py --check > /dev/null \
   || { echo "ERROR: the decks differ from make_decks.py's"; exit 2; }
echo "commit $HEAD, seas_fd_cycle_driver built at $BUILT"
for c in $CASES; do
   [ -f simple_cases/configs/$c.toml ] || { echo "ERROR: no deck for $c"; exit 2; }
   m=$(mesh_of $c); f=simple_cases/mesh/$m.msh
   SHA=$( (sha256sum "$f" 2>/dev/null || shasum -a 256 "$f") | cut -d' ' -f1)
   [ "$SHA" = "$(mesh_sha $m)" ] || { echo "ERROR: $f sha256 '$SHA', expected $(mesh_sha $m)"; exit 2; }
done

mkdir "$LOCK" 2>/dev/null || { echo "ERROR: another run holds $LOCK"; exit 2; }
trap 'rmdir "$LOCK"' EXIT
short=$(git rev-parse --short=8 HEAD)
for c in $CASES; do
   deck=simple_cases/configs/$c.toml
   out=$RES/${c}_$short
   [ ! -e "$out" ] || { echo "$c: $out exists, skipped"; continue; }
   mkdir -p "$out" || exit 2
   cp "$deck" "$out/"
   m=$(mesh_of $c)
   { echo "commit $HEAD"; echo "built $BUILT"; echo "deck $deck"
     echo "mesh simple_cases/mesh/$m.msh $(mesh_sha $m)"
     echo "ranks $NP"; echo "wall_limit $WALL s"; echo "start $(date -u +%Y-%m-%dT%H:%M:%SZ)"; } > "$out/RUN_INFO.txt"
   s=$(date +%s)
   mpirun -np "$NP" ./seas_fd_cycle_driver --deck "$deck" --output-dir "$out" --wall-limit "$WALL" 2>&1 \
      | awk '/\[FRIC-GUARD\] sigma_n_eff <= 0/ { n++; if (n > 200) next } { print; fflush() }
             END { if (n > 200) print "[run_simple_mac] " n - 200 " more [FRIC-GUARD] sigma_n_eff <= 0 lines not kept (" n " in all)" }' \
      > "$out/run.log"
   rc=${PIPESTATUS[0]}
   echo "exit $rc after $(( $(date +%s) - s )) s" >> "$out/RUN_INFO.txt"
   echo "$c: exit $rc in $(( $(date +%s) - s )) s; $(grep -m1 -E 'MFEM abort|Verification failed|STOP:|wall limit' "$out/run.log" | cut -c1-160)"
done
echo "all done; next: python3 simple_cases/scripts/analyse_simple.py $RES"
