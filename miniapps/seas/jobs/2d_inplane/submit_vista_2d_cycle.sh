#!/bin/bash
# 2D map-view cycle decks on Vista (Chunhui Zhao, 2026-11-02: "apply to UCERF
# 2D plane strain, with 10MPa cap, and rerun the problem, ... let's also
# explore cluster to run this test because 7hr is too long"): the build and,
# for each deck, a chain of run links, each job waiting on the one before
# it.  On the login node, in the clone (vista_build.sbatch's header), after
# git fetch + checkout of the pushed commit:
#   CYCLE_DECKS="config/fd_cycle/ucerf2d_saf_sjf_500m_R70_avg_tip_relax_fd.toml" \
#      bash miniapps/seas/jobs/2d_inplane/submit_vista_2d_cycle.sh
# Jobs:
#   build    vista_build.sbatch (1 node, <= 4 h): PETSc, MFEM, the driver
#            and its 2D tests (test-relaxation among them)
#   per deck, link 1   vista_2d_cycle.sbatch (1 node, <= 24 h), after the
#                      build: the run from t = 0
#            link 2..  the same run folder, each after the link before
#                      (CYCLE_LINKS in all, default 2: the expected one and a
#                      spare, which exits at once when the run is done)
# The decks' chains run side by side (each its own node).  Every dependency
# is afterok with --kill-on-invalid-dep=yes: a job that fails cancels the
# jobs waiting on it, so the runs start only when the build and every check
# pass.  Prints the job ids and the run folders.
# CYCLE_BUILD_JOB=<id>: use a build job already submitted at this code
# (pending or running) instead of submitting one; CYCLE_BUILD_JOB=none: the
# binaries are already built at this code (BUILD_COMMIT_VISTA), no build.
# CYCLE_NTASKS: ranks per run (default 144, one gg node); CYCLE_WALL_S the
# driver's wall limit (vista_2d_cycle.sbatch).
# Cost: build <= 1.3 SU, <= 8 SU per link (EAR20006, gg, 0.33 SU per
# node-hour).
# Submit only on Chunhui Zhao's typed go for these jobs.
cd "$(dirname "$0")" || exit 2                  # sbatch from jobs/2d_inplane
ROOT=$(cd ../../../.. && pwd)
LINKS=${CYCLE_LINKS:-2}
NTASKS=${CYCLE_NTASKS:-144}
DECKS=${CYCLE_DECKS:?set CYCLE_DECKS to one or more decks, relative to miniapps/seas}
git -C "$ROOT" fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git -C "$ROOT" status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git -C "$ROOT" status --short; exit 2; }
[ -n "$(git -C "$ROOT" branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git -C "$ROOT" rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "commit $(git -C "$ROOT" rev-parse HEAD)"
[[ $NTASKS =~ ^[0-9]+$ ]] && [ "$NTASKS" -ge 1 ] && [ "$NTASKS" -le 144 ] \
   || { echo "ERROR: CYCLE_NTASKS '$NTASKS' is not 1 to 144"; exit 2; }
for d in $DECKS; do
   [ -f "$ROOT/miniapps/seas/$d" ] || { echo "ERROR: no deck miniapps/seas/$d"; exit 2; }
   P=$(sed -n 's/^prefix *= *"\([^"]*\)"$/\1/p' "$ROOT/miniapps/seas/$d" | head -1)
   [ -n "$P" ] || { echo "ERROR: no [cycle.output] prefix in $d"; exit 2; }
done

sub() {   # sbatch arguments...; prints the job id
   # TACC's sbatch prints its welcome banner on stdout before the id, even
   # with --parsable: the id is the last line ("<id>" or "<id>;<cluster>").
   local out id
   out=$(sbatch --parsable "$@") || { echo "ERROR: sbatch $* failed:" >&2; echo "$out" >&2; exit 1; }
   id=$(printf '%s\n' "$out" | tail -1)
   id=${id%%;*}
   [[ $id =~ ^[0-9]+$ ]] || { echo "ERROR: no job id from sbatch $*:" >&2; echo "$out" >&2; exit 1; }
   echo "$id"
}
DEP="--kill-on-invalid-dep=yes --dependency"
if [ "${CYCLE_BUILD_JOB:-}" = none ]; then
   jb=
   echo "build  none (CYCLE_BUILD_JOB=none: BUILD_COMMIT_VISTA $(cat ../../BUILD_COMMIT_VISTA 2>/dev/null))"
elif [ -n "${CYCLE_BUILD_JOB:-}" ]; then
   [[ $CYCLE_BUILD_JOB =~ ^[0-9]+$ ]] || { echo "ERROR: CYCLE_BUILD_JOB '$CYCLE_BUILD_JOB' is not a job id"; exit 2; }
   jb=$CYCLE_BUILD_JOB
   echo "build  $jb  (CYCLE_BUILD_JOB, not resubmitted)"
else
   jb=$(sub vista_build.sbatch) || exit 1
   echo "build  $jb"
fi
for d in $DECKS; do
   EXP=--export=ALL,CYCLE_DECK=$d${CYCLE_WALL_S:+,CYCLE_WALL_S=$CYCLE_WALL_S}
   P=$(sed -n 's/^prefix *= *"\([^"]*\)"$/\1/p' "$ROOT/miniapps/seas/$d" | head -1)
   if [ -n "$jb" ]; then
      j1=$(sub $DEP=afterok:$jb -n "$NTASKS" "$EXP" vista_2d_cycle.sbatch) || exit 1
   else
      j1=$(sub -n "$NTASKS" "$EXP" vista_2d_cycle.sbatch) || exit 1
   fi
   RUN=$ROOT-runs/${P}_$j1
   echo "$d"
   echo "  link 1 $j1${jb:+  (after $jb)}, $NTASKS ranks, run folder $RUN"
   prev=$j1
   for ((k = 2; k <= LINKS; k++)); do
      j=$(sub $DEP=afterok:$prev -n "$NTASKS" "$EXP" vista_2d_cycle.sbatch "$RUN") || exit 1
      echo "  link $k $j  (after $prev)"
      prev=$j
   done
done
