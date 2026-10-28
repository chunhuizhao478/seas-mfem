#!/bin/bash
# SCEC BP3-FD on Vista: the build and the run's chain submitted in one go,
# each job waiting on the one before it (Chunhui Zhao, 2026-10-27: "cancel
# this run, add snapshots, submit it on cluster").  On the login node, after
# the clone and STAGE=fetch (vista_build.sbatch's header):
#   bash miniapps/seas/jobs/2d_inplane/submit_vista_bp3_fd.sh
# Jobs:
#   build    vista_build.sbatch (1 node, <= 4 h): PETSc, MFEM, the driver,
#            the 2D tests of the BP3-FD path and the 2D snapshot check
#   link 1   vista_bp3_fd.sbatch (1 node, <= 24 h), after the build: the run
#            from t = 0
#   link 2.. the same run folder, each after the link before (BP3FD_LINKS in
#            all, default 2: the expected one and a spare, which exits at
#            once when the run is done)
# Every dependency is afterok with --kill-on-invalid-dep=yes: a job that
# fails cancels the jobs waiting on it, so the run starts only when the
# build and every check pass.  Prints the job ids.
# BP3FD_BUILD_JOB=<id>: use a build job already submitted at this code
# (pending or running) instead of submitting one.  BP3FD_DECK is passed to
# the run links (default config/fd_cycle/bp3_fd_d60_thrust.toml).
# Cost: build <= 1.3 SU, <= 8 SU per link (EAR20006, gg, 0.33 SU per
# node-hour).
# Submit only on Chunhui Zhao's typed go for these jobs.
cd "$(dirname "$0")" || exit 2                  # sbatch from jobs/2d_inplane
ROOT=$(cd ../../../.. && pwd)
LINKS=${BP3FD_LINKS:-2}
[ -z "$(git -C "$ROOT" status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git -C "$ROOT" status --short; exit 2; }
[ -n "$(git -C "$ROOT" branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git -C "$ROOT" rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "commit $(git -C "$ROOT" rev-parse HEAD)"

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
EXP=--export=ALL${BP3FD_DECK:+,BP3FD_DECK=$BP3FD_DECK}
if [ -n "$BP3FD_BUILD_JOB" ]; then
   [[ $BP3FD_BUILD_JOB =~ ^[0-9]+$ ]] || { echo "ERROR: BP3FD_BUILD_JOB '$BP3FD_BUILD_JOB' is not a job id"; exit 2; }
   jb=$BP3FD_BUILD_JOB
   echo "build  $jb  (BP3FD_BUILD_JOB, not resubmitted)"
else
   jb=$(sub vista_build.sbatch) || exit 1
   echo "build  $jb"
fi
j1=$(sub $DEP=afterok:$jb "$EXP" vista_bp3_fd.sbatch) || exit 1
NAME=$(basename "${BP3FD_DECK:-config/fd_cycle/bp3_fd_d60_thrust.toml}" .toml)
RUN=$ROOT-runs/bp3fd_${NAME#bp3_fd_}_$j1
echo "link 1 $j1  (after $jb), run folder $RUN"
prev=$j1
for ((k = 2; k <= LINKS; k++)); do
   j=$(sub $DEP=afterok:$prev "$EXP" vista_bp3_fd.sbatch "$RUN") || exit 1
   echo "link $k $j  (after $prev)"
   prev=$j
done
