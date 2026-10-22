#!/bin/bash
# Fully dynamic BP5 cycle on Vista: the p = 2 run and its checks submitted in
# one go, each job waiting on the ones it needs (Chunhui Zhao, 2026-10-22:
# "let's maximum the test speed by submit all three jobs together").  On the
# login node, after the clone and STAGE=fetch (vista_build.sbatch's header):
#   bash miniapps/seas/jobs/fd_cycle/submit_vista_p2.sh
# Jobs:
#   build   vista_build.sbatch, FDC_PART=build (1 node, <= 3 h): PETSc,
#           MUMPS, HDF5, MFEM and cycle-all
#   tests   vista_build.sbatch, FDC_PART=tests (1 node, <= 5 h), after the
#           build: every cycle test, the new ones included
#   check   vista_p2_check.sbatch (4 nodes, <= 5 h), after the build, side
#           by side with the tests: builds the H-matrix, then the speed-ups
#           at the run's size against the plain paths
#   link 1  vista_p2_1800yr.sbatch (4 nodes, <= 48 h), after the tests AND
#           the check: the run from t = 0
#   link 2.. the same run folder, each after the link before (FDC_LINKS in
#           all, default 3: the expected two and a spare, which exits at
#           once when the run is done)
# Every dependency is afterok with --kill-on-invalid-dep=yes: a job that
# fails cancels the jobs waiting on it, so the run starts only when the
# build, every test and every check pass.  Prints the job ids.
# FDC_BUILD_JOB=<id>: use a build job already submitted (pending or
# running; a finished job soon drops out of Slurm and can no longer be a
# dependency) instead of submitting one; the jobs/ and config/ files may
# have changed since, the code not (the tests and check refuse it).
# Cost: build <= 1 SU, tests <= 1.7 SU, check <= 6.7 SU, <= 64 SU per link
# (EAR20006, gg, 0.33 SU per node-hour).
# Submit only on Chunhui Zhao's typed go for these jobs.
cd "$(dirname "$0")" || exit 2                  # sbatch from jobs/fd_cycle
ROOT=$(cd ../../../.. && pwd)
LINKS=${FDC_LINKS:-3}
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
if [ -n "$FDC_BUILD_JOB" ]; then                # a build job already submitted at HEAD
   [[ $FDC_BUILD_JOB =~ ^[0-9]+$ ]] || { echo "ERROR: FDC_BUILD_JOB '$FDC_BUILD_JOB' is not a job id"; exit 2; }
   jb=$FDC_BUILD_JOB
   echo "build  $jb  (FDC_BUILD_JOB, not resubmitted)"
else
   jb=$(sub -t 03:00:00 --export=ALL,FDC_PART=build vista_build.sbatch) || exit 1
   echo "build  $jb"
fi
jt=$(sub -J fdc_tests -o fdc_tests_%j.out -e fdc_tests_%j.err -t 05:00:00 \
         $DEP=afterok:$jb --export=ALL,FDC_PART=tests vista_build.sbatch) || exit 1
echo "tests  $jt  (after $jb)"
jc=$(sub $DEP=afterok:$jb vista_p2_check.sbatch) || exit 1
echo "check  $jc  (after $jb)"
j1=$(sub $DEP=afterok:$jt:$jc vista_p2_1800yr.sbatch) || exit 1
RUN=$ROOT-runs/p2_1800_$j1
echo "link 1 $j1  (after $jt and $jc), run folder $RUN"
prev=$j1
for ((k = 2; k <= LINKS; k++)); do
   j=$(sub $DEP=afterok:$prev vista_p2_1800yr.sbatch "$RUN") || exit 1
   echo "link $k $j  (after $prev)"
   prev=$j
done
