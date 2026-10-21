#!/bin/bash
# Shared body of the profiling jobs (speed-next, Chunhui Zhao 2026-10-20: "use some profiling tools to
# get the breakdown of each part in the code").  Sourced by map_lts_tpv205_200m_p{1,2}.sbatch, which set
# PROB, MESH, ORDER, TFINAL, GOLD and DECK.  One run of today's fastest setup (--lts rate2, the gold deck)
# under Linaro MAP (map --profile, express launch through mpirun; samples every rank for the whole run,
# setup included), from seas_nitsche_dyn_driver_g (build_map_vista.sbatch: the same code with -g).
# Outputs in $SCRATCH/nitsche_profiling/map_<prob>_200m_p<order>_<job>/:
#   *.map            the MAP profile (opens in the MAP GUI)
#   map_functions.csv  every function's time, from MAP's main-thread view (--export-functions)
#   profile.json     the whole profile as JSON (map --export), for the breakdown chart
#   map.log          MAP's own log (--log)
#   run/, driver.log the run's own outputs and log; its stations are compared with GOLD's first
#                    TFINAL s (meant to be bitwise equal: the profiler only samples)
# If MAP fails (exit code, or no .map file), the same run is repeated with perf record on four ranks
# (perf_wrap.sh: rank 0, rank 1, the middle and the last rank) and `perf report` text per rank, and the
# job then exits 1.
# Runs only from a clean pushed checkout (repo CLAUDE.md guard) whose -g driver was built at the same
# commit (BUILD_COMMIT_VISTA_G) and whose mesh sha256 is in jobs/vista_profiling/meshes.sha256.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1      # one thread per MPI rank
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)    # repository root
cd $ROOT
for i in 1 2 3 4; do                                 # retry: concurrent jobs collide on the ref lock
   git fetch -q origin && break
   [ $i = 4 ] && { echo "ERROR: git fetch failed"; exit 2; }
   sleep $(( 5 + RANDOM % 20 ))
done
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "commit $(git rev-parse HEAD)"
UNTRACKED=$(git ls-files --others --exclude-standard -- '*.cpp' '*.hpp' '*.inl' Makefile '*.toml' '*.sbatch')
[ -z "$UNTRACKED" ] || { echo "ERROR: untracked sources in $ROOT:"; echo "$UNTRACKED"; exit 2; }
HEAD=$(git rev-parse HEAD)
cd miniapps/seas
echo "host $(hostname)  nodes $SLURM_JOB_NUM_NODES  ranks $SLURM_NTASKS"
[ -x ./seas_nitsche_dyn_driver_g ] || { echo "ERROR: seas_nitsche_dyn_driver_g not built (build_map_vista.sbatch)"; exit 2; }
[ "$(cat BUILD_COMMIT_VISTA_G 2>/dev/null)" = "$HEAD" ] \
   || { echo "ERROR: -g driver built at '$(cat BUILD_COMMIT_VISTA_G 2>/dev/null)', checkout is $HEAD"; exit 2; }
[ -f "$MESH" ] || { echo "ERROR: $MESH missing; build_lts_vista.sbatch copies it"; exit 2; }
SHA=$(sha256sum "$MESH" | cut -d' ' -f1)
grep -q "^$SHA  $MESH\$" jobs/vista_profiling/meshes.sha256 \
   || { echo "ERROR: $MESH sha256 $SHA is not the gold mesh (jobs/vista_profiling/meshes.sha256)"; exit 2; }
[ -d "$GOLD" ] || { echo "ERROR: gold folder $GOLD missing"; exit 2; }
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
# map/24.1.1 puts Forge's own libstdc++ (/home1/apps/forge/24.1.1/lib, without GLIBCXX_3.4.32) on
# LD_LIBRARY_PATH ahead of gcc 14's, and the gcc 14 driver then stops at startup on every rank, under MAP
# and under perf alike (job 1085112).  gcc 14's libstdc++ goes back in front; it also serves MAP, being newer.
GCC_LIB=$(dirname "$(g++ -print-file-name=libstdc++.so.6)")
module load map/24.1.1 || echo "WARNING: module map/24.1.1 did not load"
case "$GCC_LIB" in
   /*) export LD_LIBRARY_PATH=$GCC_LIB:$LD_LIBRARY_PATH ;;
   *) echo "ERROR: gcc 14's libstdc++ not found (g++ -print-file-name=libstdc++.so.6)"; exit 2 ;;
esac
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
echo "map: $(command -v map)  mpirun: $(command -v mpirun)  perf: $(command -v perf)  perf_event_paranoid $(cat /proc/sys/kernel/perf_event_paranoid)"
LDD=$(ldd ./seas_nitsche_dyn_driver_g 2>&1)          # the driver must load in this environment
echo "libstdc++: $(echo "$LDD" | awk '/libstdc\+\+/ {print $3}')"
if echo "$LDD" | grep -q "not found"; then
   echo "ERROR: seas_nitsche_dyn_driver_g cannot load its libraries here:"; echo "$LDD" | grep "not found"; exit 2
fi
OUT=$SCRATCH/nitsche_profiling/map_${PROB}_200m_p${ORDER}_$SLURM_JOB_ID
mkdir -p $OUT/run
cp "$DECK" $OUT/
{ echo "commit $HEAD"; echo "mesh $MESH $SHA"; echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ranks $SLURM_NTASKS"; } > $OUT/RUN_INFO.txt
ARGS="--config $DECK --mesh $MESH --order $ORDER --tfinal $TFINAL --log-every 1000 --lts rate2"
N=$SLURM_NTASKS
T0=$(date +%s)
map_rc=1
if command -v map >/dev/null; then
   # --mpi=openmpi: TACC's default Forge configuration (~/.forge_24.1.1/system.config, made at the first
   # MAP run) turns MPI auto-detection off and sets the type "OpenMPI (Compatibility)", which has no express
   # launch (map --list-mpis); this job uses express launch.  The likely cause of job 1085936's MAP exit,
   # 8 s after "Waiting for licence"; MAP's own log (--log) says why if it exits again.
   echo "=== MAP: map --profile --mpi=openmpi --log=$OUT/map.log --export-functions=$OUT/map_functions.csv mpirun -np $N --map-by core --bind-to core ./seas_nitsche_dyn_driver_g $ARGS --output-dir $OUT/run"
   map --profile --mpi=openmpi --log=$OUT/map.log --export-functions=$OUT/map_functions.csv \
       mpirun -np $N --map-by core --bind-to core ./seas_nitsche_dyn_driver_g $ARGS \
       --output-dir $OUT/run > $OUT/driver.log 2>&1
   map_rc=$?
   MAPF=$(ls -t seas_nitsche_dyn_driver_g_${N}p_*.map 2>/dev/null | head -1)
   echo "map exit $map_rc after $(( $(date +%s) - T0 )) s; profile file '${MAPF}'"
   if [ -n "$MAPF" ]; then
      mv "$MAPF" $OUT/
      map --export=$OUT/profile.json $OUT/$MAPF > $OUT/map_export.log 2>&1
      echo "map --export exit $?: $(ls -l $OUT/profile.json 2>/dev/null | awk '{print $5}') bytes"
   elif [ "$map_rc" = 0 ]; then
      map_rc=3   # exit 0 but no profile file
   fi
   grep -ai "licen[cs]e\|error" $OUT/driver.log | head -10
fi
grep -a "partition:\|LTS\|fault ParaView output:\|time in ode.Step per rank\|fault work:\|halo\|=== done" $OUT/driver.log | cut -c1-400
grep -a "^ *step [0-9]* t = " $OUT/driver.log | awk '{print "   " $0}'
if [ -f $OUT/map_functions.csv ]; then
   echo "--- map_functions.csv: the first 60 lines ---"
   head -60 $OUT/map_functions.csv
fi
if [ "$map_rc" != 0 ]; then
   echo "=== MAP did not produce a profile (exit $map_rc); the end of map.log and of MAP's console (map_driver.log):"
   tail -20 $OUT/map.log 2>/dev/null | cut -c1-300
   mv $OUT/driver.log $OUT/map_driver.log; tail -20 $OUT/map_driver.log | cut -c1-300
   echo "=== perf fallback on ranks 0 1 $((N / 2)) $((N - 1))"
   rm -rf $OUT/run; mkdir -p $OUT/run $OUT/perf
   export PERF_DIR=$OUT/perf PERF_RANKS="0 1 $((N / 2)) $((N - 1))"
   s=$(date +%s)
   ibrun jobs/vista_profiling/perf_wrap.sh ./seas_nitsche_dyn_driver_g $ARGS --output-dir $OUT/run > $OUT/driver.log 2>&1
   echo "perf run exit $? in $(( $(date +%s) - s )) s"
   grep -a "time in ode.Step per rank\|=== done" $OUT/driver.log | cut -c1-300
   for f in $OUT/perf/perf.r*.data; do
      [ -f "$f" ] || continue
      perf report -i $f --stdio --no-children --sort dso,symbol --percent-limit 0.3 > ${f%.data}.txt 2> /dev/null
      echo "--- $(basename ${f%.data}): top symbols ---"; grep -v "^#" ${f%.data}.txt | grep -v "^$" | head -25 | cut -c1-200 | sed "s/ *$//"
   done
fi
echo "--- stations against the first $TFINAL s of $GOLD (meant to be bitwise equal) ---"
python3 -I jobs/vista_profiling/compare_stations.py $GOLD $OUT/run --tmax $TFINAL --time-tol 2e-9 --groups
du -sh $OUT/*
echo "total $(( $(date +%s) - T0 )) s"
# a job without a MAP profile fails in Slurm, perf fallback or not
[ "$map_rc" = 0 ] || { echo "ERROR: MAP did not produce a profile (exit $map_rc)"; exit 1; }
