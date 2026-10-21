#!/bin/bash
# Shared body of the profiling jobs (speed-next, Chunhui Zhao 2026-10-20: "use some profiling tools to
# get the breakdown of each part in the code").  Sourced by map_lts_tpv205_200m_p{1,2}.sbatch, which set
# PROB, MESH, ORDER, TFINAL, GOLD and DECK.  One run of today's fastest setup (--lts rate2, the gold deck)
# under Linaro MAP (map --profile without express launch, through ibrun or else srun; samples every
# rank for the whole run, setup included), from seas_nitsche_dyn_driver_g (build_map_vista.sbatch: the
# same code with -g).  Outputs in $SCRATCH/nitsche_profiling/map_<prob>_200m_p<order>_<job>/:
#   profile_<route>.map  the MAP profile (opens in the MAP GUI); <route> = ibrun or srun
#   map_functions.csv  every function's time, from MAP's main-thread view (--export-functions)
#   profile.json     the whole profile as JSON (map --export), for the breakdown chart
#   map_<route>.log  MAP's own log of each try (--log); driver_<route>.log its console (the run's log)
#   run/, driver.log the run's own outputs and log; its stations are compared with GOLD's first
#                    TFINAL s (meant to be bitwise equal: the profiler only samples)
# If no try gives a profile (exit code, or no .map file), the same run is repeated with perf record on
# four ranks (perf_wrap.sh: rank 0, rank 1, the middle and the last rank) and `perf report` text per
# rank, and the job then exits 1.
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
MAP_ROUTE=
if command -v map >/dev/null; then
   # Express launch (map --profile ... mpirun ...) fails here.  Job 1085936: TACC's default Forge
   # configuration (~/.forge_24.1.1/system.config) sets the MPI type "OpenMPI (Compatibility)", which has
   # no express launch.  Job 1086945, with --mpi=openmpi: the licence and express launch started, then
   # "pindirect exited before reaching MPIR_Breakpoint"; express launch waits for mpirun's MPIR interface,
   # which Open MPI 5's mpirun no longer has (docs.open-mpi.org, app-debug/mpir-tools).  So MAP starts
   # the job itself (map [OPTION...] PROGRAM), by two routes tried in turn until one gives a profile:
   #   ibrun  MPI type openmpi-compat, MAP's mpirun being TACC's ibrun (system.config: mpirun =
   #          /usr/local/bin/ibrun, override mpirun = yes), the launcher of the gold runs; ibrun still
   #          starts Open MPI 5's mpirun, so this may fail the same way
   #   srun   MPI type slurm, srun --mpi=pmix (Vista's Slurm has pmix_v3), no mpirun at all
   # --profile implies --no-queue.  A failed start costs about a minute; a try that hangs is stopped after
   # MAP_TRY_S seconds (set by the job script, twice its estimate), so the next try and perf still fit.
   for route in ibrun srun; do
      case $route in
         ibrun) MAP_MPI=(--mpi=openmpi-compat) ;;
         srun)  MAP_MPI=(--mpi=slurm "--mpiargs=--mpi=pmix --cpu-bind=cores") ;;
      esac
      rm -rf $OUT/run; mkdir -p $OUT/run
      s=$(date +%s)
      echo "=== MAP ($route): map --profile ${MAP_MPI[*]} -n $N --nodes $SLURM_JOB_NUM_NODES --procs-per-node $((N / SLURM_JOB_NUM_NODES)) --log=$OUT/map_$route.log --output=$OUT/profile_$route.map --export-functions=$OUT/map_functions.csv ./seas_nitsche_dyn_driver_g $ARGS --output-dir $OUT/run"
      timeout --kill-after=30 ${MAP_TRY_S:-1200} map --profile "${MAP_MPI[@]}" -n $N --nodes $SLURM_JOB_NUM_NODES --procs-per-node $((N / SLURM_JOB_NUM_NODES)) \
          --log=$OUT/map_$route.log --output=$OUT/profile_$route.map --export-functions=$OUT/map_functions.csv \
          ./seas_nitsche_dyn_driver_g $ARGS --output-dir $OUT/run > $OUT/driver_$route.log 2>&1
      map_rc=$?
      echo "map ($route) exit $map_rc after $(( $(date +%s) - s )) s; profile file $(ls -l $OUT/profile_$route.map 2>/dev/null | awk '{print $5}') bytes"
      grep -ai "licen[cs]e\|error" $OUT/driver_$route.log | head -10
      if [ -f $OUT/profile_$route.map ]; then MAP_ROUTE=$route; map_rc=0; break; fi
      [ "$map_rc" = 0 ] && map_rc=3   # exit 0 but no profile file
      echo "--- MAP ($route) gave no profile; the end of map_$route.log and of driver_$route.log:"
      tail -20 $OUT/map_$route.log 2>/dev/null | cut -c1-300
      tail -20 $OUT/driver_$route.log | cut -c1-300
   done
   if [ -n "$MAP_ROUTE" ]; then
      cp $OUT/driver_$MAP_ROUTE.log $OUT/driver.log
      map --export=$OUT/profile.json $OUT/profile_$MAP_ROUTE.map > $OUT/map_export.log 2>&1
      echo "MAP profile from the $MAP_ROUTE route; map --export exit $?: $(ls -l $OUT/profile.json 2>/dev/null | awk '{print $5}') bytes"
   fi
fi
if [ -f $OUT/driver.log ]; then
   grep -a "partition:\|LTS\|fault ParaView output:\|time in ode.Step per rank\|fault work:\|halo\|setup (wall s\|=== done" $OUT/driver.log | cut -c1-400
   grep -a "^ *step [0-9]* t = " $OUT/driver.log | awk '{print "   " $0}'
fi
if [ -f $OUT/map_functions.csv ]; then
   echo "--- map_functions.csv: the first 60 lines ---"
   head -60 $OUT/map_functions.csv
fi
if [ "$map_rc" != 0 ]; then
   echo "=== MAP did not produce a profile (exit $map_rc); perf fallback on ranks 0 1 $((N / 2)) $((N - 1))"
   rm -rf $OUT/run; mkdir -p $OUT/run $OUT/perf
   export PERF_DIR=$OUT/perf PERF_RANKS="0 1 $((N / 2)) $((N - 1))"
   s=$(date +%s)
   ibrun jobs/vista_profiling/perf_wrap.sh ./seas_nitsche_dyn_driver_g $ARGS --output-dir $OUT/run > $OUT/driver.log 2>&1
   echo "perf run exit $? in $(( $(date +%s) - s )) s"
   grep -a "setup (wall s\|time in ode.Step per rank\|=== done" $OUT/driver.log | cut -c1-300
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
