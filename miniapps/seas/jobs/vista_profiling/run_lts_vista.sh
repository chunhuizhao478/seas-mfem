#!/bin/bash
# Shared body of the local time stepping (plan item 7, P5) gold jobs (Chunhui Zhao, 2026-10-19).
# Sourced by lts_tpv{205,104}_200m_p{1,2}.sbatch, which set PROB, DECK, MESH, ORDER, TFINAL and
# GOLD (the global-stepping gold of the same deck, mesh, order and rank layout).  One allocation
# runs, in turn:
#   gts   the defaults (global central difference, items 1, 3, 5 and 8a), 0 -> 1.0 s: the speed
#         reference at this commit and layout
#   lts   --lts rate2 (multirate leapfrog, cap 4 clusters, fault in one cluster), 0 -> TFINAL: the
#         new gold
# then compares the stations: gts against the first 1.0 s of GOLD (round-off), lts against gts
# over 1.0 s and against GOLD over TFINAL (both the time discretisation error).
# Runs only from a clean pushed checkout (repo CLAUDE.md guard) whose driver was built at the same
# commit (build_lts_vista.sbatch writes miniapps/seas/BUILD_COMMIT_VISTA) and whose mesh sha256 is in
# jobs/vista_profiling/meshes.sha256 (the gold meshes, copied in by build_lts_vista.sbatch).
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1      # one thread per MPI rank
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)    # repository root
cd $ROOT
# jobs that start together fetch in the same checkout and can collide on the ref lock ("cannot lock
# ref", job 1075918), so retry after a random pause
for i in 1 2 3 4; do
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
[ -x ./seas_nitsche_dyn_driver ] || { echo "ERROR: seas_nitsche_dyn_driver not built"; exit 2; }
[ "$(cat BUILD_COMMIT_VISTA 2>/dev/null)" = "$HEAD" ] \
   || { echo "ERROR: driver built at '$(cat BUILD_COMMIT_VISTA 2>/dev/null)', checkout is $HEAD; rerun build_lts_vista.sbatch"; exit 2; }
[ -f "$MESH" ] || { echo "ERROR: $MESH missing; build_lts_vista.sbatch copies it"; exit 2; }
SHA=$(sha256sum "$MESH" | cut -d' ' -f1)
grep -q "^$SHA  $MESH\$" jobs/vista_profiling/meshes.sha256 \
   || { echo "ERROR: $MESH sha256 $SHA is not the gold mesh (jobs/vista_profiling/meshes.sha256)"; exit 2; }
[ -d "$GOLD" ] || { echo "ERROR: gold folder $GOLD missing"; exit 2; }
echo "deck $DECK  mesh $MESH sha256 $SHA  order $ORDER  tfinal $TFINAL  gold $GOLD"
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
BASE=$SCRATCH/nitsche_profiling/lts_${PROB}_200m_p${ORDER}_$SLURM_JOB_ID
mkdir -p $BASE
cp "$DECK" $BASE/
{ echo "commit $HEAD"; echo "mesh $MESH $SHA"; echo "deck $DECK"; echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ranks $SLURM_NTASKS"; } > $BASE/RUN_INFO.txt
T0=$(date +%s)
rc_all=0
for v in gts lts; do
   case $v in
      gts) FLAGS=""; TF=1.0 ;;
      lts) FLAGS="--lts rate2"; TF=$TFINAL ;;
   esac
   OUT=$BASE/$v; mkdir -p $OUT
   # memory in use on this (the first) node, which also holds rank 0 and the fault output
   ( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $OUT/meminfo_kB.log &
   MON=$!
   s=$(date +%s)
   ibrun ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order $ORDER --tfinal $TF \
         --output-dir $OUT --log-every 1000 $FLAGS > $OUT/driver.log 2>&1
   rc=$?; kill $MON 2>/dev/null; wait $MON 2>/dev/null
   [ $rc -eq 0 ] || rc_all=$rc
   echo "=== $v [$FLAGS] 0 -> $TF s: driver exit $rc, $(( $(date +%s) - s )) s in the job"
   grep -a "partition:\|LTS\|fault ParaView output:\|time in ode.Step per rank\|fault work:\|=== done" $OUT/driver.log || tail -20 $OUT/driver.log
   awk '{u=($2-$3)/1048576; if(u>m)m=u; t=$2/1048576} END{printf "first node memory: peak in use %.1f GB of %.1f GB (2 s samples)\n", m, t}' $OUT/meminfo_kB.log
done
echo "=== stations (per quantity: max |a - b| over max |a| of slip, V, tau, sigma_n at that station) ==="
python3 -I jobs/vista_profiling/compare_stations.py $GOLD $BASE/gts --tmax 1.0
python3 -I jobs/vista_profiling/compare_stations.py $BASE/gts $BASE/lts --tmax 1.0 --time-tol 1e-9 --groups
python3 -I jobs/vista_profiling/compare_stations.py $GOLD $BASE/lts --tmax $TFINAL --time-tol 1e-9 --groups
du -sh $BASE/*
echo "total $(( $(date +%s) - T0 )) s"
exit $rc_all
