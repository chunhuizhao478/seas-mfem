#!/bin/bash
#SBATCH -J tpv205_200m_p2_dense
#SBATCH -o tpv205_200m_p2_dense_%j.out
#SBATCH -e tpv205_200m_p2_dense_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 00:45:00
#SBATCH -A EAR26024
# Vista measurement of plan item 8a, dense element and face blocks for K u + K_eta v (Chunhui Zhao,
# 2026-10-18): TPV205, 200 m gold mesh, p = 2, 0 -> 1.0 s, the deck of gold 1075916
# (tpv205/configs/tpv205_nitsche.toml), 16 nodes x 144 ranks as the gold.  One allocation runs, in turn:
#   off      --fast-fault 0 --fused-spmv 0 --fast-fault-output 0 --dense-blocks 0  (the old calls)
#   item4    --dense-blocks 0  (items 1, 3, 4, 5: the defaults before item 8a, "on" of job 1078825)
#   dense    the defaults (items 1, 3, 5 and 8a)
# then compares the stations: off against the first 1.0 s of gold 1075916, item4 and dense against off,
# and dense against item4.
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_profiling_vista.sbatch passed at this commit.
# Estimate from job 1078825 (321 s off, 244 s with items 1-5; about 100-115 s of each is setup): about
# 12 min = 3.2 node-h = 1.1 SU.  Requested: 16 nodes x 0.75 h x 0.33 SU = 4 SU at most.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1      # one thread per MPI rank
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)    # repository root
cd $ROOT
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
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
   || { echo "ERROR: driver built at '$(cat BUILD_COMMIT_VISTA 2>/dev/null)', checkout is $HEAD; rerun build_profiling_vista.sbatch"; exit 2; }
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_200m.msh
GOLD=tpv205/gold/nitsche-kv0p1_p2_cd_pc_200m_jobid1075916_git81759da/results
TFINAL=1.0
[ -f "$MESH" ] || { echo "ERROR: $MESH missing; build_profiling_vista.sbatch copies it"; exit 2; }
SHA=$(sha256sum "$MESH" | cut -d' ' -f1)
grep -q "^$SHA  $MESH\$" jobs/vista_profiling/meshes.sha256 \
   || { echo "ERROR: $MESH sha256 $SHA is not the gold mesh (jobs/vista_profiling/meshes.sha256)"; exit 2; }
echo "deck $DECK  mesh $MESH sha256 $SHA  order 2  tfinal $TFINAL"
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
BASE=$SCRATCH/nitsche_profiling/tpv205_200m_p2_dense_$SLURM_JOB_ID
mkdir -p $BASE
cp "$DECK" $BASE/
{ echo "commit $HEAD"; echo "mesh $MESH $SHA"; echo "deck $DECK"; echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ranks $SLURM_NTASKS"; } > $BASE/RUN_INFO.txt
T0=$(date +%s)
for v in off item4 dense; do
   case $v in
      off)   FLAGS="--fast-fault 0 --fused-spmv 0 --fast-fault-output 0 --dense-blocks 0" ;;
      item4) FLAGS="--dense-blocks 0" ;;
      dense) FLAGS="" ;;
   esac
   OUT=$BASE/$v; mkdir -p $OUT
   # memory in use on this (the first) node, which also holds rank 0 and the fault output
   ( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $OUT/meminfo_kB.log &
   MON=$!
   s=$(date +%s)
   ibrun ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order 2 --tfinal $TFINAL \
         --output-dir $OUT --log-every 1000 $FLAGS > $OUT/driver.log 2>&1
   rc=$?; kill $MON 2>/dev/null; wait $MON 2>/dev/null
   echo "=== $v [$FLAGS]: driver exit $rc, $(( $(date +%s) - s )) s in the job"
   grep -a "partition:\|fault ParaView output:\|time in ode.Step per rank\|fault work:\|=== done" $OUT/driver.log || tail -20 $OUT/driver.log
   awk '{u=($2-$3)/1048576; if(u>m)m=u; t=$2/1048576} END{printf "first node memory: peak in use %.1f GB of %.1f GB (2 s samples)\n", m, t}' $OUT/meminfo_kB.log
done
echo "=== stations (per quantity: max |a - b| over max |a| of slip, V, tau, sigma_n at that station) ==="
python3 -I jobs/vista_profiling/compare_stations.py $GOLD $BASE/off --tmax $TFINAL
python3 -I jobs/vista_profiling/compare_stations.py $BASE/off $BASE/item4
python3 -I jobs/vista_profiling/compare_stations.py $BASE/off $BASE/dense
python3 -I jobs/vista_profiling/compare_stations.py $BASE/item4 $BASE/dense
du -sh $BASE/*
echo "total $(( $(date +%s) - T0 )) s"
