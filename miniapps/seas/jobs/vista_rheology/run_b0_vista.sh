#!/bin/bash
# Shared body of the B0 jobs of rheology plan rev 1 (Chunhui Zhao, 2026-10-20), sourced by
# b0_tpv{205,104}_200m_p{1,2}.sbatch, which set PROB, DECK, MESH, ORDER, GOLD (the global-stepping gold's
# results), REF (the LTS gold folder beside it) and REF_DIFF.
# B0: with every TPV26/27 key absent (the TPV205/TPV104 decks), this branch must compute what the trunk
# computes.  The LTS gold jobs (1082586, 1082384, 1082587, 1082385) ran the trunk's code (ac3adeb, the code
# of trunk 4c97e32) with the defaults over 0 -> 1.0 s, with the gold's deck, mesh and rank layout, and
# kept that run's driver.log as $REF/results/gts_driver.log; their stdout records its station difference
# from the gold.  This job repeats that run at this commit and passes when
#   1. the "=== done" line equals the one in gts_driver.log (all but the wall time), and
#   2. compare_stations.py against the gold over 0 -> 1.0 s reports the same largest difference, at the
#      same station and quantity, as the LTS job did ($REF_DIFF).
# The cloud check behind it (3 ranks, 1000 m meshes, TPV205, TPV104 and TPV102 decks, 0 -> 0.6 s) found
# every output file byte-identical to the trunk code's.
# Runs only from a clean pushed checkout whose driver was built at the same commit
# (build_tpv2627_vista.sbatch writes BUILD_COMMIT_VISTA and copies the gold meshes).
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1      # one thread per MPI rank
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)    # repository root
cd $ROOT
for i in 1 2 3 4; do                                 # jobs that start together share the ref lock
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
   || { echo "ERROR: driver built at '$(cat BUILD_COMMIT_VISTA 2>/dev/null)', checkout is $HEAD; rerun build_tpv2627_vista.sbatch"; exit 2; }
[ -f "$MESH" ] || { echo "ERROR: $MESH missing; build_tpv2627_vista.sbatch copies it"; exit 2; }
SHA=$(sha256sum "$MESH" | cut -d' ' -f1)
grep -q "^$SHA  $MESH\$" jobs/vista_profiling/meshes.sha256 \
   || { echo "ERROR: $MESH sha256 $SHA is not the gold mesh (jobs/vista_profiling/meshes.sha256)"; exit 2; }
[ -d "$GOLD" ] && [ -f "$REF/results/gts_driver.log" ] || { echo "ERROR: $GOLD or $REF/results/gts_driver.log missing"; exit 2; }
echo "deck $DECK  mesh $MESH sha256 $SHA  order $ORDER  tfinal 1.0  gold $GOLD"
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
OUT=$SCRATCH/rheology/b0_${PROB}_200m_p${ORDER}_$SLURM_JOB_ID
mkdir -p $OUT
cp "$DECK" $OUT/
{ echo "commit $HEAD"; echo "mesh $MESH $SHA"; echo "deck $DECK"; echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ranks $SLURM_NTASKS"; } > $OUT/RUN_INFO.txt
T0=$(date +%s)
# --mesh-setup every-rank: the mesh setup of the LTS gold jobs.  The default since the merge of the Nitsche
# profiling branch (bc91ac9), rank0, differs from it by round-off only (profiling jobs 1090565-66), which
# would fail the two equality tests below.
ibrun ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order $ORDER --tfinal 1.0 \
      --mesh-setup every-rank --output-dir $OUT --log-every 1000 > $OUT/driver.log 2>&1
rc=$?
echo "driver exit $rc after $(( $(date +%s) - T0 )) s"
grep -a "time in ode.Step per rank\|=== done" $OUT/driver.log || { tail -20 $OUT/driver.log; exit 1; }
[ $rc -eq 0 ] || exit $rc
pass=1
DONE_NEW=$(grep -a "=== done" $OUT/driver.log | sed 's/; wall .*//')
DONE_REF=$(grep -a "=== done" $REF/results/gts_driver.log | sed 's/; wall .*//')
echo "done line, this run: $DONE_NEW"
echo "done line, trunk:    $DONE_REF"
[ "$DONE_NEW" = "$DONE_REF" ] && echo "1. done lines equal" || { echo "1. done lines DIFFER"; pass=0; }
CMP=$(python3 -I jobs/vista_profiling/compare_stations.py $GOLD $OUT --tmax 1.0)
echo "$CMP"
GOT=$(echo "$CMP" | sed -n 's/.*largest difference //p')
echo "largest difference from the gold, this run: $GOT"
echo "largest difference from the gold, trunk:    $REF_DIFF"
[ "$GOT" = "$REF_DIFF" ] && echo "2. same difference from the gold" || { echo "2. difference from the gold DIFFERS"; pass=0; }
echo "total $(( $(date +%s) - T0 )) s"
[ $pass -eq 1 ] && { echo "B0 $PROB p$ORDER: PASS"; exit 0; }
echo "B0 $PROB p$ORDER: FAIL"; exit 1
