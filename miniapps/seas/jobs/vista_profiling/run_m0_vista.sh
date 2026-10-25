#!/bin/bash
# Shared body of the M0 measurement jobs (speed-next plan rev 1.1, §6; Chunhui Zhao, 2026-10-20) and of
# the V1 verification jobs (v1_*.sbatch, TAG=v1; the output folder starts with TAG, m0 by default).
# Sourced by m0_lts_tpv205_200m_p{1,2}.sbatch and v1_lts_tpv205_200m_p{1,2}.sbatch, which set PROB,
# MESH, ORDER, TFINAL, GOLD (the LTS gold of the same mesh, order and rank layout) and VARIANTS, one
# "name|deck|flags" per variant.  Every variant runs local time stepping (--lts rate2), today's fastest
# setup; the M0 variants use existing flags and a deck copy only.  One allocation runs the variants in
# turn, then compares each variant's stations with the first TFINAL s of GOLD (variants whose stations are
# meant to be identical) and prints the timing lines the plan needs: wall time at each log line,
# time in ode.Step per rank, halo waits, the fault-frame time and the fault file size.
# Runs only from a clean pushed checkout (repo CLAUDE.md guard) whose driver was built at the same
# commit (build_lts_vista.sbatch or build_speed_next_vista.sbatch, through build_vista.sh, writes
# miniapps/seas/BUILD_COMMIT_VISTA) and whose mesh sha256 is in
# jobs/vista_profiling/meshes.sha256.
# V2 (v2_*.sbatch, fault-writer note rev1, Chunhui Zhao 2026-10-24) sets GOLD to V1's kept run of the
# same layout, and FAULT_COMPARE, one "a|b|flags" per comparison of two fault.vtkhdf files after the
# variants (@BASE@ = this job's output folder), run on one rank with seas_test_fault_frame_spread
# --compare a b flags; a difference fails the job.
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
echo "mesh $MESH sha256 $SHA  order $ORDER  tfinal $TFINAL  gold $GOLD"
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
BASE=$SCRATCH/nitsche_profiling/${TAG:-m0}_${PROB}_200m_p${ORDER}_$SLURM_JOB_ID
mkdir -p $BASE
{ echo "commit $HEAD"; echo "mesh $MESH $SHA"; echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ranks $SLURM_NTASKS"; } > $BASE/RUN_INFO.txt
T0=$(date +%s)
rc_all=0
for spec in "${VARIANTS[@]}"; do
   IFS='|' read -r v DECK FLAGS COMPARE <<< "$spec"
   OUT=$BASE/$v; mkdir -p $OUT
   cp "$DECK" $OUT/
   # memory in use on this (the first) node, which also holds rank 0 and the fault output
   ( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $OUT/meminfo_kB.log &
   MON=$!
   s=$(date +%s)
   ibrun ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order $ORDER --tfinal $TFINAL \
         --output-dir $OUT --log-every 1000 --lts rate2 $FLAGS > $OUT/driver.log 2>&1
   rc=$?; kill $MON 2>/dev/null; wait $MON 2>/dev/null
   [ $rc -eq 0 ] || rc_all=$rc
   echo "=== $v: deck $DECK, flags [--lts rate2 $FLAGS], 0 -> $TFINAL s: driver exit $rc, $(( $(date +%s) - s )) s in the job"
   grep -a "partition:\|LTS\|fault ParaView output:\|time in ode.Step per rank\|fault work:\|rank balance\|slowest ranks\|halo\|setup (wall\|station output:\|=== done" $OUT/driver.log || tail -20 $OUT/driver.log
   grep -a "^ *step [0-9]* t = " $OUT/driver.log | awk '{print "   " $0}'
   awk '{u=($2-$3)/1048576; if(u>m)m=u; t=$2/1048576} END{printf "first node memory: peak in use %.1f GB of %.1f GB (2 s samples)\n", m, t}' $OUT/meminfo_kB.log
   ls -l $OUT/fault.vtkhdf 2>/dev/null | awk '{printf "fault.vtkhdf %.1f MB\n", $5/1e6}'
   if [ "$COMPARE" = "gold" ]; then
      echo "--- $v stations against the first $TFINAL s of $GOLD (meant to be bitwise equal) ---"
      python3 -I jobs/vista_profiling/compare_stations.py $GOLD $OUT --tmax $TFINAL --time-tol 2e-9 --groups
   fi
done
for spec in "${FAULT_COMPARE[@]}"; do
   IFS='|' read -r A B CFLAGS <<< "$spec"
   A=${A//@BASE@/$BASE}; B=${B//@BASE@/$BASE}
   echo "--- fault file $A against $B ${CFLAGS:+($CFLAGS)} ---"
   [ -x ./seas_test_fault_frame_spread ] && [ ! ./seas_test_fault_frame_spread -ot ./seas_nitsche_dyn_driver ] \
      || { echo "ERROR: seas_test_fault_frame_spread missing or older than the driver; rerun the build job"; rc_all=2; continue; }
   [ -f "$A" ] && [ -f "$B" ] || { echo "ERROR: $( [ -f "$A" ] && echo "$B" || echo "$A") missing"; rc_all=2; continue; }
   ibrun -n 1 -o 0 ./seas_test_fault_frame_spread --compare "$A" "$B" $CFLAGS || rc_all=1
done
du -sh $BASE/*
echo "total $(( $(date +%s) - T0 )) s"
exit $rc_all
