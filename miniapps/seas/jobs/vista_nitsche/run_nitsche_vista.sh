#!/bin/bash
# Shared body of the Vista benchmark jobs (Chunhui Zhao, 2026-10-16).  Sourced by
# tpv{205,104}_200m_p{1,2,3}.sbatch, tpv102_200m_p3.sbatch, tpv{205,104}_100m_p1.sbatch,
# tpv{205,104}_mirror_*.sbatch, tpv{6,7}_100m_p{1,2,3}.sbatch, tpv31_200m_p{1,2,3,4}.sbatch, tpv104_200m_p4.sbatch and
# tpv34_200m_p{1,2,3,4}.sbatch and tpv34_200m_p3_check.sbatch, which set PROB, DECK, MESH, ORDER, TFINAL and RES (default 200m).
# TPV34 also sets TABLE, the 3-D material table (not in the repository: 444 MB): the job passes it
# as --material-table and runs only when its sha256 is in jobs/vista_nitsche/tables.sha256.
# Optional: LOG_EVERY (default 1000), STALL_SAMPLE_AFTER (seconds; stall sampler, off by default) and
# STALL_ABORT=1 (with the sampler: cancel the job when it found no step line three times).
# DRIVER_ARGS: further driver options, appended as given (TPV34 p = 4: --dense-fill direct).
# Runs only from a clean pushed checkout (repo CLAUDE.md guard) whose driver was built at
# the same commit (build_verify_vista.sbatch writes miniapps/seas/BUILD_COMMIT_VISTA) and
# whose mesh checksum is in MESH_SHA256_VISTA or MESH_SHA256_MIRROR_VISTA (meshes made on Vista by
# build_verify_vista.sbatch and make_mirror_meshes_vista.sbatch), MESH_SHA256_TPV2627_VISTA (the
# TPV26/27 meshes of jobs/vista_rheology/build_tpv2627_vista.sbatch) or jobs/vista_nitsche/meshes.sha256
# (meshes copied in).  Every log records the commit and the mesh sha256.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1      # one thread per MPI rank
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)    # repository root
cd $ROOT
# Jobs that start together fetch in the same checkout and can collide on the ref lock ("cannot lock
# ref", job 1075918), so retry after a random pause.
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
   || { echo "ERROR: driver built at '$(cat BUILD_COMMIT_VISTA 2>/dev/null)', checkout is $HEAD; rerun build_verify_vista.sbatch"; exit 2; }
[ -f "$MESH" ] || { echo "ERROR: $MESH missing; build_verify_vista.sbatch generates the 200 m meshes and the TPV6 100 m one, make_mirror_meshes_vista.sbatch the mirror ones, the TPV205 and TPV104 100 m ones are copied in"; exit 2; }
SHA=$(sha256sum "$MESH" | cut -d' ' -f1)
grep -q "^$SHA  $MESH\$" MESH_SHA256_VISTA MESH_SHA256_MIRROR_VISTA MESH_SHA256_TPV34_VISTA MESH_SHA256_TPV2627_VISTA \
     jobs/vista_nitsche/meshes.sha256 2>/dev/null \
   || { echo "ERROR: $MESH sha256 $SHA not in MESH_SHA256_VISTA, MESH_SHA256_MIRROR_VISTA, MESH_SHA256_TPV34_VISTA, MESH_SHA256_TPV2627_VISTA or jobs/vista_nitsche/meshes.sha256"; exit 2; }
echo "deck $DECK  mesh $MESH sha256 $SHA  order $ORDER  tfinal $TFINAL${DRIVER_ARGS:+  driver args: $DRIVER_ARGS}"
TABLE_ARGS=""
if [ -n "$TABLE" ]; then
   [ -f "$TABLE" ] || { echo "ERROR: material table $TABLE missing"; exit 2; }
   TSHA=$(sha256sum "$TABLE" | cut -d' ' -f1)
   grep -q "^$TSHA  $(basename "$TABLE")\$" jobs/vista_nitsche/tables.sha256 \
      || { echo "ERROR: $TABLE sha256 $TSHA not in jobs/vista_nitsche/tables.sha256"; exit 2; }
   echo "material table $TABLE sha256 $TSHA"
   TABLE_ARGS="--material-table $TABLE"
fi
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load openmpi/5.0.11 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
OUT=$SCRATCH/nitsche_benchmarks/${PROB}_${RES:-200m}_p${ORDER}_$SLURM_JOB_ID
mkdir -p $OUT
cp "$DECK" $OUT/
{ echo "commit $HEAD"; echo "mesh $MESH $SHA"; echo "deck $DECK"; echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ranks $SLURM_NTASKS"; [ -z "$TABLE" ] || echo "material table $TABLE $TSHA"; [ -z "$DRIVER_ARGS" ] || echo "driver args $DRIVER_ARGS"; } > $OUT/RUN_INFO.txt
# Memory on the first node, 2 s samples, for the speed record: time MemTotal MemAvailable CommitLimit
# Committed_AS (kB).  Vista compute nodes use strict overcommit, so allocations fail once
# Committed_AS reaches CommitLimit (half the RAM), well before MemAvailable runs out.
( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable|^CommitLimit|^Committed_AS/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $OUT/meminfo_kB.log &
MON=$!
# Optional stall sampler (off unless the job sets STALL_SAMPLE_AFTER, seconds): if driver.log has no
# "  step " line that long after the start, record what the ranks on the first node are doing, three
# times 60 s apart, in $OUT/stall_sample_<n>.txt: CPU and state of every rank, and the call stack of
# three ranks (gstack or gdb if the node has one and allows it; the attempt and its error are logged
# either way).  Read-only: it never signals the driver.
SMP=""
T0=$(date +%s)
if [ -n "$STALL_SAMPLE_AFTER" ]; then
   ( sleep "$STALL_SAMPLE_AFTER"
     for n in 1 2 3; do
        grep -qa "^  step " $OUT/driver.log 2>/dev/null && exit 0
        F=$OUT/stall_sample_$n.txt
        { echo "sample $n at $(date +%s) host $(hostname), $(( $(date +%s) - T0 )) s after start; no step line yet"
          echo "--- last 5 lines of driver.log"; tail -5 $OUT/driver.log
          echo "--- ranks on this node: pid state %cpu cputime rss_kB"
          ps -C seas_nitsche_dyn_driver -o pid=,stat=,pcpu=,time=,rss= | sort -n
          echo "--- ptrace_scope $(cat /proc/sys/kernel/yama/ptrace_scope 2>&1)  gstack $(command -v gstack || echo none)  gdb $(command -v gdb || echo none)"
          PIDS=$(ps -C seas_nitsche_dyn_driver -o pid= | sort -n)
          for pid in $(echo "$PIDS" | sed -n '1p;2p;$p'); do
             echo "--- stack of pid $pid"
             if command -v gstack >/dev/null; then gstack $pid 2>&1 | head -80
             elif command -v gdb >/dev/null; then gdb -batch -p $pid -ex "bt" 2>&1 | tail -60
             else echo "no gstack or gdb on this node"; fi
             echo "--- wait channel / syscall of pid $pid: $(cat /proc/$pid/wchan 2>&1) / $(cat /proc/$pid/syscall 2>&1)"
          done; } > $F 2>&1
        sleep 60
     done
     # Optional: STALL_ABORT=1 ends the job after the third sample when there is still no step line,
     # so a stalled run does not hold its nodes to the time limit.
     if [ "$STALL_ABORT" = 1 ] && ! grep -qa "^  step " $OUT/driver.log 2>/dev/null; then
        echo "no step line $(( $(date +%s) - T0 )) s after start: STALL_ABORT=1, cancelling job $SLURM_JOB_ID" | tee -a $OUT/stall_sample_3.txt
        scancel $SLURM_JOB_ID
     fi ) &
   SMP=$!
fi
ibrun ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order $ORDER --tfinal $TFINAL \
      --output-dir $OUT --log-every ${LOG_EVERY:-1000} $TABLE_ARGS $DRIVER_ARGS > $OUT/driver.log 2>&1
rc=$?
kill $MON $SMP 2>/dev/null; wait $MON $SMP 2>/dev/null
[ -z "$STALL_SAMPLE_AFTER" ] || { ls $OUT/stall_sample_*.txt 2>/dev/null && cat $OUT/stall_sample_1.txt || echo "no stall sample taken (step lines appeared)"; grep -a "^  step " $OUT/driver.log | tail -3; }
echo "driver exit $rc after $(( $(date +%s) - T0 )) s"; grep "=== done" $OUT/driver.log || tail -20 $OUT/driver.log
grep -a "time in ode.Step per rank\|fault work:\|fault ParaView output: [0-9]* frames" $OUT/driver.log
awk '{u=($2-$3)/1048576; if(u>m)m=u; t=$2/1048576; if($5/1048576>c)c=$5/1048576; l=$4/1048576} END{printf "first node memory: peak in use %.1f GB of %.1f GB, peak committed %.1f of %.1f GiB (2 s samples)\n", m, t, c, l}' $OUT/meminfo_kB.log
du -sh $OUT
exit $rc
