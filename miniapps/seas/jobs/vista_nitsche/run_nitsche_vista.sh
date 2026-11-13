#!/bin/bash
# Shared body of the Vista benchmark jobs (Chunhui Zhao, 2026-10-16).  Sourced by
# tpv{205,104}_200m_p{1,2,3}.sbatch, tpv102_200m_p3.sbatch, tpv{205,104}_100m_p1.sbatch,
# tpv{205,104}_mirror_*.sbatch, tpv{6,7}_100m_p{1,2,3}.sbatch, tpv31_200m_p{1,2,3}.sbatch and
# tpv34_200m_p{1,2,3}.sbatch, which set PROB, DECK, MESH, ORDER, TFINAL and RES (default 200m).
# TPV34 also sets TABLE, the 3-D material table (not in the repository: 444 MB): the job passes it
# as --material-table and runs only when its sha256 is in jobs/vista_nitsche/tables.sha256.
# Runs only from a clean pushed checkout (repo CLAUDE.md guard) whose driver was built at
# the same commit (build_verify_vista.sbatch writes miniapps/seas/BUILD_COMMIT_VISTA) and
# whose mesh checksum is in MESH_SHA256_VISTA or MESH_SHA256_MIRROR_VISTA (meshes made on Vista by
# build_verify_vista.sbatch and make_mirror_meshes_vista.sbatch) or jobs/vista_nitsche/meshes.sha256
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
grep -q "^$SHA  $MESH\$" MESH_SHA256_VISTA MESH_SHA256_MIRROR_VISTA MESH_SHA256_TPV34_VISTA jobs/vista_nitsche/meshes.sha256 2>/dev/null \
   || { echo "ERROR: $MESH sha256 $SHA not in MESH_SHA256_VISTA, MESH_SHA256_MIRROR_VISTA, MESH_SHA256_TPV34_VISTA or jobs/vista_nitsche/meshes.sha256"; exit 2; }
echo "deck $DECK  mesh $MESH sha256 $SHA  order $ORDER  tfinal $TFINAL"
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
{ echo "commit $HEAD"; echo "mesh $MESH $SHA"; echo "deck $DECK"; echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ranks $SLURM_NTASKS"; [ -z "$TABLE" ] || echo "material table $TABLE $TSHA"; } > $OUT/RUN_INFO.txt
# Memory on the first node, 2 s samples, for the speed record: time MemTotal MemAvailable CommitLimit
# Committed_AS (kB).  Vista compute nodes use strict overcommit, so allocations fail once
# Committed_AS reaches CommitLimit (half the RAM), well before MemAvailable runs out.
( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable|^CommitLimit|^Committed_AS/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $OUT/meminfo_kB.log &
MON=$!
T0=$(date +%s)
ibrun ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order $ORDER --tfinal $TFINAL \
      --output-dir $OUT --log-every 1000 $TABLE_ARGS > $OUT/driver.log 2>&1
rc=$?
kill $MON 2>/dev/null; wait $MON 2>/dev/null
echo "driver exit $rc after $(( $(date +%s) - T0 )) s"; grep "=== done" $OUT/driver.log || tail -20 $OUT/driver.log
grep -a "time in ode.Step per rank\|fault work:\|fault ParaView output: [0-9]* frames" $OUT/driver.log
awk '{u=($2-$3)/1048576; if(u>m)m=u; t=$2/1048576; if($5/1048576>c)c=$5/1048576; l=$4/1048576} END{printf "first node memory: peak in use %.1f GB of %.1f GB, peak committed %.1f of %.1f GiB (2 s samples)\n", m, t, c, l}' $OUT/meminfo_kB.log
du -sh $OUT
exit $rc
