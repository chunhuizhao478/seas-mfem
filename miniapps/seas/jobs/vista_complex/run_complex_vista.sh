#!/bin/bash
# Shared body of the Vista runs of the complex fault geometries (complex-geometries plan rev 2, phase P5;
# Chunhui Zhao, 2026-10-21).  Sourced by tpv{12,22,15}_*_p*.sbatch, which set PROB, DECK, MESH, RES,
# ORDER, TFINAL and FLAGS (extra driver flags, e.g. --lts rate2).
# Runs only from a clean pushed checkout (repo CLAUDE.md guard) whose driver was built and verified at the
# same commit (build_verify_complex_vista.sbatch writes miniapps/seas/BUILD_COMMIT_VISTA, through
# build_vista.sh, and VERIFIED_COMMIT_COMPLEX_VISTA only when every check passed) and whose mesh sha256 is
# in MESH_SHA256_COMPLEX_VISTA (build_verify_complex_vista.sbatch) or MESH_SHA256_TPV12_100M_VISTA
# (make_tpv12_100m_mesh_vista.sbatch).  Every log records the commit and the mesh sha256.
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
UNTRACKED=$(git ls-files --others --exclude-standard -- '*.cpp' '*.hpp' '*.inl' Makefile '*.toml' '*.sbatch' '*.geo' '*.csv')
[ -z "$UNTRACKED" ] || { echo "ERROR: untracked sources in $ROOT:"; echo "$UNTRACKED"; exit 2; }
HEAD=$(git rev-parse HEAD)
cd miniapps/seas
echo "host $(hostname)  nodes $SLURM_JOB_NUM_NODES  ranks $SLURM_NTASKS"
[ -x ./seas_nitsche_dyn_driver ] || { echo "ERROR: seas_nitsche_dyn_driver not built"; exit 2; }
[ "$(cat BUILD_COMMIT_VISTA 2>/dev/null)" = "$HEAD" ] \
   || { echo "ERROR: driver built at '$(cat BUILD_COMMIT_VISTA 2>/dev/null)', checkout is $HEAD; rerun build_verify_complex_vista.sbatch"; exit 2; }
[ "$(cat VERIFIED_COMMIT_COMPLEX_VISTA 2>/dev/null)" = "$HEAD" ] \
   || { echo "ERROR: build_verify_complex_vista.sbatch has not passed at $HEAD (VERIFIED_COMMIT_COMPLEX_VISTA: '$(cat VERIFIED_COMMIT_COMPLEX_VISTA 2>/dev/null)')"; exit 2; }
[ -f "$MESH" ] || { echo "ERROR: $MESH missing; build_verify_complex_vista.sbatch makes the 200 m meshes, make_tpv12_100m_mesh_vista.sbatch the 100 m one"; exit 2; }
SHA=$(sha256sum "$MESH" | cut -d' ' -f1)
grep -q "^$SHA  $MESH\$" MESH_SHA256_COMPLEX_VISTA MESH_SHA256_TPV12_100M_VISTA 2>/dev/null \
   || { echo "ERROR: $MESH sha256 $SHA not in MESH_SHA256_COMPLEX_VISTA or MESH_SHA256_TPV12_100M_VISTA"; exit 2; }
echo "deck $DECK  mesh $MESH sha256 $SHA  order $ORDER  tfinal $TFINAL  flags [$FLAGS]"
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
OUT=$SCRATCH/complex_geometries/${PROB}_${RES}_p${ORDER}_$SLURM_JOB_ID
mkdir -p $OUT
cp "$DECK" $OUT/
{ echo "commit $HEAD"; echo "mesh $MESH $SHA"; echo "deck $DECK"; echo "order $ORDER tfinal $TFINAL flags $FLAGS"
  echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ranks $SLURM_NTASKS"; } > $OUT/RUN_INFO.txt
# memory in use on this (the first) node, which also holds rank 0 and the fault output
( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $OUT/meminfo_kB.log &
MON=$!
T0=$(date +%s)
ibrun ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order $ORDER --tfinal $TFINAL \
      --output-dir $OUT --log-every 1000 $FLAGS > $OUT/driver.log 2>&1
rc=$?; kill $MON 2>/dev/null; wait $MON 2>/dev/null
echo "driver exit $rc after $(( $(date +%s) - T0 )) s"
grep -a "partition:\|LTS (plan\|prestress: \[\|fault options\|fault tips:\|forced rupture:\|tension floor:\|time in ode.Step per rank\|fault work:" $OUT/driver.log
grep -a "=== done" $OUT/driver.log || tail -20 $OUT/driver.log
awk '{u=($2-$3)/1048576; if(u>m)m=u; t=$2/1048576} END{printf "first node memory: peak in use %.1f GB of %.1f GB (2 s samples)\n", m, t}' $OUT/meminfo_kB.log
du -sh $OUT
exit $rc
