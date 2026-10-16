#!/bin/bash
# Shared body of the four Vista benchmark jobs (Chunhui Zhao, 2026-10-16).  Sourced by
# tpv{205,104}_200m_p{1,2}.sbatch, which set PROB, DECK, MESH, ORDER and TFINAL.
# Runs only from a clean pushed checkout (repo CLAUDE.md guard) whose driver was built at
# the same commit (build_verify_vista.sbatch writes miniapps/seas/BUILD_COMMIT_VISTA) and
# whose mesh checksum matches MESH_SHA256_VISTA.  Every log records the commit and the mesh sha256.
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
   || { echo "ERROR: driver built at '$(cat BUILD_COMMIT_VISTA 2>/dev/null)', checkout is $HEAD; rerun build_verify_vista.sbatch"; exit 2; }
[ -f "$MESH" ] || { echo "ERROR: $MESH missing; build_verify_vista.sbatch generates it"; exit 2; }
SHA=$(sha256sum "$MESH" | cut -d' ' -f1)
grep -q "^$SHA  $MESH\$" MESH_SHA256_VISTA || { echo "ERROR: $MESH sha256 $SHA not in MESH_SHA256_VISTA"; exit 2; }
echo "deck $DECK  mesh $MESH sha256 $SHA  order $ORDER  tfinal $TFINAL"
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
OUT=$SCRATCH/nitsche_benchmarks/${PROB}_200m_p${ORDER}_$SLURM_JOB_ID
mkdir -p $OUT
cp "$DECK" $OUT/
{ echo "commit $HEAD"; echo "mesh $MESH $SHA"; echo "deck $DECK"; echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ranks $SLURM_NTASKS"; } > $OUT/RUN_INFO.txt
T0=$(date +%s)
ibrun ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order $ORDER --tfinal $TFINAL \
      --output-dir $OUT --log-every 1000 > $OUT/driver.log 2>&1
rc=$?
echo "driver exit $rc after $(( $(date +%s) - T0 )) s"; grep "=== done" $OUT/driver.log || tail -20 $OUT/driver.log
du -sh $OUT
exit $rc
