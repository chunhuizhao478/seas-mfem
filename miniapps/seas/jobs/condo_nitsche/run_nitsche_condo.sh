#!/bin/bash
# Shared body of the SCEC condo benchmark jobs (Chunhui Zhao, 2026-10-16).  Sourced by
# tpv{205,104}_100m_p1.sbatch, which set PROB, RES, DECK, MESH, ORDER and TFINAL.
# Runs only from a clean pushed checkout (repo CLAUDE.md guard) whose driver was built at
# the same commit (build_verify_condo.sbatch writes miniapps/seas/BUILD_COMMIT_DISCOVERY) and
# whose mesh checksum matches meshes.sha256.  Every log records the commit and the mesh sha256.
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
echo "host $(hostname)  nodes $SLURM_JOB_NUM_NODES  ranks $SLURM_NTASKS  per node $SLURM_NTASKS_PER_NODE"
[ -n "$SLURM_NTASKS_PER_NODE" ] || { echo "ERROR: SLURM_NTASKS_PER_NODE not set"; exit 2; }
[ -x ./seas_nitsche_dyn_driver ] || { echo "ERROR: seas_nitsche_dyn_driver not built"; exit 2; }
[ "$(cat BUILD_COMMIT_DISCOVERY 2>/dev/null)" = "$HEAD" ] \
   || { echo "ERROR: driver built at '$(cat BUILD_COMMIT_DISCOVERY 2>/dev/null)', checkout is $HEAD; rerun build_verify_condo.sbatch"; exit 2; }
[ -f "$MESH" ] || { echo "ERROR: $MESH missing; build_verify_condo.sbatch copies it in"; exit 2; }
SHA=$(sha256sum "$MESH" | cut -d' ' -f1)
grep -q "^$SHA  $MESH\$" jobs/condo_nitsche/meshes.sha256 \
   || { echo "ERROR: $MESH sha256 $SHA not in jobs/condo_nitsche/meshes.sha256"; exit 2; }
echo "deck $DECK  mesh $MESH sha256 $SHA  order $ORDER  tfinal $TFINAL"
module purge
module load gcc/13.3.0 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0 hdf5/1.14.3
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
OUT=/scratch1/$USER/nitsche_benchmarks/${PROB}_${RES}_p${ORDER}_$SLURM_JOB_ID
mkdir -p $OUT
cp "$DECK" $OUT/
{ echo "commit $HEAD"; echo "mesh $MESH $SHA"; echo "deck $DECK"; echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ranks $SLURM_NTASKS"; } > $OUT/RUN_INFO.txt
T0=$(date +%s)
# ppr:N:node, not ppr:N:package: the package form crashes prterun silently on these nodes (condo map test 6459051)
mpirun -np $SLURM_NTASKS --map-by ppr:$SLURM_NTASKS_PER_NODE:node --bind-to core \
       ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order $ORDER --tfinal $TFINAL \
       --output-dir $OUT --log-every 1000 > $OUT/driver.log 2>&1
rc=$?
echo "driver exit $rc after $(( $(date +%s) - T0 )) s"; grep "=== done" $OUT/driver.log || tail -20 $OUT/driver.log
du -sh $OUT
exit $rc
