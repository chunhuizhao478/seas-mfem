#!/bin/bash
#SBATCH -J d1_gg_build
#SBATCH -o d1_gg_build_%j.out
#SBATCH -e d1_gg_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# GPU dynamic-rupture plan D1 (device plumbing), CPU builds for the gg checks
# (jobs/vista/d1_gg_check.sbatch, which only runs binaries).  The same two
# clones as D0a, both built by the trunk's build_vista.sh:
#   this checkout   gpu-dyn-nitsche, at the commit its upstream points to
#                   ($WORK/seas-mfem-gpu-dyn-nit)
#   REF             trunk a4c9dee, detached ($WORK/seas-mfem-gpu-dyn-ref, as left
#                   by D0a)
# Login node, before submitting:
#   cd $WORK/seas-mfem-gpu-dyn-nit && git pull --ff-only origin gpu-dyn-nitsche
#   cd miniapps/seas && sbatch jobs/vista/d1_gg_build.sbatch
# Steps:
#   0. the repo CLAUDE.md guard in both clones; REF must be exactly a4c9dee and
#      an ancestor of this checkout
#   1. build_vista.sh in both clones at once (each writes BUILD_COMMIT_VISTA)
#   2. the TPV205 and TPV104 1000 m meshes in this checkout: D0a's files (job
#      1080349) when their sha256 still match, otherwise regenerated from the
#      committed .geo files (TACC gmsh, msh 2.2); the sha256 must equal D0a's
#   3. the six Nitsche test targets in both clones: exit codes and PASSED/FAILED
#      counts must match, with 0 FAILED
# PASS needs all of 0-3.  Cost: 1 node x <= 2 h x 0.33 SU = <= 0.7 SU.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REF_SHA=a4c9dee0aa96a019792e48cce864fc68054aab26
ROOT=$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)              # submitted from miniapps/seas
REF=${REF:-$WORK/seas-mfem-gpu-dyn-ref}
[ -f "$ROOT/build_vista.sh" ] || { echo "ERROR: submit from miniapps/seas ($ROOT has no build_vista.sh)"; exit 2; }
[ -f "$REF/build_vista.sh" ]  || { echo "ERROR: REF $REF has no build_vista.sh"; exit 2; }
# shellcheck source=d1_env.sh
source "$ROOT/miniapps/seas/jobs/vista/d1_env.sh"
T0=$(date +%s)

# --- 0. guards
d1_guard "$ROOT"
[ "$(git -C "$ROOT" rev-parse HEAD)" = "$(git -C "$ROOT" rev-parse '@{u}' 2>/dev/null)" ] \
   || { echo "ERROR: $ROOT HEAD is not its upstream; pull first"; exit 2; }
d1_guard "$REF"
[ "$(git -C "$REF" rev-parse HEAD)" = "$REF_SHA" ] \
   || { echo "ERROR: REF $REF is at $(git -C "$REF" rev-parse HEAD), expected $REF_SHA"; exit 2; }
git -C "$ROOT" merge-base --is-ancestor "$REF_SHA" HEAD \
   || { echo "ERROR: $REF_SHA is not an ancestor of $ROOT HEAD"; exit 2; }
NEW_SHA=$(git -C "$ROOT" rev-parse HEAD)
OUT=$SCRATCH/gpu_dyn/d1_gg_build_$SLURM_JOB_ID
mkdir -p "$OUT"
{ echo "job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ)"; echo "new $NEW_SHA ($ROOT)"; echo "ref $REF_SHA ($REF)"; } \
   | tee "$OUT/RUN_INFO.txt"

# --- 1. builds, both clones at once (72 cores each)
for d in "$ROOT" "$REF"; do
   ( cd "$d" && JOBS=70 bash build_vista.sh > "$OUT/build_$(basename "$d").log" 2>&1; echo $? > "$OUT/build_$(basename "$d").rc" ) &
done
wait
for d in "$ROOT" "$REF"; do
   rc=$(cat "$OUT/build_$(basename "$d").rc")
   echo "build_vista.sh in $d: exit $rc"; tail -3 "$OUT/build_$(basename "$d").log"
   [ "$rc" = 0 ] || { echo "D1 gg build: FAIL (build in $d)"; exit 2; }
done
echo "builds done after $(( $(date +%s) - T0 )) s"

d1_gg_modules
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')" | tee -a "$OUT/RUN_INFO.txt"
fail=0

# --- 2. meshes, in this checkout only (D0a's files are kept when they match)
cd "$ROOT/miniapps/seas" || exit 2
if ! d1_check_meshes "$ROOT/miniapps/seas" > /dev/null 2>&1; then
   for g in tpv205/mesh/tpv2053d_1000m tpv104/mesh/tpv104_1000m; do
      t=$(date +%s); gmsh -3 $g.geo -format msh22 -o $g.msh > "$OUT/gmsh_$(basename $g).log" 2>&1
      echo "$g.msh: gmsh exit $? in $(( $(date +%s) - t )) s"
   done
fi
d1_check_meshes "$ROOT/miniapps/seas" | tee "$OUT/MESH_SHA256"
[ "${PIPESTATUS[0]}" = 0 ] || fail=1

# --- 3. Nitsche test targets in both clones
mkdir -p "$OUT/tests"
for side in new ref; do
   d=$ROOT; [ $side = ref ] && d=$REF
   for t in test-nitsche-gates test-nitsche-pf test-nitsche-hook-spy test-nitsche-lsw test-nitsche-srw test-nitsche-lsw-inline; do
      L="$OUT/tests/${side}_$t.log"
      ( cd "$d/miniapps/seas" && make MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= $t ) > "$L" 2>&1
      echo "$t exit $? PASSED $(grep -v echo "$L" | grep -c PASSED) FAILED $(grep -v echo "$L" | grep -c FAILED)"
   done > "$OUT/tests/summary_$side.txt"
   sed "s/^/$side: /" "$OUT/tests/summary_$side.txt"
done
if ! cmp -s "$OUT/tests/summary_new.txt" "$OUT/tests/summary_ref.txt"; then
   echo "TESTS: FAIL (new and ref differ)"; fail=1
elif grep -qv "exit 0 .* FAILED 0$" "$OUT/tests/summary_new.txt"; then
   echo "TESTS: FAIL (a test exited nonzero or reported FAILED)"; fail=1
else
   echo "TESTS: PASS (same counts in both clones, 0 FAILED)"
fi

echo "D1 gg build: $([ $fail = 0 ] && echo PASS || echo FAIL); new $NEW_SHA, ref $REF_SHA; total $(( $(date +%s) - T0 )) s"
exit $fail
