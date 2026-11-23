#!/bin/bash
#SBATCH -J d0a_nitsche_check
#SBATCH -o d0a_nitsche_check_%j.out
#SBATCH -e d0a_nitsche_check_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 03:00:00
#SBATCH -A EAR20006

# GPU dynamic-rupture plan D0a (gpu-dyn-nitsche-plan-rev2 §4.1), Nitsche side:
# the merge of gpu-qd-bp5 into gpu-dyn-nitsche must leave the Nitsche driver's
# results bitwise unchanged against trunk a4c9dee.  (The quasi-dynamic side is
# checked by the merged branch's own cpu_clone_build.sbatch + g2_cpu_gate.sbatch
# in a second clone.)
#
# Two clones, both built by the trunk's build_vista.sh (same modules, same flags):
#   this checkout   gpu-dyn-nitsche, at the commit its upstream points to
#   REF             trunk a4c9dee, detached (default $WORK/seas-mfem-gpu-dyn-ref)
# Network work on the login node, compiling here:
#   cd $WORK
#   git clone --recurse-submodules -b gpu-dyn-nitsche \
#       https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-gpu-dyn-nit
#   git clone --recurse-submodules https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-gpu-dyn-ref
#   git -C seas-mfem-gpu-dyn-ref checkout --detach a4c9dee0aa96a019792e48cce864fc68054aab26
#   cd seas-mfem-gpu-dyn-nit/miniapps/seas && sbatch jobs/vista/d0a_nitsche_check.sbatch
# Steps:
#   0. the repo CLAUDE.md guard in both clones; REF must be exactly a4c9dee and
#      an ancestor of this checkout
#   1. build_vista.sh in both clones at once (MFEM + seas_nitsche_dyn_driver + the
#      six Nitsche test programs)
#   2. the TPV205 and TPV104 1000 m meshes from the committed .geo files (TACC
#      gmsh, msh 2.2) in this checkout, sha256 recorded; both binaries read the
#      same mesh files
#   3. the six Nitsche test targets in both clones: exit codes and PASSED/FAILED
#      counts must match, with 0 FAILED
#   4. TPV205 and TPV104 on the 1000 m meshes, p = 1 and p = 2, the committed
#      decks, tfinal 1.2 s, 16 ranks, both binaries: same exit code, the same set
#      of output files, and every station file, station-copy file, station snap
#      list and rupture-time file byte-identical.  fault.vtkhdf and the logs are
#      not compared (HDF5 stores object times; logs hold timings).
#      The TPV104 1000 m mesh is the coarse debug mesh: a bitwise check, not a
#      benchmark run.
# PASS needs all of 0-4.  Cost: 1 node x <= 3 h x 0.33 SU = <= 1 SU.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REF_SHA=a4c9dee0aa96a019792e48cce864fc68054aab26
ROOT=$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)              # submitted from miniapps/seas
REF=${REF:-$WORK/seas-mfem-gpu-dyn-ref}
[ -f "$ROOT/build_vista.sh" ] || { echo "ERROR: submit from miniapps/seas ($ROOT has no build_vista.sh)"; exit 2; }
[ -f "$REF/build_vista.sh" ]  || { echo "ERROR: REF $REF has no build_vista.sh"; exit 2; }
T0=$(date +%s)

# --- 0. guards
guard () {
   local d=$1 i U
   cd "$d" || exit 2
   for i in 1 2 3 4; do
      git fetch -q origin && break
      [ $i = 4 ] && { echo "ERROR: git fetch failed in $d"; exit 2; }
      sleep $(( 5 + RANDOM % 20 ))
   done
   [ -z "$(git status --porcelain --untracked-files=no)" ] \
      || { echo "ERROR: uncommitted changes in $d"; git status --short; exit 2; }
   [ -n "$(git branch -r --contains HEAD)" ] \
      || { echo "ERROR: HEAD $(git rev-parse --short HEAD) in $d is not on GitHub"; exit 2; }
   U=$(git ls-files --others --exclude-standard -- '*.cpp' '*.hpp' '*.inl' Makefile '*.toml' '*.sbatch')
   [ -z "$U" ] || { echo "ERROR: untracked sources in $d:"; echo "$U"; exit 2; }
   echo "commit $(git rev-parse HEAD) in $d"
}
guard "$ROOT"
[ "$(git -C "$ROOT" rev-parse HEAD)" = "$(git -C "$ROOT" rev-parse '@{u}' 2>/dev/null)" ] \
   || { echo "ERROR: $ROOT HEAD is not its upstream; pull first"; exit 2; }
guard "$REF"
[ "$(git -C "$REF" rev-parse HEAD)" = "$REF_SHA" ] \
   || { echo "ERROR: REF $REF is at $(git -C "$REF" rev-parse HEAD), expected $REF_SHA"; exit 2; }
git -C "$ROOT" merge-base --is-ancestor "$REF_SHA" HEAD \
   || { echo "ERROR: $REF_SHA is not an ancestor of $ROOT HEAD"; exit 2; }
NEW_SHA=$(git -C "$ROOT" rev-parse HEAD)
OUT=$SCRATCH/gpu_dyn/d0a_nitsche_check_$SLURM_JOB_ID
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
   [ "$rc" = 0 ] || { echo "D0a Nitsche check: FAIL (build in $d)"; exit 2; }
done
echo "builds done after $(( $(date +%s) - T0 )) s"

module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6 gmsh/4.15.0
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')" | tee -a "$OUT/RUN_INFO.txt"
fail=0

# --- 2. meshes, in this checkout only
cd "$ROOT/miniapps/seas"
for g in tpv205/mesh/tpv2053d_1000m tpv104/mesh/tpv104_1000m; do
   t=$(date +%s); gmsh -3 $g.geo -format msh22 -o $g.msh > "$OUT/gmsh_$(basename $g).log" 2>&1
   echo "$g.msh: gmsh exit $? in $(( $(date +%s) - t )) s"
done
sha256sum tpv205/mesh/tpv2053d_1000m.msh tpv104/mesh/tpv104_1000m.msh | tee "$OUT/MESH_SHA256"

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

# --- 4. bitwise runs
for c in "tpv205 tpv205/configs/tpv205_nitsche.toml tpv205/mesh/tpv2053d_1000m.msh" \
         "tpv104 tpv104/configs/tpv104_nitsche.toml tpv104/mesh/tpv104_1000m.msh"; do
   set -- $c; PROB=$1; DECK=$2; MESH=$ROOT/miniapps/seas/$3
   cmp -s "$ROOT/miniapps/seas/$DECK" "$REF/miniapps/seas/$DECK" \
      || { echo "RUN $PROB: FAIL (the deck differs between the clones)"; fail=1; continue; }
   for P in 1 2; do
      CASE=${PROB}_1000m_p$P
      for side in new ref; do
         d=$ROOT; [ $side = ref ] && d=$REF
         O="$OUT/runs/$CASE/$side"; mkdir -p "$O"
         t=$(date +%s)
         ( cd "$d/miniapps/seas" && ibrun -n 16 -o 0 ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" \
              --order $P --tfinal 1.2 --output-dir "$O" --log-every 1000 ) > "$O/driver.log" 2>&1
         echo $? > "$O/exit_code"
         echo "$CASE $side: exit $(cat "$O/exit_code") in $(( $(date +%s) - t )) s; $(grep -a "=== done" "$O/driver.log" || echo 'no done line')"
      done
      N="$OUT/runs/$CASE/new"; R="$OUT/runs/$CASE/ref"; ok=1
      cmp -s "$N/exit_code" "$R/exit_code" || { echo "  exit codes differ"; ok=0; }
      find "$N" -maxdepth 1 -type f -printf '%f\n' | grep -v -E '\.vtkhdf$|\.log$|^exit_code$' | sort > "$OUT/runs/$CASE/files_new.txt"
      find "$R" -maxdepth 1 -type f -printf '%f\n' | grep -v -E '\.vtkhdf$|\.log$|^exit_code$' | sort > "$OUT/runs/$CASE/files_ref.txt"
      cmp -s "$OUT/runs/$CASE/files_new.txt" "$OUT/runs/$CASE/files_ref.txt" || { echo "  file sets differ"; ok=0; }
      n=0
      while read -r f; do
         [ -f "$R/$f" ] || continue
         n=$((n + 1))
         cmp -s "$N/$f" "$R/$f" || { echo "  differs: $f"; ok=0; }
      done < "$OUT/runs/$CASE/files_new.txt"
      if [ $ok = 1 ] && [ $n -gt 0 ]; then
         echo "RUN $CASE: PASS ($n files byte-identical, exit $(cat "$N/exit_code"))"
      else
         echo "RUN $CASE: FAIL ($n files compared)"; fail=1
      fi
   done
done

echo "D0a Nitsche check: $([ $fail = 0 ] && echo PASS || echo FAIL); new $NEW_SHA, ref $REF_SHA; total $(( $(date +%s) - T0 )) s"
exit $fail
