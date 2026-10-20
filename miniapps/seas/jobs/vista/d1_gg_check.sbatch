#!/bin/bash
#SBATCH -J d1_gg_check
#SBATCH -o d1_gg_check_%j.out
#SBATCH -e d1_gg_check_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -A EAR20006

# GPU dynamic-rupture plan D1 (device plumbing), gg gate (plan rev2 §4.2): on
# NP ranks, the Nitsche driver of this checkout with --device cpu must give
# results bitwise equal to trunk a4c9dee (D0a's reference), and with MFEM's
# debug device too.  The debug device keeps a second "device" copy of every
# vector used on the device and protects the host copy whenever the device
# copy is taken, so a host access that skips a sync aborts the run.
#
#   sbatch -t 36:00:00 -n 1 jobs/vista/d1_gg_check.sbatch
#   sbatch -n 2 jobs/vista/d1_gg_check.sbatch
#   sbatch -n 4 jobs/vista/d1_gg_check.sbatch
#   D1_DEBUG=1 sbatch -n 8 jobs/vista/d1_gg_check.sbatch
#   (from miniapps/seas of $WORK/seas-mfem-gpu-dyn-nit, after
#   d1_gg_build.sbatch passed at the same commit)
#
# SIDES (default "ref cpu", "ref cpu debug" when D1_DEBUG=1):
#   ref    the REF clone's binary (a4c9dee), no device flag
#   cpu    this checkout, --device cpu
#   debug  this checkout, --device debug, at NP >= 8 only.  The debug device
#          maps every allocation on its own, and freeing them at exit splits
#          the maps: in the cloud checks a rank peaked at 0.86 to 0.93 Linux
#          memory maps per tet it holds (p = 1 and 2).  Past vm.max_map_count
#          (65530 by default) the exit fails with "Dealloc error!", as at NP = 2
#          on the TPV205 mesh (about 84,000 per rank); the TPV104 mesh (216,015
#          tets) needs about 25,000 per rank at NP = 8.  The debug device also
#          ran 3 to 7 times slower than cpu there.
# Each of the four D1 cases (d1_env.sh) runs once per side, one run at a time,
# to the decks' tfinal of 12 s.  Every non-ref side must pass
# jobs/horizon/check_dyn_gate.sh against ref at TOL = 1e-12 (plan rev2 §4.1:
# tol_phase for host syncs on the CPU); whether it is also bitwise equal to
# ref (d1_bitwise) is reported on its own line, and is expected at every NP
# (on the cpu device the syncs change no arithmetic).  The ref runs are the gg
# reference the gh gate compares against (d1_gh_check.sbatch,
# REF_RUNS=$OUT/runs).
# Cost (inferred from D0a's 16-rank 1.2 s runs, 226 s for the four cases, and
# 10x the steps): per side about 11 h at NP = 1, 5.5 h at 2, 3 h at 4, 1.5 h
# at 8, the debug side about 3x that; 0.33 SU per node-hour.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REF_SHA=a4c9dee0aa96a019792e48cce864fc68054aab26
TOL=${TOL:-1e-12}
ROOT=$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)              # submitted from miniapps/seas
REF=${REF:-$WORK/seas-mfem-gpu-dyn-ref}
NP=${SLURM_NTASKS:?run with sbatch -n NP}
SIDES=${SIDES:-ref cpu$([ "${D1_DEBUG:-0}" = 1 ] && echo " debug")}
[ -f "$ROOT/build_vista.sh" ] || { echo "ERROR: submit from miniapps/seas ($ROOT has no build_vista.sh)"; exit 2; }
# shellcheck source=d1_env.sh
source "$ROOT/miniapps/seas/jobs/vista/d1_env.sh"
case " $SIDES " in *" debug "*) [ "$NP" -ge 8 ] || { echo "ERROR: the debug side needs NP >= 8 (Linux memory maps per rank)"; exit 2; } ;; esac
case " $SIDES " in *" ref "*) ;; *) echo "ERROR: SIDES must include ref"; exit 2 ;; esac
T0=$(date +%s)

# --- guards: both clones clean and pushed, REF at a4c9dee, both binaries built
#     at HEAD by d1_gg_build.sbatch, D0a's meshes
d1_guard "$ROOT"
[ "$(git -C "$ROOT" rev-parse HEAD)" = "$(git -C "$ROOT" rev-parse '@{u}' 2>/dev/null)" ] \
   || { echo "ERROR: $ROOT HEAD is not its upstream; pull first"; exit 2; }
d1_guard "$REF"
[ "$(git -C "$REF" rev-parse HEAD)" = "$REF_SHA" ] \
   || { echo "ERROR: REF $REF is at $(git -C "$REF" rev-parse HEAD), expected $REF_SHA"; exit 2; }
d1_check_build "$ROOT" BUILD_COMMIT_VISTA
d1_check_build "$REF" BUILD_COMMIT_VISTA
d1_check_meshes "$ROOT/miniapps/seas" || exit 2
NEW_SHA=$(git -C "$ROOT" rev-parse HEAD)
OUT=$SCRATCH/gpu_dyn/d1_gg_check_r${NP}_$SLURM_JOB_ID
mkdir -p "$OUT"
d1_gg_modules
{ echo "job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ), NP $NP, TOL $TOL, sides: $SIDES"
  echo "new $NEW_SHA ($ROOT)"; echo "ref $REF_SHA ($REF)"
  echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
  echo "vm.max_map_count $(cat /proc/sys/vm/max_map_count)"; } | tee "$OUT/RUN_INFO.txt"

# --- runs, then the bitwise checks
fail=0
cd "$ROOT/miniapps/seas" || exit 2
for c in "${D1_CASES[@]}"; do
   read -r PROB DECK MESH <<< "$c"
   cmp -s "$ROOT/miniapps/seas/$DECK" "$REF/miniapps/seas/$DECK" \
      || { echo "RUN $PROB: FAIL (the deck differs between the clones)"; fail=1; continue; }
   for P in "${D1_ORDERS[@]}"; do
      CASE=${PROB}_1000m_p$P
      for side in $SIDES; do
         d=$ROOT; dev=(--device "$side"); [ "$side" = ref ] && { d=$REF; dev=(); }
         O="$OUT/runs/$CASE/$side"; mkdir -p "$O"
         t=$(date +%s)
         ( cd "$d/miniapps/seas" && ibrun -n "$NP" -o 0 ./seas_nitsche_dyn_driver --config "$DECK" \
              --mesh "$ROOT/miniapps/seas/$MESH" --order "$P" --output-dir "$O" --log-every 1000 "${dev[@]}" ) \
            > "$O/driver.log" 2>&1
         echo $? > "$O/exit_code"
         echo "$CASE $side: exit $(cat "$O/exit_code") in $(( $(date +%s) - t )) s; $(grep -a "=== done" "$O/driver.log" || echo 'no done line')"
      done
      for side in $SIDES; do
         [ "$side" = ref ] && continue
         N="$OUT/runs/$CASE/$side"; R="$OUT/runs/$CASE/ref"
         echo "bitwise $CASE $side vs ref, $NP ranks: $(d1_bitwise "$N" "$R")"
         if bash jobs/horizon/check_dyn_gate.sh "$N" "$PROB" "$R" "$PROB" "$TOL" "$N/driver.log" \
               "$(cat "$N/exit_code")" > "$N/gate.log" 2>&1; then
            echo "RUN $CASE $side vs ref, $NP ranks: PASS"
         else
            echo "RUN $CASE $side vs ref, $NP ranks: FAIL"; fail=1
         fi
         grep -E 'max relative|byte-identical files|FAIL' "$N/gate.log" | sed 's/^/   /'
      done
   done
done

echo "D1 gg check, $NP ranks: $([ $fail = 0 ] && echo PASS || echo FAIL); new $NEW_SHA, ref $REF_SHA; runs in $OUT/runs; total $(( $(date +%s) - T0 )) s"
exit $fail
