#!/bin/bash
#SBATCH -J d2_gg_mf_gate
#SBATCH -o d2_gg_mf_gate_%j.out
#SBATCH -e d2_gg_mf_gate_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -t 08:00:00
#SBATCH -A EAR20006

# GPU dynamic-rupture plan D2 (plan rev2 §4.3), CPU gate of the matrix-free
# product: the Nitsche driver of this checkout with --matrix-free 1 (K u +
# K_eta v from NitscheSipgOperator; K and K_eta never assembled) on the gg
# CPU must pass jobs/horizon/check_dyn_gate.sh against the reference runs of
# a D1 gg check at the same NP (trunk a4c9dee, dense blocks) at TOL = 1e-10
# (plan rev2 §4.3 acceptance), for the four D1 cases (d1_env.sh) to 12 s.
#
#   REF_RUNS=$SCRATCH/gpu_dyn/d1_gg_check_r4_1082890/runs \
#      sbatch -n 4 jobs/vista/d2_gg_mf_gate.sbatch
#   (from miniapps/seas of $WORK/seas-mfem-gpu-dyn-nit, after
#   d1_gg_build.sbatch passed at the same commit; REF_RUNS is the gg check at
#   the same NP, whose ref runs are the reference)
#
# MF_LAYOUT (default auto: the element layout on the host) is passed as
# --mf-layout.  Each run's peak memory is not measured; the driver's summary
# line reports the matrix-free data in MB.
# Cost (inferred from the D1 gg check at NP = 4, about 2 to 3 h for the four
# cases): 1 node x <= 8 h x 0.33 SU = <= 2.7 SU.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
TOL=${TOL:-1e-10}
MF_LAYOUT=${MF_LAYOUT:-auto}
ROOT=$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)              # submitted from miniapps/seas
NP=${SLURM_NTASKS:?run with sbatch -n NP}
[ -f "$ROOT/build_vista.sh" ] || { echo "ERROR: submit from miniapps/seas ($ROOT has no build_vista.sh)"; exit 2; }
[ -n "${REF_RUNS:-}" ] && [ -f "$REF_RUNS/../RUN_INFO.txt" ] \
   || { echo "ERROR: set REF_RUNS to a d1_gg_check job's runs directory"; exit 2; }
grep -q ", NP $NP," "$REF_RUNS/../RUN_INFO.txt" \
   || { echo "ERROR: $REF_RUNS is not a gg check at NP = $NP:"; head -1 "$REF_RUNS/../RUN_INFO.txt"; exit 2; }
# shellcheck source=d1_env.sh
source "$ROOT/miniapps/seas/jobs/vista/d1_env.sh"
T0=$(date +%s)

d1_guard "$ROOT"
[ "$(git -C "$ROOT" rev-parse HEAD)" = "$(git -C "$ROOT" rev-parse '@{u}' 2>/dev/null)" ] \
   || { echo "ERROR: $ROOT HEAD is not its upstream; pull first"; exit 2; }
d1_check_build "$ROOT" BUILD_COMMIT_VISTA
d1_check_meshes "$ROOT/miniapps/seas" || exit 2
NEW_SHA=$(git -C "$ROOT" rev-parse HEAD)
OUT=$SCRATCH/gpu_dyn/d2_gg_mf_gate_r${NP}_$SLURM_JOB_ID
mkdir -p "$OUT"
d1_gg_modules
{ echo "job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ), NP $NP, TOL $TOL, --mf-layout $MF_LAYOUT"
  echo "new $NEW_SHA ($ROOT)"; echo "reference $REF_RUNS"; head -3 "$REF_RUNS/../RUN_INFO.txt"
  echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"; } | tee "$OUT/RUN_INFO.txt"

fail=0
cd "$ROOT/miniapps/seas" || exit 2
for c in "${D1_CASES[@]}"; do
   read -r PROB DECK MESH <<< "$c"
   for P in "${D1_ORDERS[@]}"; do
      CASE=${PROB}_1000m_p$P
      O="$OUT/runs/$CASE/mf"; G="$REF_RUNS/$CASE/ref"
      mkdir -p "$O"
      [ -f "$G/exit_code" ] || { echo "RUN $CASE: FAIL (no reference in $G)"; fail=1; continue; }
      t=$(date +%s)
      ibrun -n "$NP" -o 0 ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order "$P" \
            --output-dir "$O" --log-every 1000 --device cpu --matrix-free 1 --mf-layout "$MF_LAYOUT" \
            > "$O/driver.log" 2>&1
      rc=$?
      echo $rc > "$O/exit_code"
      echo "$CASE mf: exit $rc in $(( $(date +%s) - t )) s; $(grep -a "=== done" "$O/driver.log" || echo 'no done line')"
      grep -a "time in ode.Step" "$O/driver.log" | sed 's/^/   /'
      grep -a "time in ode.Step" "$G/driver.log" | sed 's/^ */   ref: /'
      if ! grep -a -q "matrix-free on" "$O/driver.log"; then
         echo "RUN $CASE mf: FAIL (the summary does not say the matrix-free product was used)"; fail=1
      elif bash jobs/horizon/check_dyn_gate.sh "$O" "$PROB" "$G" "$PROB" "$TOL" "$O/driver.log" "$rc" \
            > "$O/gate.log" 2>&1; then
         echo "RUN $CASE mf vs gg ref, $NP ranks: PASS"
      else
         echo "RUN $CASE mf vs gg ref, $NP ranks: FAIL"; fail=1
      fi
      grep -E 'max relative|byte-identical files|FAIL' "$O/gate.log" | sed 's/^/   /'
   done
done

echo "D2 gg matrix-free gate, $NP ranks: $([ $fail = 0 ] && echo PASS || echo FAIL); new $NEW_SHA; runs in $OUT/runs; total $(( $(date +%s) - T0 )) s"
exit $fail
