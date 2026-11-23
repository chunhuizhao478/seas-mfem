#!/bin/bash
#SBATCH -J d1_gh_check
#SBATCH -o d1_gh_check_%j.out
#SBATCH -e d1_gh_check_%j.err
#SBATCH -p gh
#SBATCH -t 24:00:00
#SBATCH -A EAR20006

# GPU dynamic-rupture plan D1 (plan rev2 §4.2), gh gate: the unported run.
# seas_nitsche_dyn_driver from the CUDA build with --device cuda and nothing
# ported, so the ODE solver's vector updates run on the GPU and everything
# else on the host through the syncs (hypre stays on the host).  A stale host
# or device copy shows up as an error of order one, not as rounding.
#
#   REF_RUNS=$SCRATCH/gpu_dyn/d1_gg_check_r<NP>_<jobid>/runs \
#      sbatch -N NP -n NP jobs/vista/d1_gh_check.sbatch
#   (NP = 1, 2, 4 GPUs, one rank per GH200 node; from miniapps/seas of
#   $WORK/seas-mfem-gpu-dyn-cuda after d1_cuda_build.sbatch at the same commit;
#   REF_RUNS is the gg check at the same NP, whose ref runs are the reference)
#
# Each of the four D1 cases (d1_env.sh) runs to 12 s and passes when
# jobs/horizon/check_dyn_gate.sh passes against REF_RUNS/<case>/ref at
# TOL = 1e-8 (plan rev2 §4.1: tol_phase for any GPU run; no D0 floors exist,
# D0c was dropped).  The errors it reports are the gg-to-gh floor for the
# later GPU gates.
# Cost (inferred: the volume work stays on one Grace core per rank): per run
# about the gg run at the same NP, so about 11 h at NP = 1, 5.5 h at 2, 3 h at
# 4; 1 SU per node-hour.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
TOL=${TOL:-1e-8}
ROOT=$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)              # submitted from miniapps/seas
NP=${SLURM_NTASKS:?run with sbatch -N NP -n NP}
[ "$NP" = "${SLURM_NNODES:-0}" ] || { echo "ERROR: one rank per GPU: sbatch -N NP -n NP (got $SLURM_NNODES nodes, $NP tasks)"; exit 2; }
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
grep -q '^MFEM_USE_CUDA *= YES' "$ROOT/config/config.mk" \
   || { echo "ERROR: $ROOT is not a CUDA build (d1_cuda_build.sbatch)"; exit 2; }
d1_check_build "$ROOT" BUILD_COMMIT_CUDA
d1_check_meshes "$ROOT/miniapps/seas" || exit 2
NEW_SHA=$(git -C "$ROOT" rev-parse HEAD)
OUT=$SCRATCH/gpu_dyn/d1_gh_check_r${NP}_$SLURM_JOB_ID
mkdir -p "$OUT"
d1_gh_modules
{ echo "job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ), NP $NP, TOL $TOL"
  echo "new $NEW_SHA ($ROOT)"; echo "reference $REF_RUNS"; head -3 "$REF_RUNS/../RUN_INFO.txt"
  echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
  nvidia-smi --query-gpu=name,compute_mode,driver_version --format=csv,noheader; } | tee "$OUT/RUN_INFO.txt"

fail=0
cd "$ROOT/miniapps/seas" || exit 2
for c in "${D1_CASES[@]}"; do
   read -r PROB DECK MESH <<< "$c"
   for P in "${D1_ORDERS[@]}"; do
      CASE=${PROB}_1000m_p$P
      O="$OUT/runs/$CASE/cuda"; G="$REF_RUNS/$CASE/ref"
      mkdir -p "$O"
      [ -f "$G/exit_code" ] || { echo "RUN $CASE: FAIL (no reference in $G)"; fail=1; continue; }
      t=$(date +%s)
      ibrun -n "$NP" ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order "$P" \
            --output-dir "$O" --log-every 1000 --device cuda > "$O/driver.log" 2>&1
      rc=$?
      echo $rc > "$O/exit_code"
      echo "$CASE cuda: exit $rc in $(( $(date +%s) - t )) s; $(grep -a "=== done" "$O/driver.log" || echo 'no done line')"
      if bash jobs/horizon/check_dyn_gate.sh "$O" "$PROB" "$G" "$PROB" "$TOL" "$O/driver.log" "$rc" \
            > "$O/gate.log" 2>&1; then
         echo "RUN $CASE cuda vs gg ref, $NP GPU(s): PASS"
      else
         echo "RUN $CASE cuda vs gg ref, $NP GPU(s): FAIL"; fail=1
      fi
      grep -E 'max|error|PASS|FAIL' "$O/gate.log" | tail -6
   done
done

echo "D1 gh check, $NP GPU(s): $([ $fail = 0 ] && echo PASS || echo FAIL); new $NEW_SHA; runs in $OUT/runs; total $(( $(date +%s) - T0 )) s"
exit $fail
