#!/bin/bash
#SBATCH -J d2_gh_gate
#SBATCH -o d2_gh_gate_%j.out
#SBATCH -e d2_gh_gate_%j.err
#SBATCH -p gh
#SBATCH -t 04:00:00
#SBATCH -A EAR20006

# GPU dynamic-rupture plan D2 (plan rev2 §4.3), GPU gate of the device step:
# seas_nitsche_dyn_driver of this clone's CUDA build with --device cuda
# --matrix-free 1, so K u + K_eta v, the slip loads, the dashpot and the mass
# inverse run on the GPU and the fault on the host, fed by the gather.  It
# must pass jobs/horizon/check_dyn_gate.sh against the reference runs of a
# D1 gg check at the same NP (trunk a4c9dee, dense blocks) at TOL = 1e-8
# (plan rev2 §4.3 acceptance), for the four D1 cases (d1_env.sh) to 12 s.
#
#   REF_RUNS=$SCRATCH/gpu_dyn/d1_gg_check_r<NP>_<jobid>/runs \
#      sbatch -N NP -n NP jobs/vista/d2_gh_gate.sbatch
#   (NP = 1, 2, 4 GPUs, one rank per GH200 node; from miniapps/seas of
#   $WORK/seas-mfem-gpu-dyn-d2 after d2_gh_build.sbatch passed at the same
#   commit; REF_RUNS is the gg check at the same NP: 1082888, 1082889, 1082890)
# Each rank needs its own GPU: at most as many ranks per node as the node has
# GPUs.  On Horizon (4 GB200 per node) all three fit on one node,
# jobs/horizon/h_d2_gate.sbatch with sbatch -n NP.
#
# GPU dynamic-rupture plan D3 (plan rev2 §4.4 acceptance, "gate at 1e-8 on
# 1, 2 and 4 GPUs"): DEVICE_FAULT=1 adds --device-fault 1, so the fault runs
# on the GPU too, and a case then also fails when the summary does not say
# so ("fault on the device: N owned fault DOFs"); the output folder is
# d2_gh_gate_r<NP>_df1_<jobid>.
#
# MF_LAYOUT (default auto: element at p = 1, team2 at p = 2 on the GPU) is
# passed as --mf-layout.  A case fails when the driver's summary does not say
# the device step ran.  Both step times are printed (the reference is the gg
# CPU at the same NP).
# Cost (inferred; no GPU timing of the whole step exists yet): the fault on
# one Grace core per rank sets the step time; at most 4 h x NP nodes at 1 SU
# per node-hour.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
TOL=${TOL:-1e-8}
MF_LAYOUT=${MF_LAYOUT:-auto}
DEVICE_FAULT=${DEVICE_FAULT:-0}
case "$DEVICE_FAULT" in 0|1) ;; *) echo "ERROR: DEVICE_FAULT is 0 or 1, not '$DEVICE_FAULT'"; exit 2 ;; esac
DF_TAG=""; [ "$DEVICE_FAULT" = 0 ] || DF_TAG="_df1"
ROOT=$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)              # submitted from miniapps/seas
NP=${SLURM_NTASKS:?run with sbatch -N NODES -n NP}
NGPU=$(nvidia-smi -L | grep -c '^GPU ')
RPN=$(( (NP + ${SLURM_NNODES:-1} - 1) / ${SLURM_NNODES:-1} ))
[ "$RPN" -le "$NGPU" ] || { echo "ERROR: one rank per GPU: $NP tasks on $SLURM_NNODES node(s) of $NGPU GPU(s)"; exit 2; }
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
   || { echo "ERROR: $ROOT is not a CUDA build (d2_gh_build.sbatch)"; exit 2; }
d1_check_build "$ROOT" BUILD_COMMIT_CUDA
d1_check_meshes "$ROOT/miniapps/seas" || exit 2
NEW_SHA=$(git -C "$ROOT" rev-parse HEAD)
OUT=$SCRATCH/gpu_dyn/d2_gh_gate_r${NP}${DF_TAG}_$SLURM_JOB_ID
mkdir -p "$OUT"
d1_gh_modules
{ echo "job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ), NP $NP, TOL $TOL, --mf-layout $MF_LAYOUT, --device-fault $DEVICE_FAULT"
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
            --output-dir "$O" --log-every 1000 --device cuda --matrix-free 1 --mf-layout "$MF_LAYOUT" \
            --device-fault "$DEVICE_FAULT" > "$O/driver.log" 2>&1
      rc=$?
      echo $rc > "$O/exit_code"
      echo "$CASE cuda: exit $rc in $(( $(date +%s) - t )) s; $(grep -a "=== done" "$O/driver.log" || echo 'no done line')"
      grep -a "time in ode.Step" "$O/driver.log" | sed 's/^/   /'
      grep -a "time in ode.Step" "$G/driver.log" | sed 's/^ */   ref: /'
      if ! grep -a -q "device step on" "$O/driver.log"; then
         echo "RUN $CASE cuda: FAIL (the summary does not say the device step ran)"; fail=1
      elif [ "$DEVICE_FAULT" = 1 ] && ! grep -a -q "fault on the device: " "$O/driver.log"; then
         echo "RUN $CASE cuda: FAIL (the summary does not say the fault ran on the device)"; fail=1
      elif bash jobs/horizon/check_dyn_gate.sh "$O" "$PROB" "$G" "$PROB" "$TOL" "$O/driver.log" "$rc" \
            > "$O/gate.log" 2>&1; then
         echo "RUN $CASE cuda vs gg ref, $NP GPU(s): PASS"
      else
         echo "RUN $CASE cuda vs gg ref, $NP GPU(s): FAIL"; fail=1
      fi
      grep -E 'max relative|byte-identical files|FAIL' "$O/gate.log" | sed 's/^/   /'
   done
done

echo "D2 gh gate, $NP GPU(s), --device-fault $DEVICE_FAULT: $([ $fail = 0 ] && echo PASS || echo FAIL); new $NEW_SHA; runs in $OUT/runs; total $(( $(date +%s) - T0 )) s"
exit $fail
