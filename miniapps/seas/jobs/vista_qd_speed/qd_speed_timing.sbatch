#!/bin/bash
#SBATCH -J qd_speed_timing
#SBATCH -o qd_speed_timing_%j.out
#SBATCH -e qd_speed_timing_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 01:30:00
#SBATCH -A EAR20006

# QD speed items 1 and 2 (plan rev 2): time the production BP5 QD deck at
# one rank layout.  The layout comes from the sbatch line (-N, -n, -c,
# --ntasks-per-node); jobs/vista_qd_speed/qd_speed_submit.sh holds every
# layout of the sweep, so submit through it:
#
#   cd miniapps/seas
#   bash jobs/vista_qd_speed/qd_speed_submit.sh          # dry run
#   SUBMIT=1 bash jobs/vista_qd_speed/qd_speed_submit.sh
#
# Per job:
#   1. gemv_bw (jobs/vista_qd_speed/gemv_bw.c, compiled here): where the
#      ranks sit (ranks per node and socket, bound or not) and the memory
#      bandwidth this layout gets.  An unbalanced layout stops the job
#      before any timing run;
#   2. for each mode in MODES (default "loops default loops default" on one
#      full node, else "default"), 1,000 accepted steps of the v63 deck at
#      p=ORDER (default 1) with --stage-timers 100, as G0's
#      jobs/vista/g0_timing.sbatch, without checkpoints:
#        loops    --qd-operators loops (the face loops, the old default);
#        default  no --qd-operators flag: the new default, which must say
#                 "qd-operators: matrix (default)" in its log;
#        matrix   --qd-operators matrix;
#        green    --qd-green --qd-green-cache GREEN_CACHE (QD speed item
#                 3): reads G from GREEN_CACHE, or builds it and writes it
#                 there when the file does not exist.  Two jobs must not
#                 build the same file at once: qd_speed_submit.sh makes the
#                 other green jobs wait for the first (afterok);
#      repeating a mode measures the run-to-run spread;
#   3. timing_summary.txt: one line per run (s/step, solve ms/call and the
#      other parts) from the end-of-run [STAGE-TIMERS] table, then the
#      [QD-GREEN] build or read lines.  In a green run "solve" is the
#      G product.
#
# Ranks per node below 144 leave cores idle: the BLAS (NVPL *_seq) and
# MUMPS of the G0 build are sequential, so threads per rank would not be
# used.  The sweep measures what fewer ranks per node do to the solve with
# the same build; threaded libraries would need a separate deps build.
#
# GitHub is the version control: the checkout must be clean and at the
# commit its upstream branch points to, or the job stops.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"
g0_check_binary

cd "${SEAS_MFEM_ROOT}"
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
   echo "ERROR: checkout has uncommitted changes; run only pushed commits." >&2
   git status --short --untracked-files=no >&2
   exit 1
fi
if [ "$(git rev-parse HEAD)" != "$(git rev-parse '@{u}' 2>/dev/null)" ]; then
   echo "ERROR: HEAD $(git rev-parse --short HEAD) is not its upstream" \
        "($(git rev-parse --short '@{u}' 2>/dev/null || echo none)); pull first." >&2
   exit 1
fi
echo "Checkout: ${SEAS_MFEM_ROOT} at $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD))"

ORDER="${ORDER:-1}"
STEPS="${STEPS:-1000}"
NNODES="${SLURM_JOB_NUM_NODES:-1}"
NRANKS="${SLURM_NTASKS:-144}"
CPT="${SLURM_CPUS_PER_TASK:-1}"
if [ -z "${MODES:-}" ]; then
   if [ "${NNODES}" = 1 ] && [ "${NRANKS}" = 144 ]; then
      MODES="loops default loops default"
   else
      MODES="default"
   fi
fi
GREEN_CACHE="${GREEN_CACHE:-${SCRATCH:-${WORK}}/seas-qd-green/bp5_prod_p${ORDER}.bin}"
GOLD_README_KNOBS+=(ORDER STEPS MODES GREEN_CACHE)

cd "${SEAS_DIR}"
OUT_ROOT="bp5/results_qd_speed/timing_n${NNODES}_r${NRANKS}_c${CPT}_p${ORDER}_job${SLURM_JOB_ID}"
mkdir -p "${OUT_ROOT}"
g0_provenance "${OUT_ROOT}"
{
   echo "layout: ${NNODES} node(s), ${NRANKS} ranks, ${CPT} CPUs per rank," \
        "SLURM_TASKS_PER_NODE=${SLURM_TASKS_PER_NODE:-?}"
   echo "modes: ${MODES}; order ${ORDER}; ${STEPS} steps"
   echo "green cache: ${GREEN_CACHE}"
} >> "${OUT_ROOT}/provenance.txt"

# ---- 1. Layout and bandwidth ----
mpicc -O3 -mcpu=native jobs/vista_qd_speed/gemv_bw.c -o "${OUT_ROOT}/gemv_bw" || exit 1
ibrun "${OUT_ROOT}/gemv_bw" 8 > "${OUT_ROOT}/gemv_bw.txt" 2>&1
rc=$?
cat "${OUT_ROOT}/gemv_bw.txt"
if [ "${rc}" != 0 ]; then
   echo "ERROR: gemv_bw exit ${rc} (3 = unbalanced layout); no timing runs." >&2
   exit 1
fi

# ---- 2. Timing runs ----
overall=0
k=0
for MODE in ${MODES}; do
   k=$((k + 1))
   TAG="run${k}_${MODE}"
   OUT="${OUT_ROOT}/${TAG}"
   PFX="qd_speed_${TAG}"
   mkdir -p "${OUT}"
   case "${MODE}" in
      loops|matrix) OPS=(--qd-operators "${MODE}") ;;
      default)      OPS=() ;;
      green)        mkdir -p "$(dirname "${GREEN_CACHE}")"
                    OPS=(--qd-green --qd-green-cache "${GREEN_CACHE}") ;;
      *) echo "ERROR: unknown mode ${MODE}" >&2; overall=1; continue ;;
   esac
   echo "=== ${TAG}: start $(date -u +%T) ==="
   ibrun ./seas_bp5_full --mesh "${G0_MESH_PROD}" "${G0_V63_FLAGS[@]}" \
         --order "${ORDER}" "${OPS[@]}" --max-steps "${STEPS}" \
         --stage-timers 100 --checkpoint-interval 0 \
         --output-dir "${OUT}" --output-prefix "${PFX}" > "${OUT}/run.log" 2>&1
   rc=$?
   echo "${rc}" > "${OUT}/exit_code"
   echo "=== ${TAG}: exit ${rc} at $(date -u +%T) ==="
   [ "${rc}" = 0 ] || overall=1
   if [ "${MODE}" = default ] && ! grep -q "qd-operators: matrix (default)" "${OUT}/run.log"; then
      echo "ERROR: ${TAG}: the log does not say 'qd-operators: matrix (default)'." >&2
      overall=1
   fi
   sed -n '/^\[STAGE-TIMERS\] end of run/,$p' "${OUT}/run.log" \
      | grep '^\[STAGE-TIMERS\]' > "${OUT}/stage_timers_final.txt"
   grep '\[QD-OPERATORS\]' "${OUT}/run.log" > "${OUT}/qd_operators.txt"
   if [ "${MODE}" = green ]; then
      grep '\[QD-GREEN\]' "${OUT}/run.log" > "${OUT}/qd_green.txt"
      if ! grep -q '\[QD-GREEN\] check against the solve path' "${OUT}/qd_green.txt"; then
         echo "ERROR: ${TAG}: no '[QD-GREEN] check against the solve path' line." >&2
         overall=1
      fi
   fi
done

# ---- 3. Summary ----
python3 jobs/vista_qd_speed/qd_speed_summary.py "${OUT_ROOT}" \
   | tee "${OUT_ROOT}/timing_summary.txt"
echo "QD speed timing overall: $([ "${overall}" = 0 ] && echo PASS || echo FAIL)"
exit "${overall}"
