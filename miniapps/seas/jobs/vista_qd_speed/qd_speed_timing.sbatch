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
#      full node, else "default"), RUN_STEPS (default 1,000) accepted steps
#      of the v63 deck at
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
#        hgreen   G kept compressed (QD speed item 4): --qd-green
#                 --qd-green-cache GREEN_CACHE --qd-green-compress HEPS
#                 --qd-green-hcache HCACHE --qd-green-compress-mem HMEM.
#                 Reads the H-matrix from HCACHE; when that file does not
#                 exist, reads dense G from GREEN_CACHE, compresses it and
#                 writes HCACHE.  Again only one job may write it: the
#                 submit script makes the others wait (afterok);
#        hcount, hbytes  hgreen with --qd-green-layout count or bytes (QD
#                 speed item C), reading HCACHE only: the job stops before
#                 its first run when HCACHE does not exist;
#        hrank, hnode  hbytes with --qd-green-gather rank or node (QD speed
#                 item B: the slip all-gather to every rank, or one copy
#                 per node in shared memory), reading HCACHE only;
#        hnode8, hnode16  hnode with SEAS_QD_HMAT_PIECE_DIV=8 or 16 (V^T cut
#                 into smaller pieces than the default divisor 4);
#      repeating a mode measures the run-to-run spread;
#   3. timing_summary.txt: one line per run (s/step, solve ms/call and the
#      other parts) from the end-of-run [STAGE-TIMERS] table, then the
#      [QD-GREEN] and [QD-HMAT] build or read lines.  In a green or hgreen
#      run "solve" is the G product.
#   With EVENT_REF="DIR PREFIX" (the accuracy gate's event runs,
#   qd_speed_submit.sh STEPS=events), each run is also compared with that
#   reference run through its earthquakes by jobs/vista_qd_speed/
#   qd_event_compare.py, into <run>/event_compare.txt.
#   With LOG_VIEW=1 each run also writes PETSc's -log_view table to
#   <run>/petsc_log.txt (QD speed item D1).
#   4. layout_compare.txt (QD speed item C), for each hbytes run that
#      follows an hcount run: every output file of the hbytes run must be
#      byte-identical (cmp) to the same file of the hcount run, both runs
#      must write the same files, and the bytes layout's largest rank must
#      hold at most 1.2x the mean ([QD-HMAT] layout line); otherwise the
#      job fails.  With LAYOUT_REF_JOB=<job id> the first hcount run is
#      also compared with that earlier job's hgreen run, over the lines
#      both files have (reported, not gated).
#   5. gather_compare.txt (QD speed item B), for each hnode* run that
#      follows an hrank run: every output file must be byte-identical to
#      the hrank run's, as in 4; otherwise the job fails.  Each run's log
#      must also hold its gather line ("[QD-GREEN] slip gather node: N
#      node(s)" with N this job's nodes, or "[QD-GREEN] slip gather rank:").
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
# Not STEPS: qd_speed_submit.sh's STEPS (its list of job kinds) reaches
# this script through sbatch's environment, and was taken for the step
# count (jobs 1080401-1080409 ran 0 steps with --max-steps "build timing
# green").
RUN_STEPS="${RUN_STEPS:-1000}"
if ! [[ "${RUN_STEPS}" =~ ^[1-9][0-9]*$ ]]; then
   echo "ERROR: RUN_STEPS='${RUN_STEPS}' is not a positive integer." >&2
   exit 1
fi
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
# hgreen: EPS 1e-10 is design rev2's choice (round-off-level product error,
# production G 2.85x smaller); HMEM is the GiB per node the compression may
# hold at once in gathered blocks (the driver's default, 64).
HEPS="${HEPS:-1e-10}"
HMEM="${HMEM:-64}"
HCACHE="${HCACHE:-${SCRATCH:-${WORK}}/seas-qd-green/bp5_prod_p${ORDER}_h${HEPS}.bin}"
if [[ " ${MODES} " =~ \ (hcount|hbytes|hrank|hnode|hnode8|hnode16)\  ]] \
   && [ ! -f "${HCACHE}" ]; then
   echo "ERROR: MODES='${MODES}': no H file ${HCACHE}; hcount, hbytes," \
        "hrank and hnode* only read it." >&2
   exit 1
fi
LOG_VIEW="${LOG_VIEW:-0}"
LAYOUT_REF_JOB="${LAYOUT_REF_JOB:-}"
EVENT_REF="${EVENT_REF:-}"
if [ -n "${EVENT_REF}" ]; then
   read -r EVENT_REF_DIR EVENT_REF_PFX <<< "${EVENT_REF}"
   if [ ! -f "${SEAS_DIR}/${EVENT_REF_DIR}/${EVENT_REF_PFX}_global.txt" ]; then
      echo "ERROR: EVENT_REF='${EVENT_REF}': no ${EVENT_REF_DIR}/${EVENT_REF_PFX}_global.txt" \
           "under ${SEAS_DIR}." >&2
      exit 1
   fi
fi
GOLD_README_KNOBS+=(ORDER RUN_STEPS MODES GREEN_CACHE HEPS HMEM HCACHE EVENT_REF
                    LOG_VIEW LAYOUT_REF_JOB)

cd "${SEAS_DIR}"
OUT_ROOT="bp5/results_qd_speed/timing_n${NNODES}_r${NRANKS}_c${CPT}_p${ORDER}_job${SLURM_JOB_ID}"
mkdir -p "${OUT_ROOT}"
g0_provenance "${OUT_ROOT}"
{
   echo "layout: ${NNODES} node(s), ${NRANKS} ranks, ${CPT} CPUs per rank," \
        "SLURM_TASKS_PER_NODE=${SLURM_TASKS_PER_NODE:-?}"
   echo "modes: ${MODES}; order ${ORDER}; ${RUN_STEPS} steps"
   echo "green cache: ${GREEN_CACHE}"
   echo "H-matrix cache: ${HCACHE} (eps ${HEPS}, ${HMEM} GiB per node while compressing)"
   echo "event reference: ${EVENT_REF:-none}"
   echo "PETSc -log_view: ${LOG_VIEW}; layout reference job: ${LAYOUT_REF_JOB:-none}"
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
   PIECE_DIV=""
   case "${MODE}" in
      loops|matrix) OPS=(--qd-operators "${MODE}") ;;
      default)      OPS=() ;;
      green)        mkdir -p "$(dirname "${GREEN_CACHE}")"
                    OPS=(--qd-green --qd-green-cache "${GREEN_CACHE}") ;;
      hgreen)       mkdir -p "$(dirname "${HCACHE}")"
                    OPS=(--qd-green --qd-green-cache "${GREEN_CACHE}"
                         --qd-green-compress "${HEPS}" --qd-green-hcache "${HCACHE}"
                         --qd-green-compress-mem "${HMEM}") ;;
      hcount|hbytes)
                    OPS=(--qd-green --qd-green-cache "${GREEN_CACHE}"
                         --qd-green-compress "${HEPS}" --qd-green-hcache "${HCACHE}"
                         --qd-green-compress-mem "${HMEM}"
                         --qd-green-layout "${MODE#h}") ;;
      hrank|hnode|hnode8|hnode16)
                    gather="${MODE#h}"; gather="${gather%%[0-9]*}"
                    PIECE_DIV="${MODE#hnode}"; [ "${PIECE_DIV}" != "${MODE}" ] || PIECE_DIV=""
                    OPS=(--qd-green --qd-green-cache "${GREEN_CACHE}"
                         --qd-green-compress "${HEPS}" --qd-green-hcache "${HCACHE}"
                         --qd-green-compress-mem "${HMEM}"
                         --qd-green-layout bytes --qd-green-gather "${gather}") ;;
      *) echo "ERROR: unknown mode ${MODE}" >&2; overall=1; continue ;;
   esac
   echo "=== ${TAG}: start $(date -u +%T) ==="
   (
      # PETSc reads PETSC_OPTIONS in PetscInitialize; only this run gets it
      if [ "${LOG_VIEW}" = 1 ]; then
         export PETSC_OPTIONS="-log_view ascii:${OUT}/petsc_log.txt"
      fi
      # The H-matrix test hooks: only hnode8 and hnode16 set one
      unset SEAS_QD_HMAT_NODE_SIZE SEAS_QD_HMAT_PIECE_DIV
      if [ -n "${PIECE_DIV}" ]; then export SEAS_QD_HMAT_PIECE_DIV="${PIECE_DIV}"; fi
      ibrun ./seas_bp5_full --mesh "${G0_MESH_PROD}" "${G0_V63_FLAGS[@]}" \
            --order "${ORDER}" "${OPS[@]}" --max-steps "${RUN_STEPS}" \
            --stage-timers 100 --checkpoint-interval 0 \
            --output-dir "${OUT}" --output-prefix "${PFX}" > "${OUT}/run.log" 2>&1
   )
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
   if ! grep -q "end of run: ${RUN_STEPS} accepted steps" "${OUT}/stage_timers_final.txt"; then
      echo "ERROR: ${TAG}: the run did not take ${RUN_STEPS} accepted steps:" \
           "$(grep 'end of run' "${OUT}/stage_timers_final.txt")" >&2
      overall=1
   fi
   grep '\[QD-OPERATORS\]' "${OUT}/run.log" > "${OUT}/qd_operators.txt"
   if [ -n "${EVENT_REF}" ]; then
      if python3 jobs/vista_qd_speed/qd_event_compare.py \
            --ref "${EVENT_REF_DIR}" "${EVENT_REF_PFX}" --test "${OUT}" "${PFX}" \
            > "${OUT}/event_compare.txt" 2>&1; then
         cat "${OUT}/event_compare.txt"
      else
         echo "ERROR: ${TAG}: qd_event_compare.py failed:" >&2
         cat "${OUT}/event_compare.txt" >&2
         overall=1
      fi
   fi
   if [ "${MODE}" != loops ] && [ "${MODE}" != matrix ] && [ "${MODE}" != default ]; then
      grep -E '\[QD-(GREEN|HMAT)\]' "${OUT}/run.log" > "${OUT}/qd_green.txt"
      if ! grep -q '\[QD-GREEN\] check against the solve path' "${OUT}/qd_green.txt"; then
         echo "ERROR: ${TAG}: no '[QD-GREEN] check against the solve path' line." >&2
         overall=1
      fi
      if [ "${MODE}" != green ] && ! grep -q '\[QD-HMAT\] stored' "${OUT}/qd_green.txt"; then
         echo "ERROR: ${TAG}: no '[QD-HMAT] stored' line: G was not compressed." >&2
         overall=1
      fi
   fi
   if [ "${MODE}" = hcount ] || [ "${MODE}" = hbytes ]; then
      if ! grep -q "\[QD-HMAT\] read from " "${OUT}/qd_green.txt"; then
         echo "ERROR: ${TAG}: no '[QD-HMAT] read from' line: the H file was not read." >&2
         overall=1
      fi
      if ! grep -q "\[QD-HMAT\] layout ${MODE#h} " "${OUT}/qd_green.txt"; then
         echo "ERROR: ${TAG}: no '[QD-HMAT] layout ${MODE#h}' line." >&2
         overall=1
      fi
   fi
   case "${MODE}" in
      hrank|hnode|hnode8|hnode16)
         if ! grep -q "\[QD-HMAT\] read from " "${OUT}/qd_green.txt"; then
            echo "ERROR: ${TAG}: no '[QD-HMAT] read from' line: the H file was not read." >&2
            overall=1
         fi
         if ! grep -q "\[QD-HMAT\] layout bytes " "${OUT}/qd_green.txt"; then
            echo "ERROR: ${TAG}: no '[QD-HMAT] layout bytes' line." >&2
            overall=1
         fi
         if [ "${MODE}" = hrank ]; then
            want="[QD-GREEN] slip gather rank: "
         else
            want="[QD-GREEN] slip gather node: ${NNODES} node(s), "
         fi
         if ! grep -qF "${want}" "${OUT}/qd_green.txt"; then
            echo "ERROR: ${TAG}: no '${want}' line." >&2
            overall=1
         fi ;;
   esac
done

# ---- 3. Summary ----
python3 jobs/vista_qd_speed/qd_speed_summary.py "${OUT_ROOT}" \
   | tee "${OUT_ROOT}/timing_summary.txt"

# ---- 4. Layouts: bytes against count (QD speed item C) ----
# same_files A_DIR A_PFX B_DIR B_PFX: cmp every output file (not the logs
# and tables, which hold times) of run A with run B's; 0 when all match.
same_files() {
   local a_dir="$1" a_pfx="$2" b_dir="$3" b_pfx="$4" bad=0 n=0 f sfx
   for f in "${a_dir}/${a_pfx}"_*; do
      [ -f "${f}" ] || continue
      sfx="${f#"${a_dir}/${a_pfx}"}"
      case "${sfx}" in _matrix_stats.json) continue ;; esac
      n=$((n + 1))
      if [ ! -f "${b_dir}/${b_pfx}${sfx}" ]; then
         echo "  MISSING ${b_pfx}${sfx}"; bad=1
      elif ! cmp -s "${f}" "${b_dir}/${b_pfx}${sfx}"; then
         echo "  DIFFER  ${sfx#_} ($(cmp "${f}" "${b_dir}/${b_pfx}${sfx}" 2>&1 | head -n 1))"
         bad=1
      fi
   done
   for f in "${b_dir}/${b_pfx}"_*; do
      [ -f "${f}" ] || continue
      sfx="${f#"${b_dir}/${b_pfx}"}"
      [ -f "${a_dir}/${a_pfx}${sfx}" ] || { echo "  EXTRA   ${b_pfx}${sfx}"; bad=1; }
   done
   echo "  ${n} files compared"
   [ "${n}" -gt 0 ] || bad=1
   return "${bad}"
}
CMP="${OUT_ROOT}/layout_compare.txt"
: > "${CMP}"
k=0
last_count=""
first_count=""
for MODE in ${MODES}; do
   k=$((k + 1))
   TAG="run${k}_${MODE}"
   if [ "${MODE}" = hcount ]; then
      last_count="${TAG}"
      [ -n "${first_count}" ] || first_count="${TAG}"
   elif [ "${MODE}" = hbytes ] && [ -n "${last_count}" ]; then
      echo "${TAG} against ${last_count}:" >> "${CMP}"
      if same_files "${OUT_ROOT}/${last_count}" "qd_speed_${last_count}" \
                    "${OUT_ROOT}/${TAG}" "qd_speed_${TAG}" >> "${CMP}"; then
         echo "  PASS: byte-identical" >> "${CMP}"
      else
         echo "  FAIL" >> "${CMP}"
         overall=1
      fi
      line="$(grep '\[QD-HMAT\] layout bytes' "${OUT_ROOT}/${TAG}/qd_green.txt" | tail -n 1)"
      ratio="$(sed -n 's/.*, \([0-9.]*\)x the mean.*/\1/p' <<< "${line}")"
      if [ -n "${ratio}" ] && awk -v r="${ratio}" 'BEGIN { exit !(r <= 1.2) }'; then
         echo "  PASS: largest rank ${ratio}x the mean (at most 1.2)" >> "${CMP}"
      else
         echo "  FAIL: largest rank '${ratio:-?}'x the mean (at most 1.2)" >> "${CMP}"
         overall=1
      fi
   fi
done
if [ -n "${LAYOUT_REF_JOB}" ] && [ -n "${first_count}" ]; then
   ref="$(ls -d bp5/results_qd_speed/timing_*_job"${LAYOUT_REF_JOB}"/run*_hgreen 2>/dev/null \
          | head -n 1)"
   echo "${first_count} against job ${LAYOUT_REF_JOB} (${ref:-not found}), over the" \
        "lines both files have; reported, not gated:" >> "${CMP}"
   if [ -n "${ref}" ]; then
      rpfx="qd_speed_$(basename "${ref}")"
      for f in "${ref}/${rpfx}"_*.txt; do
         [ -f "${f}" ] || continue
         sfx="${f#"${ref}/${rpfx}"}"
         g="${OUT_ROOT}/${first_count}/qd_speed_${first_count}${sfx}"
         if [ ! -f "${g}" ]; then echo "  ${sfx#_}: not in ${first_count}"; continue; fi
         nl="$(( $(wc -l < "${f}") < $(wc -l < "${g}") ? $(wc -l < "${f}") : $(wc -l < "${g}") ))"
         if cmp -s <(head -n "${nl}" "${f}") <(head -n "${nl}" "${g}"); then
            echo "  ${sfx#_}: same over ${nl} lines"
         else
            echo "  ${sfx#_}: differs within ${nl} lines" \
                 "($(cmp <(head -n "${nl}" "${f}") <(head -n "${nl}" "${g}") 2>&1 | head -n 1))"
         fi
      done >> "${CMP}"
   fi
fi
if [ -s "${CMP}" ]; then cat "${CMP}"; else rm -f "${CMP}"; fi

# ---- 5. Slip gather: node against rank (QD speed item B) ----
CMP="${OUT_ROOT}/gather_compare.txt"
: > "${CMP}"
k=0
last_rank=""
for MODE in ${MODES}; do
   k=$((k + 1))
   TAG="run${k}_${MODE}"
   case "${MODE}" in
      hrank) last_rank="${TAG}" ;;
      hnode|hnode8|hnode16)
         [ -n "${last_rank}" ] || continue
         echo "${TAG} against ${last_rank}:" >> "${CMP}"
         if same_files "${OUT_ROOT}/${last_rank}" "qd_speed_${last_rank}" \
                       "${OUT_ROOT}/${TAG}" "qd_speed_${TAG}" >> "${CMP}"; then
            echo "  PASS: byte-identical" >> "${CMP}"
         else
            echo "  FAIL" >> "${CMP}"
            overall=1
         fi
         # Reported: the gather line and the V^T pieces of this run
         grep -F "[QD-GREEN] slip gather" "${OUT_ROOT}/${TAG}/qd_green.txt" \
            | tail -n 1 | sed 's/^/  /' >> "${CMP}"
         sed -n 's/.*\(V^T in [0-9]* pieces\).*/  \1/p' \
            "${OUT_ROOT}/${TAG}/qd_green.txt" | tail -n 1 >> "${CMP}" ;;
   esac
done
if [ -s "${CMP}" ]; then cat "${CMP}"; else rm -f "${CMP}"; fi
echo "QD speed timing overall: $([ "${overall}" = 0 ] && echo PASS || echo FAIL)"
exit "${overall}"
