#!/bin/bash
#SBATCH -J qd_integrator_compare
#SBATCH -o qd_integrator_compare_%j.out
#SBATCH -e qd_integrator_compare_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# QD speed-up plan v2, item A: compare the integrator-study runs
# (jobs/vista_qd_speed/qd_integrator.sbatch) with the reference run, after
# all of them (qd_speed_submit.sh STEPS=integrator submits this with
# afterany, so a failed run is reported, not hidden).
#
#   RUNS  "tag:jobid ..." of every run; the first is the reference
#   LINE  "tag tag": today's setting at two tolerances (default "a0 a3")
#
# Writes bp5/results_qd_speed/integrator/compare_job<id>.txt.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"
cd "${SEAS_DIR}"
RUNS="${RUNS:-}"
LINE="${LINE:-a0 a3}"
if [ -z "${RUNS}" ]; then
   echo "ERROR: RUNS is empty." >&2
   exit 1
fi
ARGS=()
first=1
for spec in ${RUNS}; do
   tag="${spec%%:*}"
   jid="${spec##*:}"
   d="bp5/results_qd_speed/integrator/${tag}_job${jid}"
   echo "${tag}: ${d}, exit code $(cat "${d}/exit_code" 2>/dev/null || echo none)"
   if [ ! -d "${d}" ]; then
      echo "ERROR: no run directory ${d}." >&2
      exit 1
   fi
   if [ "${first}" = 1 ]; then
      ARGS+=(--ref "${tag}" "${d}")
      first=0
   else
      ARGS+=(--run "${tag}" "${d}")
   fi
done
# shellcheck disable=SC2086
python3 jobs/vista_qd_speed/qd_integrator_compare.py "${ARGS[@]}" --line ${LINE} \
   | tee "bp5/results_qd_speed/integrator/compare_job${SLURM_JOB_ID}.txt"
exit "${PIPESTATUS[0]}"
