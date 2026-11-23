#!/bin/bash
# h_submit.sh — submit the GPU plan chain on Horizon, each job after the
# previous one succeeded (afterok), so nothing waits on a person between
# jobs:
#   h_cuda_build -> h_gpu_checks -> h_green_cache -> h_gpu_gates -> h_full_cycle
#
#   cd miniapps/seas                     (of the Horizon clone)
#   QD_PATH=green bash jobs/horizon/h_submit.sh            # dry run
#   QD_PATH=green SUBMIT=1 bash jobs/horizon/h_submit.sh   # submit
#
# QD_PATH (cudss or green, no default) is the full cycle's stage solve.
# H_CHAIN (default: all five, in order) submits part of the chain, e.g.
# H_CHAIN="gates cycle"; AFTER=<job id> makes its first job wait on that one.
# Job ids go to h_submit_<date>.txt.

set -euo pipefail
if [ ! -f jobs/horizon/h_env.sh ]; then
   echo "ERROR: run from miniapps/seas of the Horizon clone." >&2
   exit 1
fi
QD_PATH="${QD_PATH:?set QD_PATH=cudss or green for the full cycle}"
case "${QD_PATH}" in
   cudss|green) ;;
   *) echo "ERROR: QD_PATH must be cudss or green, not '${QD_PATH}'." >&2; exit 1 ;;
esac
H_CHAIN="${H_CHAIN:-build checks green gates cycle}"
dep="${AFTER:-}"
log="h_submit_$(date +%Y%m%d_%H%M%S).txt"

for step in ${H_CHAIN}; do
   extra=()
   case "${step}" in
      build)  script=h_cuda_build.sbatch ;;
      checks) script=h_gpu_checks.sbatch ;;
      green)  script=h_green_cache.sbatch ;;
      gates)  script=h_gpu_gates.sbatch ;;
      cycle)  script=h_full_cycle.sbatch; extra=(--export=ALL,QD_PATH="${QD_PATH}") ;;
      *) echo "ERROR: unknown step '${step}' in H_CHAIN." >&2; exit 1 ;;
   esac
   cmd=(sbatch --parsable)
   [ -n "${dep}" ] && cmd+=(--dependency=afterok:"${dep}" --kill-on-invalid-dep=yes)
   cmd+=("${extra[@]}" "jobs/horizon/${script}")
   if [ "${SUBMIT:-0}" = 1 ]; then
      # TACC's sbatch prints a banner on stdout even with --parsable; the
      # job id ("id" or "id;cluster") is the last line.
      out="$("${cmd[@]}")"
      id="$(printf '%s\n' "${out}" | tail -n 1)"
      id="${id%%;*}"
      if ! [[ "${id}" =~ ^[0-9]+$ ]]; then
         echo "ERROR: ${step}: no job id in sbatch output:" >&2
         printf '%s\n' "${out}" >&2
         exit 1
      fi
      echo "${step}: ${id}  (${cmd[*]})" | tee -a "${log}"
      dep="${id}"
   else
      echo "dry run, ${step}: ${cmd[*]}"
      dep="<${step}>"
   fi
done
