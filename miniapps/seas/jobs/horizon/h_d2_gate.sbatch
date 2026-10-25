#!/bin/bash
#SBATCH -J h_d2_gate
#SBATCH -o h_d2_gate_%j.out
#SBATCH -e h_d2_gate_%j.err
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 10:00:00
#SBATCH -A EAR26024

# GPU dynamic-rupture plan D2 on Horizon: jobs/vista/d2_gh_gate.sbatch (the
# four D1 cases to 12 s with the device step, against the D1 gg references
# at the same rank count, TOL = 1e-8) on NP GPUs of one GB200 node, one rank
# each, with the settings of jobs/horizon/h_dyn_env.sh.  REF_RUNS, unless
# set, is the one d1_gg_check_r<NP>_*/runs under D1_REF_ROOT.
#
#   sbatch -n NP jobs/horizon/h_d2_gate.sbatch   (NP = 1, 2 or 4, from
#                                  miniapps/seas of the Horizon clone after
#                                  h_d2_build passed at the same commit; or
#                                  h_d2_submit.sh)

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/horizon/h_dyn_env.sh" || exit 1
NP=${SLURM_NTASKS:?run with sbatch -n NP}
if [ -z "${REF_RUNS:-}" ]; then
   mapfile -t _refs < <(ls -d "${D1_REF_ROOT}/d1_gg_check_r${NP}_"*/runs 2>/dev/null)
   [ "${#_refs[@]}" = 1 ] \
      || { echo "ERROR: want one d1_gg_check_r${NP}_*/runs in ${D1_REF_ROOT}, found ${#_refs[@]}"; exit 2; }
   export REF_RUNS="${_refs[0]}"
fi
exec bash "${SLURM_SUBMIT_DIR}/jobs/vista/d2_gh_gate.sbatch"
