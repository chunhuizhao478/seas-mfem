#!/bin/bash
# h_d2_submit.sh — submit the GPU dynamic-rupture plan's jobs on Horizon:
# h_d2_build, then h_d2_gate on 1, 2 and 4 GPUs and/or the 200 m runs
# h_d4_200m once the build job has ended (afterany: each checks for itself
# that the binary was built at its commit, so a failed trace in the build
# job does not hold them, as it held Vista's gates after job 1096115).
#
#   cd miniapps/seas                          (of the Horizon clone)
#   bash jobs/horizon/h_d2_submit.sh          # checks, then a dry run
#   SUBMIT=1 bash jobs/horizon/h_d2_submit.sh # checks, then submits
#
# H_CHAIN (default "build gates") submits part of it, in order: build,
# gates, d4.  GATE_NPS (default "1 2 4") are the gates' GPU counts, D4_CASES
# (default "tpv205:1 tpv205:2 tpv104:1 tpv104:2") the 200 m cases, PROB:P,
# one h_d4_200m job each.  AFTER=<job id> makes the first job wait on that
# one.  DEVICE_FAULT=1 (GPU plan D3) runs the gates and the 200 m runs with
# the fault on the GPU (sbatch passes it on; the build job always traces
# both).  Job ids go to h_d2_submit_<date>.txt.
#
# Before anything is submitted it checks, on this login node, that the clone
# is clean at its upstream commit, that the meshes, the D1 references (for
# gates) and the 200 m meshes (for d4) of h_dyn_env.sh are in place, and that
# HDF5 is installed or its tarball fetched; it stops with what is missing.

set -euo pipefail
if [ ! -f jobs/horizon/h_dyn_env.sh ]; then
   echo "ERROR: run from miniapps/seas of the Horizon clone." >&2
   exit 1
fi
# shellcheck source=jobs/horizon/h_dyn_env.sh
source jobs/horizon/h_dyn_env.sh
# shellcheck source=jobs/vista/d1_env.sh
source jobs/vista/d1_env.sh
H_CHAIN="${H_CHAIN:-build gates}"
GATE_NPS="${GATE_NPS:-1 2 4}"
D4_CASES="${D4_CASES:-tpv205:1 tpv205:2 tpv104:1 tpv104:2}"
export DEVICE_FAULT="${DEVICE_FAULT:-0}"
case "${DEVICE_FAULT}" in 0|1) ;; *) echo "ERROR: DEVICE_FAULT is 0 or 1, not '${DEVICE_FAULT}'." >&2; exit 1 ;; esac

missing=0
if [ -n "$(git status --porcelain --untracked-files=no)" ] ||
   [ "$(git rev-parse HEAD)" != "$(git rev-parse '@{u}' 2>/dev/null)" ]; then
   echo "MISSING: the clone is not clean at its upstream commit (git pull --ff-only)." >&2
   missing=1
fi
s=$(cd "${MESH_SRC}" 2>/dev/null && sha256sum tpv205/mesh/tpv2053d_1000m.msh tpv104/mesh/tpv104_1000m.msh 2>/dev/null || true)
if [ "${s}" != "${D1_MESH_SHA256}" ]; then
   echo "MISSING: D0a's meshes (sha256 of d1_env.sh) in ${MESH_SRC}/tpv205/mesh and tpv104/mesh." >&2
   missing=1
fi
if [[ " ${H_CHAIN} " == *" gates "* ]]; then
   for np in ${GATE_NPS}; do
      n=$( (ls -d "${D1_REF_ROOT}/d1_gg_check_r${np}_"*/runs 2>/dev/null || true) | wc -l)
      if [ "${n}" != 1 ] || ! grep -q ", NP ${np}," "${D1_REF_ROOT}/d1_gg_check_r${np}_"*/RUN_INFO.txt 2>/dev/null; then
         echo "MISSING: one D1 gg check at NP ${np} (d1_gg_check_r${np}_<jobid>/RUN_INFO.txt and runs/) in ${D1_REF_ROOT}." >&2
         missing=1
      fi
   done
fi
if [[ " ${H_CHAIN} " == *" d4 "* ]]; then
   for c in ${D4_CASES}; do
      case "${c}" in
         tpv205:1|tpv205:2) m=tpv205/mesh/tpv2053d_200m.msh ;;
         tpv104:1|tpv104:2) m=tpv104/mesh/tpv104_200m.msh ;;
         *) echo "ERROR: D4_CASES entry '${c}' is not tpv205:P or tpv104:P (P 1 or 2)." >&2; exit 1 ;;
      esac
      d=.; [ -f "${m}" ] || d="${MESH_SRC}"
      if ! d4_check_mesh "${d}" "${m}" > /dev/null; then
         echo "MISSING: the gold's 200 m mesh ${m} (sha256 of h_dyn_env.sh) here or in ${MESH_SRC}." >&2
         missing=1
      fi
   done
fi
HDF5_VERSION="${HDF5_VERSION:-1.14.6}"
if [ ! -f "${HDF5_ROOT}/hdf5-${HDF5_VERSION}-cuda/include/hdf5.h" ] &&
   [ ! -f "${HDF5_ROOT}/downloads/hdf5-${HDF5_VERSION}.tar.gz" ]; then
   echo "MISSING: HDF5 in ${HDF5_ROOT}.  Fetch its tarball here (login node, internet):" >&2
   echo "   source jobs/horizon/h_dyn_env.sh && MODE=cuda STAGE=hdf5-fetch WITH_HDF5=1 \\" >&2
   echo "      PROVENANCE_FILE=\$HDF5_ROOT/fetch_provenance.txt bash ../../build_horizon.sh" >&2
   missing=1
fi
[ "${missing}" = 0 ] || exit 1
echo "Checks passed for '${H_CHAIN}' (DEVICE_FAULT=${DEVICE_FAULT}): clone at $(git rev-parse HEAD), meshes, references, HDF5."

log="h_d2_submit_$(date +%Y%m%d_%H%M%S).txt"
dep="${AFTER:-}"
# submit NAME DEPENDENCY SBATCH-ARGS...: prints the job id
submit () {
   local name=$1 d=$2 out id
   shift 2
   local cmd=(sbatch --parsable)
   [ -n "${d}" ] && cmd+=(--dependency="${d}" --kill-on-invalid-dep=yes)
   cmd+=("$@")
   if [ "${SUBMIT:-0}" = 1 ]; then
      # TACC's sbatch prints a banner on stdout even with --parsable; the
      # job id ("id" or "id;cluster") is the last line.
      out="$("${cmd[@]}")"
      id="$(printf '%s\n' "${out}" | tail -n 1)"
      id="${id%%;*}"
      if ! [[ "${id}" =~ ^[0-9]+$ ]]; then
         echo "ERROR: ${name}: no job id in sbatch output:" >&2
         printf '%s\n' "${out}" >&2
         exit 1
      fi
      echo "${name}: ${id}  (DEVICE_FAULT=${DEVICE_FAULT} ${cmd[*]})" | tee -a "${log}" >&2
      echo "${id}"
   else
      echo "dry run, ${name}: DEVICE_FAULT=${DEVICE_FAULT} ${cmd[*]}" >&2
      echo "<${name}>"
   fi
}

for step in ${H_CHAIN}; do
   case "${step}" in
      build)
         id=$(submit build "${dep:+afterok:${dep}}" jobs/horizon/h_d2_build.sbatch)
         dep="${id}"
         ;;
      gates)
         for np in ${GATE_NPS}; do
            submit "gate_np${np}" "${dep:+afterany:${dep}}" -n "${np}" jobs/horizon/h_d2_gate.sbatch > /dev/null
         done
         ;;
      d4)
         for c in ${D4_CASES}; do
            submit "d4_${c%%:*}_p${c##*:}" "${dep:+afterany:${dep}}" jobs/horizon/h_d4_200m.sbatch "${c%%:*}" "${c##*:}" > /dev/null
         done
         ;;
      *) echo "ERROR: unknown step '${step}' in H_CHAIN." >&2; exit 1 ;;
   esac
done
