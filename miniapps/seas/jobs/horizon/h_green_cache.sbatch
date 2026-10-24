#!/bin/bash
#SBATCH -J h_green_cache
#SBATCH -o h_green_cache_%j.out
#SBATCH -e h_green_cache_%j.err
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 02:00:00
#SBATCH -A EAR26024

# GPU plan on Horizon: write the dense Green's-function cache GREEN_CACHE
# (production mesh, p=1, 37.5 GB) that the G3c gate and the dense-G full
# cycle read.  The Vista file cannot be copied here directly, and on Vista
# the QD speed job built it in about 4 minutes on one 144-rank node.
#
# Runs jobs/vista_qd_speed/qd_speed_timing.sbatch with MODES=green and
# RUN_STEPS=10 on the 144 Grace cores of one node, on the CPU (no
# --device): it builds G with the solve path, writes it with its read-back
# check, checks it against the solve path, then takes 10 steps.  If the
# file already exists, the run reads and checks it instead.
#
# Horizon's $SCRATCH is NFS (VAST), not Lustre as on Vista.  There the
# 144-rank MPI-IO write of G read back wrong (job 68408: "G rows differ
# from what was written (first on rank 12)"; the file was deleted).  So G
# is written to the node's local /tmp, copied to GREEN_CACHE by one
# process, and the copy is compared byte for byte with the local file.
# The jobs that read GREEN_CACHE check G against the solve path at
# startup (elasticity_operator_qd.inl), with a test slip on every DOF.
#
#   sbatch jobs/horizon/h_green_cache.sbatch   (from miniapps/seas of the Horizon clone)

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/horizon/h_env.sh" || exit 1
export MODES=green RUN_STEPS=10 ORDER=1
if [ -s "${GREEN_CACHE}" ]; then
   exec bash "${SLURM_SUBMIT_DIR}/jobs/vista_qd_speed/qd_speed_timing.sbatch"
fi

FINAL="${GREEN_CACHE}"
LOCAL_DIR="/tmp/seas-qd-green-${SLURM_JOB_ID:-$$}"
mkdir -p "${LOCAL_DIR}" "$(dirname "${FINAL}")" || exit 1
export GREEN_CACHE="${LOCAL_DIR}/$(basename "${FINAL}")"
bash "${SLURM_SUBMIT_DIR}/jobs/vista_qd_speed/qd_speed_timing.sbatch"
rc=$?
if [ "${rc}" = 0 ]; then
   echo "copying ${GREEN_CACHE} ($(stat -c %s "${GREEN_CACHE}") bytes) to ${FINAL}"
   if cp "${GREEN_CACHE}" "${FINAL}.part" && sync "${FINAL}.part" &&
      cmp "${GREEN_CACHE}" "${FINAL}.part" && mv "${FINAL}.part" "${FINAL}"; then
      echo "G file: ${FINAL} ($(stat -c '%s bytes, modified %y' "${FINAL}"))," \
           "identical to the local copy"
   else
      echo "ERROR: copying G to ${FINAL} failed or the copy differs." >&2
      rm -f "${FINAL}.part"
      rc=1
   fi
fi
rm -rf "${LOCAL_DIR}"
exit "${rc}"
