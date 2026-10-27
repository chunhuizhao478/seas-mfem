#!/bin/bash
#SBATCH -J qd_green_compress
#SBATCH -o qd_green_compress_%j.out
#SBATCH -e qd_green_compress_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 144
#SBATCH -t 03:00:00
#SBATCH -A EAR20006

# QD speed item 4, before any code: how far the production Green's
# function compresses with the stage product kept at the round-off level.
# Runs jobs/vista_qd_speed/qd_green_compress.py on the cache that job
# 1081359 wrote (GREEN_CACHE, p = 1, 1000 m mesh, 34.9 GiB), for each
# (leaf, eta) in TREES (default "64:2 128:2"), at eps 1e-8, 1e-10, 1e-12
# and 1e-13.  Reads the file only; no solver run.
#
#   sbatch jobs/vista_qd_speed/qd_green_compress.sbatch   (from miniapps/seas;
#   jobs/vista_qd_speed/qd_speed_submit.sh STEPS=compress does this)
#
# numpy: Vista's default python3 (miniforge) has none; the python3/3.11.8
# module (gcc/13.2.0 toolchain) has numpy 2.2 on OpenBLAS, so its
# interpreter is called by path (QD_PY).  The SVDs run in PROCS processes
# (default 36) with OPENBLAS_NUM_THREADS (default 4) threads each.
#
# GitHub is the version control: the checkout must be clean and at the
# commit its upstream branch points to, or the job stops.

set -uo pipefail
cd "${SLURM_SUBMIT_DIR}" || exit 1
if [ ! -f jobs/vista_qd_speed/qd_green_compress.py ]; then
   echo "ERROR: submit from miniapps/seas" >&2
   exit 1
fi
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

module load gcc/13.2.0 python3/3.11.8 || exit 1
QD_PY="${QD_PY:-/opt/apps/gcc14/cuda12/python3/3.11.8/bin/python3}"
"${QD_PY}" -c "import numpy; print('numpy', numpy.__version__)" || {
   echo "ERROR: ${QD_PY} has no numpy; set QD_PY." >&2
   exit 1
}
GREEN_CACHE="${GREEN_CACHE:-${SCRATCH:-${WORK}}/seas-qd-green/bp5_prod_p1.bin}"
if [ ! -f "${GREEN_CACHE}" ]; then
   echo "ERROR: no Green's cache ${GREEN_CACHE}" >&2
   exit 1
fi
TREES="${TREES:-64:2 128:2}"
PROCS="${PROCS:-36}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-4}"
export OMP_NUM_THREADS="${OPENBLAS_NUM_THREADS}"

OUT_ROOT="bp5/results_qd_speed/compress_job${SLURM_JOB_ID}"
mkdir -p "${OUT_ROOT}"
{
   echo "job: ${SLURM_JOB_ID} on $(hostname), $(date -u +%FT%TZ)"
   echo "source: $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD))"
   echo "cache: ${GREEN_CACHE}, $(stat -c '%s bytes, modified %y' "${GREEN_CACHE}")"
   echo "python: ${QD_PY}; procs ${PROCS} x ${OPENBLAS_NUM_THREADS} threads"
   echo "trees (leaf:eta): ${TREES}"
} | tee "${OUT_ROOT}/provenance.txt"

overall=0
for tree in ${TREES}; do
   IFS=: read -r LEAF ETA <<< "${tree}"
   OUT="${OUT_ROOT}/compress_leaf${LEAF}_eta${ETA}.txt"
   echo "=== leaf ${LEAF}, eta ${ETA}: start $(date -u +%T) ==="
   "${QD_PY}" jobs/vista_qd_speed/qd_green_compress.py "${GREEN_CACHE}" \
      --leaf "${LEAF}" --eta "${ETA}" --procs "${PROCS}" > "${OUT}" 2>&1
   rc=$?
   cat "${OUT}"
   echo "=== leaf ${LEAF}, eta ${ETA}: exit ${rc} at $(date -u +%T) ==="
   [ "${rc}" = 0 ] || overall=1
done
echo "QD Green compression overall: $([ "${overall}" = 0 ] && echo PASS || echo FAIL)"
exit "${overall}"
