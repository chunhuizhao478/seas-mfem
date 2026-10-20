#!/bin/bash
#SBATCH -J qd_integrator
#SBATCH -o qd_integrator_%j.out
#SBATCH -e qd_integrator_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 04:00:00
#SBATCH -A EAR20006

# QD speed-up plan v2, item A (integrator study): one run of the production
# BP5 QD deck (v63 flags, p=1, compressed G read from the H file) from t = 0
# to TFINAL (default 7.2e9 s, about 2.9 yr after gold 1074099's second
# event starts) with the PETSc TS options file TSCFG.  One run per job, so
# the runs go side by side; jobs/vista_qd_speed/qd_speed_submit.sh
# STEPS=integrator submits them all, then qd_integrator_compare.sbatch.
#
#   TSCFG   options file, relative to miniapps/seas (required)
#   TAG     run name (required): output in
#           bp5/results_qd_speed/integrator/<TAG>_job<id>
#   TFINAL  end time in s (default 7.2e9)
#   HCACHE  the H file (default the one job 1084131 wrote, eps 1e-10).  It
#           is only read: the job stops if it does not exist, so it never
#           compresses or writes there.
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

cd "${SEAS_DIR}"
TSCFG="${TSCFG:-}"
TAG="${TAG:-}"
TFINAL="${TFINAL:-7.2e9}"
HEPS="${HEPS:-1e-10}"
HCACHE="${HCACHE:-${SCRATCH:-${WORK}}/seas-qd-green/bp5_prod_p1_h${HEPS}.bin}"
if [ -z "${TSCFG}" ] || [ ! -f "${TSCFG}" ]; then
   echo "ERROR: TSCFG='${TSCFG}' is not an options file under ${SEAS_DIR}." >&2
   exit 1
fi
if ! [[ "${TAG}" =~ ^[A-Za-z0-9_]+$ ]]; then
   echo "ERROR: TAG='${TAG}' must be letters, digits and _." >&2
   exit 1
fi
if [ ! -f "${HCACHE}" ]; then
   echo "ERROR: no H file ${HCACHE}; this job only reads it." >&2
   exit 1
fi
GOLD_README_KNOBS+=(TSCFG TAG TFINAL HEPS HCACHE)

# The v63 flags with this run's options file in place of
# petsc_ts_rk45_tandem.cfg.
FLAGS=()
k=0
while [ "${k}" -lt "${#G0_V63_FLAGS[@]}" ]; do
   FLAGS+=("${G0_V63_FLAGS[k]}")
   if [ "${G0_V63_FLAGS[k]}" = --petsc-ts-options ]; then
      FLAGS+=("${TSCFG}")
      k=$((k + 1))
   fi
   k=$((k + 1))
done

OUT="bp5/results_qd_speed/integrator/${TAG}_job${SLURM_JOB_ID}"
PFX="qd_int_${TAG}"
mkdir -p "${OUT}"
g0_provenance "${OUT}"
{
   echo "integrator run ${TAG}: options ${TSCFG}, t_final ${TFINAL} s"
   echo "layout: ${SLURM_JOB_NUM_NODES:-1} node(s), ${SLURM_NTASKS:-144} ranks"
   echo "H file (read only): ${HCACHE}"
   echo "options file:"
   sed 's/^/   /' "${TSCFG}"
} >> "${OUT}/provenance.txt"

echo "=== ${TAG}: start $(date -u +%T) ==="
ibrun ./seas_bp5_full --mesh "${G0_MESH_PROD}" "${FLAGS[@]}" --order 1 \
      --qd-green --qd-green-compress "${HEPS}" --qd-green-hcache "${HCACHE}" \
      --tfinal "${TFINAL}" --stage-timers 5000 --checkpoint-interval 0 \
      --output-dir "${OUT}" --output-prefix "${PFX}" > "${OUT}/run.log" 2>&1
rc=$?
echo "${rc}" > "${OUT}/exit_code"
echo "=== ${TAG}: exit ${rc} at $(date -u +%T) ==="

fail=0
[ "${rc}" = 0 ] || fail=1
grep -F '[TS]' "${OUT}/run.log" | head -n 1 > "${OUT}/ts.txt"
if [ ! -s "${OUT}/ts.txt" ]; then
   echo "ERROR: ${TAG}: no '[TS]' line in run.log." >&2
   fail=1
fi
if ! grep -q '\[QD-HMAT\]' "${OUT}/run.log" \
   || ! grep -q "H-matrix read from ${HCACHE}" "${OUT}/run.log"; then
   echo "ERROR: ${TAG}: the log does not say the H-matrix was read from ${HCACHE}." >&2
   fail=1
fi
sed -n '/^\[STAGE-TIMERS\] end of run/,$p' "${OUT}/run.log" \
   | grep '^\[STAGE-TIMERS\]' > "${OUT}/stage_timers_final.txt"
grep -E '^  (Final time|Total steps|Step rejections|Seismic events):' \
   "${OUT}/run.log" > "${OUT}/summary.txt"
cat "${OUT}/ts.txt" "${OUT}/summary.txt"
head -n 2 "${OUT}/stage_timers_final.txt"
# The run must end at TFINAL: the last global row's time.
last_t="$(tail -n 1 "${OUT}/${PFX}_global.txt" 2>/dev/null | awk '{print $1}')"
if ! python3 -c "import sys; t=float(sys.argv[1]); T=float(sys.argv[2]); sys.exit(0 if abs(t-T) <= 1e-9*T else 1)" \
      "${last_t:-nan}" "${TFINAL}"; then
   echo "ERROR: ${TAG}: last global row at t = ${last_t:-none}, not ${TFINAL}." >&2
   fail=1
fi
echo "QD integrator run ${TAG}: $([ "${fail}" = 0 ] && echo PASS || echo FAIL)"
exit "${fail}"
