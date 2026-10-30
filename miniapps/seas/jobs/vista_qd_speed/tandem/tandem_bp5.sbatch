#!/bin/bash
#SBATCH -J tandem_bp5
#SBATCH -o tandem_bp5_%j.out
#SBATCH -e tandem_bp5_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# Tandem vs MFEM-SEAS speed test: one Tandem BP5 QD run.
#
#   sbatch --export=ALL,TAG=t1chk,DEGREE=1,MODE=QD,MAXSTEPS=200,PROBE_TMAX=0 \
#          jobs/vista_qd_speed/tandem/tandem_bp5.sbatch
#   (from miniapps/seas of the qd-speed checkout; -N/-n as for the MFEM run
#   it is compared with)
#
#   TAG         run name (required): output in
#               bp5/results_qd_speed/tandem/<TAG>_job<id>
#   DEGREE      1, 2 or 3: uses ${TANDEM_ROOT}/build-d3-p<DEGREE>/app/tandem
#   MODE        QD (an elasticity solve per stage) or QDGreen (Tandem's
#               discrete Green's function; its G is written under
#               GF_DIR and reused by a later run with the same DEGREE and
#               rank count)
#   TFINAL      end time in s (default 3.15576e8 = 10 yr)
#   MAXSTEPS    stop after this many steps (default none: run to TFINAL)
#   TSCFG       PETSc TS options, relative to miniapps/seas (default
#               tests/verification/petsc_ts_rk45_tandem.cfg, the file the
#               MFEM production deck reads: rk5dp, atol 1e-7, rtol 1e-50,
#               infinity norm, dt_max 0.1 yr)
#   PROBE_TMAX  Tandem fault probe t_max in s (default 9460800 = upstream's
#               deck; 0 writes every step, for the 200-step check)
#
# The deck is upstream examples/tandem/3d/bp5.toml with three changes, all
# so both codes solve the same discrete problem the same way:
#   - mesh_file: our bp5/mesh/reference/bp5_tandem_exact.msh (the MFEM
#     production mesh, sha256 checked; km, fault y = 0, tags 1/3/5);
#   - matrix_free = false, no mg_strategy: both codes factor the assembled
#     operator with MUMPS LU (examples/options/lu_mumps.cfg);
#   - final_time and probe t_max from TFINAL and PROBE_TMAX.
# The Lua scenario is upstream's bp5.lua (bp5_outside), copied unchanged.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista_qd_speed/tandem/tandem_env.sh"
seas_guard
tandem_guard

TAG="${TAG:-}"; DEGREE="${DEGREE:-}"; MODE="${MODE:-}"
TFINAL="${TFINAL:-3.15576e8}"; MAXSTEPS="${MAXSTEPS:-}"
TSCFG="${TSCFG:-tests/verification/petsc_ts_rk45_tandem.cfg}"
PROBE_TMAX="${PROBE_TMAX:-9460800}"
[[ "${TAG}" =~ ^[A-Za-z0-9_]+$ ]] || { echo "ERROR: TAG='${TAG}' must be letters, digits and _." >&2; exit 1; }
[[ "${DEGREE}" =~ ^[1-3]$ ]] || { echo "ERROR: DEGREE='${DEGREE}' must be 1, 2 or 3." >&2; exit 1; }
[[ "${MODE}" =~ ^(QD|QDGreen)$ ]] || { echo "ERROR: MODE='${MODE}' must be QD or QDGreen." >&2; exit 1; }
[ -z "${MAXSTEPS}" ] || [[ "${MAXSTEPS}" =~ ^[0-9]+$ ]] || { echo "ERROR: MAXSTEPS='${MAXSTEPS}'" >&2; exit 1; }
[ -f "${TSCFG}" ] || { echo "ERROR: TSCFG='${TSCFG}' is not a file under ${SEAS_DIR}" >&2; exit 1; }
BIN="$(tandem_bin "${DEGREE}")"
[ -x "${BIN}" ] || { echo "ERROR: ${BIN} missing; run tandem_build.sbatch DEGREE=${DEGREE} first." >&2; exit 1; }
MESH="${SEAS_DIR}/${G0_MESH_PROD}"
echo "${TANDEM_MESH_SHA256}  ${MESH}" | sha256sum -c --quiet \
   || { echo "ERROR: ${MESH} is not the production mesh (sha256)" >&2; exit 1; }

OUT="${SEAS_DIR}/bp5/results_qd_speed/tandem/${TAG}_job${SLURM_JOB_ID}"
mkdir -p "${OUT}"
cp "${TANDEM_SRC}/examples/tandem/3d/bp5.lua" "${OUT}/bp5.lua"
cp "${TANDEM_SRC}/examples/options/lu_mumps.cfg" "${OUT}/lu_mumps.cfg"
cp "${TSCFG}" "${OUT}/ts.cfg"
GF_DIR="${SCRATCH:-${WORK}}/seas-qd-green/tandem_p${DEGREE}_n${SLURM_NTASKS:-144}"
mkdir -p "${GF_DIR}"

# Probe list and the rest of the deck from upstream, with the changes above.
python3 - "${TANDEM_SRC}/examples/tandem/3d/bp5.toml" "${OUT}/bp5.toml" \
   "${MESH}" "${TFINAL}" "${MODE}" "${PROBE_TMAX}" "${GF_DIR}/gf" <<'PY' || exit 1
import re, sys
src, dst, mesh, tfinal, mode, tmax, gf = sys.argv[1:]
s = open(src).read()
def sub(pat, rep):
    global s
    s, n = re.subn(pat, rep, s, count=1, flags=re.M)
    if n != 1:
        sys.exit(f"ERROR: deck pattern {pat!r} not found in {src}")
sub(r'^final_time = .*$', f'final_time = {float(tfinal):.17g}')
sub(r'^mesh_file = .*$', f'mesh_file = "{mesh}"')
sub(r'^mode = .*$', f'mode = "{mode}"')
sub(r'^mg_strategy = .*\n', '')
sub(r'^matrix_free = .*$', 'matrix_free = false')
sub(r'^t_max = .*$', f't_max = {float(tmax):.17g}')
if mode == "QDGreen":
    s = s.replace("\n[fault_probe_output]",
                  f'\n[gf_checkpoint]\nprefix = "{gf}"\n\n[fault_probe_output]', 1)
open(dst, "w").write(s)
PY

{
   echo "job: ${SLURM_JOB_ID:-none} on $(hostname), $(date -u +%FT%TZ)"
   echo "nodes: ${SLURM_JOB_NUM_NODES:-?}, tasks: ${SLURM_NTASKS:-?}"
   echo "seas-mfem: $(git -C "${SEAS_MFEM_ROOT}" rev-parse HEAD)"
   echo "tandem: ${TANDEM_COMMIT}, binary ${BIN} sha256 $(sha256sum "${BIN}" | cut -c1-16)"
   cat "${TANDEM_ROOT}/build-d3-p${DEGREE}/provenance.txt" 2>/dev/null | sed 's/^/   build: /'
   echo "run ${TAG}: mode ${MODE}, degree ${DEGREE}, t_final ${TFINAL}, max steps ${MAXSTEPS:-none}, probe t_max ${PROBE_TMAX}"
   echo "mesh: ${MESH} (sha256 ${TANDEM_MESH_SHA256})"
   echo "lua: upstream bp5.lua sha256 $(sha256sum "${OUT}/bp5.lua" | cut -c1-16)"
   echo "Green's function checkpoint: ${GF_DIR}"
   echo "threads: OMP=${OMP_NUM_THREADS}"
   echo "--- ts.cfg"; cat "${OUT}/ts.cfg"
   echo "--- lu_mumps.cfg"; cat "${OUT}/lu_mumps.cfg"
   echo "--- bp5.toml"; cat "${OUT}/bp5.toml"
} > "${OUT}/provenance.txt"

PETSC_ARGS=(-options_file ts.cfg -options_file lu_mumps.cfg -ts_monitor -log_view :petsc_log.txt)
[ -z "${MAXSTEPS}" ] || PETSC_ARGS+=(-ts_max_steps "${MAXSTEPS}")

cd "${OUT}" || exit 1
echo "=== ${TAG}: start $(date -u +%T) ==="
start=$(date +%s.%N)
ibrun "${BIN}" bp5.toml --petsc "${PETSC_ARGS[@]}" > run.log 2>&1
rc=$?
end=$(date +%s.%N)
echo "${rc}" > exit_code
steps=$(grep -cE '^ *[0-9]+ TS dt ' run.log)
wall=$(python3 -c "print(f'{${end}-${start}:.1f}')")
echo "wall ${wall} s, ts_monitor lines ${steps}" | tee wall.txt
echo "=== ${TAG}: exit ${rc} at $(date -u +%T) ==="
exit "${rc}"
