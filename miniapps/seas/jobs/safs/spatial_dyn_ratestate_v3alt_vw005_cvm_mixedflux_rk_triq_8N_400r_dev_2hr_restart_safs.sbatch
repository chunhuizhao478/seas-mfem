#!/bin/bash
#SBATCH -J safs_rs_v3alt_rst
#SBATCH -o safs_rs_v3alt_rst_%j.out
#SBATCH -e safs_rs_v3alt_rst_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# =========================================================================
# SAFS RATESTATE v3.0.0 ALT — MFEM port with CVM material + MIXED FLUX
# (adjacent) + RK4.  Config: config/spatial_friction_safs_seisol_v3_0_0_RSSRW_ALT_vw005.toml
# Plan: safs/project_7.0_alternative/document/PLAN_mixed_flux_seissol_port_2026-06-19.md
# =========================================================================
# RS-SRW vw005 (f_0=0.6, f_w=0.10, Dc=0.10, V_w=0.05->3.0), multiscale_statewise CVM material (the PRODUCTION material), spatially-
# varying CSM d12 effective stress, mixed_flux=adjacent + time_integrator=rk4
# (central flux is UNSTABLE under ADER; RK is required).  CFL 0.5 is a STARTING
# calibration for the RK central-flux de-rate (adjacent 0.6x) — sweep via SAFS_CFL.
#
# RESTART/CONTINUATION (8N/400r/2 h): resumes a 400-rank checkpoint (RESTART_FROM) toward tfinal 150 s. 400-rank checkpoints can ONLY resume at 400 ranks (driver hard-checks num_ranks).
#
# PRODUCTION CHECKLIST (enforced below):
#   * NO --no-sidecar-material  (production loads the CVM sidecar)
#   * [velocity].use_sidecar=true + velocity_safs.h5 present under dataset_root
#   * the CSM stress .h5 present (gitignored/large; scp to Frontera)
#
# Overrides: SAFS_CONFIG, SAFS_MESH, SAFS_MIXED_FLUX (default adjacent),
# SAFS_TIME_INTEGRATOR (default rk4), SAFS_TFINAL (default 150s), SAFS_CFL,
# RESTART_FROM.

export LC_ALL=C
export LANG=C
export SEAS_DIAG_BLOWUP="${SEAS_DIAG_BLOWUP:-1}"   # V_max / sigma_n drift monitor ON

# === USER: RESTART SOURCE (set via --export=ALL,RESTART_FROM=/path/to/run_dir) ===
# Output directory of the run to CONTINUE (holds cp_checkpoint_r*.txt, 400 files).
# The driver reads "<dir>/cp_checkpoint_r<rank>.txt"; the restart OUTPUT dir is
# distinct (job-ID-keyed) so the driver does not abort (restart != output parent).
RESTART_FROM="${RESTART_FROM:-}"
# ================================================================================

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

# Resolve seas-mfem root by walking UP until miniapps/seas/Makefile is found.
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR}}"
_walk_up_dir="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5; do
    if [[ -f "${_walk_up_dir}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"; break
    fi
    _walk_up_dir="${_walk_up_dir}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}."
    echo "  Set SEAS_MFEM_ROOT via sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi
cd "${SEAS_MFEM_ROOT}"

grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
# No in-job build (avoid the concurrent relink race); the binary must be pre-built.
if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not found.  Build it FIRST:"
    echo "   sbatch miniapps/seas/jobs/safs/spatial_dyn_build_safs.sbatch   (or make on a login/idev node)"
    exit 1
fi

SAFS_ROOT=safs/project_7.0_alternative
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_safs_seisol_v3_0_0_RSSRW_ALT_vw005.toml}"

# Mesh: DEFAULT = the sliver-repaired (triq) 500 m mesh; override with SAFS_MESH.
TRIQ_MESH="${SAFS_ROOT}/experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh"
MESH="${SAFS_MESH:-${TRIQ_MESH}}"
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: triq mesh not found: $(pwd)/${MESH}  (gitignored — scp it first)"; exit 1
fi

# PRODUCTION SIDECAR PRE-FLIGHT — this run uses CVM material + the CSM stress
# field; BOTH .h5 sidecars are gitignored/large and MUST be present (scp'd).
VEL_H5="${SAFS_ROOT}/velocity/results/multiscale_statewise_cvm/velocity_safs.h5"
STR_H5="${SAFS_ROOT}/stress/results/csm_d12_mfem/safs_stress_csm_d12_mfem.h5"
for _h5 in "${VEL_H5}" "${STR_H5}"; do
    if [[ ! -f "${_h5}" ]]; then
        echo "ERROR: required sidecar missing: $(pwd)/${_h5}"
        echo "       It is gitignored/large; scp it to Frontera before submitting."
        echo "       (Do NOT work around this with --no-sidecar-material — that forces the"
        echo "        constant fallback, which aborts on the matrix path.)"
        exit 1
    fi
done

# Floor pre-flight (v3 ratestate keeps the 10 MPa MFEM-only sigma_n floor).
grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null     || { echo "ERROR: ${CONFIG_TOML} has no sigma_n_strength_floor_pa"; exit 1; }

# Flux: MIXED (adjacent) by DEFAULT — the production choice for v3 ALT (central-on-
# fault-adjacent, the sigma_n-leak cure).  Requires RK (set below).
SAFS_MIXED_FLUX="${SAFS_MIXED_FLUX:-adjacent}"
case "${SAFS_MIXED_FLUX}" in
    none|adjacent|all_continuous) : ;;
    *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid"; exit 1 ;;
esac
MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"

# Time integrator: RK4 by DEFAULT — mixed/central flux is UNSTABLE under ADER
# (MatrixMixedFluxUnderAder hard guard); RK is mandatory with mixed flux.
SAFS_TIME_INTEGRATOR="${SAFS_TIME_INTEGRATOR:-rk4}"
case "${SAFS_TIME_INTEGRATOR}" in
    ader|rk4|rk45) : ;;
    *) echo "ERROR: SAFS_TIME_INTEGRATOR='${SAFS_TIME_INTEGRATOR}' invalid (ader|rk4|rk45)"; exit 1 ;;
esac
if [[ "${SAFS_MIXED_FLUX}" != "none" && "${SAFS_TIME_INTEGRATOR}" == "ader" ]]; then
    echo "ERROR: mixed_flux=${SAFS_MIXED_FLUX} + ADER is UNSTABLE (driver aborts). Use rk4/rk45."; exit 1
fi
TI_ARG="--time-integrator ${SAFS_TIME_INTEGRATOR}"

TFINAL="${SAFS_TFINAL:-150s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Restart (MANDATORY for this job): continue RESTART_FROM's 400-rank checkpoint.
if [[ -z "${RESTART_FROM}" ]]; then
    echo "ERROR: RESTART_FROM is empty.  Set it to the run dir holding cp_checkpoint_r*.txt:"
    echo "         sbatch --export=ALL,RESTART_FROM=/scratch/.../run_dir <this script>"
    exit 1
fi
RESTART_FROM="${RESTART_FROM%/}"
if [[ ! -f "${RESTART_FROM}/cp_checkpoint_r0.txt" ]]; then
    echo "ERROR: no checkpoint at ${RESTART_FROM}/cp_checkpoint_r0.txt"; exit 1
fi
_n_cp=$(ls "${RESTART_FROM}"/cp_checkpoint_r*.txt 2>/dev/null | wc -l)
if [[ "${_n_cp}" -ne 400 ]]; then
    echo "ERROR: found ${_n_cp} checkpoint files in ${RESTART_FROM} (expected 400; -n 400)."; exit 1
fi
RESTART_ARG="--restart ${RESTART_FROM}/cp"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_rs_v3alt_rst_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_rs_v3alt_rst_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Queue/size:  development, 8N x 400r, 02:00:00"
echo "EXPERIMENT:  ratestate v3 ALT, CVM material + CSM stress, MIXED FLUX + RK"
echo "Config:      ${CONFIG_TOML}"
echo "Mesh:        ${MESH}"
echo "Sidecars:    vel=${VEL_H5}  stress=${STR_H5}"
echo "Mixed flux:  ${SAFS_MIXED_FLUX}    Integrator: ${SAFS_TIME_INTEGRATOR}"
echo "tfinal:      ${TFINAL}    CFL: ${CFL_ARG:-<config default 0.5>}"
echo "Restart:     ${RESTART_ARG:-<none — fresh start from t=0>}"
echo "Output dir:  ${OUT}"

# NOTE: production invocation deliberately OMITS --no-sidecar-material.
ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    ${MIXFLUX_ARG} \
    ${TI_ARG} \
    --tfinal "${TFINAL}" \
    ${CFL_ARG} \
    ${RESTART_ARG} \
    --paraview \
    --paraview-fault-vtu \
    --paraview-free-surface vtu \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

# Banner check (the parser silently ignores unknown flags, so verify the EFFECT):
if grep -qE 'FreeSurf:' "${LOG}"; then
    echo "FreeSurf output: $(grep -m1 -E 'FreeSurf:' "${LOG}")"
else
    echo "WARNING: no 'FreeSurf:' banner line in ${LOG} — free-surface output may be OFF."
fi

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG} (check [derived] lines + any guard ABORT)."; exit ${RC}
fi
echo "driver_rc=${RC} — done.  Inspect ${LOG}: V_max(t) bounded + on-fault sigma_n stable in time."
exit ${RC}
