#!/bin/bash
#SBATCH -J safs_sw_prf_rst
#SBATCH -o safs_sw_prf_rst_%j.out
#SBATCH -e safs_sw_prf_rst_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# =========================================================================
# SAFS SLIPWEAKENING v3.0.0 PREFERRED — MFEM port with CVM material + MIXED FLUX (adjacent) + RK45.
# *** FRONTERA (TACC) job — 10N x 50r = 500 ranks / 48 h / normal queue. ***
# Config: ${SAFS_ROOT}/config/spatial_friction_safs_seisol_v3_0_0_LSW_PREFERRED.toml
# Plan:   safs/project_7.0_alternative/document/PLAN_mixed_flux_seissol_port_2026-06-19.md
# =========================================================================
# LSW (mu_s=0.318, mu_d=0.10, d_c=1.2),  multiscale_statewise CVM material (the PRODUCTION material),
# spatially-varying CSM effective stress, mixed_flux=adjacent + time_integrator=rk45
# (central flux needs an RK integrator — UNSTABLE under ADER; rk45 is FSAL so it
# avoids the np>1 endpoint-Mult issue rk4 hit).  CFL 0.5 = STARTING calibration for
# the RK central-flux de-rate (adjacent 0.6x) — sweep via SAFS_CFL.
#
# ---- FRONTERA contract (matches spatial_dyn_build_safs_frontera.sbatch) ------
#   * Launcher : ibrun   (NOT srun; there is no --mpi=pmi2 on Frontera).
#   * Toolchain: intel/19.1.1 + impi/19.0.9; hypre/2.31.0, mumps/5.3, parmetis,
#     petsc/3.15, fftw3 from TACC modules (NOT from-source like Expanse).
#   * HYPRE pin: TACC_HYPRE_LIB FIRST on LD_LIBRARY_PATH so hypre/2.31.0 is not
#     shadowed by the pre-2.31 libHYPRE that petsc/3.15+mumps/5.3 drag in
#     (the run-time "undefined symbol: HYPRE_Initialize" / MPI exit 127 abort).
#   * Sizing : normal queue, 10 nodes x 50 ranks = 500 ranks, 48 h.  Output goes
#     UNDER the repo: <root>/miniapps/seas/safs (override via SAFS_OUTDIR).
#   * Account : EDIT the -A line above, or submit with
#       sbatch --account=<your_tacc_alloc> <this file>
#
# RESTART/CONTINUATION (normal, 10N x 50r = 500 ranks, 48 h): resumes a 500-rank
# checkpoint (RESTART_FROM) toward tfinal=150 s.  Chain as many of these as needed.
#
# Scale by OVERRIDING the sbatch headers at submit time, e.g.
#   sbatch -N 20 -n 1000 <this file>
# (the restart sibling resumes a checkpoint at the SAME rank count it was written).
#
# PRODUCTION CHECKLIST (enforced below):
#   * NO --no-sidecar-material  (production loads the CVM sidecar)
#   * [velocity].use_sidecar=true + the velocity .h5 present under dataset_root
#   * the CSM stress .h5 present (gitignored/large; scp to Frontera scratch)
#
# Overrides: SAFS_CONFIG, SAFS_MESH, SAFS_MIXED_FLUX (default adjacent),
# SAFS_TIME_INTEGRATOR (default rk45), SAFS_TFINAL (default 150s), SAFS_CFL,
# SAFS_DIAG (0/1, default 0 = quiet log), RESTART_FROM (MANDATORY), SAFS_OUTDIR (repo base).

export LC_ALL=C
export LANG=C
# Per-step [DIAG]/[DIAG-ONSET] monitor (V_max / sigma_n drift).  OFF by default —
# once V_max>10 it prints EVERY step (+ per-DOF onset dumps) and floods the log.
# Set SAFS_DIAG=1 to enable.  (getenv() treats SEAS_DIAG_BLOWUP=0 as "set", so we
# must UNSET it to disable, not set it to 0.)
if [[ "${SAFS_DIAG:-0}" == "1" ]]; then export SEAS_DIAG_BLOWUP=1; else unset SEAS_DIAG_BLOWUP; fi
# R-101 shared-fault cross-rank consistency tripwire — on the CURVILINEAR SAFS
# fault the per-QP strike frame can flip between the two ranks owning a shared
# face, so this check can abort during nucleation (only ~0.02% of fault DOFs sit
# on shared faces).  SKIPPED by default for SAFS; SAFS_R101=1 re-enables it.
# (Needs the SEAS_R101_SKIP gate in the driver — rebuild to pick it up.)
if [[ "${SAFS_R101:-0}" == "1" ]]; then unset SEAS_R101_SKIP; else export SEAS_R101_SKIP=1; fi

# --- Frontera modules (MUST match spatial_dyn_build_safs_frontera.sbatch) ---
module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

# HYPRE 2.31.0 FIRST on the runtime search path (see header) — avoids the
# pre-2.31 libHYPRE shadowing -> "undefined symbol: HYPRE_Initialize" abort.
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"

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
    echo "   sbatch miniapps/seas/jobs/safs/safs_frontera/spatial_dyn_build_safs_frontera.sbatch"
    exit 1
fi

SAFS_ROOT=safs/project_7.0_alternative
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_safs_seisol_v3_0_0_LSW_PREFERRED.toml}"

# Mesh: DEFAULT = the safv4_deep 500 m mesh (safstags tags 101/102/103/104, non-flat topo); override with SAFS_MESH.
DEFAULT_MESH="safs/project_7.0_preferred/meshing/results/safv4_deep_500m_opt_safstags.msh"
MESH="${SAFS_MESH:-${DEFAULT_MESH}}"
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: mesh not found: $(pwd)/${MESH}  (gitignored/large — scp it first)"; exit 1
fi

# PRODUCTION SIDECAR PRE-FLIGHT — this run uses CVM material + the CSM stress
# field; BOTH .h5 sidecars are gitignored/large and MUST be present (scp'd).
VEL_H5="safs/project_7.0_preferred/velocity/results/multiscale_statewise_cvm/velocity_safs.h5"
STR_H5="safs/project_7.0_preferred/stress/results/csm_d12_lsw_mfem/safs_stress_csm_d12_mfem.h5"
for _h5 in "${VEL_H5}" "${STR_H5}"; do
    if [[ ! -f "${_h5}" ]]; then
        echo "ERROR: required sidecar missing: $(pwd)/${_h5}"
        echo "       It is gitignored/large; scp it to Frontera before submitting."
        echo "       (Do NOT work around this with --no-sidecar-material — that forces the"
        echo "        constant fallback, which aborts on the matrix path.)"
        exit 1
    fi
done

# Floor pre-flight (v3 PREFERRED slipweakening keeps the 5 MPa MFEM-only sigma_n floor).
grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null     || { echo "ERROR: ${CONFIG_TOML} has no sigma_n_strength_floor_pa"; exit 1; }

# Flux: MIXED (adjacent) by DEFAULT — the production choice (central-on-fault-
# adjacent, the sigma_n-leak cure).  Requires RK (set below).
SAFS_MIXED_FLUX="${SAFS_MIXED_FLUX:-adjacent}"
case "${SAFS_MIXED_FLUX}" in
    none|adjacent|all_continuous) : ;;
    *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid"; exit 1 ;;
esac
MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"

# Time integrator: RK45 by DEFAULT — mixed/central flux is UNSTABLE under ADER
# (MatrixMixedFluxUnderAder hard guard); RK is mandatory with mixed flux.  rk45
# (Dormand-Prince) is FSAL, so it skips the endpoint re-eval Mult that rk4 needed
# (the np>1 deadlock fix); use rk45 for mixed flux.
SAFS_TIME_INTEGRATOR="${SAFS_TIME_INTEGRATOR:-rk45}"
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

# Restart (MANDATORY for this job): continue RESTART_FROM's checkpoint set.  Set
# RESTART_FROM to the run dir holding cp_checkpoint_r*.txt (one per rank):
#   sbatch --export=ALL,RESTART_FROM=/scratch1/<acct>/<user>/.../run_dir <this file>
# A checkpoint can ONLY resume at the SAME rank count it was written at (the driver
# hard-checks num_ranks; this job is 500 ranks).  The restart OUTPUT dir is distinct
# (job-ID-keyed), so the driver does not abort (restart != output parent).
RESTART_FROM="${RESTART_FROM:-}"
if [[ -z "${RESTART_FROM}" ]]; then
    echo "ERROR: RESTART_FROM is empty.  Set it to the run dir holding cp_checkpoint_r*.txt:"
    echo "         sbatch --export=ALL,RESTART_FROM=/scratch1/.../run_dir <this file>"
    exit 1
fi
RESTART_FROM="${RESTART_FROM%/}"
if [[ ! -f "${RESTART_FROM}/cp_checkpoint_r0.txt" ]]; then
    echo "ERROR: no checkpoint at ${RESTART_FROM}/cp_checkpoint_r0.txt"; exit 1
fi
_n_cp=$(ls "${RESTART_FROM}"/cp_checkpoint_r*.txt 2>/dev/null | wc -l)
if [[ "${_n_cp}" -ne "${SLURM_NTASKS}" ]]; then
    echo "ERROR: found ${_n_cp} checkpoint files in ${RESTART_FROM} (expected ${SLURM_NTASKS} = this job's rank count)."; exit 1
fi
RESTART_ARG="--restart ${RESTART_FROM}/cp"

# Output UNDER the repo; the per-run dir is job-ID-keyed.  Override base with SAFS_OUTDIR.
OUT_BASE="${SAFS_OUTDIR:-${SEAS_MFEM_ROOT}/miniapps/seas/safs}"
OUT="${OUT_BASE}/safs_dyn_sw_prf_rst_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs/safs_frontera"
LOG="${LOG_DIR}/spatial_dyn_sw_prf_rst_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Cluster:     Frontera (TACC)   Launcher: ibrun"
echo "Queue/size:  normal, 10N x 50r = ${SLURM_NTASKS} ranks, 48:00:00"
echo "EXPERIMENT:  slipweakening v3 PREFERRED, CVM material + CSM stress, MIXED FLUX + RK"
echo "Config:      ${CONFIG_TOML}"
echo "Mesh:        ${MESH}"
echo "Sidecars:    vel=${VEL_H5}  stress=${STR_H5}"
echo "Mixed flux:  ${SAFS_MIXED_FLUX}    Integrator: ${SAFS_TIME_INTEGRATOR}"
echo "tfinal:      ${TFINAL}    CFL: ${CFL_ARG:-<config default 0.5>}"
echo "Restart:     ${RESTART_ARG:-<none — fresh start from t=0>}"
echo "Output dir:  ${OUT}"

# NOTE: production invocation deliberately OMITS --no-sidecar-material.
# stdbuf -oL -eL line-buffers the driver stdout/stderr so ${LOG} updates LIVE
# (output-init messages + the per-100-step [DIAG-SIGN] progress) instead of only
# flushing in ~4 KB blocks — so a stalled run is visible immediately.
ibrun stdbuf -oL -eL ./seas_spatial_dyn_driver \
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
