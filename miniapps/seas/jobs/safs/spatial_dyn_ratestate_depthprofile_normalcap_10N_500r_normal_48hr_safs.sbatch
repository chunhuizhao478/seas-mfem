#!/bin/bash
#SBATCH -J safs_rsDP_cap
#SBATCH -o safs_rsDP_cap_%j.out
#SBATCH -e safs_rsDP_cap_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# SAFS + RATE-AND-STATE (aging, depth-varying a(z)/b(z)) WITH sigma_n CAP
# — NORMAL-QUEUE single-allocation variant (10N/500r/48h).
# =========================================================================
# Same physics/config as the dev-queue sibling
#   spatial_dyn_ratestate_depthprofile_normalcap_8N_400r_dev_2hr_safs.sbatch
# but on the NORMAL queue at 10N/500r/48h.  WHY normal queue:
#   tfinal is 100 s (~285 k macro-steps at dt~3.5e-4 s) — far beyond what the
#   2 h development wall can reach.  Normal/48h + 500 ranks gets as far as
#   possible per allocation; the run checkpoints and is resumed via
#   SAFS_RESTART until it reaches 100 s.
#
# Config (floor enabled):
#   config/spatial_friction_rate_state_safs_projected_stress_depthprofile.toml
#   ([friction].sigma_n_strength_floor_pa = 10.0e6; a(z) <- param_a.csv,
#    (a-b)(z) <- param_a_minus_b.csv, b = a - (a-b); V_0=1e-6, Dc=0.02, f_0=0.6,
#    V_init=1e-12, eta=auto; gradual_overstress nucleation).
#
# CAP CAVEAT (REVIEW R-002): on the RATE-AND-STATE path the floor is applied
# only to the post-solve slip-rate DECOMPOSITION strength; the Brent solve that
# sets the slip-rate MAGNITUDE |V| still uses the unfloored |sigma_n| (plan
# §Phase 3 decision 2, intentional v1).  Tensile |sigma_n| makes RS friction
# STRONGER (the solver locks |V|~0 and the floored decomposition strength is
# never reached), so the cap is largely INERT on RS — RS does not have the LSW
# tensile free-slip failure mode.  This run is for completeness / comparison to
# the no-floor RS depthprofile run; do NOT treat the floor as an active
# tensile-slip safeguard here.
#
# RISK GATE (R-004 / R-020): shared-fault RS psi is 1st-order in time at
# ader_order>=2 and the np>1 RS path is only PARTIALLY validated.  The driver
# prints "Treat parallel RS results as PRELIMINARY" at np>1 — do NOT report
# multi-rank SAFS-RS results as validated until the interior-vs-shared temporal
# drift bound passes.
#
# RUN WINDOW (R-025): the SAFS box is reflection-free only for ~6.94 s
# (min_box_dim/cp) and PML is NOT wired.  tfinal defaults to 100 s (per the
# all-normal-queue=100s directive), so RS seismic fields after ~7 s are
# reflection-CONTAMINATED — treat this as a long cap/stability monitor (does
# V_max stay bounded?), NOT a physical long rupture.  For physical RS fields
# use the dev sibling's 6 s window (or wire a real PML layer).
#
# Flux: uses the config's mixed_flux ("adjacent") by NOT passing --mixed-flux;
# override with SAFS_MIXED_FLUX={none|adjacent|all_continuous}.

export LC_ALL=C
export LANG=C

# Lean production diagnostics: keep the low-overhead V_max trace
# (SEAS_DIAG_BLOWUP) as the cap safety monitor; the verbose per-DOF frame/macro
# and slip traces default OFF (re-enable via env for a debug rerun).
export SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-0}"
export SEAS_DIAG_BLOWUP="${SEAS_DIAG_BLOWUP:-1}"
export SEAS_DIAG_SLIP_VTHR="${SEAS_DIAG_SLIP_VTHR:-10.0}"
export SEAS_NOOPENING="${SEAS_NOOPENING:-0}"
export SEAS_DIAG_FRAME="${SEAS_DIAG_FRAME:-0}"
export SEAS_DIAG_FRAME_XYZ="${SEAS_DIAG_FRAME_XYZ:-607518,3706359,-4543}"
export SEAS_DIAG_FRAME_R="${SEAS_DIAG_FRAME_R:-300}"
export SEAS_DIAG_MACRO="${SEAS_DIAG_MACRO:-0}"

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
        SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"
        break
    fi
    _walk_up_dir="${_walk_up_dir}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR="
    echo "  ${SLURM_SUBMIT_DIR}.  Set SEAS_MFEM_ROOT explicitly via"
    echo "  sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi
cd "${SEAS_MFEM_ROOT}"

# ZFP HDF5 filter plugin (libh5zzfp.so), required at runtime for
# --paraview-*-zfp-tol.
H5Z_ZFP_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install"
if [[ -e "${H5Z_ZFP_BUILT_DIR}/plugin/libh5zzfp.so" ]]; then
    export HDF5_PLUGIN_PATH="${H5Z_ZFP_BUILT_DIR}/plugin:${HDF5_PLUGIN_PATH:-}"
else
    echo "ERROR: ZFP plugin not found at"
    echo "         ${H5Z_ZFP_BUILT_DIR}/plugin/libh5zzfp.so"
    echo "       Build it on Frontera with:  bash build_frontera.sh"
    echo "       (or omit --paraview-*-zfp-tol from the driver invocation)."
    exit 1
fi

grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'   config/config.mk || { echo "ERROR: no HDF5"; exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: no H5Z_ZFP"; exit 1; }

cd miniapps/seas
# This run job does NOT build — the binary must be pre-built (by
# spatial_dyn_build_safs.sbatch, or a login/idev-node `make`).  Building here
# would let concurrently-running jobs race on the same relink (ETXTBSY /
# truncated binary).  Fail fast if it is missing.
if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not found/executable in $(pwd)."
    echo "       Build it FIRST (no in-job build, to avoid the relink race):"
    echo "         sbatch miniapps/seas/jobs/safs/spatial_dyn_build_safs.sbatch"
    echo "       or on a login/idev node:  cd miniapps/seas && make seas_spatial_dyn_driver"
    exit 1
fi

SAFS_ROOT=safs/project_7.0_alternative
# Default = the depth-varying rate-and-state variant (two-CSV a(z)/b(z)) WITH
# [friction].sigma_n_strength_floor_pa = 10.0e6.
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_rate_state_safs_projected_stress_depthprofile.toml}"

# PRE-FLIGHT: this is the CAP run — fail loudly if the config does not enable the
# floor (catches an accidentally-reverted/wrong config).
if ! grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null; then
    echo "ERROR: ${CONFIG_TOML} does not set sigma_n_strength_floor_pa — this is"
    echo "       the NORMAL-STRESS-CAP run and REQUIRES the floor enabled."
    echo "       (It must be a TOP-LEVEL [friction] key, e.g. 10.0e6.)"
    exit 1
fi

# Mesh: read [mesh].path from the config (500 m z0embed); override with SAFS_MESH.
if [[ -n "${SAFS_MESH:-}" ]]; then
    MESH="${SAFS_MESH}"
else
    MESH=$(sed -n '/^\[mesh\]/,/^\[/p' "${CONFIG_TOML}" \
            | grep -E '^[[:space:]]*path[[:space:]]*=' | head -1 \
            | sed -E 's/^[^=]*=[[:space:]]*"?([^"#]*)"?.*/\1/; s/[[:space:]]+$//')
fi
if [[ -z "${MESH}" ]]; then
    echo "ERROR: could not resolve [mesh].path from ${CONFIG_TOML}"
    exit 1
fi
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: SAFS mesh not found:"
    echo "         $(pwd)/${MESH}"
    echo "       The .msh is gitignored and is NOT pulled with the repo."
    echo "       scp it to the path above before submitting this job."
    exit 1
fi

# Depth-profile CSVs (read at config-load time by the parser).  Fail early here
# with a clearer message than a mid-parse abort.
for _csv in \
    "${SAFS_ROOT}/friction/rate-and-state/param_a.csv" \
    "${SAFS_ROOT}/friction/rate-and-state/param_a_minus_b.csv"; do
    if [[ ! -f "${_csv}" ]]; then
        echo "ERROR: depth-profile CSV not found: $(pwd)/${_csv}"
        exit 1
    fi
done

# Flux: default = config's mixed_flux (adjacent) by NOT passing --mixed-flux.
MIXFLUX_ARG=""
if [[ -n "${SAFS_MIXED_FLUX:-}" ]]; then
    case "${SAFS_MIXED_FLUX}" in
        none|adjacent|all_continuous) : ;;
        *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid (none|adjacent|all_continuous)"; exit 1 ;;
    esac
    MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"
fi

# tfinal default 100 s (all-normal-queue directive).  This is past the ~6.94 s
# reflection-free window (R-025) — see the RUN WINDOW note above.  The 48 h wall
# will not reach 100 s in one allocation; the job checkpoints and is resumed via
# SAFS_RESTART.
TFINAL="${SAFS_TFINAL:-100s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Restart hook — resume if the 48 h wall is still hit.  Set SAFS_RESTART to the
# PREVIOUS job's checkpoint PREFIX (driver writes "<output-dir>/cp").
RESTART_ARG=""
[[ -n "${SAFS_RESTART:-}" ]] && RESTART_ARG="--restart ${SAFS_RESTART}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_rsDP_cap_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_rsDP_cap_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Queue/size:  normal, 10N x 500r, 48 h"
echo "EXPERIMENT:  RS depth-profile WITH sigma_n strength floor (10 MPa)"
echo "             (cap largely INERT on RS — see R-002 caveat above)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}"
echo "Config:      ${CONFIG_TOML}   (RS aging; depth-varying a(z)/b(z); floor=10 MPa)"
echo "Diag:        SEAS_DIAG_SLIP=${SEAS_DIAG_SLIP} VTHR=${SEAS_DIAG_SLIP_VTHR} SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP} SEAS_NOOPENING=${SEAS_NOOPENING}"
echo "Diag(frame): SEAS_DIAG_FRAME=${SEAS_DIAG_FRAME} XYZ=${SEAS_DIAG_FRAME_XYZ} R=${SEAS_DIAG_FRAME_R} SEAS_DIAG_MACRO=${SEAS_DIAG_MACRO}"
echo "Mixed flux:  ${SAFS_MIXED_FLUX:-<config default = adjacent>}"
echo "tfinal:      ${TFINAL}  (PAST the R-025 ~6.94 s reflection-free window; long stability monitor)"
echo "CFL arg:     ${CFL_ARG:-<config default = 0.5>}"
echo "Restart:     ${SAFS_RESTART:-<none — fresh start from t=0>}"

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --no-sidecar-material \
    ${MIXFLUX_ARG} \
    --tfinal "${TFINAL}" \
    ${CFL_ARG} \
    ${RESTART_ARG} \
    --paraview \
    --paraview-fault-hdf5 \
    --paraview-fault-zfp-tol 1.0e-12 \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

# Confirm the floor was actually active (banner echo).
if grep -qE 'sigma_n strength floor:[[:space:]]*10[[:space:]]*MPa' "${LOG}"; then
    echo "floor: ACTIVE (banner reports 10 MPa)"
else
    echo "WARNING: did not find 'sigma_n strength floor: 10 MPa' in the banner —"
    echo "         the floor may be DISABLED.  Inspect ${LOG}."
fi

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG} (check the [derived] RS gate lines)."
    exit ${RC}
fi

# Post-run verifier — OPTIONAL (needs h5py, absent on Frontera python3).
module load python3
VERIFIER="${SEAS_MFEM_ROOT}/miniapps/seas/safs/project_7.0_alternative/spatial/code/scripts/verify_spatial_dyn_smoke_safs.py"
CONFIG_TOML_ABS="${SEAS_MFEM_ROOT}/miniapps/seas/${CONFIG_TOML}"
if python3 -c "import h5py" 2>/dev/null; then
    python3 "${VERIFIER}" \
        --fault-vtkhdf "${OUT}/fault.vtkhdf" \
        --toml-config "${CONFIG_TOML_ABS}" \
        --log "${LOG}" \
        > "${LOG_DIR}/verify_summary_rsDP_cap_${SLURM_JOB_ID}.txt"
    CHECK_RC=$?
    echo "driver_rc=${RC}, verify_rc=${CHECK_RC}"
    exit $(( RC | CHECK_RC ))
else
    echo "driver_rc=${RC}, verify=SKIPPED (no h5py on this host)."
    echo "  Run the verifier on a host that has h5py:"
    echo "    python3 ${VERIFIER} \\"
    echo "        --fault-vtkhdf ${OUT}/fault.vtkhdf \\"
    echo "        --toml-config ${CONFIG_TOML_ABS} \\"
    echo "        --log ${LOG}"
    exit ${RC}
fi
