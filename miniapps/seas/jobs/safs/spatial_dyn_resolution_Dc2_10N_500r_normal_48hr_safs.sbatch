#!/bin/bash
#SBATCH -J safs_resDc2_norm
#SBATCH -o safs_resDc2_norm_%j.out
#SBATCH -e safs_resDc2_norm_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# NORMAL-queue production long-run — post-R-004-fix (commit 550c640+).
# =========================================================================
# Same physics/mesh/config as the dev resolution-test
# (spatial_dyn_resolution_Dc2_8N_400r_dev_2hr_safs.sbatch):
#   config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml
#   (D_c=2.0 m, radius=4000 m, cfl=0.5, mu_s=0.85, mu_d=0.30, delta_tau=20 MPa,
#    500 m z0embed mesh), 50 ranks/node (10N x 500r), 48 h wall.
#
# Purpose: now that R-004 (shared-fault byNODES ghost read) is fixed and the
# runaway is gone (dev job 7748716: V_substep==V_max, sigma_n~40 MPa,
# max_slip physical), carry the SAFS rupture far in time on the normal queue.
# 48 h reaches MUCH further than the 2 h dev window.
#
# DIAGNOSTICS FULLY QUIET here (production run): SEAS_DIAG_FRAME/_MACRO/_SLIP
# AND SEAS_DIAG_BLOWUP ([DIAG]/[DIAG-ONSET]) all default to 0, so a 48 h run
# emits neither the per-sub-step seed traces nor the per-100-step [DIAG] line.
# The plain "step N/M  t = ... s  V_max = ..." progress line (NOT gated by any
# of these; every 100 steps) still prints t and V_max for monitoring.  Re-enable
# via --export (SEAS_DIAG_BLOWUP=1 -> V_substep_max / max_slip / [DIAG-ONSET];
# SEAS_DIAG_FRAME/_MACRO/_SLIP=1 -> the seed traces).
#
# REFLECTION CAVEAT (physical, NOT a code issue): tfinal=100 s exceeds
# min_box_dim/cp_max (~6.94 s), so the first-order absorbing BC (bottom+sides;
# top is the free surface) leaks reflected waves back into the rupture after
# ~7 s.  NB: --pml is NOT wired into the spatial driver (the flag only silences
# the in-driver warning — it builds NO absorbing layer), so do NOT rely on it.
# For physically-clean late-time use SAFS_TFINAL<=~6 s or a larger box; a real
# PML / higher-order ABC is deferred.  The early rupture (t<~7 s) is clean.
#
# Restart: --checkpoint-every 2000 + the SAFS_RESTART hook let a follow-up job
# resume past the 48 h wall (or resume a dev checkpoint here, and vice versa —
# the OUT dir is job-id-keyed and interchangeable).
#
# Overrides: SAFS_CONFIG, SAFS_MESH, SAFS_MIXED_FLUX{none|adjacent|all_continuous},
# SAFS_TFINAL, SAFS_CFL, SAFS_RESTART, SEAS_DIAG_*.

export LC_ALL=C
export LANG=C

# Per-sub-step seed diagnostics DEFAULT OFF on this production run (see header).
# Zero overhead / byte-exact when 0.  Override with --export=ALL,SEAS_DIAG_FRAME=1,...
export SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-0}"
export SEAS_DIAG_SLIP_VTHR="${SEAS_DIAG_SLIP_VTHR:-10.0}"
export SEAS_NOOPENING="${SEAS_NOOPENING:-0}"
export SEAS_DIAG_FRAME="${SEAS_DIAG_FRAME:-0}"
export SEAS_DIAG_FRAME_XYZ="${SEAS_DIAG_FRAME_XYZ:-607518,3706359,-4543}"
export SEAS_DIAG_FRAME_R="${SEAS_DIAG_FRAME_R:-300}"
export SEAS_DIAG_MACRO="${SEAS_DIAG_MACRO:-0}"
# [DIAG]/[DIAG-ONSET] CLOSED by default (per request).  The ungated
# "step N/M t=.. V_max=.." line still prints t/V_max every 100 steps.
# Set SEAS_DIAG_BLOWUP=1 to restore V_substep_max / max_slip / [DIAG-ONSET].
export SEAS_DIAG_BLOWUP="${SEAS_DIAG_BLOWUP:-0}"

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

JOBS=16

grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'   config/config.mk || { echo "ERROR: no HDF5"; exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: no H5Z_ZFP"; exit 1; }

cd miniapps/seas
# NB: build once on a login/idev node BEFORE submitting (ETXTBSY-safe); this
# in-job make should then be a no-op ("Nothing to be done").
make seas_spatial_dyn_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

SAFS_ROOT=safs/project_7.0_alternative
# Default = the resolution-test variant (D_c=2.0 / radius=4000 / cfl=0.5).
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml}"

# Mesh: read [mesh].path from the config (single source of truth = 500 m
# z0embed); override with SAFS_MESH.
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

# Flux: default = config's mixed_flux (adjacent) by NOT passing --mixed-flux.
MIXFLUX_ARG=""
if [[ -n "${SAFS_MIXED_FLUX:-}" ]]; then
    case "${SAFS_MIXED_FLUX}" in
        none|adjacent|all_continuous) : ;;
        *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid (none|adjacent|all_continuous)"; exit 1 ;;
    esac
    MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"
fi

# tfinal default 100 s (see header REFLECTION CAVEAT); override with SAFS_TFINAL.
TFINAL="${SAFS_TFINAL:-100s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Restart hook — resume a previous run's checkpoint (dev or normal).  Set
# SAFS_RESTART to the PREVIOUS job's checkpoint PREFIX "<output-dir>/cp", e.g.
#   SAFS_RESTART=${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_resDc2_<prevjobid>/cp
# This job writes to its OWN job-id OUT, so the driver's "--restart parent dir
# must differ from --output-dir" guard holds automatically.
RESTART_ARG=""
[[ -n "${SAFS_RESTART:-}" ]] && RESTART_ARG="--restart ${SAFS_RESTART}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_resDc2_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_resDc2_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Queue/size:  normal, 10N x 500r, 48 h"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}"
echo "Config:      ${CONFIG_TOML}   (D_c=2.0 / radius=4000, production long-run)"
echo "Diag:        SEAS_DIAG_SLIP=${SEAS_DIAG_SLIP} VTHR=${SEAS_DIAG_SLIP_VTHR} SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP} SEAS_NOOPENING=${SEAS_NOOPENING}"
echo "Diag(frame): SEAS_DIAG_FRAME=${SEAS_DIAG_FRAME} XYZ=${SEAS_DIAG_FRAME_XYZ} R=${SEAS_DIAG_FRAME_R} SEAS_DIAG_MACRO=${SEAS_DIAG_MACRO}"
echo "Mixed flux:  ${SAFS_MIXED_FLUX:-<config default = adjacent>}"
echo "tfinal:      ${TFINAL}"
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
    --paraview-max-snapshots 200 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG} (probably --print-derived ABORTed; "
    echo "  check the [derived] lines)."
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
        > "${LOG_DIR}/verify_summary_resDc2_${SLURM_JOB_ID}.txt"
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
