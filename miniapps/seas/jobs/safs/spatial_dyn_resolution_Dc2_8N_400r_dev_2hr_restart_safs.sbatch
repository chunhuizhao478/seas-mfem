#!/bin/bash
#SBATCH -J safs_resDc2_restart
#SBATCH -o safs_resDc2_restart_%j.out
#SBATCH -e safs_resDc2_restart_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# DEV-queue RESTART job — resume a previous SAFS Dc2 run from its checkpoint.
# =========================================================================
# Identical setup to spatial_dyn_resolution_Dc2_8N_400r_dev_2hr_safs.sbatch
# (D_c=2.0 / radius=4000 / cfl=0.5, 500 m mesh, 8N x 400r, 2 h) but REQUIRES
# SAFS_RESTART — it continues an existing run past the 2 h wall toward tfinal
# instead of starting from t=0.
#
# Usage:
#   # use the in-file default checkpoint (SAFS_RESTART_DEFAULT below):
#   sbatch jobs/safs/spatial_dyn_resolution_Dc2_8N_400r_dev_2hr_restart_safs.sbatch
#   # or override the checkpoint at submit time:
#   sbatch --export=ALL,SAFS_RESTART=<prevOUT>/cp \
#          jobs/safs/spatial_dyn_resolution_Dc2_8N_400r_dev_2hr_restart_safs.sbatch
# where <prevOUT> is the previous job's output dir, e.g.
#   ${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_resDc2_<prevjobid>
# (the driver writes/overwrites "<OUT>/cp" every --checkpoint-every steps).
# The checkpoint PREFIX is a USER INPUT — edit SAFS_RESTART_DEFAULT below (or
# pass SAFS_RESTART via --export to override it).  If the resolved prefix has no
# checkpoint files, this job lists the available checkpoints and exits
# (fail-fast, before the build — no wasted node time).  This job writes to its
# OWN job-id OUT, so the driver's "--restart parent dir must differ from
# --output-dir" guard holds automatically; chain job N -> N+1 by pointing
# SAFS_RESTART at job N's <OUT>/cp.
#
# DIAGNOSTICS FULLY QUIET (continuation run): SEAS_DIAG_FRAME/_MACRO/_SLIP AND
# SEAS_DIAG_BLOWUP ([DIAG]/[DIAG-ONSET]) all default to 0.  The ungated
# "step N/M t=.. V_max=.." line still prints t/V_max every 100 steps.  Override
# via --export (SEAS_DIAG_BLOWUP=1 -> V_substep_max/max_slip/[DIAG-ONSET]).
#
# Overrides: SAFS_CONFIG, SAFS_MESH, SAFS_MIXED_FLUX, SAFS_TFINAL, SAFS_CFL,
# SEAS_DIAG_*.

export LC_ALL=C
export LANG=C

# =========================  USER INPUT  =========================
# Checkpoint PREFIX ('<OUT>/cp' of a previous run) to resume from.
# Edit this to point at the run you want to continue, or override at
# submit time with:  sbatch --export=ALL,SAFS_RESTART=<prefix> ...
SAFS_RESTART_DEFAULT="/scratch2/10024/zhaochun/seas-project/seas-mfem-safs/miniapps/seas/safs/safs_dyn_resDc2_7748716/cp"
# ================================================================

# Per-sub-step seed diagnostics DEFAULT OFF on this continuation run.
export SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-0}"
export SEAS_DIAG_SLIP_VTHR="${SEAS_DIAG_SLIP_VTHR:-10.0}"
export SEAS_NOOPENING="${SEAS_NOOPENING:-0}"
export SEAS_DIAG_FRAME="${SEAS_DIAG_FRAME:-0}"
export SEAS_DIAG_FRAME_XYZ="${SEAS_DIAG_FRAME_XYZ:-607518,3706359,-4543}"
export SEAS_DIAG_FRAME_R="${SEAS_DIAG_FRAME_R:-300}"
export SEAS_DIAG_MACRO="${SEAS_DIAG_MACRO:-0}"
# [DIAG]/[DIAG-ONSET] CLOSED by default (per request); the ungated
# "step N/M t=.. V_max=.." line still prints.  SEAS_DIAG_BLOWUP=1 restores it.
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

# Resolve the checkpoint prefix: --export SAFS_RESTART wins, else the in-file
# USER INPUT default (SAFS_RESTART_DEFAULT).  Fail fast (before the build) if
# the resolved prefix has NO checkpoint files, listing the ones that ARE
# available so the resubmit is one copy-paste.
SAFS_RESTART="${SAFS_RESTART:-${SAFS_RESTART_DEFAULT}}"
if [[ -z "${SAFS_RESTART}" ]] || ! ls "${SAFS_RESTART}"* >/dev/null 2>&1; then
    if [[ -z "${SAFS_RESTART}" ]]; then
        echo "ERROR: no checkpoint prefix — SAFS_RESTART is empty and the in-file"
        echo "       SAFS_RESTART_DEFAULT is unset."
    else
        echo "ERROR: no checkpoint files found for prefix:"
        echo "         ${SAFS_RESTART}"
        echo "       (expected files like '${SAFS_RESTART}...')."
    fi
    echo "       Set it to a previous run's checkpoint PREFIX '<OUT>/cp'."
    echo
    echo "Available checkpoints (newest last):"
    _found=0
    for _d in $(ls -dtr "${SEAS_MFEM_ROOT}"/miniapps/seas/safs/safs_dyn_resDc2_*/ 2>/dev/null); do
        if ls "${_d}"cp* >/dev/null 2>&1; then
            echo "    SAFS_RESTART=${_d}cp"
            _found=1
        fi
    done
    [[ ${_found} -eq 0 ]] && echo "    (none found under ${SEAS_MFEM_ROOT}/miniapps/seas/safs)"
    echo
    echo "Resubmit, e.g. (or edit SAFS_RESTART_DEFAULT in the sbatch):"
    echo "    sbatch --export=ALL,SAFS_RESTART=<prefix above> \\"
    echo "           jobs/safs/spatial_dyn_resolution_Dc2_8N_400r_dev_2hr_restart_safs.sbatch"
    exit 1
fi

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

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'   config/config.mk || { echo "ERROR: no HDF5"; exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: no H5Z_ZFP"; exit 1; }

cd miniapps/seas
# Build once on a login/idev node BEFORE submitting (ETXTBSY-safe); this in-job
# make should then be a no-op.
make seas_spatial_dyn_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

SAFS_ROOT=safs/project_7.0_alternative
# Default = the resolution-test variant (D_c=2.0 / radius=4000 / cfl=0.5).
# MUST match the config the checkpoint was written with.
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml}"

# Mesh: read [mesh].path from the config; override with SAFS_MESH.
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

# tfinal default 100 s (the resume target); override with SAFS_TFINAL.
TFINAL="${SAFS_TFINAL:-100s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# SAFS_RESTART is guaranteed set here (checked above).
RESTART_ARG="--restart ${SAFS_RESTART}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_resDc2_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_resDc2_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Queue/size:  development, 8N x 400r, 2 h (RESTART)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}"
echo "Config:      ${CONFIG_TOML}   (D_c=2.0 / radius=4000)"
echo "Diag:        SEAS_DIAG_SLIP=${SEAS_DIAG_SLIP} VTHR=${SEAS_DIAG_SLIP_VTHR} SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP} SEAS_NOOPENING=${SEAS_NOOPENING}"
echo "Diag(frame): SEAS_DIAG_FRAME=${SEAS_DIAG_FRAME} XYZ=${SEAS_DIAG_FRAME_XYZ} R=${SEAS_DIAG_FRAME_R} SEAS_DIAG_MACRO=${SEAS_DIAG_MACRO}"
echo "Mixed flux:  ${SAFS_MIXED_FLUX:-<config default = adjacent>}"
echo "tfinal:      ${TFINAL}"
echo "CFL arg:     ${CFL_ARG:-<config default = 0.5>}"
echo "Restart:     ${SAFS_RESTART}"

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
    echo "driver_rc=${RC} — see ${LOG} (restart read failure shows 'ERROR: restart"
    echo "  prefix ... ' near the top; check SAFS_RESTART points at a valid <OUT>/cp)."
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
    exit ${RC}
fi
