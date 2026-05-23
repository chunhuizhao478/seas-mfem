#!/bin/bash
#SBATCH -J safs_resDc2diag
#SBATCH -o safs_resDc2diag_%j.out
#SBATCH -e safs_resDc2diag_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# CROSS-RANK DIP-CHANNEL DIAGNOSTIC run (R-DIP2).
# =========================================================================
# IDENTICAL simulation to spatial_dyn_resolution_Dc2_8N_400r_dev_2hr_safs.sbatch
# (same config / mesh / driver flags) EXCEPT the runtime env-gated cross-rank
# trace SEAS_DIAG_XRANK is turned ON (see the export block below).  Purpose:
# diagnose WHY the R-101 guard trips at t~0.455 s on a cross-rank divergence of
# V1 (dip slip-rate) at the shared fault QP (607518, 3706359, -4543).
#
# The run is EXPECTED to abort at the R-101 guard (~t=0.455 s) — that is fine;
# the [XRANK] trace fires on every sub-step up to the abort, so the FULL
# divergence history is captured.  The script extracts the trace into a
# separate file and splits it by the (two) ranks that own the diverging QP, so
# you can diff them directly:
#
#   diff <per-rank-A>.log <per-rank-B>.log   (or paste/awk step-by-step)
#
# How to read the diff — find the first sub-step where the two ranks differ and
# which field goes first:
#   can_t1 differs              => H-C  (canonical frame not bit-identical on the
#                                        curvilinear face)
#   frame same, I*_can differ   => H-B  (bulk-Q FP seed; self vs nbr = local vs
#                                        ghost) — then V1 growth/step = the gain
#   tau1_0 / tau1_nuc differ    => H-A  (prestress/nucleation gap)
#
# Requires commit a9bd4d2 or later on branch `safs` (the SEAS_DIAG_XRANK trace +
# the non-aborting VerifySharedFaultDOFDataConsistency mode).

export LC_ALL=C
export LANG=C

# ---- R-DIP2 cross-rank diagnostic flags (ALL baked in here) -------------
# Enable the env-gated [XRANK] trace.  Default-off in the code => only this job
# (and any that export SEAS_DIAG_XRANK) pays any cost; production jobs are
# byte-exact.
export SEAS_DIAG_XRANK=1
# Target the diverging shared fault QP reported by the R-101 abort (job
# 7746601).  Change if a future run reports a different centroid.
export SEAS_DIAG_XRANK_QP="607518,3706359,-4543"
# Capture radius [m] around the target (catches the face's ~3 QPs; the 500 m
# mesh's nearest QPs are well inside this).
export SEAS_DIAG_XRANK_R=200
# Earliest sim-time [s] to print (0 => full history, incl. the bit-identical
# pre-divergence phase so the FIRST divergence step is visible).  Raise to e.g.
# 0.3 if the log volume is too large.
export SEAS_DIAG_XRANK_TMIN=0.0
# -------------------------------------------------------------------------

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

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'   config/config.mk || { echo "ERROR: no HDF5"; exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: no H5Z_ZFP"; exit 1; }

cd miniapps/seas
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

# tfinal default 2 s (covers the divergence + the t~1.0 s blow-up window);
# cfl default = config (0.5).
TFINAL="${SAFS_TFINAL:-2s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_resDc2diag_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_resDc2diag_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}"
echo "Config:      ${CONFIG_TOML}   (RESOLUTION test: D_c=2.0 / radius=4000)"
echo "XRANK diag:  SEAS_DIAG_XRANK=${SEAS_DIAG_XRANK} QP=${SEAS_DIAG_XRANK_QP} R=${SEAS_DIAG_XRANK_R} TMIN=${SEAS_DIAG_XRANK_TMIN}"
echo "tfinal:      ${TFINAL}"
echo "CFL arg:     ${CFL_ARG:-<config default = 0.5>}"

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --no-sidecar-material \
    ${MIXFLUX_ARG} \
    --tfinal "${TFINAL}" \
    ${CFL_ARG} \
    --paraview \
    --paraview-fault-hdf5 \
    --paraview-fault-zfp-tol 1.0e-12 \
    --paraview-max-snapshots 200 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

# --- Extract the [XRANK] cross-rank trace for offline diffing -------------
# Runs regardless of RC (the R-101 abort is EXPECTED; the trace is the point).
XRANK_LOG="${LOG_DIR}/xrank_resDc2diag_${SLURM_JOB_ID}.log"
grep -h '\[XRANK\]' "${LOG}" > "${XRANK_LOG}" 2>/dev/null || true
NX=$(wc -l < "${XRANK_LOG}" 2>/dev/null || echo 0)
echo "XRANK trace: ${XRANK_LOG}  (${NX} lines)"
if [[ "${NX}" -gt 0 ]]; then
    # Split by the (exactly two) ranks owning the diverging shared QP so the
    # two streams can be diffed directly.
    for r in $(grep -ho 'rank=[0-9][0-9]*' "${XRANK_LOG}" | sort -u); do
        rn=${r#rank=}
        grep -h "${r} " "${XRANK_LOG}" \
            > "${LOG_DIR}/xrank_resDc2diag_${SLURM_JOB_ID}_r${rn}.log"
        echo "  per-rank: ${LOG_DIR}/xrank_resDc2diag_${SLURM_JOB_ID}_r${rn}.log "\
             "($(wc -l < "${LOG_DIR}/xrank_resDc2diag_${SLURM_JOB_ID}_r${rn}.log") lines)"
    done
    echo "  Diff the two per-rank files step-by-step to localise the seed"
    echo "  (can_t1 => H-C; I*_can => bulk seed/H-B; tau1_0/nuc => H-A)."
else
    echo "  WARNING: no [XRANK] lines captured.  Either SEAS_DIAG_XRANK was not"
    echo "  exported into the driver env, the run aborted before the target QP"
    echo "  was processed, or the QP centroid (${SEAS_DIAG_XRANK_QP}) is wrong"
    echo "  for this run — check the R-101 abort centroid in ${LOG}."
fi

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — EXPECTED (R-101 guard abort at the divergence, or"
    echo "  --print-derived).  The [XRANK] trace above is the deliverable; see"
    echo "  ${LOG} for the R-101 abort centroid + the [derived] lines."
    exit 0   # diagnostic job: a guard-abort is success, not failure.
fi

echo "driver_rc=${RC} (ran to completion without tripping the R-101 guard)."
exit 0
