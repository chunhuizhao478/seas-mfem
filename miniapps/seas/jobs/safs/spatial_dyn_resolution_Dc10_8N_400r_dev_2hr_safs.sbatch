#!/bin/bash
#SBATCH -J safs_resDc10
#SBATCH -o safs_resDc10_%j.out
#SBATCH -e safs_resDc10_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# Cohesive-zone EXTREME OVER-RESOLUTION probe for the job-7744292/7745103/
# 7745138 (and zerodip) speckle / opening blow-up.
# =========================================================================
# Runs the D_c=10 over-resolution config variant on the 500 m mesh:
#   config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc10.toml
#   (D_c=10.0 m, radius=4000 m, cfl=0.5, mu_s=0.85, mu_d=0.30, delta_tau=20 MPa)
#
# Cloned from the D_c=2.0 resolution job; ONLY D_c differs (2.0 -> 10.0 m).
# It pushes the patch L_nuc/h_min from ~12 to ~60 -- grossly over-resolving
# the cohesive zone on the SAME 500 m mesh.
#
# *** NOT a spontaneous-rupture test — it CANNOT propagate at D_c=10 ***
# r_c (supercritical core) is fixed by radius/delta_tau (~2.70 km) and does
# NOT scale with D_c, while L_nuc grows to ~11.8 km, so r_c/L_nuc ~= 0.23 << 1
# and the patch is far subcritical: it slips while the gradual_overstress ramp
# drives it (T_nuc=1 s) then ARRESTS.  Restoring spontaneous rupture would need
# radius > ~17.5 km (> the fault), so it is intentionally NOT attempted.
# This is a DRIVEN-PATCH probe of the normal-channel (fault-opening) blow-up:
#   bounded   -> the opening blow-up needs an under-resolved cohesive zone OR
#                a fast propagating front (both absent here); NOT intrinsic to
#                a smooth driven patch.
#   blows up  -> the normal-channel opening runaway is intrinsic and
#                resolution-independent (does not need front grid-noise) ->
#                the fix is in the fault normal-stress / opening numerics, not
#                the mesh.  Compare against the D_c=2.0 (resDc2) leg.
# For a clean SPONTANEOUS-rupture resolution test, use resDc2 or the 250 m mesh.
#
# Flux: uses the config's mixed_flux ("adjacent", matching production + the
# baseline 7745103) by NOT passing --mixed-flux.  Override with
# SAFS_MIXED_FLUX={none|adjacent|all_continuous} for an A/B leg.
# tfinal: 2 s (covers the t~1.0 s T_nuc ramp + margin); override with SAFS_TFINAL.

export LC_ALL=C
export LANG=C

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
# Default = the D_c=10 over-resolution variant (D_c=10.0 / radius=4000 / cfl=0.5).
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc10.toml}"

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
# Set SAFS_MIXED_FLUX to force a specific mode for an A/B leg.
MIXFLUX_ARG=""
if [[ -n "${SAFS_MIXED_FLUX:-}" ]]; then
    case "${SAFS_MIXED_FLUX}" in
        none|adjacent|all_continuous) : ;;
        *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid (none|adjacent|all_continuous)"; exit 1 ;;
    esac
    MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"
fi

# tfinal default 2 s (covers the t~1.0 s T_nuc ramp); cfl default = config (0.5).
TFINAL="${SAFS_TFINAL:-2s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_resDc10_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_resDc10_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}"
echo "Config:      ${CONFIG_TOML}   (OVER-RESOLUTION probe: D_c=10.0 / radius=4000; NON-PROPAGATING, r_c/L_nuc~0.23)"
echo "Mixed flux:  ${SAFS_MIXED_FLUX:-<config default = adjacent>}"
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
        > "${LOG_DIR}/verify_summary_resDc10_${SLURM_JOB_ID}.txt"
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
