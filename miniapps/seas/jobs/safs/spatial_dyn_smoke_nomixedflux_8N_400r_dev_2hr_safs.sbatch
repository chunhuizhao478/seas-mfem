#!/bin/bash
#SBATCH -J safs_nomixflux
#SBATCH -o safs_nomixflux_%j.out
#SBATCH -e safs_nomixflux_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# Mixed-flux A/B harness for the job-7744292 salt-and-pepper blow-up.
# ===================================================================
# Clone of spatial_dyn_smoke_8N_400r_dev_2hr_safs.sbatch with ONE knob:
# the interior non-fault flux mode (Zhang et al. 2023 mixed flux).
#
# WHY: the checkerboard speckle is a grid-scale HIGH-FREQUENCY oscillation
# (HFO).  The production config uses mixed_flux = "adjacent" -> a
# NON-DISSIPATIVE central flux on faces adjacent to the fault, exactly
# where the speckle lives; central flux "allows minor HFOs"
# (wave_operator.hpp:51-62).  --mixed-flux none uses the DISSIPATIVE
# upwind/Godunov flux everywhere, which damps grid-scale modes.  This run
# tests whether that dissipation suppresses the speckle.
#
# A/B usage (vary ONE knob at a time, same config + mesh + CFL):
#   default                                  -> --mixed-flux none      (this job)
#   sbatch --export=ALL,SAFS_MIXED_FLUX=adjacent ...   -> the A/B partner
#   sbatch --export=ALL,SAFS_MIXED_FLUX=all_continuous ...
#
# CLEAN ISOLATION NOTE: the committed config is the mu_s=0.85 subcritical
# recalibration (commit 404bdd1), so by default this stacks "no mixed flux"
# ON TOP of the criticality fix.  To isolate mixed flux's role against the
# ORIGINAL blow-up, point SAFS_CONFIG at the pre-fix mu_s=0.65 config and
# run none-vs-adjacent there:
#   sbatch --export=ALL,SAFS_CONFIG=<mu_s=0.65 baseline>.toml ...
# (ask the assistant to regenerate that baseline TOML from git if needed.)
#
# --tfinal 6s: the full onset(t~0.84s)->growth->bend-arrest/recovery(t~4.6s)
# cycle completes well within 6 s, and [0,6s] is byte-identical to the 20 s
# run (tfinal only sets the stop time), so this is a cheaper, complete A/B.
# Override with SAFS_TFINAL if you want the reflection regime (t>7s) too.

export LC_ALL=C
export LANG=C

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

# Resolve seas-mfem root by walking UP until miniapps/seas/Makefile is found
# (same logic as the baseline smoke sbatch).  Set SEAS_MFEM_ROOT to override.
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
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_slip_weakening_safs_projected_stress.toml}"

# Mesh: read [mesh].path from the config (single source of truth); override
# with SAFS_MESH for mesh-refinement A/B.
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

# THE knob this job exists for.  Default = none (upwind everywhere).
MIXED_FLUX="${SAFS_MIXED_FLUX:-none}"
case "${MIXED_FLUX}" in
    none|adjacent|all_continuous) : ;;
    *) echo "ERROR: SAFS_MIXED_FLUX='${MIXED_FLUX}' invalid (none|adjacent|all_continuous)"; exit 1 ;;
esac

# tfinal: default 6 s (covers the full blow-up/recovery cycle); override
# with SAFS_TFINAL=20s for the reflection regime.
TFINAL="${SAFS_TFINAL:-6s}"

# CFL override (default = config cfl=0.5).
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Output + log keyed by job id AND flux mode so A/B runs never collide.
OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_${MIXED_FLUX}_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_${MIXED_FLUX}_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}"
echo "Config:      ${CONFIG_TOML}"
echo "Mixed flux:  ${MIXED_FLUX}   (THIS is the A/B knob; SAFS_MIXED_FLUX)"
echo "tfinal:      ${TFINAL}"
echo "CFL arg:     ${CFL_ARG:-<config default>}"

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --no-sidecar-material \
    --mixed-flux "${MIXED_FLUX}" \
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
        > "${LOG_DIR}/verify_summary_${MIXED_FLUX}_${SLURM_JOB_ID}.txt"
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
