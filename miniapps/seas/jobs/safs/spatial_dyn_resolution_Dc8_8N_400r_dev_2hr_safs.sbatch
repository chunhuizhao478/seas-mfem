#!/bin/bash
#SBATCH -J safs_resDc8
#SBATCH -o safs_resDc8_%j.out
#SBATCH -e safs_resDc8_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# SPECKLE-RUNAWAY DISCRIMINATOR — large-D_c leg (paired with the D_c=2.0 baseline
# spatial_dyn_resolution_Dc2_8N_400r_dev_2hr_safs.sbatch).
# =========================================================================
# Runs the D_c=8 discriminator config on the 500 m mesh:
#   config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc8.toml
#   (D_c=8.0 m, radius=4000 m, cfl=0.5, mu_s=0.85, mu_d=0.30, delta_tau=20 MPa)
#
# Context: after the cross-rank reconcile fix (job 7747525) the max_slip runaway
# (-> 2.7e6 m, macro V_max ~ 4 m/s) LOCALIZED to ~11 DOFs within ~7.7 km of the
# nucleation patch, all COMPRESSIVE (sigma_n 46-52 MPa => tension R-005d refuted).
# The remaining cause is a predictor overshoot at the slip-onset kink: either
# cohesive-zone under-resolution (R-007d) or within-step non-relaxation (R-006d).
# D_c 2.0->8.0 raises patch L_nuc/h_min ~12 -> ~48 (gentler kink); if the nucleus
# overshoot is a resolution artifact it should shrink/vanish here.
#
# *** REQUIRES the Phase-1 [SLIP] trace (SEAS_DIAG_SLIP, set below). ***
# This run is the discriminator ONLY with the per-sub-step trace.  Compare the
# [SLIP] lines in the DRIVEN ramp window t in [0, T_nuc=1 s] against the D_c=2.0
# baseline leg:
#   nucleus tau_abs overshoot present at D_c=2.0, GONE at D_c=8 -> R-007d
#       (kink/under-resolution): production fix = finer mesh / resolved zone.
#   overshoot PERSISTS across sub-steps at both D_c              -> R-006d
#       (within-step non-relaxation): fix = smaller macro dt / corrected predictor.
#
# *** CAVEAT — D_c=8 IS SUBCRITICAL (r_c/L_nuc ~0.29, vs ~1.14 at D_c=2.0). ***
# The patch slips during the T_nuc=1 s ramp and ARRESTS after (no spontaneous
# propagation; r_c is set by radius/delta_tau, not D_c).  This is intentional and
# does NOT confound the kink probe (the nucleus crosses the kink during the driven
# ramp regardless).  Do NOT read post-ramp (t > 1 s) max_slip as a resolution
# signal -- the propagation regimes differ between the legs.  See the config
# header and document/PLAN_speckle_slip_runaway_2026-05-23.md.
#
# Flux: config default ("adjacent") by NOT passing --mixed-flux; override with
# SAFS_MIXED_FLUX={none|adjacent|all_continuous}.
# tfinal: 2 s (onset ~t=0.6 s + the runaway window); override with SAFS_TFINAL.

export LC_ALL=C
export LANG=C

# Phase-1 speckle-runaway diagnostics (default ON for this discriminator run;
# override by exporting 0).  SEAS_DIAG_SLIP emits the per-sub-step [SLIP] trace
# at nodes with |V| > SEAS_DIAG_SLIP_VTHR (default 10 m/s); SEAS_DIAG_BLOWUP
# emits the [DIAG] line (now incl. V_substep_max, the honest sub-step peak).
export SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-1}"
export SEAS_DIAG_BLOWUP="${SEAS_DIAG_BLOWUP:-1}"
export SEAS_DIAG_SLIP_VTHR="${SEAS_DIAG_SLIP_VTHR:-10.0}"
# No-opening cap EXPERIMENT (default OFF): SEAS_NOOPENING=1 zeros frictional
# slip under tension (sigma_n_total<=0) in SolveLSW_TPV205 -- the decisive test
# of whether the tension free-slide path is the runaway amplifier
# (REVIEW_speckle_tension_analysis test 1).  OFF => byte-exact.
export SEAS_NOOPENING="${SEAS_NOOPENING:-0}"

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
# Default = the D_c=8 discriminator variant (radius=4000 / cfl=0.5).
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc8.toml}"

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

# tfinal default 2 s (covers the t~1.05 s onset); cfl default = config (0.5).
TFINAL="${SAFS_TFINAL:-2s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_resDc8_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_resDc8_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}"
echo "Config:      ${CONFIG_TOML}   (DISCRIMINATOR: D_c=8.0 / radius=4000, SUBCRITICAL)"
echo "Diag:        SEAS_DIAG_SLIP=${SEAS_DIAG_SLIP} VTHR=${SEAS_DIAG_SLIP_VTHR} SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP} SEAS_NOOPENING=${SEAS_NOOPENING}"
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
        > "${LOG_DIR}/verify_summary_resDc8_${SLURM_JOB_ID}.txt"
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
