#!/bin/bash
#SBATCH -J safs_resDc2_upwind
#SBATCH -o safs_resDc2_upwind_%j.out
#SBATCH -e safs_resDc2_upwind_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# CENTRAL-FLUX HYPOTHESIS TEST (pure upwind) — debug 2026-05-25.
# =========================================================================
# Question: job 7748818 (LSW Dc2, mixed_flux=adjacent) ran a clean ~10-12 m/s
# rupture for ~25 s, then V_max climbed off the plateau at t~28 s and ran away
# super-exponentially to 5e64 by t~64.6 s, localizing on mesh-sliver fault
# elements (x~365 km z~-2.5 km, and the x~535 km shallow trace).
#
# Established (spatial_dynamic_rupture_sliver_blowup_2026-05-25.md):
#   - dt IS already set by the worst sliver tet (insphere 42 m -> dt 3.5e-4 s),
#     so the blow-up is NOT a CFL/time-stepping problem.  No dt cut / no LTS
#     fixes a positive-eigenvalue SPATIAL mode.
#   - Leading suspect: the NON-DISSIPATIVE central flux.  job 7748818 used
#     mixed_flux=adjacent -> 200 308 interior faces on central flux.  Central
#     flux removes the dissipation that would damp the sigma_n/[[v_n]] leak on
#     the badly-conditioned sliver elements.  SeisSol uses upwind (dissipative)
#     interior flux + insphere-CFL and is robust.
#
# THIS RUN: identical mesh / config / stress / nucleation / ranks as 7748818,
# the ONLY change is **mixed_flux = none (pure upwind)**.
#   - If V_max STAYS ~10-15 m/s through t~28-35 s (no plateau departure) ->
#     central flux is CONFIRMED as the instability source.
#   - If it still departs at t~28 s -> central flux is NOT it; next look at the
#     ADER predictor/ghost [[v_n]] path on slivers (R-1303/R-1601).
# A from-t=0 run (NOT a restart) on purpose: the only surviving 7748818
# checkpoint is the blown-up one (prefix "cp" has no step in the filename and is
# overwritten every 2000 steps), AND a restart would carry central-flux history.
# From t=0 with upwind is the controlled experiment.
#
# tfinal default 35 s: a few seconds past the original t~28 s onset (and past
# where 7748818 hit V~100 at t~42 s).  If stable it simply runs to 35 s and
# exits early.  SEAS_DIAG_BLOWUP=1 so the [DIAG]/[DIAG-ONSET] V_max trace is
# directly comparable to the 7748818 log.
#
# REFLECTION CAVEAT unchanged: tfinal 35 s > min_box_dim/cp_max (~6.94 s); the
# 1st-order ABC leaks after ~7 s.  Reflections are present in BOTH this run and
# 7748818, so they cancel as a variable — this isolates the flux effect.  --pml
# is NOT wired (flag only silences the warning).  A separate PML/tfinal=6 s run
# quantifies the reflection amplifier.
#
# Overrides: SAFS_MIXED_FLUX (default none here), SAFS_TFINAL (default 35s),
# SAFS_CONFIG, SAFS_MESH, SAFS_CFL, SAFS_RESTART, SEAS_DIAG_*.

export LC_ALL=C
export LANG=C

# DIAGNOSTICS ON for this test (compare V_max(t) to the 7748818 blow-up log).
export SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-0}"
export SEAS_DIAG_SLIP_VTHR="${SEAS_DIAG_SLIP_VTHR:-10.0}"
export SEAS_NOOPENING="${SEAS_NOOPENING:-0}"
export SEAS_DIAG_FRAME="${SEAS_DIAG_FRAME:-0}"
export SEAS_DIAG_FRAME_XYZ="${SEAS_DIAG_FRAME_XYZ:-365100,3838920,-2488}"
export SEAS_DIAG_FRAME_R="${SEAS_DIAG_FRAME_R:-300}"
export SEAS_DIAG_MACRO="${SEAS_DIAG_MACRO:-0}"
# [DIAG]/[DIAG-ONSET] ON (default 1 here) -> V_substep_max / max_slip /
# per-step argmax DOF, the trace used to pinpoint the 7748818 onset.
export SEAS_DIAG_BLOWUP="${SEAS_DIAG_BLOWUP:-1}"

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
# Same config as the 7748818 blow-up (D_c=2.0 / radius=4000 / cfl=0.5).
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml}"

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

# Flux: THIS TEST defaults to pure upwind (mixed_flux=none).  Override with
# SAFS_MIXED_FLUX=adjacent to reproduce the 7748818 blow-up as a control.
SAFS_MIXED_FLUX="${SAFS_MIXED_FLUX:-none}"
case "${SAFS_MIXED_FLUX}" in
    none|adjacent|all_continuous) : ;;
    *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid (none|adjacent|all_continuous)"; exit 1 ;;
esac
MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"

# tfinal default 35 s: a few s past the t~28 s onset; exits early if stable.
TFINAL="${SAFS_TFINAL:-35s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Restart hook (unused for this controlled from-t=0 test; here for parity).
RESTART_ARG=""
[[ -n "${SAFS_RESTART:-}" ]] && RESTART_ARG="--restart ${SAFS_RESTART}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_resDc2_upwind_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_resDc2_upwind_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Queue/size:  normal, 10N x 500r, 48 h"
echo "EXPERIMENT:  central-flux hypothesis test — mixed_flux=${SAFS_MIXED_FLUX} (pure upwind)"
echo "Compare to:  job 7748818 (mixed_flux=adjacent) blow-up at t~28-64 s"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}"
echo "Config:      ${CONFIG_TOML}   (D_c=2.0 / radius=4000)"
echo "Diag:        SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP} (V_max/[DIAG-ONSET] trace ON)"
echo "Mixed flux:  ${SAFS_MIXED_FLUX}"
echo "tfinal:      ${TFINAL}"
echo "CFL arg:     ${CFL_ARG:-<config default = 0.5>}"

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

echo "driver_rc=${RC} — done.  Inspect V_max(t) in ${LOG}:"
echo "  grep -E '^\\[DIAG\\] step' ${LOG} | awk 'NR%200==1{print \$4,\$5}'"
echo "  STABLE (V_max ~10-15 through t>=30 s)  -> central flux CONFIRMED as the source."
echo "  RUNAWAY (departs plateau at t~28 s)    -> central flux NOT it; check ADER/ghost [[v_n]]."
exit ${RC}
