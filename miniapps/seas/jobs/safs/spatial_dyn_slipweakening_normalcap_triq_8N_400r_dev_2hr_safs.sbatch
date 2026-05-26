#!/bin/bash
#SBATCH -J safs_sw_triq
#SBATCH -o safs_sw_triq_%j.out
#SBATCH -e safs_sw_triq_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# LSW Dc2 + sigma_n CAP, on the SLIVER-REPAIRED (triq) MESH — dev smoke.
# =========================================================================
# Plan: spatial_dynamic_rupture_sliver_blowup_2026-05-25.md §6 ("rerun on a
# sliver-repaired mesh ... If stable -> the sliver conditioning was the
# amplifier") + §Phase 4 (the sigma_n strength floor).
#
# Two fixes for the job-7748818 LSW Dc2 blow-up are active TOGETHER here:
#   (1) the sigma_n strength floor (10 MPa), baked into
#       config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml;
#   (2) the re-triangulated (TMOP/retriangulated) mesh
#       experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh,
#       which repairs the x~365 km + shallow-trace sliver patches that localized
#       the runaway.
# => This run CONFLATES the two fixes (cap + repaired mesh).  To isolate the
#    MESH effect alone (plan §6), rerun with the floor DISABLED (a config whose
#    [friction].sigma_n_strength_floor_pa is removed / < 0, via SAFS_CONFIG).
#
# DEV-QUEUE SCOPE: at ~1.7 s wall/step-batch (see the normal-queue sibling), a
# 2 h dev allocation reaches only ~1.5-2 s of simulation — past nucleation
# (T_nuc=1 s) and into breakout, but NOT to the t~28 s no-floor onset.  This is
# therefore a "does the refined mesh start cleanly (nucleation + breakout) with
# the cap?" SMOKE.  For the full arrest test to t~35 s use the normal-queue
# sibling spatial_dyn_slipweakening_normalcap_10N_500r_normal_48hr_safs.sbatch
# with SAFS_MESH set to this triq mesh.
#
# REFLECTION CAVEAT: physical only to ~6.94 s (1st-order ABC; PML not wired) —
# moot here since dev reaches only ~1.5-2 s.
#
# Overrides: SAFS_CONFIG, SAFS_MESH (default = triq mesh), SAFS_MIXED_FLUX
# (default = config adjacent), SAFS_TFINAL (default 35s), SAFS_CFL, SAFS_RESTART,
# SEAS_DIAG_*.

export LC_ALL=C
export LANG=C

# Lean diag: V_max trace ON (the cap/stability monitor); verbose frame/macro/slip
# OFF.  Frame point preset to the 7748818 cluster-A blow-up DOF in case the
# refined mesh still destabilizes there.
export SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-0}"
export SEAS_DIAG_SLIP_VTHR="${SEAS_DIAG_SLIP_VTHR:-10.0}"
export SEAS_NOOPENING="${SEAS_NOOPENING:-0}"
export SEAS_DIAG_FRAME="${SEAS_DIAG_FRAME:-0}"
export SEAS_DIAG_FRAME_XYZ="${SEAS_DIAG_FRAME_XYZ:-365100,3838920,-2488}"
export SEAS_DIAG_FRAME_R="${SEAS_DIAG_FRAME_R:-300}"
export SEAS_DIAG_MACRO="${SEAS_DIAG_MACRO:-0}"
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
# Same config as the 7748818 blow-up (D_c=2.0 / radius=4000 / cfl=0.5) — with
# [friction].sigma_n_strength_floor_pa = 10.0e6.
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml}"

# PRE-FLIGHT: this is the CAP run — fail loudly if the config does not enable the
# floor.  (To run cap-OFF for the mesh-isolation experiment, point SAFS_CONFIG at
# a config without sigma_n_strength_floor_pa.)
if ! grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null; then
    echo "ERROR: ${CONFIG_TOML} does not set sigma_n_strength_floor_pa — this CAP"
    echo "       run REQUIRES the floor.  (For a cap-OFF mesh-isolation run, set"
    echo "       SAFS_CONFIG to a config without the floor.)"
    exit 1
fi

# Mesh: DEFAULT = the sliver-repaired (triq) mesh; override with SAFS_MESH.  This
# overrides the config's [mesh].path via the driver's --mesh flag.
TRIQ_MESH="${SAFS_ROOT}/experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh"
MESH="${SAFS_MESH:-${TRIQ_MESH}}"
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: triq mesh not found:"
    echo "         $(pwd)/${MESH}"
    echo "       The .msh is gitignored and is NOT pulled with the repo."
    echo "       scp the re-triangulated mesh to the path above before submitting"
    echo "       (it lives under ${SAFS_ROOT}/experimental_mesh_refinement/)."
    exit 1
fi

# Flux: default = config's mixed_flux (adjacent) by NOT passing --mixed-flux —
# the FLOOR-ONLY (+ repaired-mesh) experiment.  Override with SAFS_MIXED_FLUX.
MIXFLUX_ARG=""
if [[ -n "${SAFS_MIXED_FLUX:-}" ]]; then
    case "${SAFS_MIXED_FLUX}" in
        none|adjacent|all_continuous) : ;;
        *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid (none|adjacent|all_continuous)"; exit 1 ;;
    esac
    MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"
fi

# tfinal default 35 s (the arrest-test target); dev/2 h reaches only ~1.5-2 s and
# checkpoints — use the normal-queue sibling (SAFS_MESH=triq) for the full run.
TFINAL="${SAFS_TFINAL:-35s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Restart hook — resume past the 2 h dev wall toward tfinal.  Set SAFS_RESTART to
# the PREVIOUS job's checkpoint PREFIX (driver writes "<output-dir>/cp").
RESTART_ARG=""
[[ -n "${SAFS_RESTART:-}" ]] && RESTART_ARG="--restart ${SAFS_RESTART}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_sw_triq_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_sw_triq_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Queue/size:  development, 8N x 400r, 2 h"
echo "EXPERIMENT:  LSW Dc2 + sigma_n floor (10 MPa) on the SLIVER-REPAIRED (triq) mesh"
echo "             (CONFLATES cap + mesh repair; dev smoke = nucleation+breakout)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}  (triq / sliver-repaired)"
echo "Config:      ${CONFIG_TOML}   (D_c=2.0 / radius=4000 / floor=10 MPa)"
echo "Diag:        SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP} (V_max/[DIAG-ONSET] trace ON)"
echo "Mixed flux:  ${SAFS_MIXED_FLUX:-<config default = adjacent>}"
echo "tfinal:      ${TFINAL}  (dev reaches ~1.5-2 s then checkpoints)"
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

# Confirm the floor was actually active (banner echo).
if grep -qE 'sigma_n strength floor:[[:space:]]*10[[:space:]]*MPa' "${LOG}"; then
    echo "floor: ACTIVE (banner reports 10 MPa)"
else
    echo "WARNING: did not find 'sigma_n strength floor: 10 MPa' in the banner —"
    echo "         the floor may be DISABLED.  Inspect ${LOG} before trusting the result."
fi

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG} (probably --print-derived ABORTed; "
    echo "  check the [derived] lines, incl. the L_nuc/h_min resolution on the triq mesh)."
    exit ${RC}
fi

echo "driver_rc=${RC} — done.  Inspect V_max(t) in ${LOG}:"
echo "  grep -E '^\\[DIAG\\] step' ${LOG} | awk 'NR%200==1{print \$4,\$5}'"
echo "  Clean ~10-12 m/s through breakout (~1.5-2 s) -> refined mesh starts cleanly with the cap."
echo "  Early V_max departure                        -> the refined mesh did NOT remove the"
echo "                                                   instability in the reachable window."
exit ${RC}
