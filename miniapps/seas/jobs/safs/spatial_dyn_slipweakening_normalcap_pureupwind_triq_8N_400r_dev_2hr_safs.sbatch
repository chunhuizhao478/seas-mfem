#!/bin/bash
#SBATCH -J safs_sw_triqN_pu_dev
#SBATCH -o safs_sw_triqN_pu_dev_%j.out
#SBATCH -e safs_sw_triqN_pu_dev_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# LSW Dc2 + sigma_n CAP, PURE UPWIND (mixed_flux=none), on the SLIVER-REPAIRED
# (triq) MESH — DEV SMOKE (development queue, 8N/400r/2 h).
# =========================================================================
# NO MIXED FLUX: the central/mixed flux is confirmed as the SAFS LSW blow-up
# driver — the clean controlled experiment (triq nocap) had pure-upwind
# job 7751974 STABLE to t~15 s while the IDENTICAL mixed_flux=adjacent job
# 7751975 ran away to V_max~9e6 m/s (453 km slip).  This run forces
# --mixed-flux none.  Analysis:
#   safs/project_7.0_alternative/debug_document/spatial_dynamic_rupture_mixed_flux_analysis_2026-05-26.md
#
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
# DEV-QUEUE SCOPE: the 2 h dev wall reaches only ~1.5-2 s of sim (past
# nucleation into breakout), then checkpoints (--checkpoint-every) — a "does it
# start cleanly on pure upwind?" smoke.  Use the 10N/500r/48 h normal-queue
# parent (spatial_dyn_slipweakening_normalcap_pureupwind_triq_10N_500r_normal_48hr_safs.sbatch)
# for the full t>=28 s arrest test.  ACCEPTANCE here is just: clean start + V_max
# bounded through breakout.  tfinal default stays 100 s (the run checkpoints
# early on the dev wall).
#
# REFLECTION CAVEAT: tfinal 100 s > min_box_dim/cp_max (~6.94 s); the 1st-order
# ABC leaks after ~7 s.  Reflections are present in both this run and 7748818,
# so they cancel as a variable for the V_max(t) comparison.  --pml is NOT wired
# (flag only silences the warning).
#
# Overrides: SAFS_CONFIG, SAFS_MESH (default = triq mesh), SAFS_MIXED_FLUX
# (default = none / pure upwind), SAFS_TFINAL (default 100s), SAFS_CFL,
# SAFS_RESTART, SEAS_DIAG_*.

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

# Flux: PURE UPWIND (--mixed-flux none) by DEFAULT — central/mixed flux is the
# confirmed SAFS LSW blow-up driver (header note).  Pure upwind is the most
# dissipative DG flux and carries CFL factor 1.0.  Still overridable via
# SAFS_MIXED_FLUX={none|adjacent|all_continuous}.
SAFS_MIXED_FLUX="${SAFS_MIXED_FLUX:-none}"
case "${SAFS_MIXED_FLUX}" in
    none|adjacent|all_continuous) : ;;
    *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid (none|adjacent|all_continuous)"; exit 1 ;;
esac
MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"

# tfinal default 100 s (matches the normal-queue parent).  On the 2 h dev wall
# the run reaches only ~1.5-2 s, then checkpoints (--checkpoint-every) and is
# resumed via SAFS_RESTART (or continue on the normal-queue parent).
TFINAL="${SAFS_TFINAL:-100s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Restart hook — resume past the 2 h dev wall toward tfinal.  Set SAFS_RESTART to
# the PREVIOUS job's checkpoint PREFIX (driver writes "<output-dir>/cp").
RESTART_ARG=""
[[ -n "${SAFS_RESTART:-}" ]] && RESTART_ARG="--restart ${SAFS_RESTART}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_sw_triqN_pu_dev_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_sw_triqN_pu_dev_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Queue/size:  development, 8N x 400r, 2 h"
echo "EXPERIMENT:  LSW Dc2 + sigma_n floor (10 MPa), PURE UPWIND (mixed_flux=none),"
echo "             on the SLIVER-REPAIRED (triq) mesh (cap + mesh repair + no central flux; tfinal 100 s)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}  (triq / sliver-repaired)"
echo "Config:      ${CONFIG_TOML}   (D_c=2.0 / radius=4000 / floor=10 MPa)"
echo "Diag:        SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP} (V_max/[DIAG-ONSET] trace ON)"
echo "Mixed flux:  ${SAFS_MIXED_FLUX} (pure upwind; CFL factor 1.0)"
echo "tfinal:      ${TFINAL}  (2 h dev wall reaches ~1.5-2 s; checkpoint + SAFS_RESTART for the rest)"
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
    --paraview-fault-vtu \
    --paraview-fault-dt 1.0 \
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
echo "  PASS (V_max bounded, e.g. <1e3, through t>=28 s) -> cap + refined mesh arrest the runaway."
echo "  FAIL (departs plateau ~t=28 s)                   -> neither the cap nor the mesh repair"
echo "                                                       (together) removed the instability."
exit ${RC}
