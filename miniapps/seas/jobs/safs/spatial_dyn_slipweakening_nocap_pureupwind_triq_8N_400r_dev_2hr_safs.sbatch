#!/bin/bash
#SBATCH -J safs_sw_triq_nocap_pu
#SBATCH -o safs_sw_triq_nocap_pu_%j.out
#SBATCH -e safs_sw_triq_nocap_pu_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# LSW Dc2, NO sigma_n CAP, PURE UPWIND (mixed_flux=none), triq MESH — dev smoke.
# =========================================================================
# MESH + FLUX isolation (sliver_blowup plan 2026-05-25 §6).  This is the
# pure-upwind sibling of spatial_dyn_slipweakening_nocap_triq_* : same nocap
# triq mesh-isolation setup, but additionally forces --mixed-flux none (NO
# central-flux faces).  Central flux is non-dissipative and the prime suspect
# for the per-sub-step predictor/ghost speckle; pure upwind is the most
# dissipative DG flux (and carries the largest CFL factor, 1.0).  So this run
# answers: does the speckle survive even with the cap OFF AND zero central
# faces?  If stable here but not under adjacent flux -> central flux is the
# amplifier; if still speckles -> it is intrinsic to the predictor/ghost path.
#
#   * Config: ..._resolution_Dc2_nocap.toml  (identical to ..._Dc2.toml EXCEPT
#     [friction].sigma_n_strength_floor_pa is OMITTED => cap OFF).
#   * Mesh:   experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh
#
# Expected reads:
#   - Clean V_max (~10-12 m/s) through breakout and to the t~6.94 s reflection
#     time  => the repaired mesh ALONE tamed the Dc2 instability (cap not
#     needed).  (Compare to 7751832 cap+triq, which was also clean to ~6.9 s.)
#   - Early V_max runaway in the physical window (t < 6.94 s) like 7748818
#     => the mesh repair did NOT remove the instability; the cap was doing
#     the work.
#   NB: reflected-wave contamination after t~6.94 s (1st-order ABC, no PML)
#   blows V_max up regardless of mesh/cap — physical only to ~6.94 s.
#
# DEV-QUEUE SCOPE: 2 h dev reaches ~1.5-2 s of sim (past nucleation T_nuc=1 s,
# into breakout) then checkpoints; this is a "does the refined mesh start
# cleanly WITHOUT the cap?" smoke.  Use the normal-queue sibling for the full
# arrest test.
#
# Overrides: SAFS_CONFIG, SAFS_MESH (default = triq mesh), SAFS_MIXED_FLUX
# (default = config adjacent), SAFS_TFINAL (default 35s), SAFS_CFL, SAFS_RESTART,
# SEAS_DIAG_*.

export LC_ALL=C
export LANG=C

# Diag: V_max trace ON (stability monitor) + [DIAG-SIGN] tensile/floor tracker
# (auto, via SEAS_DIAG_BLOWUP).  Slip trace ON (REVIEW R-002): [DIAG-SIGN]
# localizes WHERE sigma_n goes tensile, but only the per-sub-step slip trace
# emits the sn_vjump/sn_sterm split that says WHY (predictor/ghost OPENING vs
# tensile bulk stress).  VTHR=15 m/s fires it on a speckle DOF (which spikes
# locally above the ~13 m/s legitimate rupture front) without tracing the
# whole front.  Frame point preset to the 7748818 cluster-A blow-up DOF.
export SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-1}"
export SEAS_DIAG_SLIP_VTHR="${SEAS_DIAG_SLIP_VTHR:-15.0}"
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
# NO-CAP config: same D_c=2.0 / radius=4000 / cfl=0.5 as the 7748818 blow-up,
# but WITHOUT [friction].sigma_n_strength_floor_pa (cap OFF).
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2_nocap.toml}"

# PRE-FLIGHT (INVERTED vs the cap run): this is the NO-CAP isolation run — fail
# loudly if the config DOES enable the floor.  (To run cap-ON, use the
# normalcap sibling sbatch instead.)
if grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null; then
    echo "ERROR: ${CONFIG_TOML} sets sigma_n_strength_floor_pa — this is the"
    echo "       NO-CAP mesh-isolation run; the config must NOT enable the floor."
    echo "       Use spatial_dyn_slipweakening_normalcap_triq_*.sbatch for cap-ON,"
    echo "       or point SAFS_CONFIG at a config without the floor key."
    exit 1
fi

# Mesh: DEFAULT = the sliver-repaired (triq) mesh; override with SAFS_MESH.
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

# Flux: this variant forces PURE UPWIND (--mixed-flux none) — NO central-flux
# faces.  Rationale: central flux is non-dissipative and the prime suspect for
# the per-sub-step predictor/ghost speckle; pure upwind is the most dissipative
# (most stable) DG flux and also carries the largest CFL factor (1.0 vs 0.9 for
# adjacent), so this isolates whether the speckle survives without ANY central
# face.  Default = none; still overridable via SAFS_MIXED_FLUX.
SAFS_MIXED_FLUX="${SAFS_MIXED_FLUX:-none}"
case "${SAFS_MIXED_FLUX}" in
    none|adjacent|all_continuous) : ;;
    *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid (none|adjacent|all_continuous)"; exit 1 ;;
esac
MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"

# tfinal default 35 s (the arrest-test target); dev/2 h reaches only ~1.5-2 s.
TFINAL="${SAFS_TFINAL:-35s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Restart hook — resume past the 2 h dev wall toward tfinal.  Set SAFS_RESTART to
# the PREVIOUS job's checkpoint PREFIX (driver writes "<output-dir>/cp").
RESTART_ARG=""
[[ -n "${SAFS_RESTART:-}" ]] && RESTART_ARG="--restart ${SAFS_RESTART}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_sw_triq_nocap_pu_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_sw_triq_nocap_pu_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Queue/size:  development, 8N x 400r, 2 h"
echo "EXPERIMENT:  LSW Dc2 NO cap, PURE UPWIND (mixed_flux=none), triq mesh"
echo "             (MESH-ISOLATION + flux-isolation: cap off AND no central-flux"
echo "              faces — tests if the speckle survives the most dissipative"
echo "              flux; plan 2026-05-25 §6)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}  (triq / sliver-repaired)"
echo "Config:      ${CONFIG_TOML}   (D_c=2.0 / radius=4000 / cap OFF)"
echo "Diag:        SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP} (V_max/[DIAG-ONSET] trace ON)"
echo "Mixed flux:  ${SAFS_MIXED_FLUX} (pure upwind; CFL factor 1.0)"
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

# Confirm the floor was actually OFF (banner echo).  The cap run greps for
# "sigma_n strength floor: 10 MPa"; here we expect it ABSENT/disabled.
if grep -qE 'sigma_n strength floor:[[:space:]]*10[[:space:]]*MPa' "${LOG}"; then
    echo "WARNING: banner reports 'sigma_n strength floor: 10 MPa' — the cap is"
    echo "         ON, but this is the NO-CAP run.  Inspect ${LOG}; the config"
    echo "         override may have leaked the floor."
else
    echo "cap: OFF (no '10 MPa floor' banner — mesh-isolation as intended)"
fi

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG} (probably --print-derived ABORTed; "
    echo "  check the [derived] lines, incl. the L_nuc/h_min resolution on the triq mesh)."
    exit ${RC}
fi

echo "driver_rc=${RC} — done.  Inspect V_max(t) in ${LOG}:"
echo "  grep -E '^\\[DIAG\\] step' ${LOG} | awk 'NR%200==1{print \$4,\$5}'"
echo "  Clean ~10-12 m/s through breakout (to t~6.94 s) -> repaired mesh ALONE"
echo "                                                     tamed Dc2 (cap not needed)."
echo "  Early V_max runaway (t < 6.94 s, like 7748818)   -> mesh did NOT fix it;"
echo "                                                     the cap was doing the work."
exit ${RC}
