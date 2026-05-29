#!/bin/bash
#SBATCH -J safs_rsDP_triqN_pu
#SBATCH -o safs_rsDP_triqN_pu_%j.out
#SBATCH -e safs_rsDP_triqN_pu_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# SAFS RATE-AND-STATE — SLIP LAW + STRONG RATE WEAKENING (TPV104 FL=103,
# depth-varying a(z)/b(z)) + sigma_n CAP, PURE UPWIND (mixed_flux=none), on the
# 500 m _triq MESH — normal-queue full run (10N/500r/48 h; tfinal 100 s).
# =========================================================================
# NO MIXED FLUX: central/mixed flux is the confirmed SAFS LSW blow-up driver
# (clean triq-nocap isolation: pure-upwind 7751974 STABLE, identical
# mixed_flux=adjacent 7751975 ran away to V_max~9e6 m/s).  RS regularizes
# friction so it is less prone, but we hold the flux fixed across SW/RS for an
# apples-to-apples comparison.  Analysis:
#   safs/project_7.0_alternative/debug_document/spatial_dynamic_rupture_mixed_flux_analysis_2026-05-26.md
#
# This NORMAL-QUEUE job runs the SLIP-LAW + STRONG-RATE-WEAKENING (TPV104 FL=103)
# config on the 500 m _triq mesh: state_evolution="slip_srw", f_w=0.2, V_w=0.1
# m/s, Dc=0.10, a=0.0127, b=0.0261, 8 km-DIAMETER nucleation — see CONFIG_TOML /
# TRIQ_MESH below.  The proven Dc/a/b set the nucleation forcing; f_w/V_w add the
# TPV104-style fast dynamic weakening.  The DEV sibling (8N/400r/2 h) runs the
# IDENTICAL physics; the aging Dc010 normal-queue file is UNCHANGED.  NORMAL
# queue at 10N/500r/48 h, on the shape-preserving _triq mesh
#   experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh
# (overrides the config's [mesh].path via the driver's --mesh flag).
#
# Config (floor enabled):
#   config/spatial_friction_rate_state_safs_projected_stress_srw_Dc010_nuc8km_500m.toml
#   ([friction].sigma_n_strength_floor_pa = 10.0e6; state_evolution="slip_srw",
#    f_w=0.2, V_w=0.1 m/s; Dc=0.10; a(z) <- param_a_500m_strongvw.csv,
#    (a-b)(z) <- param_a_minus_b_500m_strongvw.csv (VW->VS transition ~11 km;
#    CAVEAT R-002: under SRW f_ss->f_w at V>>V_w regardless of (b-a) sign, so the
#    a>b region does NOT guarantee aging-style down-dip arrest — MONITOR extent);
#    gradual_overstress nucleation, 8 km-DIAMETER patch / 4 km e-fold radius).
# => CONFLATES the cap + the SRW law.  To isolate, rerun cap-OFF via SAFS_CONFIG
#    (a config without sigma_n_strength_floor_pa).
#
# CAP CAVEAT (REVIEW R-002): on the RS path the floor touches only the post-solve
# slip-rate DECOMPOSITION strength; the Brent solve setting |V| still uses the
# unfloored |sigma_n|.  Tensile |sigma_n| strengthens RS friction (|V|~0), so the
# cap is largely INERT on RS — RS lacks the LSW tensile free-slip mode.  This run
# is for completeness / mesh comparison, not an active tensile-slip safeguard.
#
# RISK GATE (R-004 / R-020): shared-fault RS psi is 1st-order in time at
# ader_order>=2 and the np>1 RS path is only PARTIALLY validated.  The driver
# prints "Treat parallel RS results as PRELIMINARY" at np>1.
#
# RUN WINDOW (R-025): the SAFS box is reflection-free only for ~6.94 s
# (min_box_dim/cp) and PML is NOT wired.  tfinal defaults to 100 s (per the
# all-normal-queue=100s directive), so RS seismic fields after ~7 s are
# reflection-CONTAMINATED — treat this as a long cap/stability monitor (does
# V_max stay bounded?), NOT a physical long rupture.  For physical RS fields use
# a 6 s window (or wire a real PML layer).  The 48 h wall will not reach 100 s in
# one allocation; the job checkpoints and resumes via SAFS_RESTART.
#
# Flux: PURE UPWIND (--mixed-flux none) by DEFAULT (header note);
# override with SAFS_MIXED_FLUX={none|adjacent|all_continuous}.

export LC_ALL=C
export LANG=C

# Lean diag (match the slip-weakening jobs): the V_max trace (SEAS_DIAG_BLOWUP)
# is the cap/stability monitor; the verbose per-step [SLIP]/[FRAME]/[MACRO]
# traces (frame-orthonormality + predictor/macro investigations, both already
# characterized — NOT live errors) default OFF.  Re-enable via env for an RS
# speckle debug rerun.
export SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-0}"
export SEAS_DIAG_BLOWUP="${SEAS_DIAG_BLOWUP:-1}"
export SEAS_DIAG_SLIP_VTHR="${SEAS_DIAG_SLIP_VTHR:-10.0}"
export SEAS_NOOPENING="${SEAS_NOOPENING:-0}"
export SEAS_DIAG_FRAME="${SEAS_DIAG_FRAME:-0}"
export SEAS_DIAG_FRAME_XYZ="${SEAS_DIAG_FRAME_XYZ:-607518,3706359,-4543}"
export SEAS_DIAG_FRAME_R="${SEAS_DIAG_FRAME_R:-300}"
export SEAS_DIAG_MACRO="${SEAS_DIAG_MACRO:-0}"

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
# Default = the SLIP-LAW + STRONG-RATE-WEAKENING (TPV104 FL=103) variant on the
# 500 m _triq mesh: state_evolution="slip_srw", f_w=0.2, V_w=0.1 m/s, Dc=0.10,
# a=0.0127, b=0.0261 (depth-varying a(z)/b(z), two-CSV), WITH
# [friction].sigma_n_strength_floor_pa = 10.0e6.  DEDICATED config (the aging
# Dc010 sibling and the dev-queue sibling share it).  Plan:
#   debug_document/spatial_dynamic_rupture_slip_law_srw_plan_2026-05-28.md
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_rate_state_safs_projected_stress_srw_Dc010_nuc8km_500m.toml}"

# PRE-FLIGHT: this is the CAP run — fail loudly if the config does not enable the
# floor.  (For a cap-OFF mesh-isolation run, set SAFS_CONFIG to a no-floor config.)
if ! grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null; then
    echo "ERROR: ${CONFIG_TOML} does not set sigma_n_strength_floor_pa — this CAP"
    echo "       run REQUIRES the floor.  (For a cap-OFF mesh-isolation run, set"
    echo "       SAFS_CONFIG to a config without the floor.)"
    exit 1
fi

# Mesh: DEFAULT = the 500 m shape-preserving _triq mesh (1.16 M tets); override
# with SAFS_MESH.  Resolves L_b(dip)=2.51 km at ~5 cells @ 500 m for Dc=0.10.
TRIQ_MESH="${SAFS_ROOT}/experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh"

# Mesh INPUT MODE (R-OOM).  The 500 m _triq mesh has ~1.16 M elements; a serial
# read replicates it on every rank but fits comfortably at 10N/500r, so this job
# DEFAULTS to serial (no presplit needed — unlike the 16 M-tet 250 m mesh):
#   SAFS_MESH_MODE=serial   (DEFAULT) --mesh <.msh>; read + in-memory METIS on
#                           all ranks.  OOM-safe for the 500 m mesh.
#   SAFS_MESH_MODE=presplit           --mesh <prefix>; each rank reads ONLY its
#                           own seas_partition_mesh file (no replication).  Only
#                           needed for much larger meshes; requires pre-built
#                           per-rank files at THIS task count.
SAFS_MESH_MODE="${SAFS_MESH_MODE:-serial}"
case "${SAFS_MESH_MODE}" in
    serial|presplit) : ;;
    *) echo "ERROR: SAFS_MESH_MODE='${SAFS_MESH_MODE}' invalid (serial|presplit)"; exit 1 ;;
esac

# presplit prefix default is keyed by task count so a stale partition set built
# for a different np cannot be silently consumed (the driver re-checks too).
PRESPLIT_PREFIX="${SAFS_MESH_PREFIX:-${SAFS_ROOT}/experimental_mesh_refinement/presplit/safs_500m_triq_np${SLURM_NTASKS}}"

if [[ "${SAFS_MESH_MODE}" == "presplit" ]]; then
    MESH="${SAFS_MESH:-${PRESPLIT_PREFIX}}"
    if [[ ! -f "${MESH}.np" ]]; then
        echo "ERROR: pre-split manifest not found: $(pwd)/${MESH}.np"
        echo "       The 500 m mesh is OOM-safe in SERIAL mode (the default) — just unset"
        echo "       SAFS_MESH_MODE, or run: sbatch --export=ALL,SAFS_MESH_MODE=serial ..."
        echo "       (presplit is only needed for much larger meshes; build per-rank files"
        echo "        at -n ${SLURM_NTASKS} and set SAFS_MESH_PREFIX if you really want it)."
        exit 1
    fi
else
    MESH="${SAFS_MESH:-${TRIQ_MESH}}"
    if [[ ! -f "${MESH}" ]]; then
        echo "ERROR: triq mesh not found:"
        echo "         $(pwd)/${MESH}"
        echo "       The .msh is gitignored and is NOT pulled with the repo."
        echo "       scp the re-triangulated mesh to the path above before submitting"
        echo "       (it lives under ${SAFS_ROOT}/experimental_mesh_refinement/)."
        exit 1
    fi
fi
MESHMODE_ARG="--mesh-mode ${SAFS_MESH_MODE}"

# Depth-profile CSVs (read at config-load time by the parser).  Fail early here
# with a clearer message than a mid-parse abort.
for _csv in \
    "${SAFS_ROOT}/friction/rate-and-state/param_a_500m_strongvw.csv" \
    "${SAFS_ROOT}/friction/rate-and-state/param_a_minus_b_500m_strongvw.csv"; do
    if [[ ! -f "${_csv}" ]]; then
        echo "ERROR: depth-profile CSV not found: $(pwd)/${_csv}"
        exit 1
    fi
done

# Flux: PURE UPWIND (--mixed-flux none) by DEFAULT — central/mixed flux is the
# confirmed SAFS LSW blow-up driver (header note); pure upwind carries CFL
# factor 1.0.  Still overridable via SAFS_MIXED_FLUX={none|adjacent|all_continuous}.
SAFS_MIXED_FLUX="${SAFS_MIXED_FLUX:-none}"
case "${SAFS_MIXED_FLUX}" in
    none|adjacent|all_continuous) : ;;
    *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid (none|adjacent|all_continuous)"; exit 1 ;;
esac
MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"

# tfinal default 100 s (all-normal-queue directive).  This is past the ~6.94 s
# reflection-free window (R-025) — see the RUN WINDOW note above.  The 48 h wall
# will not reach 100 s in one allocation; the job checkpoints and is resumed via
# SAFS_RESTART.
TFINAL="${SAFS_TFINAL:-100s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Restart hook — resume past the 48 h wall toward tfinal.  Set SAFS_RESTART to
# the PREVIOUS job's checkpoint PREFIX (driver writes "<output-dir>/cp").
RESTART_ARG=""
[[ -n "${SAFS_RESTART:-}" ]] && RESTART_ARG="--restart ${SAFS_RESTART}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_rsDP_triqN_pu_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_rsDP_triqN_pu_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Queue/size:  normal, 10N x 500r, 48 h"
echo "EXPERIMENT:  RS SLIP-LAW + STRONG RATE WEAKENING (TPV104 FL=103: f_w=0.2, V_w=0.1 m/s) + sigma_n floor (10 MPa), PURE UPWIND (mixed_flux=none),"
echo "             on the 500 m _triq mesh, Dc=0.10, a=0.0127, b=0.0261, 8 km-dia nucleation (no central flux; full run, tfinal 100 s; cap largely INERT on RS — see R-002)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}  (500 m _triq, 1.16 M tets)"
echo "Mesh mode:   ${SAFS_MESH_MODE}  (serial=read+partition on every rank; presplit=per-rank files, no replication)"
echo "Config:      ${CONFIG_TOML}   (RS slip-law SRW; depth-varying a(z)/b(z); Dc=0.10; 8 km-dia nuc; f_w=0.2, V_w=0.1; floor=10 MPa)"
echo "Diag:        SEAS_DIAG_SLIP=${SEAS_DIAG_SLIP} VTHR=${SEAS_DIAG_SLIP_VTHR} SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP} SEAS_NOOPENING=${SEAS_NOOPENING}"
echo "Diag(frame): SEAS_DIAG_FRAME=${SEAS_DIAG_FRAME} XYZ=${SEAS_DIAG_FRAME_XYZ} R=${SEAS_DIAG_FRAME_R} SEAS_DIAG_MACRO=${SEAS_DIAG_MACRO}"
echo "Mixed flux:  ${SAFS_MIXED_FLUX} (pure upwind; CFL factor 1.0)"
echo "tfinal:      ${TFINAL}  (PAST the R-025 ~6.94 s reflection-free window; long stability monitor; 48 h wall + SAFS_RESTART)"
echo "CFL arg:     ${CFL_ARG:-<config default = 0.5>}"
echo "Restart:     ${SAFS_RESTART:-<none — fresh start from t=0>}"

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    ${MESHMODE_ARG} \
    --no-sidecar-material \
    ${MIXFLUX_ARG} \
    --tfinal "${TFINAL}" \
    ${CFL_ARG} \
    ${RESTART_ARG} \
    --paraview \
    --paraview-fault-vtu \
    --paraview-fault-dt 0.5 \
    --paraview-max-snapshots 300 \
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
    echo "         the floor may be DISABLED.  Inspect ${LOG}."
fi

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG} (check the [derived] RS gate lines)."
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
        > "${LOG_DIR}/verify_summary_rsDP_triqN_pu_${SLURM_JOB_ID}.txt"
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
