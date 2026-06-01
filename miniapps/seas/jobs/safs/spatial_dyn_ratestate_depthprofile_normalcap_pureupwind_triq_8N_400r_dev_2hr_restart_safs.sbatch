#!/bin/bash
#SBATCH -J safs_rsDP_triqN_pu_rst
#SBATCH -o safs_rsDP_triqN_pu_rst_%j.out
#SBATCH -e safs_rsDP_triqN_pu_rst_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# SAFS RATE-AND-STATE — SLIP LAW + STRONG RATE WEAKENING (TPV104 FL=103,
# depth-varying a(z)/b(z), f_w=0.3 / V_w=0.1 m/s) + sigma_n CAP, PURE UPWIND
# (mixed_flux=none), on the 500 m _triq MESH — RESTART/CONTINUATION
# (development queue, 8N/400r/2 h, tfinal 100 s).
# [Ported from seas-mfem-safs 2026-06-01; MUST use the SAME SRW config as the run
#  it resumes — restarting an SRW checkpoint with the aging config corrupts physics.]
# =========================================================================
# NO MIXED FLUX: central/mixed flux is the confirmed SAFS LSW blow-up driver
# (clean triq-nocap isolation: pure-upwind 7751974 STABLE, identical
# mixed_flux=adjacent 7751975 ran away to V_max~9e6 m/s).  RS regularizes
# friction so it is less prone, but we hold the flux fixed across SW/RS for an
# apples-to-apples comparison.  Analysis:
#   safs/project_7.0_alternative/debug_document/spatial_dynamic_rupture_mixed_flux_analysis_2026-05-26.md
#
# Same physics/config as the 10N/500r normal-queue parent, but RESTARTS a
# 400-rank checkpoint at 8N/400r/2 h dev (the only rank count that matches the
# dev run's per-rank checkpoints; num_ranks is hard-checked in the driver), on the
# re-triangulated (TMOP/retriangulated) mesh
#   experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh
# (overrides the config's [mesh].path via the driver's --mesh flag).
#
# Config (floor enabled):
#   config/spatial_friction_rate_state_safs_projected_stress_srw_Dc010_nuc8km_500m.toml
#   ([friction].sigma_n_strength_floor_pa = 10.0e6; state_evolution="slip_law_srw",
#    f_w_default=0.3, V_w_default=0.1 m/s; Dc=0.10; a(z) <- param_a_500m_strongvw.csv,
#    (a-b)(z) <- param_a_minus_b_500m_strongvw.csv (VW->VS transition ~11 km);
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
# a 6 s window (or wire a real PML layer).  Each 2 h dev wall advances the run a
# little; it checkpoints and is re-resumed by pointing RESTART_FROM at THIS job's
# output dir.
#
# Flux: PURE UPWIND (--mixed-flux none) by DEFAULT (header note);
# override with SAFS_MIXED_FLUX={none|adjacent|all_continuous}.

export LC_ALL=C
export LANG=C

# === USER: RESTART SOURCE (edit this line) ======================================
# Output directory of the run to CONTINUE — the one holding cp_checkpoint_r*.txt
# (400 per-rank files).  The driver reads "<dir>/cp_checkpoint_r<rank>.txt".
# MUST be a 400-rank checkpoint (this job runs -n 400 to match; a 500-rank
# checkpoint aborts with a num_ranks mismatch).  Override at submit time:
#   sbatch --export=ALL,RESTART_FROM=/path/to/run_dir <this script>
RESTART_FROM="${RESTART_FROM:-/scratch2/10024/zhaochun/seas-project/seas-mfem-safs/miniapps/seas/safs/safs_dyn_rsDP_triqN_pu_dev_7752381}"
# ================================================================================

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

grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

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
# Default = the SLIP-LAW + STRONG-RATE-WEAKENING (TPV104 FL=103) depth-varying
# rate-and-state variant (two-CSV a(z)/b(z); state_evolution="slip_law_srw",
# f_w_default=0.3, V_w_default=0.1, Dc=0.10) WITH
# [friction].sigma_n_strength_floor_pa = 10.0e6.
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_rate_state_safs_projected_stress_srw_Dc010_nuc8km_500m.toml}"

# PRE-FLIGHT: this is the CAP run — fail loudly if the config does not enable the
# floor.  (For a cap-OFF mesh-isolation run, set SAFS_CONFIG to a no-floor config.)
if ! grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null; then
    echo "ERROR: ${CONFIG_TOML} does not set sigma_n_strength_floor_pa — this CAP"
    echo "       run REQUIRES the floor.  (For a cap-OFF mesh-isolation run, set"
    echo "       SAFS_CONFIG to a config without the floor.)"
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

# tfinal default 100 s.  This is past the ~6.94 s reflection-free window (R-025)
# — see the RUN WINDOW note above.  Each 2 h dev wall advances the run a little;
# it checkpoints and is re-resumed by pointing RESTART_FROM at THIS job's output
# dir next time.
TFINAL="${SAFS_TFINAL:-100s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Restart (MANDATORY for this job): continue RESTART_FROM's 400-rank checkpoint.
# The driver checkpoint prefix is "<run-dir>/cp" (files cp_checkpoint_r<rank>.txt).
RESTART_FROM="${RESTART_FROM%/}"          # strip any trailing slash
SAFS_RESTART="${RESTART_FROM}/cp"
if [[ ! -f "${RESTART_FROM}/cp_checkpoint_r0.txt" ]]; then
    echo "ERROR: no checkpoint at ${RESTART_FROM}/cp_checkpoint_r0.txt"
    echo "       Set RESTART_FROM (top of this script, or --export=ALL,RESTART_FROM=...)"
    echo "       to the run dir that holds cp_checkpoint_r*.txt."
    exit 1
fi
_n_cp=$(ls "${RESTART_FROM}"/cp_checkpoint_r*.txt 2>/dev/null | wc -l)
if [[ "${_n_cp}" -ne 400 ]]; then
    echo "ERROR: found ${_n_cp} checkpoint files in ${RESTART_FROM} (expected 400)."
    echo "       This job runs -n 400; the checkpoint num_ranks MUST match."
    exit 1
fi
RESTART_ARG="--restart ${SAFS_RESTART}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_rsDP_triqN_pu_rst_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_rsDP_triqN_pu_rst_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Queue/size:  development, 8N x 400r, 2 h (RESTART/continuation)"
echo "Restart from: ${RESTART_FROM}  (${_n_cp} cp files; prefix ${SAFS_RESTART})"
echo "EXPERIMENT:  RS SLIP-LAW + STRONG RATE WEAKENING (TPV104 FL=103: f_w=0.3, V_w=0.1 m/s) + sigma_n floor (10 MPa), PURE UPWIND (mixed_flux=none),"
echo "             on the 500 m _triq mesh, Dc=0.10, 8 km-dia nucleation (no central flux; DEV restart, tfinal 100 s; cap largely INERT on RS — see R-002)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}  (triq / sliver-repaired)"
echo "Config:      ${CONFIG_TOML}   (RS slip-law SRW; depth-varying a(z)/b(z); Dc=0.10; f_w=0.3, V_w=0.1; floor=10 MPa)"
echo "Diag:        SEAS_DIAG_SLIP=${SEAS_DIAG_SLIP} VTHR=${SEAS_DIAG_SLIP_VTHR} SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP} SEAS_NOOPENING=${SEAS_NOOPENING}"
echo "Diag(frame): SEAS_DIAG_FRAME=${SEAS_DIAG_FRAME} XYZ=${SEAS_DIAG_FRAME_XYZ} R=${SEAS_DIAG_FRAME_R} SEAS_DIAG_MACRO=${SEAS_DIAG_MACRO}"
echo "Mixed flux:  ${SAFS_MIXED_FLUX} (pure upwind; CFL factor 1.0)"
echo "tfinal:      ${TFINAL}  (PAST the R-025 ~6.94 s reflection-free window; 2 h dev wall; re-resume via RESTART_FROM for the rest)"
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
    --paraview-volume-vtu \
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
        > "${LOG_DIR}/verify_summary_rsDP_triqN_pu_rst_${SLURM_JOB_ID}.txt"
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
