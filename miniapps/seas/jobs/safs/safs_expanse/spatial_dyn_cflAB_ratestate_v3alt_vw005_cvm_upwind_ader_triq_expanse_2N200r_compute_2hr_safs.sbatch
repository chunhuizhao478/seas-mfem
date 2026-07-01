#!/bin/bash
#SBATCH -J safs_cflAB_rs_alt_exp
#SBATCH -o safs_cflAB_rs_alt_exp_%j.out
#SBATCH -e safs_cflAB_rs_alt_exp_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=100
#SBATCH --time=02:00:00
#SBATCH --mem=200000M
#SBATCH -A lbl107

# =========================================================================
# SAFS RATESTATE v3.0.0 ALT (vw005) — A/B "SEISSOL-CFL" FINGERPRINT (Lever A).
# *** EXPANSE (SDSC) — 2N x 100r = 200 ranks / 2 h / 200 GB-per-node. ***
#
# WHAT THIS A/B TESTS -----------------------------------------------------
#   SAFS carries an EXTRA x3 DG safety factor on top of the mandatory Cockburn-
#   Shu order factor 1/(2N+1), so its explicit ADER step is ~3x SMALLER than the
#   validated SeisSol default (GlobalTimestep.cpp:46: dt = cfl*2r/((2N+1)*cp)).
#   Lever A exposes that extra factor as [numerics].cfl_dg_safety (default 3.0 =
#   today, byte-exact; 1.0 = SeisSol-equivalent, dt x3).  The Courant number
#   `cfl` is a SEPARATE knob and STAYS 0.5 (<= 1) in BOTH arms — we do NOT raise
#   it.  The x3 comes ONLY from cfl_dg_safety, via the new --cfl-dg-safety flag.
#
#   * Arm A (baseline):  SAFS_CFL_DG_SAFETY=3.0  -> dt = today's dt (reference).
#   * Arm B (SeisSol):   SAFS_CFL_DG_SAFETY=1.0  -> dt x3 (fewer, larger steps).
#   Everything else is IDENTICAL (same config, mesh, ranks, flux, integrator,
#   opt stack, tfinal, snapshot schedule) -> the ONLY difference is the step.
#
#   SUBMIT BOTH ARMS (outputs are job-ID + cfl_dg_safety keyed, never collide):
#       SAFS_CFL_DG_SAFETY=3.0 sbatch spatial_dyn_cflAB_..._safs.sbatch   # A
#       SAFS_CFL_DG_SAFETY=1.0 sbatch spatial_dyn_cflAB_..._safs.sbatch   # B
#
# WHAT TO COMPARE ---------------------------------------------------------
#   (1) SPEEDUP: the driver logs "[time] nsteps = N" and "[time] dt_cfl = ...".
#       Expect nsteps_A / nsteps_B ~= 3 and dt_cfl_B ~= 3 x dt_cfl_A.  Since per-
#       macro-step cost is dt-independent, wall-to-tfinal (or sim-time reached in
#       a fixed 2 h wall) should scale the same ~3x.
#   (2) FINGERPRINT (C9, the SOLE science arbiter): diff the two runs' ParaView
#       fault output at MATCHED sim-times (both arms share tfinal + snapshot
#       schedule, so snapshots land at identical sim-times over the overlap).
#       Compare peak slip-rate, rupture-front arrival, and dip-slip.  Bank
#       cfl_dg_safety=1.0 ONLY if the fingerprint is clean; else the largest
#       clean value is the ceiling (sweep 3.0 -> 2.0 -> 1.5 -> 1.0).
#
# ****  MUST REBUILD THE BINARY ON EXPANSE FROM THIS BRANCH  ****
#   --cfl-dg-safety is NEW (2026-07-01).  A binary built before it SILENTLY
#   ignores the flag (argv-scan) and runs cfl_dg_safety=3.0 for BOTH arms ->
#   arm B == arm A with NO error.  DETECTOR: the driver logs a
#   "[time] cfl = 0.5 ... cfl_dg_safety = <v>" line; if that line is ABSENT
#   from the log the binary is stale — rebuild:  bash <root>/build_expanse.sh.
#   The post-run check below greps for it and WARNS if missing.
#
# WALL vs SIM: tfinal=100 s will NOT complete in 2 h at this scale (ADER is CFL-
#   bound).  Arm A advances ~40-55 s then CHECKPOINTS; arm B (dt x3) advances
#   ~3x further in the same wall.  For a CLEAN completed-wall speedup instead,
#   submit BOTH arms with SAFS_TFINAL=<small, e.g. 10s> so arm A finishes; then
#   compare total wall.  Otherwise compare sim-time-reached-in-2h (ratio ~= 3)
#   plus the fingerprint on the [0, t_A] overlap (covers nucleation + breakout).
#
# Config: ${SAFS_ROOT}/config/spatial_friction_safs_seisol_v3_0_0_RSSRW_ALT_vw005.toml
# NOTE (ALT nucleation): the v3 ALT deck's Tnuc center was re-snapped on-fault
#   (project memory project_safs_alt_nucleation_offfault_recenter); confirm the
#   log shows nucleation (V_max rises) — a non-nucleating run has no fingerprint.
# =========================================================================
# RS-SRW vw005 (f_0=0.6, f_w=0.10, Dc=0.10, V_w=0.05->3.0), multiscale_statewise CVM material,
# spatially-varying CSM effective stress, mixed_flux=none (Godunov upwind) + time_integrator=ader
# (the SeisSol-faithful, central-flux-free combo).  cfl=0.5 (Courant number) is UNCHANGED.
#
# ---- EXPANSE contract (matches build_expanse.sh) ----------------------------
#   * Launcher : srun --mpi=pmi2 -n $SLURM_NTASKS   (there is NO ibrun on Expanse)
#   * Toolchain: cpu/0.17.3b + gcc/10.2.0 + openmpi/4.1.x  (NOT intel/impi).
#   * Runtime libs: libhdf5.so + libzfp.so live under extern/*/install (below).
#   * Sizing : compute (128-core/256 GB nodes), 2N x 100r = 200 ranks, 200 GB/node,
#     2 h.  SAME 200-rank topology as the caliper baseline 51443560.
#   * Account : EDIT the -A line above, or submit with --account=<your_alloc>.
#
# Overrides: SAFS_CFL_DG_SAFETY (default 3.0; THE A/B knob), SAFS_CONFIG,
# SAFS_MESH, SAFS_MIXED_FLUX (default none), SAFS_TIME_INTEGRATOR (default ader),
# SAFS_TFINAL (default 100s), SAFS_CFL (Courant number; default config 0.5),
# SAFS_SHARED_CK / SAFS_DERIV_CACHE / SAFS_FACE_CACHE (1=on,0=off; all default ON),
# SAFS_DIAG (0/1, default 0), SAFS_RESTART, SAFS_OUTDIR (qstore base).

export LC_ALL=C
export LANG=C
if [[ "${SAFS_DIAG:-0}" == "1" ]]; then export SEAS_DIAG_BLOWUP=1; else unset SEAS_DIAG_BLOWUP; fi
if [[ "${SAFS_R101:-0}" == "1" ]]; then unset SEAS_R101_SKIP; else export SEAS_R101_SKIP=1; fi

# --- Expanse modules (MUST match build_expanse.sh: cpu/0.17.3b gcc openmpi) ---
module purge
module load slurm 2>/dev/null || true
module load cpu/0.17.3b
module load gcc/10.2.0
for _mpi in openmpi/4.1.5 openmpi/4.1.3 openmpi/4.1.1 openmpi; do
    if module load "${_mpi}" 2>/dev/null && module is-loaded "${_mpi}" 2>/dev/null; then
        echo "OpenMPI module: ${_mpi}"; break
    fi
done

# Resolve seas-mfem root (walk UP to miniapps/seas/Makefile).
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR}}"
_walk_up_dir="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5; do
    if [[ -f "${_walk_up_dir}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"; break
    fi
    _walk_up_dir="${_walk_up_dir}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}."
    echo "  Set SEAS_MFEM_ROOT via sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi

export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"

cd "${SEAS_MFEM_ROOT}"
grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not found.  Build it FIRST on Expanse from"
    echo "  THIS branch (carries --cfl-dg-safety):  bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    exit 1
fi

SAFS_ROOT=safs/project_7.0_alternative
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_safs_seisol_v3_0_0_RSSRW_ALT_vw005.toml}"

DEFAULT_MESH="${SAFS_ROOT}/experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh"
MESH="${SAFS_MESH:-${DEFAULT_MESH}}"
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: mesh not found: $(pwd)/${MESH}  (gitignored — scp it first)"; exit 1
fi

# PRODUCTION SIDECAR PRE-FLIGHT — CVM material + CSM stress; both .h5 required.
VEL_H5="${SAFS_ROOT}/velocity/results/multiscale_statewise_cvm/velocity_safs.h5"
STR_H5="${SAFS_ROOT}/stress/results/csm_d12_mfem/safs_stress_csm_d12_mfem.h5"
for _h5 in "${VEL_H5}" "${STR_H5}"; do
    if [[ ! -f "${_h5}" ]]; then
        echo "ERROR: required sidecar missing: $(pwd)/${_h5}  (gitignored — scp it first)."
        exit 1
    fi
done

grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null \
    || { echo "ERROR: ${CONFIG_TOML} has no sigma_n_strength_floor_pa"; exit 1; }

# ===================== THE A/B KNOB: cfl_dg_safety ========================
# Default 3.0 (arm A, byte-exact today).  Arm B: SAFS_CFL_DG_SAFETY=1.0 (dt x3).
# Guard >= 1.0 (1.0 = SeisSol's order-factor floor; the driver also re-checks).
SAFS_CFL_DG_SAFETY="${SAFS_CFL_DG_SAFETY:-3.0}"
if ! awk -v v="${SAFS_CFL_DG_SAFETY}" 'BEGIN{exit !(v+0 >= 1.0)}'; then
    echo "ERROR: SAFS_CFL_DG_SAFETY='${SAFS_CFL_DG_SAFETY}' must be >= 1.0 (1.0 = SeisSol floor)."; exit 1
fi
CFL_DG_ARG="--cfl-dg-safety ${SAFS_CFL_DG_SAFETY}"
CFLDG_TAG="cfldg${SAFS_CFL_DG_SAFETY//./p}"   # e.g. 3.0 -> cfldg3p0, 1.0 -> cfldg1p0
# =========================================================================

SAFS_MIXED_FLUX="${SAFS_MIXED_FLUX:-none}"
case "${SAFS_MIXED_FLUX}" in
    none|adjacent|all_continuous) : ;;
    *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid"; exit 1 ;;
esac
MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"

SAFS_TIME_INTEGRATOR="${SAFS_TIME_INTEGRATOR:-ader}"
case "${SAFS_TIME_INTEGRATOR}" in
    ader|rk4|rk45) : ;;
    *) echo "ERROR: SAFS_TIME_INTEGRATOR='${SAFS_TIME_INTEGRATOR}' invalid (ader|rk4|rk45)"; exit 1 ;;
esac
if [[ "${SAFS_MIXED_FLUX}" != "none" && "${SAFS_TIME_INTEGRATOR}" == "ader" ]]; then
    echo "ERROR: mixed_flux=${SAFS_MIXED_FLUX} + ADER is UNSTABLE (driver aborts). Use rk4/rk45."; exit 1
fi
TI_ARG="--time-integrator ${SAFS_TIME_INTEGRATOR}"
# cfl_dg_safety is the ADER Dg-path lever; on rk4/rk45 it is INERT (RkCflFactor
# governs).  Keep the A/B on ADER (the default + the baseline scheme).
if [[ "${SAFS_TIME_INTEGRATOR}" != "ader" ]]; then
    echo "WARNING: cfl_dg_safety only affects the ADER Dg path; SAFS_TIME_INTEGRATOR=${SAFS_TIME_INTEGRATOR} makes arm B == arm A."
fi

OPT_ARGS=""
if [[ "${SAFS_TIME_INTEGRATOR}" == "ader" && "${SAFS_SHARED_CK:-1}" == "1" ]]; then
    OPT_ARGS="${OPT_ARGS} --shared-ck-recursion"
fi
[[ "${SAFS_DERIV_CACHE:-1}" == "1" ]] && OPT_ARGS="${OPT_ARGS} --deriv-cache"
if [[ "${SAFS_TIME_INTEGRATOR}" == "ader" && "${SAFS_FACE_CACHE:-1}" == "1" ]]; then
    OPT_ARGS="${OPT_ARGS} --face-cache"
fi

TFINAL="${SAFS_TFINAL:-100s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"   # Courant number (default config 0.5)

RESTART_ARG=""
[[ -n "${SAFS_RESTART:-}" ]] && RESTART_ARG="--restart ${SAFS_RESTART}"

OUT_BASE="${SAFS_OUTDIR:-${SEAS_MFEM_ROOT}/miniapps/seas/safs}"
OUT="${OUT_BASE}/safs_dyn_cflAB_rs_alt_${CFLDG_TAG}_exp_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs/safs_expanse"
LOG="${LOG_DIR}/spatial_dyn_cflAB_rs_alt_${CFLDG_TAG}_exp_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Queue/size:  compute, 2N x 100r = ${SLURM_NTASKS} ranks, 200 GB/node, 02:00:00"
echo "EXPERIMENT:  A/B SEISSOL-CFL (Lever A) — cfl_dg_safety=${SAFS_CFL_DG_SAFETY} (${CFLDG_TAG})"
echo "             arm A=3.0 (baseline), arm B=1.0 (SeisSol-equiv, dt x3); cfl stays 0.5."
echo "Config:      ${CONFIG_TOML}"
echo "Mesh:        ${MESH}"
echo "Mixed flux:  ${SAFS_MIXED_FLUX}    Integrator: ${SAFS_TIME_INTEGRATOR}"
echo "Opt levers:  ${OPT_ARGS:-<none>}"
echo "tfinal:      ${TFINAL}    cfl(Courant): ${CFL_ARG:-<config default 0.5>}    ${CFL_DG_ARG}"
echo "Restart:     ${RESTART_ARG:-<none — fresh start from t=0>}"
echo "Output dir:  ${OUT}"

srun --mpi=pmi2 -n "${SLURM_NTASKS}" stdbuf -oL -eL ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    ${MIXFLUX_ARG} \
    ${TI_ARG} \
    ${OPT_ARGS} \
    --tfinal "${TFINAL}" \
    ${CFL_ARG} \
    ${CFL_DG_ARG} \
    ${RESTART_ARG} \
    --paraview \
    --paraview-fault-vtu \
    --paraview-free-surface vtu \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

# --- STALE-BINARY / flag-honored check: the [time] cfl_dg_safety line must exist.
if grep -qE '\[time\][[:space:]]+cfl[[:space:]]*=.*cfl_dg_safety[[:space:]]*=' "${LOG}"; then
    echo "cfl_dg_safety line: $(grep -m1 -E '\[time\][[:space:]]+cfl[[:space:]]*=.*cfl_dg_safety' "${LOG}")"
    echo "dt_cfl / nsteps:    $(grep -m1 -E '\[time\][[:space:]]+dt_cfl' "${LOG}")  |  $(grep -m1 -E '\[time\][[:space:]]+nsteps' "${LOG}")"
else
    echo "WARNING: no '[time] cfl = ... cfl_dg_safety = ...' line in ${LOG}."
    echo "         The binary may be STALE (pre --cfl-dg-safety) — arm B may equal arm A."
    echo "         Rebuild on Expanse: bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
fi

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG} (check [derived] lines + any guard ABORT)."; exit ${RC}
fi
echo "driver_rc=${RC} — done.  A/B: diff arm A (cfldg3p0) vs arm B (cfldg1p0) ParaView output"
echo "  at matched sim-times (fingerprint), and compare nsteps / dt_cfl / wall (speedup ~3x)."
exit ${RC}
