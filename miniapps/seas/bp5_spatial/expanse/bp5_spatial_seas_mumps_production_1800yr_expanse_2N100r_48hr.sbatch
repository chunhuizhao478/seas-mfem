#!/bin/bash
#SBATCH -J bp5seas_mumps_prod
#SBATCH -o bp5seas_mumps_prod_%j.out
#SBATCH -e bp5seas_mumps_prod_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=100
#SBATCH --time=48:00:00
#SBATCH --mem=243G
#SBATCH -A ddp408

# =========================================================================
# BP5 SPATIAL quasi-dynamic SEAS — PRODUCTION run, DIRECT (MUMPS) solver.
# *** EXPANSE (SDSC) job — 2N x 100r = 200 ranks / 48 h / compute / ddp408. ***
# Driver:  seas_spatial_seas_driver   (the generic QD SEAS driver, --config flag)
# Config:  config/safs_qd/bp5_spatial_seas_verification.toml  ([solver].type="mumps")
# Run:     FULL ~1800-yr BP5-QD seismic cycle, to t_final (NO --max-steps cap),
#          WITH checkpoint/restart so the cycle spans multiple 48 h segments.
# =========================================================================
# This is the Expanse, DIRECT-solver, production sibling of:
#   * bp5_spatial/expanse/bp5_spatial_seas_cgamg_caliper_expanse_2N80r_100steps.sbatch
#       (Expanse, iterative CG-AMG, Caliper PROFILE, first 100 steps) — same
#       modules / launcher / root-resolution scaffold, inverted solver gate.
#   * bp5_spatial/bp5_production_gmres_restart_normal_48hr.sbatch
#       (FRONTERA, iterative GMRES, seas_driver) — same chained-segment restart
#       idea, but that job uses a DIFFERENT driver + DIFFERENT checkpoint flag
#       (--checkpoint-interval / io/checkpoint.hpp).  This job uses the QD
#       driver's own checkpoint path (--checkpoint-every / io/seas_qd_checkpoint.hpp);
#       the two checkpoint formats are NOT interchangeable.
#
# === DIRECT SOLVER ("direct solver for full 1800 yr run") ==================
# "[solver].type = mumps" -> mfem::MUMPSSolver, an exact sparse LU factorization
# (the BP5 golden / parity solver; ARCHITECTURE.md names "MUMPS (or MUMPS_BLR)"
# as the BP5 direct default).  On Expanse only MUMPS is compiled into the QD
# driver (build_expanse.sh builds it from source via PETSc --download-mumps;
# SuperLU_DIST and STRUMPACK are NOT built), so this job HARD-CHECKS:
#   (1) config.mk has MFEM_USE_MUMPS=YES *and* MFEM_USE_PETSC=YES (the on-Expanse
#       MUMPS is delivered through PETSc; a MUMPS=YES/PETSC=NO config.mk is the
#       stale committed file, not a real build_expanse.sh output), and
#   (2) the chosen config's [solver].type is a DIRECT solver, and specifically a
#       MUMPS variant (mumps / mumps_blr) — superlu/strumpack would MFEM_ABORT
#       at solver construction on this build.
# Override the config with BP5_CONFIG (must still be a MUMPS-variant direct
# config).  For a memory-hungry factorization, BP5_CONFIG can point at a
# mumps_blr (Block Low-Rank) variant to cut factorization memory/time.
#
# === FULL 1800-yr RUN ======================================================
# The verification config ships tfinal = 3.15e11 s (~10 kyr) as a large cap that
# the parity gate stops early with --max-steps 50.  This PRODUCTION job removes
# the step cap and OVERRIDES the terminal time with --tfinal so the run goes to
# the SCEC BP5 final time of ~1800 yr.  Default BP5_TFINAL = 56844000000 s — the
# exact value shipped in config/bp5_production_gmres.toml ("# ~1800 yr"), so the
# Expanse direct run terminates at the SAME t_final as the Frontera production
# GMRES runs (clean cross-solver comparison).  For reference, an EXACT 1800 yr is
# 56803680000 s using the code constant seconds_per_year = 365.25*24*3600 =
# 31557600 (config/bp5_params.hpp:172,175).  Override with BP5_TFINAL.
#
# === OPEN CHECKPOINTS (chained-segment restart) ============================
# 1800 yr will NOT finish in one 48 h walltime, so checkpointing is REQUIRED to
# resume.  --checkpoint-every CKPT writes a QD checkpoint every CKPT accepted
# steps (plus one guaranteed final checkpoint at clean exit) to
#   <output_dir>/spatial_seas_checkpoint_r<rank>.txt   (one ASCII file per rank).
# Each submission is ONE segment writing into BASE/segment_NNN.  Submit segment 1
# (cold start: runs the 4-phase init), let it run/timeout, then submit segment 2
# (SEGMENT=2 -> --restart BASE/segment_001/spatial_seas, which restores
# state+(t,dt,step) and SKIPS the 4-phase init), etc., until t_final is reached:
#
#   # 0) build the QD driver ONCE on Expanse (full build, then fast rebuilds):
#   #      USE_MUMPS=YES bash build_expanse.sh        # interactive shared session
#   #      sbatch miniapps/seas/bp5_spatial/expanse/bp5_spatial_seas_build_expanse.sbatch
#   # 1) cold-start segment 1:
#   sbatch miniapps/seas/bp5_spatial/expanse/bp5_spatial_seas_mumps_production_1800yr_expanse_2N100r_48hr.sbatch
#   # 2) after it finishes/timeouts, chain the next segment(s):
#   sbatch --export=ALL,SEGMENT=2 miniapps/seas/bp5_spatial/expanse/bp5_spatial_seas_mumps_production_1800yr_expanse_2N100r_48hr.sbatch
#   sbatch --export=ALL,SEGMENT=3 ...                # until "SEAS Simulation Complete"
#
# The restart aborts if the mesh/partition/rank-count changed or if the friction
# params drift > 1e-9 vs the checkpoint (config-drift guard) — keep the same
# config, node/rank layout, and mesh across segments.  A run that reaches t_final
# exits 0; a 48 h walltime kill is a SLURM TIMEOUT (resume via the next segment).
# Station traces are written per-segment under BASE/segment_NNN — stitch them in
# post-processing.
#
# *** RESUME CAVEATS (read before relying on the chain) ***
#  (1) Checkpoint cadence is STEP-based, and the driver has NO SIGTERM handler:
#      a 48 h TIMEOUT is a SIGKILL, so the guaranteed post-loop "final" checkpoint
#      (spatial_seas_driver.cpp:1330) does NOT fire on a kill.  A segment therefore
#      leaves a resumable checkpoint ONLY if it completed >= CKPT accepted steps
#      before the kill.  Pick CKPT comfortably below the worst-case step count a
#      segment can reach in 48 h.  This matters MOST for SEGMENT=1: a segment-1
#      timeout that did < CKPT steps leaves NOTHING to resume and the chain must
#      restart from scratch.  The default below (250) is deliberately small (the
#      QD checkpoint is a tiny per-rank ASCII file; the MFEM stiffness matrix is
#      factorized ONCE and reused, so accepted steps are cheap back-solves and
#      thousands typically fit in 48 h — frequent checkpoints cost ~nothing).
#      Lower CKPT further if a segment risks a slow/short window.
#  (2) Each segment keeps only ONE checkpoint set (fixed prefix, overwritten in
#      place; the write is not atomic — io/seas_qd_checkpoint.hpp).  A SIGKILL
#      landing mid-write can truncate that segment's checkpoint.  The fallback is
#      the PREVIOUS segment's (intact) checkpoint — so always keep the prior
#      segment_NNN dir until the next segment has written a fresh checkpoint.
#      SEGMENT=1 has no previous segment, so re-cold-start if its checkpoint is bad.
#
# Overrides: BASE (results base dir), SEGMENT (1..), CKPT (checkpoint interval in
# accepted steps, default 250), BP5_CONFIG (MUMPS-variant direct TOML),
# BP5_MESH (mesh override), BP5_TFINAL (terminal time, s).

export LC_ALL=C
export LANG=C

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

# Resolve seas-mfem root (walk UP to miniapps/seas/Makefile) BEFORE the lib paths.
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR}}"
_walk_up_dir="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5 6 7; do
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

# Runtime shared libs (hdf5/zfp from-source; petsc lib dir; caliper paths are
# harmless if absent on a non-Caliper production build).  PETSc/hypre/MUMPS/metis
# are STATIC from-source (linked INTO the binary), so no LD_LIBRARY_PATH ordering
# can shadow them — unlike Frontera there is no HYPRE-2.31 dance.
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"

cd "${SEAS_MFEM_ROOT}"

# --- Build pre-flight: this is a DIRECT (MUMPS) run.  Both MUMPS and its PETSc
# delivery vehicle MUST be in config.mk, else the driver MFEM_ABORTs at solver
# construction.  (The CG-AMG caliper job needs only PETSc; this one needs MUMPS.)
grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || {
    echo "ERROR: config.mk has no MFEM_USE_PETSC=YES."
    echo "  On Expanse, MUMPS is built from source through PETSc; rebuild with:"
    echo "      USE_MUMPS=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || {
    echo "ERROR: config.mk has no MFEM_USE_MUMPS=YES — the [solver].type=\"mumps\""
    echo "  DIRECT solve would MFEM_ABORT at construction.  Rebuild MUMPS:"
    echo "      USE_MUMPS=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    echo "  (build_expanse.sh delivers MUMPS via PETSc --download-mumps --download-scalapack.)"
    exit 1; }

cd miniapps/seas
if [[ ! -x ./seas_spatial_seas_driver ]]; then
    echo "ERROR: ./seas_spatial_seas_driver not found.  Build it FIRST on Expanse:"
    echo "  (full, with MUMPS)  USE_MUMPS=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    echo "  (fast incremental)  sbatch bp5_spatial/expanse/bp5_spatial_seas_build_expanse.sbatch"
    exit 1
fi

# --- Config: the MUMPS direct-solver BP5 production config by DEFAULT. ---------
CONFIG_TOML="${BP5_CONFIG:-config/safs_qd/bp5_spatial_seas_verification.toml}"
[[ -f "${CONFIG_TOML}" ]] || { echo "ERROR: config not found: $(pwd)/${CONFIG_TOML}"; exit 1; }

# Enforce the "direct solver" requirement: refuse an iterative config, and refuse
# a direct solver that is NOT compiled on Expanse (only MUMPS is built here).
SOLVER_TYPE="$(grep -E '^[[:space:]]*type[[:space:]]*=' "${CONFIG_TOML}" | head -1 | sed -E 's/.*=[[:space:]]*"?([a-zA-Z_]+)"?.*/\1/')"
case "${SOLVER_TYPE}" in
    mumps|mumps_blr)
        echo "Solver:      ${SOLVER_TYPE}  (direct, MUMPS variant — OK on Expanse)" ;;
    superlu|strumpack)
        echo "ERROR: ${CONFIG_TOML} uses [solver].type=\"${SOLVER_TYPE}\" (a direct solver"
        echo "  that is NOT compiled on Expanse — only MUMPS is built).  Use a MUMPS"
        echo "  config (config/safs_qd/bp5_spatial_seas_verification.toml) or pass a"
        echo "  mumps/mumps_blr BP5_CONFIG."
        exit 1 ;;
    cg_amg|gmres_amg|gmres_ilu)
        echo "ERROR: ${CONFIG_TOML} uses [solver].type=\"${SOLVER_TYPE}\" (an ITERATIVE solver)."
        echo "  This is the DIRECT-solver production job.  Use the MUMPS config"
        echo "  (config/safs_qd/bp5_spatial_seas_verification.toml) or pass a"
        echo "  mumps/mumps_blr BP5_CONFIG.  (For the iterative profile, use the"
        echo "  bp5_spatial_seas_cgamg_caliper_expanse job instead.)"
        exit 1 ;;
    *)
        echo "WARNING: could not determine [solver].type from ${CONFIG_TOML}"
        echo "         (parsed '${SOLVER_TYPE:-<none>}'); proceeding — verify it is mumps/mumps_blr." ;;
esac

# --- Mesh pre-flight: the mesh path lives in the config ([mesh].path); resolve it
# (BP5_MESH overrides) and confirm it exists (the .msh is gitignored/large). -----
CONFIG_MESH="$(grep -E '^[[:space:]]*path[[:space:]]*=' "${CONFIG_TOML}" | head -1 | sed -E 's/.*=[[:space:]]*"([^"]+)".*/\1/')"
MESH="${BP5_MESH:-${CONFIG_MESH}}"
MESH_ARG=""
[[ -n "${BP5_MESH:-}" ]] && MESH_ARG="--mesh ${BP5_MESH}"
if [[ -n "${MESH}" && ! -f "${MESH}" ]]; then
    echo "ERROR: mesh not found: $(pwd)/${MESH}"
    echo "  It is referenced by ${CONFIG_TOML} ([mesh].path); the .msh is gitignored/large."
    echo "  scp it to Expanse, or pass an existing one with BP5_MESH=<path>."
    exit 1
fi

# --- Run length + checkpoint/segment bookkeeping --------------------------------
# 1800 yr full run.  Default = the value shipped in config/bp5_production_gmres.toml
# ("# ~1800 yr"), so this matches the Frontera production GMRES terminal time.
TFINAL="${BP5_TFINAL:-56844000000}"
BASE="${BASE:-${SEAS_MFEM_ROOT}/miniapps/seas/bp5/spatial_seas_expanse/bp5_mumps_prod_1800yr}"
SEGMENT="${SEGMENT:-1}"
CKPT="${CKPT:-250}"                        # checkpoint every CKPT accepted steps
                                           # (small on purpose: a TIMEOUT kill
                                           # only resumes from the last STEP-based
                                           # checkpoint — see RESUME CAVEATS above)
seg() { printf "segment_%03d" "$1"; }

OUT="${BASE}/$(seg "${SEGMENT}")"
mkdir -p "${OUT}"

RESTART_ARG=""
if [[ "${SEGMENT}" -gt 1 ]]; then
    PREV="${BASE}/$(seg "$((SEGMENT-1))")"
    # The QD driver writes its checkpoint with prefix "<output_dir>/spatial_seas"
    # (spatial_seas_driver.cpp:1137); the per-rank files are
    # "<prefix>_checkpoint_r<rank>.txt" (io/seas_qd_checkpoint.hpp:43).
    if ! ls "${PREV}/spatial_seas"_checkpoint_r*.txt >/dev/null 2>&1; then
        echo "WARNING: no checkpoint found under ${PREV}/ — segment ${SEGMENT} may have nothing to restart from."
    fi
    RESTART_ARG="--restart ${PREV}/spatial_seas"
    echo "Restarting segment ${SEGMENT} from ${PREV}/spatial_seas"
else
    echo "Initial segment ${SEGMENT} (cold start: runs the 4-phase init)."
fi

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Queue/size:  compute, 2N x 100r = ${SLURM_NTASKS} ranks, 48:00:00  (-A ddp408)"
echo "EXPERIMENT:  BP5 spatial QUASI-DYNAMIC, MUMPS DIRECT solver, FULL ~1800 yr"
echo "Config:      ${CONFIG_TOML}   (solver=${SOLVER_TYPE})"
echo "Mesh:        ${MESH:-<from config>} ${MESH_ARG:+(CLI override)}"
echo "t_final:     ${TFINAL} s  (~$(awk "BEGIN{printf \"%.1f\", ${TFINAL}/31557600}") yr)"
echo "Checkpoint:  every ${CKPT} accepted steps -> ${OUT}/spatial_seas_checkpoint_r<rank>.txt"
echo "BASE=${BASE}  SEGMENT=${SEGMENT}  OUT=${OUT}"

# Full run to t_final (NO --max-steps).  --checkpoint-every keeps periodic
# checkpoints so a 48 h TIMEOUT can resume from the next segment.  ParaView is
# forced off defensively (the QD driver does not emit volume/fault PV in its time
# loop today; only checkpoints + BP5 station traces are written).
srun --mpi=pmi2 -n "${SLURM_NTASKS}" stdbuf -oL -eL ./seas_spatial_seas_driver \
    --config "${CONFIG_TOML}" \
    ${MESH_ARG} \
    --tfinal "${TFINAL}" \
    --checkpoint-every "${CKPT}" \
    --paraview-volume off \
    --paraview-fault  off \
    --output-dir "${OUT}" \
    ${RESTART_ARG}
RC=$?

echo
echo "=== outputs (segment ${SEGMENT}) ==="
echo "checkpoints  -> ${OUT}/spatial_seas_checkpoint_r*.txt   (resume: SEGMENT=$((SEGMENT+1)))"
echo "BP5 stations -> ${OUT}/spatial_seas_*                    (10 SCEC traces + global V_max)"
if [[ ${RC} -eq 0 ]]; then
    echo "driver_rc=0 — segment ${SEGMENT} reached t_final (SEAS Simulation Complete) OR clean exit."
    echo "  If the log does NOT say 'Complete', it stopped early — submit the next segment."
else
    echo "driver_rc=${RC} — segment ${SEGMENT} ended non-zero (TIMEOUT or error); see the .out/.err."
    echo "  On a 48 h TIMEOUT, resume:  sbatch --export=ALL,SEGMENT=$((SEGMENT+1)) <this file>"
fi
exit ${RC}
