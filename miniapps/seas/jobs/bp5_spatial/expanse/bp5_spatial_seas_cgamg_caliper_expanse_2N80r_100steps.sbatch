#!/bin/bash
#SBATCH -J bp5seas_cg_cali
#SBATCH -o bp5seas_cg_cali_%j.out
#SBATCH -e bp5seas_cg_cali_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=80
#SBATCH --time=02:00:00
#SBATCH --mem=100000M
#SBATCH -A ddp408

# =========================================================================
# BP5 SPATIAL quasi-dynamic SEAS — CALIPER PROFILING run (CG-AMG iterative solver).
# *** EXPANSE (SDSC) job — 2N x 80r = 160 ranks / FIRST 100 STEPS ONLY. ***
# Driver:  seas_spatial_seas_driver   (the generic QD SEAS driver, --config flag)
# Config:  config/safs_qd/bp5_spatial_seas_verification_cgamg.toml  ([solver].type = "cg_amg")
# =========================================================================
# A SHORT run (--max-steps 100, NOT to t_final) whose ONLY purpose is a
# PER-COMPONENT + WHOLE-PROGRAM SPEED BREAKDOWN of the BP5-spatial QD solver via
# Caliper.  Same physics as the production parity gate (same config / mesh /
# CG-AMG) — just capped at 100 accepted steps so the run EXITS CLEANLY (Caliper
# only flushes its report at clean program exit; a walltime KILL writes nothing).
# This is the bp5_spatial sibling of the SAFS dynamic Caliper job
#   jobs/safs/safs_expanse/spatial_dyn_ratestate_..._caliper_..._expanse_2N100r_compute_4hr_safs.sbatch
# adapted for the QUASI-DYNAMIC (seas_spatial_seas_driver) BP5 problem on 160 ranks.
#
# === ITERATIVE SOLVER ======================================================
# "Use the iterative solver for BP5" -> [solver].type = "cg_amg" (Conjugate
# Gradient + BoomerAMG) in the config above.  The driver also supports
# "gmres_amg" / "gmres_ilu" (iterative) and "mumps"/"mumps_blr"/"superlu"/
# "strumpack" (direct).  This job HARD-CHECKS that the chosen config's solver is
# iterative (refuses a direct solver) so the profile is of the scalable KSP path.
# Override the config with BP5_CONFIG (must still be an iterative-solver config).
#
# === HOW THE PROFILE IS PRODUCED ==========================================
# Caliper is driven ENTIRELY by the CALI_CONFIG env var (read at init) — there is
# NO driver flag.  The QD driver has FEW instrumented MFEM_PERF_SCOPE regions
# (only the elasticity SetupSolver::CG_AMG / GMRES_AMG region fires, once/rank at
# construction — there is NO per-time-step region in the QD path), so unlike the
# dynamic driver the runtime-report TABLE is SPARSE.  Therefore the WHOLE-PROGRAM
# SAMPLING channel (sample-report) is the PRIMARY signal here: it samples every
# function — init, the 4-phase equilibrium solve, per-step CG-AMG iterations,
# HYPRE BoomerAMG setup/solve, friction Brent solves, MPI, I/O — instrumented or
# not.  (Same lesson as the BP5 Frontera Caliper job.)
#
# *** REQUIRES A CALIPER-ENABLED BINARY *** (MFEM_USE_CALIPER=YES).  On Expanse:
#     USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh        # full build
#   then (after later source edits) the fast incremental rebuild:
#     sbatch miniapps/seas/jobs/bp5_spatial/expanse/bp5_spatial_seas_build_expanse.sbatch
# build_expanse.sh ships Caliper v2.11.0 (WITH_SAMPLER always; WITH_LIBDW when
# elfutils is found -> sample-report works).  USE_CALIPER defaults to NO, so the
# NORMAL binary has the regions compiled out (empty report); this job hard-checks
# the libcaliper link before running (override with BP5_SKIP_CALI_CHECK=1).
#
# Overrides: BP5_CONFIG (iterative-solver TOML), BP5_MESH, BP5_MAXSTEPS (default 100),
# BP5_CALI_CONFIG (full CALI_CONFIG override), BP5_CALI_HZ (sampling rate/rank,
# default 200), BP5_CALI_SAMPLE (1=add whole-program sample-report, default 1;
# needs a WITH_LIBDW Caliper build), BP5_CALI_CALLPATH (1=use callpath-sample-report
# for sample CALL-TREES instead of flat, default 0; needs a WITH_LIBUNWIND build),
# BP5_PV_VOLUME / BP5_PV_FAULT (ParaView modes, default "off" to keep the profile
# clean), BP5_SKIP_CALI_CHECK (1=run even if the binary lacks Caliper), BP5_OUTDIR.

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

# Runtime shared libs.  libcaliper.so FIRST (this job profiles); then the usual
# hdf5/zfp/petsc from build_expanse.sh's from-source deps.
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"

cd "${SEAS_MFEM_ROOT}"
grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_PETSC=YES"; exit 1; }

cd miniapps/seas
if [[ ! -x ./seas_spatial_seas_driver ]]; then
    echo "ERROR: ./seas_spatial_seas_driver not found.  Build it FIRST on Expanse:"
    echo "  (full, with Caliper)  USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    echo "  (fast incremental)    sbatch jobs/bp5_spatial/expanse/bp5_spatial_seas_build_expanse.sbatch"
    exit 1
fi

# --- Caliper pre-flight: the binary MUST be linked against libcaliper, else the
# MFEM_PERF_SCOPE regions are compiled out and the report is EMPTY. -------------
if ! ldd ./seas_spatial_seas_driver 2>/dev/null | grep -qi 'libcaliper'; then
    echo "ERROR: ./seas_spatial_seas_driver is NOT linked against libcaliper"
    echo "       (MFEM_USE_CALIPER=NO -> the MFEM_PERF_SCOPE regions are no-ops and"
    echo "        the whole-program sampler still runs, but the region channel is empty)."
    echo "  Rebuild WITH Caliper:  USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    echo "  To run anyway (e.g. just to time the bulk):"
    echo "        sbatch --export=ALL,BP5_SKIP_CALI_CHECK=1 <this file>"
    [[ "${BP5_SKIP_CALI_CHECK:-0}" == "1" ]] || exit 1
    echo "  BP5_SKIP_CALI_CHECK=1 set — continuing without Caliper."
fi

# --- Config: the CG-AMG iterative-solver BP5 parity config by DEFAULT. ---------
CONFIG_TOML="${BP5_CONFIG:-config/safs_qd/bp5_spatial_seas_verification_cgamg.toml}"
[[ -f "${CONFIG_TOML}" ]] || { echo "ERROR: config not found: $(pwd)/${CONFIG_TOML}"; exit 1; }

# Enforce the "iterative solver" requirement: refuse a direct-solver config.  The
# [solver].type must be one of the iterative kinds (cg_amg / gmres_amg / gmres_ilu).
SOLVER_TYPE="$(grep -E '^[[:space:]]*type[[:space:]]*=' "${CONFIG_TOML}" | head -1 | sed -E 's/.*=[[:space:]]*"?([a-zA-Z_]+)"?.*/\1/')"
case "${SOLVER_TYPE}" in
    cg_amg|gmres_amg|gmres_ilu)
        echo "Solver:      ${SOLVER_TYPE}  (iterative — OK)" ;;
    mumps|mumps_blr|superlu|strumpack)
        echo "ERROR: ${CONFIG_TOML} uses [solver].type=\"${SOLVER_TYPE}\" (a DIRECT solver)."
        echo "  This job profiles the ITERATIVE BP5 solve.  Use the cg_amg config"
        echo "  (config/safs_qd/bp5_spatial_seas_verification_cgamg.toml) or pass an"
        echo "  iterative-solver BP5_CONFIG."
        exit 1 ;;
    *)
        echo "WARNING: could not determine [solver].type from ${CONFIG_TOML}"
        echo "         (parsed '${SOLVER_TYPE:-<none>}'); proceeding — verify it is iterative." ;;
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

MAXSTEPS="${BP5_MAXSTEPS:-100}"      # FIRST 100 STEPS ONLY (clean exit -> Caliper flushes)
PV_VOLUME="${BP5_PV_VOLUME:-off}"    # ParaView off by default: keep the profile clean
PV_FAULT="${BP5_PV_FAULT:-off}"      #   (BP5 station traces still written via bp5_stations=true)

OUT_BASE="${BP5_OUTDIR:-${SEAS_MFEM_ROOT}/miniapps/seas/bp5/spatial_seas_expanse}"
OUT="${OUT_BASE}/bp5_seas_cgamg_cali_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/bp5_spatial/expanse"
LOG="${LOG_DIR}/bp5_spatial_seas_cgamg_cali_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

# === Caliper config — 3 complementary channels =============================
CALI_REGION_TXT="${LOG%.log}.cali-region-report.txt"   # instrumented region table (SPARSE for QD)
CALI_SPOT="${LOG%.log}.cali"                           # Spot .cali (Hatchet/Thicket, per-rank)
CALI_SAMPLE_TXT="${LOG%.log}.cali-sample-report.txt"   # SAMPLING: WHICH FUNCTION IS SLOWEST (PRIMARY)
#
# Channel 1  runtime-report : time of each INSTRUMENTED region (incl+excl, call
#                             counts, per-region MPI; Min/Avg/Max across ranks —
#                             Max>>Avg flags LOAD IMBALANCE).  Sparse here: only
#                             the elasticity SetupSolver::CG_AMG region fires.
# Channel 2  spot           : .cali for the LLNL Spot UI + Python Hatchet/Thicket
#                             (per-rank, comparable, deep analysis).
# Channel 3  sample-report  : timer SAMPLING of the WHOLE program — every function
#                             (init, 4-phase solve, per-step CG-AMG/BoomerAMG,
#                             friction Brent, MPI, I/O).  THE primary signal for
#                             the QD driver.  Needs Caliper's symbollookup service
#                             (WITH_LIBDW); ConfigManager aborts the ENTIRE
#                             CALI_CONFIG if a channel's service is missing, so it
#                             is a toggle (default ON — the Expanse build is
#                             WITH_LIBDW; set BP5_CALI_SAMPLE=0 if yours is not).
export CALI_SAMPLER_FREQUENCY="${BP5_CALI_HZ:-200}"    # samples/s/rank
DEFAULT_CALI_CONFIG="runtime-report(calc.inclusive,region.count,profile.mpi,max_column_width=60,output=${CALI_REGION_TXT}),spot(output=${CALI_SPOT})"
if [[ "${BP5_CALI_SAMPLE:-1}" == "1" ]]; then
    if [[ "${BP5_CALI_CALLPATH:-0}" == "1" ]]; then
        # callpath-sample-report: sample CALL-TREES (needs WITH_LIBUNWIND too).
        DEFAULT_CALI_CONFIG="${DEFAULT_CALI_CONFIG},callpath-sample-report(output=${CALI_SAMPLE_TXT})"
    else
        # flat sample-report: needs only symbollookup (WITH_LIBDW).
        DEFAULT_CALI_CONFIG="${DEFAULT_CALI_CONFIG},sample-report(output=${CALI_SAMPLE_TXT})"
    fi
fi
export CALI_CONFIG="${BP5_CALI_CONFIG:-${DEFAULT_CALI_CONFIG}}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Queue/size:  compute, 2N x 80r = ${SLURM_NTASKS} ranks, 02:00:00  (CALIPER PROFILE)"
echo "EXPERIMENT:  BP5 spatial QUASI-DYNAMIC, CG-AMG iterative solver, first ${MAXSTEPS} steps"
echo "Config:      ${CONFIG_TOML}   (solver=${SOLVER_TYPE})"
echo "Mesh:        ${MESH:-<from config>} ${MESH_ARG:+(CLI override)}"
echo "ParaView:    volume=${PV_VOLUME}  fault=${PV_FAULT}   (BP5 station traces always written)"
echo "Caliper:     ldd $(ldd ./seas_spatial_seas_driver 2>/dev/null | grep -m1 -i libcaliper | awk '{print $1" -> "$3}')"
echo "CALI_CONFIG: ${CALI_CONFIG}"
echo "Output dir:  ${OUT}"

# --max-steps ${MAXSTEPS}: stop after 100 accepted steps -> CLEAN exit so Caliper
# flushes.  ParaView off keeps I/O out of the profile; bp5_stations=true (in the
# config) still writes the SCEC station traces into --output-dir.
srun --mpi=pmi2 -n "${SLURM_NTASKS}" stdbuf -oL -eL ./seas_spatial_seas_driver \
    --config "${CONFIG_TOML}" \
    ${MESH_ARG} \
    --max-steps "${MAXSTEPS}" \
    --paraview-volume "${PV_VOLUME}" \
    --paraview-fault  "${PV_FAULT}" \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

echo
echo "=== Caliper outputs ==="
if [[ -f "${CALI_REGION_TXT}" ]]; then
    echo "region report -> ${CALI_REGION_TXT}   (instrumented regions; SPARSE for QD — SetupSolver only)"
    echo "--- instrumented regions (slowest first) ---"
    head -40 "${CALI_REGION_TXT}"
else
    echo "NOTE: ${CALI_REGION_TXT} not written — region channel failed or no clean exit; CALI messages:"
    grep -iE 'cali|caliper' "${LOG}" | head
fi
echo
if [[ -f "${CALI_SAMPLE_TXT}" ]]; then
    echo "sample report -> ${CALI_SAMPLE_TXT}   (WHICH FUNCTION IS SLOWEST — whole program; PRIMARY)"
    echo "--- top sampled functions (self time) ---"
    head -40 "${CALI_SAMPLE_TXT}"
else
    echo "NOTE: ${CALI_SAMPLE_TXT} not written — the sampler/symbollookup service is"
    echo "      missing from this Caliper build (set BP5_CALI_SAMPLE=0 to drop it, or"
    echo "      rebuild Caliper WITH_LIBDW), or the run did not exit cleanly."
    grep -iE 'sampl|callpath|symbollookup|cali' "${LOG}" | head
fi
echo
[[ -f "${CALI_SPOT}" ]] && echo "spot/.cali    -> ${CALI_SPOT}   (Hatchet/Thicket or Spot UI: PER-RANK detail)"

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG}."; exit ${RC}
fi
echo "driver_rc=${RC} — 100-step profiling run done."
exit ${RC}
