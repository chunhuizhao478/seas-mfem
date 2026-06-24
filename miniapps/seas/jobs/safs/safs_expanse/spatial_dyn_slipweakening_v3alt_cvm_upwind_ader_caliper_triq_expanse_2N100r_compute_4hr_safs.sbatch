#!/bin/bash
#SBATCH -J safs_sw_alt_uw_cali
#SBATCH -o safs_sw_alt_uw_cali_%j.out
#SBATCH -e safs_sw_alt_uw_cali_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=50
#SBATCH --time=04:00:00
#SBATCH --mem=100000M
#SBATCH -A ddp408

# =========================================================================
# SAFS SLIPWEAKENING v3.0.0 ALT — CALIPER PROFILING run (UPWIND Godunov + ADER).
# *** EXPANSE (SDSC) job — 2N x 50r = 100 ranks / 4 h / 100 GB / tfinal=1 s. ***
# Config: ${SAFS_ROOT}/config/spatial_friction_safs_seisol_v3_0_0_LSW_ALT.toml
# Plan:   safs/project_7.0_alternative/document/PLAN_mixed_flux_seissol_port_2026-06-19.md
# =========================================================================
# A SHORT run (tfinal=1 s, 4 h wall) whose ONLY purpose is a PER-COMPONENT SPEED
# BREAKDOWN of the solver via Caliper.  Identical physics to the production
# job (same config / mesh / sidecars / UPWIND Godunov + ADER) — just short + profiled.
#
# === HOW THE PROFILE IS PRODUCED ==========================================
# The solver is already instrumented with MFEM_PERF_SCOPE / MFEM_PERF_BEGIN
# regions (general/annotation.hpp -> Caliper CALI_CXX_MARK_SCOPE when
# MFEM_USE_CALIPER=YES, else a NO-OP).  Per-component regions include:
#   seas::spatial_dyn::step                 (one whole time step)
#   seas::WaveOperator::Mult                (RHS eval; called per RK stage)
#     -> seas::WaveOperator::ComputeVolumeRHS         (volume/interior derivative)
#     -> seas::WaveOperator::ComputeFaceFluxRHS       (interior face flux)
#     -> seas::WaveOperator::ComputeSharedFaceFluxRHS (rank-seam face flux)
#   seas::spatial_dyn::friction_substep     (rate-state / LSW friction solve)
#   seas::WaveOperator::ComputeMaxDt        (CFL dt reduction)
#   (ADER IS the hot integrator here; the RK45 regions are unused)
# Caliper is driven ENTIRELY by the CALI_CONFIG env var (read at init) — there is
# NO driver flag.  This job sets a detailed CALI_CONFIG below.
#
# *** REQUIRES A CALIPER-ENABLED BINARY *** (MFEM_USE_CALIPER=YES).  On Expanse:
#     USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh
# build_expanse.sh ships Caliper v2.11.0 but USE_CALIPER defaults to NO, so the
# normal binary has the regions COMPILED OUT (empty report).  This job hard-checks
# the libcaliper link before running (override with SAFS_SKIP_CALI_CHECK=1).
#
# === CALIPER CONFIG RECIPES (what Caliper can offer) =======================
# CALI_CONFIG takes a comma-separated list of built-in "config recipes"; each
# becomes an output channel.  The most useful for this solver:
#   runtime-report            per-region time TABLE (this run's default output).
#                             options (key[=val], comma-separated inside the (...)):
#                               calc.inclusive     also show inclusive (subtree) time
#                               region.count       show per-region call counts
#                               profile.mpi        break out MPI time per region
#                               mem.highwatermark  per-region peak resident memory
#                               output=<file>      where to write (else stdout)
#                               max_column_width=N widen the region-name column
#   hatchet-region-profile    .json for Python Hatchet (scriptable region tree).
#   spot                      .cali for the LLNL Spot web UI + Hatchet (richest).
#   mpi-report                standalone MPI-function breakdown (calls + time).
#   sample-report    SAMPLING profile (low overhead, no instrumentation).
#   loop-report               per-iteration timing of marked loops.
# Discover the full list + options at runtime:  CALI_CONFIG=help ./<any cali app>
# (or see Caliper's "Builtin configurations" docs).  Tuning knobs:
#   CALI_LOG_VERBOSITY=1   CALI_CONFIG_PROFILE   CALI_AGGREGATE_KEY=...
#
# Overrides: SAFS_CONFIG, SAFS_MESH, SAFS_MIXED_FLUX (default none),
# SAFS_TIME_INTEGRATOR (default ader), SAFS_TFINAL (default 1s — SHORT profile),
# SAFS_CFL, SAFS_DIAG (0/1, default 0), SAFS_CALI_CONFIG (full CALI_CONFIG override),
# SAFS_CALI_HZ (sampling rate/rank, default 500), SAFS_CALI_SAMPLE (1 = add the
# whole-program sample-report channel; needs a WITH_LIBDW Caliper build),
# SAFS_CALI_PAPI (1 = HW counters,
# needs a WITH_PAPI build) + SAFS_CALI_PAPI_COUNTERS, SAFS_SKIP_CALI_CHECK
# (1 = run even if the binary lacks Caliper), SAFS_OUTDIR.

export LC_ALL=C
export LANG=C
# Per-step [DIAG]/[DIAG-ONSET] monitor: OFF by default (it floods once V_max>10).
# Leave it OFF for profiling — the DIAG prints would distort the I/O timing.
if [[ "${SAFS_DIAG:-0}" == "1" ]]; then export SEAS_DIAG_BLOWUP=1; else unset SEAS_DIAG_BLOWUP; fi
# R-101 shared-fault tripwire SKIPPED by default for SAFS; SAFS_R101=1 re-enables.
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

# Resolve seas-mfem root (walk UP to miniapps/seas/Makefile) BEFORE the lib paths.
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

# Runtime shared libs.  libcaliper.so FIRST (this job profiles); then the usual
# hdf5/zfp/petsc from build_expanse.sh's from-source deps.
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"

cd "${SEAS_MFEM_ROOT}"
grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not found.  Build it FIRST on Expanse"
    echo "  (with Caliper, for this job):  USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    exit 1
fi

# --- Caliper pre-flight: the binary MUST be linked against libcaliper, else the
# MFEM_PERF_SCOPE regions are compiled out and the report is EMPTY. -------------
if ! ldd ./seas_spatial_dyn_driver 2>/dev/null | grep -qi 'libcaliper'; then
    echo "ERROR: ./seas_spatial_dyn_driver is NOT linked against libcaliper"
    echo "       (MFEM_USE_CALIPER=NO -> the MFEM_PERF_SCOPE regions are no-ops,"
    echo "        so Caliper would produce an EMPTY breakdown)."
    echo "  Rebuild WITH Caliper:  USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    echo "  (build_expanse.sh ships Caliper v2.11.0; USE_CALIPER just defaults to NO.)"
    echo "  To run anyway (e.g. to time the bulk without per-region data):"
    echo "        sbatch --export=ALL,SAFS_SKIP_CALI_CHECK=1 <this file>"
    [[ "${SAFS_SKIP_CALI_CHECK:-0}" == "1" ]] || exit 1
    echo "  SAFS_SKIP_CALI_CHECK=1 set — continuing without Caliper (empty report)."
fi

SAFS_ROOT=safs/project_7.0_alternative
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_safs_seisol_v3_0_0_LSW_ALT.toml}"

# Mesh: DEFAULT = the sliver-repaired (triq) 500 m mesh (flat z=0 top); override with SAFS_MESH.
DEFAULT_MESH="${SAFS_ROOT}/experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh"
MESH="${SAFS_MESH:-${DEFAULT_MESH}}"
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: mesh not found: $(pwd)/${MESH}  (gitignored — scp it first)"; exit 1
fi

# PRODUCTION SIDECAR PRE-FLIGHT (CVM material + CSM stress; both .h5 required).
VEL_H5="${SAFS_ROOT}/velocity/results/multiscale_statewise_cvm/velocity_safs.h5"
STR_H5="${SAFS_ROOT}/stress/results/csm_d12_mfem/safs_stress_csm_d12_mfem.h5"
for _h5 in "${VEL_H5}" "${STR_H5}"; do
    if [[ ! -f "${_h5}" ]]; then
        echo "ERROR: required sidecar missing: $(pwd)/${_h5}"
        echo "       It is gitignored/large; scp it to Expanse before submitting."
        exit 1
    fi
done

# Floor pre-flight (v3 slipweakening keeps the 5 MPa MFEM-only sigma_n floor).
grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null     || { echo "ERROR: ${CONFIG_TOML} has no sigma_n_strength_floor_pa"; exit 1; }

# Flux: UPWIND/none (Godunov) by DEFAULT — central-flux-free.
SAFS_MIXED_FLUX="${SAFS_MIXED_FLUX:-none}"
case "${SAFS_MIXED_FLUX}" in
    none|adjacent|all_continuous) : ;;
    *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid"; exit 1 ;;
esac
MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"

# Time integrator: ADER by DEFAULT (Godunov upwind is dissipative + stable under ADER).
SAFS_TIME_INTEGRATOR="${SAFS_TIME_INTEGRATOR:-ader}"
case "${SAFS_TIME_INTEGRATOR}" in
    ader|rk4|rk45) : ;;
    *) echo "ERROR: SAFS_TIME_INTEGRATOR='${SAFS_TIME_INTEGRATOR}' invalid (ader|rk4|rk45)"; exit 1 ;;
esac
if [[ "${SAFS_MIXED_FLUX}" != "none" && "${SAFS_TIME_INTEGRATOR}" == "ader" ]]; then
    echo "ERROR: mixed_flux=${SAFS_MIXED_FLUX} + ADER is UNSTABLE (driver aborts). Use rk4/rk45."; exit 1
fi
TI_ARG="--time-integrator ${SAFS_TIME_INTEGRATOR}"

# === ADER hot-path optimization levers (opt 2026-06-24) =====================
# CALIPER MEASUREMENT run: levers default ON to profile the OPTIMIZED path
# (document/code_optimization_dev/action_plan_2026-06-24.md, Phase 0).  Set
# SAFS_SHARED_CK=0 / SAFS_DERIV_CACHE=0 to A/B against the un-optimized baseline.
#
#   --shared-ck-recursion : ONE Cauchy-Kovalewski recursion/step instead of two
#                           (BIT-EXACT; removes the 2nd ApplySpatialDerivative
#                           subtree).  ADER-ONLY — a no-op under rk/rk45, so it
#                           is omitted automatically when integrator != ader.
#   --deriv-cache         : precomputed D_d^e / S_d^e dense mat-vec instead of
#                           per-QP quadrature for ApplySpatialDerivative (ADER)
#                           and ComputeVolumeRHS (ADER + RK) (<=1e-12, NOT
#                           bit-exact).  Matrix path gated by
#                           `make test-bimaterial-deriv-cache-parity` (also under
#                           `make test-parallel`); PRODUCTION SCIENCE runs must
#                           ALSO pass a SAFS output fingerprint check on Expanse.
OPT_ARGS=""
if [[ "${SAFS_TIME_INTEGRATOR}" == "ader" && "${SAFS_SHARED_CK:-1}" == "1" ]]; then
    OPT_ARGS="${OPT_ARGS} --shared-ck-recursion"
fi
[[ "${SAFS_DERIV_CACHE:-1}" == "1" ]] && OPT_ARGS="${OPT_ARGS} --deriv-cache"

# SHORT profiling horizon: tfinal=1 s by default (just enough to cover
# init -> nucleation -> early breakout so every component is exercised).
TFINAL="${SAFS_TFINAL:-1s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

OUT_BASE="${SAFS_OUTDIR:-${SEAS_MFEM_ROOT}/miniapps/seas/safs}"
OUT="${OUT_BASE}/safs_dyn_sw_alt_uw_cali_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs/safs_expanse"
LOG="${LOG_DIR}/spatial_dyn_sw_alt_uw_cali_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

# === Caliper config — MAXIMAL detail (3 complementary channels) =============
CALI_REGION_TXT="${LOG%.log}.cali-region-report.txt"   # instrumented region table
CALI_SPOT="${LOG%.log}.cali"                           # Spot .cali (Hatchet/Thicket, per-rank)
CALI_SAMPLE_TXT="${LOG%.log}.cali-sample-report.txt"   # SAMPLING: every function, all workflow
#
# Channel 1  runtime-report          : time of each INSTRUMENTED region (inclusive
#                                       + exclusive), call counts, per-region MPI
#                                       time, aggregated across ranks (Min/Avg/Max
#                                       — Max>>Avg flags LOAD IMBALANCE, expected
#                                       here since the fault is on a subset of ranks).
# Channel 2  spot                    : .cali for the LLNL Spot UI + Python
#                                       Hatchet/Thicket — the PER-RANK, comparable,
#                                       deep-analysis artifact (which ranks are slow).
# Channel 3  sample-report  : timer SAMPLING of the WHOLE program — every
#                                       function, instrumented or NOT (init,
#                                       PETSc/MUMPS, HDF5/VTU I/O, MPI internals,
#                                       partitioning).  THIS answers "which function
#                                       is the slowest" across the FULL workflow; the
#                                       region channels only see the marked regions.
#                                       Sampler is built-in; call-trees need the
#                                       Caliper build's callpath/libunwind (degrades
#                                       to a flat warning if absent — see sample-report).
export CALI_SAMPLER_FREQUENCY="${SAFS_CALI_HZ:-500}"   # samples/s/rank (higher = finer, more overhead)
# Channel 3 (sample-report) needs Caliper's symbollookup service (sample PC ->
# function NAME) = a Caliper built WITH_LIBDW.  ConfigManager aborts the ENTIRE
# CALI_CONFIG if any channel's required service is missing (taking runtime-report
# + spot DOWN with it), so the sampler channel is OFF by default; set
# SAFS_CALI_SAMPLE=1 once you have a WITH_LIBDW Caliper build (build_expanse.sh).
DEFAULT_CALI_CONFIG="runtime-report(calc.inclusive,region.count,profile.mpi,max_column_width=60,output=${CALI_REGION_TXT}),spot(output=${CALI_SPOT})"
if [[ "${SAFS_CALI_SAMPLE:-0}" == "1" ]]; then
    DEFAULT_CALI_CONFIG="${DEFAULT_CALI_CONFIG},sample-report(output=${CALI_SAMPLE_TXT})"
fi
export CALI_CONFIG="${SAFS_CALI_CONFIG:-${DEFAULT_CALI_CONFIG}}"
#
# The "WHY a function is slow" layer — hardware counters (memory-bound vs
# compute-bound).  SAFS_CALI_PAPI=1 turns on PAPI: IPC (PAPI_TOT_INS/PAPI_TOT_CYC)
# + cache misses (L2/L3) attributed per region, written to a separate Hatchet
# .json.  REQUIRES a Caliper build WITH_PAPI (USE_CALIPER=YES + the papi module;
# build_expanse.sh now enables it when papi is found).  Expanse is AMD EPYC, so
# these are RAW PAPI counters, NOT Intel 'topdown'; tune via SAFS_CALI_PAPI_COUNTERS
# (run 'papi_avail' on a compute node for the supported set).
if [[ "${SAFS_CALI_PAPI:-0}" == "1" ]]; then
    export CALI_PAPI_COUNTERS="${SAFS_CALI_PAPI_COUNTERS:-PAPI_TOT_CYC,PAPI_TOT_INS,PAPI_L2_DCM,PAPI_L3_TCM}"
    export CALI_CONFIG="${CALI_CONFIG},hatchet-region-profile(output=${LOG%.log}.cali-papi-profile.json,use.mpi)"
    echo "Caliper PAPI: ${CALI_PAPI_COUNTERS}  (needs a WITH_PAPI Caliper build)"
fi

echo "Git commit:  $(git rev-parse HEAD)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Queue/size:  compute, 2N x 50r = ${SLURM_NTASKS} ranks, 100 GB/node, 04:00:00  (CALIPER PROFILE)"
echo "EXPERIMENT:  slipweakening v3 ALT, CVM + CSM, UPWIND (Godunov) + ADER — Caliper per-component profile"
echo "Config:      ${CONFIG_TOML}"
echo "Mesh:        ${MESH}"
echo "Mixed flux:  ${SAFS_MIXED_FLUX}    Integrator: ${SAFS_TIME_INTEGRATOR}"
echo "Opt levers:  ${OPT_ARGS:-<none>}  (SAFS_SHARED_CK=${SAFS_SHARED_CK:-1} SAFS_DERIV_CACHE=${SAFS_DERIV_CACHE:-1}; shared-ck is ADER-only)"
echo "tfinal:      ${TFINAL}  (SHORT profile)   CFL: ${CFL_ARG:-<config default 0.5>}"
echo "Caliper:     ldd $(ldd ./seas_spatial_dyn_driver 2>/dev/null | grep -m1 -i libcaliper | awk '{print $1" -> "$3}')"
echo "CALI_CONFIG: ${CALI_CONFIG}"
echo "Output dir:  ${OUT}"

# NOTE: --paraview-max-snapshots 5 keeps I/O tiny so it does not dominate the
# 1 s profile (production uses 300).  Checkpointing is disabled (-every 0) for
# the same reason — this is a throwaway profiling run.
srun --mpi=pmi2 -n "${SLURM_NTASKS}" stdbuf -oL -eL ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    ${MIXFLUX_ARG} \
    ${TI_ARG} \
    ${OPT_ARGS} \
    --tfinal "${TFINAL}" \
    ${CFL_ARG} \
    --paraview \
    --paraview-fault-vtu \
    --paraview-free-surface vtu \
    --paraview-max-snapshots 5 \
    --print-derived \
    --checkpoint-every 0 \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

echo
echo "=== Caliper outputs ==="
if [[ -f "${CALI_REGION_TXT}" ]]; then
    echo "region report -> ${CALI_REGION_TXT}   (instrumented components; Max>>Avg = imbalance)"
    echo "--- instrumented regions (slowest first) ---"
    head -40 "${CALI_REGION_TXT}"
else
    echo "WARNING: ${CALI_REGION_TXT} not written — region channel failed; CALI messages:"
    grep -iE 'cali|caliper' "${LOG}" | head
fi
echo
if [[ -f "${CALI_SAMPLE_TXT}" ]]; then
    echo "sample report -> ${CALI_SAMPLE_TXT}   (WHICH FUNCTION IS SLOWEST — whole program)"
    echo "--- top sampled functions (self time) ---"
    head -40 "${CALI_SAMPLE_TXT}"
else
    echo "NOTE: ${CALI_SAMPLE_TXT} not written — the sampler/callpath service is missing"
    echo "      from this Caliper build (sampler is built-in; call-trees need libunwind)."
    grep -iE 'sampl|callpath|cali' "${LOG}" | head
fi
echo
[[ -f "${CALI_SPOT}" ]] && echo "spot/.cali    -> ${CALI_SPOT}   (Hatchet/Thicket or Spot UI: PER-RANK detail)"
[[ -f "${LOG%.log}.cali-papi-profile.json" ]] && echo "papi profile  -> ${LOG%.log}.cali-papi-profile.json   (hardware counters: the 'why')"

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG}."; exit ${RC}
fi
echo "driver_rc=${RC} — profiling run done."
exit ${RC}
