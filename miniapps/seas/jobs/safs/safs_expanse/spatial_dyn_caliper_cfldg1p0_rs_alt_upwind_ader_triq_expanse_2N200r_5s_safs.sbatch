#!/bin/bash
#SBATCH -J safs_cali_cfldg1p0
#SBATCH -o safs_cali_cfldg1p0_%j.out
#SBATCH -e safs_cali_cfldg1p0_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=100
#SBATCH --time=02:00:00
#SBATCH --mem=200000M
#SBATCH -A lbl107

# =========================================================================
# SAFS RATESTATE v3.0.0 ALT (vw005) — CALIPER PROFILE of the SeisSol-CFL arm.
# *** EXPANSE (SDSC) — 2N x 100r = 200 ranks / 2 h / 200 GB / tfinal=5 s. ***
#
# PURPOSE: a SHORT (tfinal=5 s) per-component speed breakdown of the run WITH
#   Lever A engaged at cfl_dg_safety = 1.0 (SeisSol-equivalent step, dt x3).  Same
#   config / mesh / sidecars / scheme / opt-stack / rank-count as the CFL A/B and
#   the caliper baseline 51443560 — so its region tree is directly comparable.
#   This profiles the SeisSol arm ONLY (arm B); default SAFS_CFL_DG_SAFETY=1.0.
#
# **** SEPARATE Caliper-enabled binary REQUIRED — do NOT clobber the A/B binary ****
#   Caliper needs a binary linked against libcaliper (MFEM_USE_CALIPER=YES); the
#   non-Caliper A/B binary is NOT.  Build it in a SEPARATE checkout so the running
#   A/B jobs' binary is untouched (a queued A/B job would otherwise pick up this
#   Caliper binary on start):
#       cd <separate-checkout-root>
#       git checkout system/spatial_dyn_driver && git pull --ff-only   # has --cfl-dg-safety
#       USE_CALIPER=YES bash build_expanse.sh
#   Then submit THIS script from that separate checkout (SEAS_MFEM_ROOT resolves
#   to whichever checkout you launch from).  The ldd pre-flight below hard-checks
#   the libcaliper link.
#
# Overrides: SAFS_CFL_DG_SAFETY (default 1.0 = the SeisSol arm; the A/B knob),
# SAFS_TFINAL (default 5s), SAFS_CONFIG, SAFS_MESH, SAFS_MIXED_FLUX (default none),
# SAFS_TIME_INTEGRATOR (default ader), SAFS_CFL (Courant number; default 0.5),
# SAFS_SHARED_CK / SAFS_DERIV_CACHE / SAFS_FACE_CACHE (1=on,0=off; all default ON),
# SAFS_CALI_CONFIG (full CALI_CONFIG override), SAFS_CALI_HZ (sample rate, 500),
# SAFS_CALI_SAMPLE (1 = add whole-program sample-report; needs WITH_LIBDW),
# SAFS_CALI_PAPI (1 = HW counters; needs WITH_PAPI) + SAFS_CALI_PAPI_COUNTERS,
# SAFS_SKIP_CALI_CHECK (1 = run even without libcaliper), SAFS_DIAG (0/1), SAFS_OUTDIR.
# =========================================================================

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

# Runtime shared libs — libcaliper.so FIRST (this job profiles).
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"

cd "${SEAS_MFEM_ROOT}"
grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not found.  Build WITH Caliper in a SEPARATE"
    echo "  checkout:  USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    exit 1
fi

# --- Caliper pre-flight: binary MUST be linked against libcaliper, else the
# MFEM_PERF_SCOPE regions are compiled out and the report is EMPTY. ------------
if ! ldd ./seas_spatial_dyn_driver 2>/dev/null | grep -qi 'libcaliper'; then
    echo "ERROR: ./seas_spatial_dyn_driver is NOT linked against libcaliper"
    echo "       (MFEM_USE_CALIPER=NO -> Caliper produces an EMPTY breakdown)."
    echo "  Rebuild WITH Caliper in a separate checkout:"
    echo "       USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    echo "  Or run anyway (bulk timing only): sbatch --export=ALL,SAFS_SKIP_CALI_CHECK=1 <this>"
    [[ "${SAFS_SKIP_CALI_CHECK:-0}" == "1" ]] || exit 1
    echo "  SAFS_SKIP_CALI_CHECK=1 set — continuing without Caliper (empty report)."
fi

SAFS_ROOT=safs/project_7.0_alternative
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_safs_seisol_v3_0_0_RSSRW_ALT_vw005.toml}"

DEFAULT_MESH="${SAFS_ROOT}/experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh"
MESH="${SAFS_MESH:-${DEFAULT_MESH}}"
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: mesh not found: $(pwd)/${MESH}  (gitignored — scp it first)"; exit 1
fi

VEL_H5="${SAFS_ROOT}/velocity/results/multiscale_statewise_cvm/velocity_safs.h5"
STR_H5="${SAFS_ROOT}/stress/results/csm_d12_mfem/safs_stress_csm_d12_mfem.h5"
for _h5 in "${VEL_H5}" "${STR_H5}"; do
    if [[ ! -f "${_h5}" ]]; then
        echo "ERROR: required sidecar missing: $(pwd)/${_h5}  (gitignored — scp it first)."; exit 1
    fi
done

grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null \
    || { echo "ERROR: ${CONFIG_TOML} has no sigma_n_strength_floor_pa"; exit 1; }

# === THE A/B KNOB — this profile is the SeisSol arm: cfl_dg_safety = 1.0 ======
SAFS_CFL_DG_SAFETY="${SAFS_CFL_DG_SAFETY:-1.0}"
if ! awk -v v="${SAFS_CFL_DG_SAFETY}" 'BEGIN{exit !(v+0 >= 1.0)}'; then
    echo "ERROR: SAFS_CFL_DG_SAFETY='${SAFS_CFL_DG_SAFETY}' must be >= 1.0 (1.0 = SeisSol floor)."; exit 1
fi
CFL_DG_ARG="--cfl-dg-safety ${SAFS_CFL_DG_SAFETY}"
CFLDG_TAG="cfldg${SAFS_CFL_DG_SAFETY//./p}"
# =============================================================================

SAFS_MIXED_FLUX="${SAFS_MIXED_FLUX:-none}"
case "${SAFS_MIXED_FLUX}" in none|adjacent|all_continuous) : ;; *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid"; exit 1 ;; esac
MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"

SAFS_TIME_INTEGRATOR="${SAFS_TIME_INTEGRATOR:-ader}"
case "${SAFS_TIME_INTEGRATOR}" in ader|rk4|rk45) : ;; *) echo "ERROR: SAFS_TIME_INTEGRATOR='${SAFS_TIME_INTEGRATOR}' invalid"; exit 1 ;; esac
if [[ "${SAFS_MIXED_FLUX}" != "none" && "${SAFS_TIME_INTEGRATOR}" == "ader" ]]; then
    echo "ERROR: mixed_flux=${SAFS_MIXED_FLUX} + ADER is UNSTABLE (driver aborts). Use rk4/rk45."; exit 1
fi
TI_ARG="--time-integrator ${SAFS_TIME_INTEGRATOR}"
if [[ "${SAFS_TIME_INTEGRATOR}" != "ader" ]]; then
    echo "WARNING: cfl_dg_safety only affects the ADER Dg path; SAFS_TIME_INTEGRATOR=${SAFS_TIME_INTEGRATOR} makes it inert."
fi

OPT_ARGS=""
if [[ "${SAFS_TIME_INTEGRATOR}" == "ader" && "${SAFS_SHARED_CK:-1}" == "1" ]]; then OPT_ARGS="${OPT_ARGS} --shared-ck-recursion"; fi
[[ "${SAFS_DERIV_CACHE:-1}" == "1" ]] && OPT_ARGS="${OPT_ARGS} --deriv-cache"
if [[ "${SAFS_TIME_INTEGRATOR}" == "ader" && "${SAFS_FACE_CACHE:-1}" == "1" ]]; then OPT_ARGS="${OPT_ARGS} --face-cache"; fi

TFINAL="${SAFS_TFINAL:-5s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

OUT_BASE="${SAFS_OUTDIR:-${SEAS_MFEM_ROOT}/miniapps/seas/safs}"
OUT="${OUT_BASE}/safs_dyn_cali_${CFLDG_TAG}_rs_alt_exp_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs/safs_expanse"
LOG="${LOG_DIR}/spatial_dyn_cali_${CFLDG_TAG}_rs_alt_exp_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

# === Caliper config — MAXIMAL detail (region table + Spot .cali) =============
CALI_REGION_TXT="${LOG%.log}.cali-region-report.txt"
CALI_SPOT="${LOG%.log}.cali"
CALI_SAMPLE_TXT="${LOG%.log}.cali-sample-report.txt"
export CALI_SAMPLER_FREQUENCY="${SAFS_CALI_HZ:-500}"
DEFAULT_CALI_CONFIG="runtime-report(calc.inclusive,region.count,profile.mpi,max_column_width=60,output=${CALI_REGION_TXT}),spot(output=${CALI_SPOT})"
if [[ "${SAFS_CALI_SAMPLE:-0}" == "1" ]]; then
    DEFAULT_CALI_CONFIG="${DEFAULT_CALI_CONFIG},sample-report(output=${CALI_SAMPLE_TXT})"
fi
export CALI_CONFIG="${SAFS_CALI_CONFIG:-${DEFAULT_CALI_CONFIG}}"
if [[ "${SAFS_CALI_PAPI:-0}" == "1" ]]; then
    export CALI_PAPI_COUNTERS="${SAFS_CALI_PAPI_COUNTERS:-PAPI_TOT_CYC,PAPI_TOT_INS,PAPI_L2_DCM,PAPI_L3_TCM}"
    export CALI_CONFIG="${CALI_CONFIG},hatchet-region-profile(output=${LOG%.log}.cali-papi-profile.json,use.mpi)"
    echo "Caliper PAPI: ${CALI_PAPI_COUNTERS}  (needs a WITH_PAPI Caliper build)"
fi

echo "Git commit:  $(git rev-parse HEAD)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Queue/size:  compute, 2N x 100r = ${SLURM_NTASKS} ranks, 200 GB/node, 02:00:00  (CALIPER PROFILE)"
echo "EXPERIMENT:  ratestate v3 ALT, UPWIND + ADER — Caliper profile @ cfl_dg_safety=${SAFS_CFL_DG_SAFETY} (${CFLDG_TAG}, SeisSol arm)"
echo "Config:      ${CONFIG_TOML}"
echo "Mesh:        ${MESH}"
echo "Integrator:  ${SAFS_TIME_INTEGRATOR}   Mixed flux: ${SAFS_MIXED_FLUX}   ${CFL_DG_ARG}"
echo "Opt levers:  ${OPT_ARGS:-<none>}"
echo "tfinal:      ${TFINAL}  (SHORT profile)   cfl(Courant): ${CFL_ARG:-<config default 0.5>}"
echo "Caliper:     ldd $(ldd ./seas_spatial_dyn_driver 2>/dev/null | grep -m1 -i libcaliper | awk '{print $1" -> "$3}')"
echo "CALI_CONFIG: ${CALI_CONFIG}"
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
    --paraview \
    --paraview-fault-vtu \
    --paraview-free-surface vtu \
    --paraview-max-snapshots 5 \
    --print-derived \
    --checkpoint-every 0 \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

# Confirm Lever A actually engaged (stale-binary detector) + show the step count.
echo
grep -E '\[time\] +(cfl *=.*cfl_dg_safety|dt_cfl|nsteps)' "${LOG}" || \
    echo "WARNING: no '[time] ... cfl_dg_safety' line — binary may predate --cfl-dg-safety (rebuild)."

echo
echo "=== Caliper outputs ==="
if [[ -f "${CALI_REGION_TXT}" ]]; then
    echo "region report -> ${CALI_REGION_TXT}   (Max>>Avg column = load imbalance)"
    echo "--- instrumented regions (slowest first) ---"
    head -40 "${CALI_REGION_TXT}"
else
    echo "WARNING: ${CALI_REGION_TXT} not written; CALI messages:"; grep -iE 'cali|caliper' "${LOG}" | head
fi
[[ -f "${CALI_SPOT}" ]] && echo "spot/.cali    -> ${CALI_SPOT}   (Hatchet/Thicket per-rank detail)"
[[ -f "${CALI_SAMPLE_TXT}" ]] && { echo "sample report -> ${CALI_SAMPLE_TXT}"; head -40 "${CALI_SAMPLE_TXT}"; }

if [[ ${RC} -ne 0 ]]; then echo "driver_rc=${RC} — see ${LOG}."; exit ${RC}; fi
echo "driver_rc=${RC} — Caliper profile (cfl_dg_safety=${SAFS_CFL_DG_SAFETY}) done."
exit ${RC}
