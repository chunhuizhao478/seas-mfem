#!/bin/bash
#SBATCH -J bp5_cali_200
#SBATCH -o bp5_caliper_200steps_%j.out
#SBATCH -e bp5_caliper_200steps_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=50
#SBATCH --time=02:00:00
#SBATCH --mem=100000M
#SBATCH -A ddp408

# =========================================================================
# BP5-QD CALIPER PROFILING run — iterative GMRES+AMG, FIRST 200 STEPS ONLY.
# *** EXPANSE (SDSC) job — 8N x 50r = 400 ranks / 2 h / compute / ddp408. ***
# Driver: seas_driver (BP5 TOML driver).  Config: config/bp5_caliper_500yr.toml
# (Expanse sibling of the Frontera bp5_caliper_gmres_500yr_normal_48hr.sbatch;
#  Expanse modules + srun --mpi=pmi2 + extern/ from-source libs from build_expanse.sh.)
# =========================================================================
# The full 500 yr run is TOO SLOW to profile end-to-end, so this run is capped
# at --max-steps 200 (seas_driver -> petsc_ode->SetMaxIter(200) on the PETSc-TS
# path this config uses) — just enough to see WHERE the per-step cost goes
# (init/assembly/factorization vs GMRES iters vs friction vs I/O/MPI).
#
# Caliper is driven ENTIRELY by CALI_CONFIG (read at init) — no driver flag.
# Requires a CALIPER-ENABLED binary built WITH_LIBDW (symbollookup) so the
# sample-report channel works.  Build this (separate) Expanse checkout with:
#     USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh   # builds MFEM + Caliper
#     ( cd ${SEAS_MFEM_ROOT}/miniapps/seas && make seas_driver )
# build_expanse.sh adds -DWITH_LIBDW=ON when libdw/elfutils is found (Expanse /usr
# has it), so symbollookup (PC -> function name) is compiled in.  We hard-check the
# libcaliper link below (override BP5_SKIP_CALI_CHECK=1).
#
# CHANNELS (3, complementary):
#   runtime-report : time of each INSTRUMENTED MFEM_PERF_SCOPE region.  NOTE:
#                    seas_driver (QD) has FEW marked regions -> this table is
#                    SPARSE; the sampling channel below is the primary signal.
#   spot           : .cali for the LLNL Spot UI / Python Hatchet (per-rank).
#   sample-report  : timer SAMPLING of the WHOLE program (every function: init,
#                    PETSc/HYPRE/MUMPS, GMRES iters, HDF5 I/O, MPI).  THIS answers
#                    "which function is slowest" end-to-end.  Needs WITH_LIBDW;
#                    default ON, drop with BP5_CALI_SAMPLE=0.
#
# Overrides: BP5_CONFIG, BP5_MAX_STEPS (default 200), BP5_CALI_CONFIG (full
# CALI_CONFIG override), BP5_CALI_HZ (sampling rate/rank, default 200),
# BP5_CALI_SAMPLE (0 = drop sample-report; needs a WITH_LIBDW build, default ON),
# BP5_SKIP_CALI_CHECK (1 = run even without libcaliper), BP5_OUTDIR.

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
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
_walk_up_dir="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5 6 7; do
    if [[ -f "${_walk_up_dir}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"; break
    fi
    _walk_up_dir="${_walk_up_dir}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-<unset>}."
    echo "  Set it explicitly: sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi

# Runtime shared libs.  libcaliper.so FIRST (this job profiles); then the usual
# hypre/hdf5/zfp/petsc from build_expanse.sh's from-source deps.
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hypre/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"

cd "${SEAS_MFEM_ROOT}"
grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
if [[ ! -x ./seas_driver ]]; then
    echo "ERROR: ./seas_driver not found.  Build it FIRST on Expanse (with Caliper):"
    echo "  USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh && make seas_driver"
    exit 1
fi

# Caliper pre-flight: the binary MUST be linked against libcaliper, else the
# MFEM_PERF_SCOPE regions are compiled out and the report is EMPTY.
if ! ldd ./seas_driver 2>/dev/null | grep -qi 'libcaliper'; then
    echo "ERROR: ./seas_driver is NOT linked against libcaliper"
    echo "       (MFEM_USE_CALIPER=NO -> the per-region breakdown would be empty)."
    echo "  Rebuild WITH Caliper:  USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh && make seas_driver"
    echo "  To run anyway (no per-region/sample data):"
    echo "        sbatch --export=ALL,BP5_SKIP_CALI_CHECK=1 <this file>"
    [[ "${BP5_SKIP_CALI_CHECK:-0}" == "1" ]] || exit 1
    echo "  BP5_SKIP_CALI_CHECK=1 — continuing without Caliper (empty report)."
fi

CONFIG_TOML="${BP5_CONFIG:-config/bp5_caliper_500yr.toml}"
[[ -f "${CONFIG_TOML}" ]] || { echo "ERROR: config not found: $(pwd)/${CONFIG_TOML}"; exit 1; }

MAX_STEPS="${BP5_MAX_STEPS:-200}"

OUT_BASE="${BP5_OUTDIR:-bp5}"
OUT="${OUT_BASE}/results_caliper_200steps_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/bp5_spatial"
LOG="${LOG_DIR}/bp5_caliper_200steps_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

# === Caliper config — 3 complementary channels (built-in recipes) ===========
CALI_REGION_TXT="${LOG%.log}.cali-region-report.txt"   # instrumented region table (sparse for QD)
CALI_SPOT="${LOG%.log}.cali"                           # Spot .cali (Hatchet/Thicket, per-rank)
CALI_SAMPLE_TXT="${LOG%.log}.cali-sample-report.txt"   # SAMPLING: every function, whole program
export CALI_SAMPLER_FREQUENCY="${BP5_CALI_HZ:-200}"    # samples/s/rank
# sample-report needs Caliper's symbollookup service (sample PC -> function name) =
# a Caliper built WITH_LIBDW.  ConfigManager aborts the ENTIRE CALI_CONFIG if a
# channel's required service is missing (taking runtime-report + spot down too).
# Default ON because the QD region table is sparse — the sampler is the primary
# signal; set BP5_CALI_SAMPLE=0 on a Caliper build that lacks libdw.
DEFAULT_CALI_CONFIG="runtime-report(calc.inclusive,region.count,profile.mpi,max_column_width=60,output=${CALI_REGION_TXT}),spot(output=${CALI_SPOT})"
if [[ "${BP5_CALI_SAMPLE:-1}" != "0" ]]; then
    DEFAULT_CALI_CONFIG="${DEFAULT_CALI_CONFIG},sample-report(output=${CALI_SAMPLE_TXT})"
fi
export CALI_CONFIG="${BP5_CALI_CONFIG:-${DEFAULT_CALI_CONFIG}}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Queue/size:  compute, 8N x 50r = ${SLURM_NTASKS} ranks, 100 GB/node, 02:00:00  (CALIPER PROFILE)"
echo "Config:      ${CONFIG_TOML}  (GMRES+AMG; --max-steps ${MAX_STEPS} caps the run)"
echo "Caliper:     ldd $(ldd ./seas_driver 2>/dev/null | grep -m1 -i libcaliper | awk '{print $1" -> "$3}')"
echo "CALI_CONFIG: ${CALI_CONFIG}"
echo "Output dir:  ${OUT}"

srun --mpi=pmi2 -n "${SLURM_NTASKS}" stdbuf -oL -eL ./seas_driver "${CONFIG_TOML}" \
     --max-steps "${MAX_STEPS}" \
     --override output.output_dir="${OUT}" \
     --override output.output_prefix="bp5_caliper_200steps_${SLURM_JOB_ID}" \
     > "${LOG}" 2>&1
RC=$?

echo
echo "=== Caliper outputs ==="
if [[ -f "${CALI_REGION_TXT}" ]]; then
    echo "region report -> ${CALI_REGION_TXT}   (instrumented regions; sparse for QD driver)"
    echo "--- instrumented regions (slowest first) ---"
    head -40 "${CALI_REGION_TXT}"
else
    echo "WARNING: ${CALI_REGION_TXT} not written; CALI messages:"
    grep -iE 'cali|caliper' "${LOG}" | head
fi
echo
if [[ -f "${CALI_SAMPLE_TXT}" ]]; then
    echo "sample report -> ${CALI_SAMPLE_TXT}   (WHICH FUNCTION IS SLOWEST — whole program)"
    echo "--- top sampled functions (self time) ---"
    head -40 "${CALI_SAMPLE_TXT}"
else
    echo "NOTE: ${CALI_SAMPLE_TXT} not written — symbollookup missing (needs WITH_LIBDW),"
    echo "      BP5_CALI_SAMPLE=0, or no clean exit.  CALI messages:"
    grep -iE 'sampl|symbollookup|callpath|cali' "${LOG}" | head
fi
echo
[[ -f "${CALI_SPOT}" ]] && echo "spot/.cali    -> ${CALI_SPOT}   (Hatchet/Thicket or Spot UI: PER-RANK detail)"

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG}."; exit ${RC}
fi
echo "driver_rc=${RC} — 200-step profiling run done."
exit ${RC}
