#!/bin/bash
#SBATCH -J bp5_cali
# Bare filenames: logs land in the dir you run `sbatch` from.
#SBATCH -o bp5_caliper_gmres_500yr_%j.out
#SBATCH -e bp5_caliper_gmres_500yr_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# =========================================================================
# BP5-QD CALIPER PROFILING run — iterative GMRES+AMG, t_final = 500 yr.
# *** FRONTERA (TACC) job — 8N x 50r = 400 ranks / 48 h / normal queue. ***
# Driver: seas_driver (BP5 TOML driver).  Config: config/bp5_caliper_500yr.toml
# (modeled on safs/safs_frontera ... adapted from the Expanse Caliper job
#  spatial_dyn_ratestate_..._caliper_..._expanse_2N100r_compute_4hr_safs.sbatch).
# =========================================================================
# Same physics as the production run; Caliper gives a PER-COMPONENT speed
# breakdown.  Caliper is driven ENTIRELY by CALI_CONFIG (read at init) — no
# driver flag.  Requires a CALIPER-ENABLED binary (the Frontera MFEM build links
# libcaliper; we hard-check below, override with BP5_SKIP_CALI_CHECK=1).
#
# CHANNELS (3, complementary):
#   runtime-report         : time of each INSTRUMENTED MFEM_PERF_SCOPE region
#                            (inclusive+exclusive, call counts, per-region MPI).
#                            NOTE: seas_driver (QD) has fewer marked regions than
#                            the dynamic driver, so this table is sparse — the
#                            sampling channel below is the primary signal here.
#   spot                   : .cali for the LLNL Spot UI / Python Hatchet (per-rank).
#   callpath-sample-report : timer SAMPLING of the WHOLE program (every function:
#                            init, PETSc/HYPRE/MUMPS, GMRES iters, HDF5 I/O, MPI).
#                            THIS answers "which function is slowest" end-to-end.
#
# *** Caliper writes its report at clean program exit.  A 48 h walltime KILL
# (SLURM TIMEOUT) may produce NO report — lower [time].t_final in the config if
# you need a guaranteed profile (the per-region breakdown is similar at 50 yr). ***
#
# Overrides: BP5_CONFIG, BP5_CALI_CONFIG (full CALI_CONFIG override),
# BP5_CALI_HZ (sampling rate/rank, default 200), BP5_SKIP_CALI_CHECK (1 = run even
# without libcaliper), BP5_OUTDIR.

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

# Resolve seas-mfem root (walk UP to miniapps/seas/Makefile) BEFORE the lib paths.
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
_walk_up_dir="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5 6 7; do
    if [[ -f "${_walk_up_dir}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"
        break
    fi
    _walk_up_dir="${_walk_up_dir}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-<unset>}."
    echo "  Set it explicitly: sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi

# libcaliper FIRST (this job profiles), then hypre/2.31.0 (avoid pre-2.31 petsc
# libHYPRE shadowing), then the rest of the TACC module libs.
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

cd "${SEAS_MFEM_ROOT}/miniapps/seas"

if [[ ! -x ./seas_driver ]]; then
    echo "ERROR: ./seas_driver not found — build it first:"
    echo "  sbatch jobs/bp5_spatial/rebuild_seas_driver_frontera.sbatch"
    exit 1
fi

# Caliper pre-flight: the binary MUST be linked against libcaliper, else the
# MFEM_PERF_SCOPE regions are compiled out and the report is EMPTY.
if ! ldd ./seas_driver 2>/dev/null | grep -qi 'libcaliper'; then
    echo "ERROR: ./seas_driver is NOT linked against libcaliper"
    echo "       (MFEM_USE_CALIPER=NO -> the per-region breakdown would be empty)."
    echo "  Rebuild MFEM/seas_driver with Caliper, or run anyway:"
    echo "        sbatch --export=ALL,BP5_SKIP_CALI_CHECK=1 <this file>"
    [[ "${BP5_SKIP_CALI_CHECK:-0}" == "1" ]] || exit 1
    echo "  BP5_SKIP_CALI_CHECK=1 — continuing without Caliper (empty report)."
fi

CONFIG_TOML="${BP5_CONFIG:-config/bp5_caliper_500yr.toml}"
[[ -f "${CONFIG_TOML}" ]] || { echo "ERROR: config not found: $(pwd)/${CONFIG_TOML}"; exit 1; }

OUT_BASE="${BP5_OUTDIR:-bp5}"
OUT="${OUT_BASE}/results_caliper_500yr_${SLURM_JOB_ID}"
mkdir -p "${OUT}"

# === Caliper config — 3 complementary channels (built-in recipes) ===========
CALI_REGION_TXT="bp5_caliper_500yr_${SLURM_JOB_ID}.cali-region-report.txt"
CALI_SPOT="bp5_caliper_500yr_${SLURM_JOB_ID}.cali"
CALI_SAMPLE_TXT="bp5_caliper_500yr_${SLURM_JOB_ID}.cali-sample-report.txt"
export CALI_SAMPLER_FREQUENCY="${BP5_CALI_HZ:-200}"      # samples/s/rank
DEFAULT_CALI_CONFIG="runtime-report(calc.inclusive,region.count,profile.mpi,max_column_width=60,output=${CALI_REGION_TXT}),spot(output=${CALI_SPOT}),callpath-sample-report(output=${CALI_SAMPLE_TXT})"
export CALI_CONFIG="${BP5_CALI_CONFIG:-${DEFAULT_CALI_CONFIG}}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Cluster:     Frontera (TACC)   Launcher: ibrun"
echo "Queue/size:  normal, 8N x 50r = ${SLURM_NTASKS} ranks, 48:00:00   (CALIPER PROFILE)"
echo "Config:      ${CONFIG_TOML}  (GMRES+AMG, t_final=500 yr)"
echo "Caliper:     $(ldd ./seas_driver 2>/dev/null | grep -m1 -i libcaliper | awk '{print $1" -> "$3}')"
echo "CALI_CONFIG: ${CALI_CONFIG}"
echo "Output dir:  ${OUT}"

ibrun ./seas_driver "${CONFIG_TOML}" \
      --override output.output_dir="${OUT}" \
      --override output.output_prefix="bp5_caliper_500yr_${SLURM_JOB_ID}"
RC=$?

echo
echo "=== Caliper outputs ==="
if [[ -f "${CALI_REGION_TXT}" ]]; then
    echo "region report -> ${CALI_REGION_TXT}   (instrumented regions; sparse for QD driver)"
    head -40 "${CALI_REGION_TXT}"
else
    echo "NOTE: ${CALI_REGION_TXT} not written (no clean exit yet, or no marked regions)."
fi
echo
if [[ -f "${CALI_SAMPLE_TXT}" ]]; then
    echo "sample report -> ${CALI_SAMPLE_TXT}   (WHICH FUNCTION IS SLOWEST — whole program)"
    head -40 "${CALI_SAMPLE_TXT}"
else
    echo "NOTE: ${CALI_SAMPLE_TXT} not written (sampler/callpath service missing, or no clean exit)."
fi
[[ -f "${CALI_SPOT}" ]] && echo "spot/.cali    -> ${CALI_SPOT}   (Hatchet/Thicket / Spot UI: per-rank)"

echo "driver_rc=${RC}"
exit ${RC}
