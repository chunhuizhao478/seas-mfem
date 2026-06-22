#!/bin/bash
#SBATCH -J bp5_prod_gmres
# Bare filenames: logs land in the dir you run `sbatch` from (no nested
# jobs/bp5_spatial/ duplication).
#SBATCH -o bp5_production_gmres_%j.out
#SBATCH -e bp5_production_gmres_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5-QD production run with the iterative GMRES+AMG solver (seas_driver).
# Runs to t_final (~1800 yr) or until the 48 h walltime; output (station traces
# + global Vmax) is flushed every step, so partial results are kept on TIMEOUT.
#
# This job does NOT build — build the driver ONCE first, then submit this:
#   sbatch jobs/bp5_spatial/rebuild_seas_driver_frontera.sbatch   # wait for it
#   sbatch jobs/bp5_spatial/bp5_production_gmres_normal_48hr.sbatch
#
# NOTE: seas_driver has no working checkpoint/restart, so a TIMEOUT cannot be
# resumed — it is a single walltime-bounded run.
#
# Acceptance / what to watch:
#   - every step's GMRES solve CONVERGES (no "did NOT converge" abort)
#   - the nucleation patch (V_max ~ V_nuc = 0.01) develops into an event (V_max
#     spikes toward ~1 m/s), then arrests and V_max drops; dt ramps to 0.1 yr
#     interseismically (the -ts_adapt_dt_max cap)
#   - no NaN, traction ~25 MPa scale (never 2.99e12)

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

# hypre/2.31.0 FIRST on LD_LIBRARY_PATH (avoid pre-2.31 petsc/mumps libHYPRE
# shadowing -> "undefined symbol: HYPRE_Initialize").
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

# Resolve seas-mfem root by walking UP until miniapps/seas/Makefile is found.
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
cd "${SEAS_MFEM_ROOT}/miniapps/seas"

# Use the binary built by rebuild_seas_driver_frontera.sbatch — this job does
# not build.  If it is missing, submit the rebuild job first.
if [[ ! -x ./seas_driver ]]; then
    echo "ERROR: ./seas_driver not found — build it first:"
    echo "  sbatch jobs/bp5_spatial/rebuild_seas_driver_frontera.sbatch"
    exit 1
fi

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"

RESULT_DIR="bp5/results_prod_gmres_${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

ibrun ./seas_driver config/bp5_production_gmres.toml \
      --override output.output_dir="${RESULT_DIR}" \
      --override output.output_prefix="bp5_prod_gmres_${SLURM_JOB_ID}"
