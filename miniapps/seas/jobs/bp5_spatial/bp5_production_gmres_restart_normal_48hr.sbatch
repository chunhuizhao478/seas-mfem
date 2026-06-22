#!/bin/bash
#SBATCH -J bp5_prod_rst
# Bare filenames: logs land in the dir you run `sbatch` from.
#SBATCH -o bp5_production_gmres_restart_%j.out
#SBATCH -e bp5_production_gmres_restart_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5-QD PRODUCTION run to t_final (~1800 yr), iterative GMRES+AMG, WITH
# checkpoint/restart so the cycle can span multiple 48 h walltime segments.
# Driver: seas_driver (the BP5 TOML driver) — NOT seas_spatial_dyn_driver.
#
# === CHAINED-SEGMENT RESTART ==============================================
# Each submission is one segment writing into BASE/segment_NNN.  Checkpoints
# are written every CKPT steps; the latest overwrites, so a segment can resume
# from the previous one's last checkpoint.  Submit segment 1, let it run/timeout,
# then submit segment 2 (which restarts from segment 1), etc.:
#
#   sbatch jobs/bp5_spatial/rebuild_seas_driver_frontera.sbatch          # build once
#   sbatch jobs/bp5_spatial/bp5_production_gmres_restart_normal_48hr.sbatch              # SEGMENT=1
#   # after it finishes/timeouts:
#   sbatch --export=ALL,SEGMENT=2 jobs/bp5_spatial/bp5_production_gmres_restart_normal_48hr.sbatch
#   sbatch --export=ALL,SEGMENT=3 jobs/bp5_spatial/bp5_production_gmres_restart_normal_48hr.sbatch
#   ...                                                                  # until t_final reached
#
# The restart loads (state, t, dt, event counters, PETSc dt-next) and SKIPS the
# 4-phase SetInitialCondition.  Output (station traces) is per-segment; stitch
# the segment_NNN dirs in post-processing.  A run that reaches t_final exits 0
# ("SEAS Simulation Complete"); a walltime kill is a SLURM TIMEOUT (resume via
# the next segment).
#
# Overrides: BASE (results base dir), SEGMENT (1..), CKPT (checkpoint interval steps).

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

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

if [[ ! -x ./seas_driver ]]; then
    echo "ERROR: ./seas_driver not found — build it first:"
    echo "  sbatch jobs/bp5_spatial/rebuild_seas_driver_frontera.sbatch"
    exit 1
fi

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"

# --- segment / restart bookkeeping ---------------------------------------
BASE="${BASE:-bp5/results_prod_gmres_restart}"
SEGMENT="${SEGMENT:-1}"
CKPT="${CKPT:-1000}"                       # checkpoint every CKPT accepted steps
PREFIX="bp5_prod_gmres"
seg() { printf "segment_%03d" "$1"; }

OUT="${BASE}/$(seg "${SEGMENT}")"
mkdir -p "${OUT}"

RESTART_ARG=""
if [[ "${SEGMENT}" -gt 1 ]]; then
    PREV="${BASE}/$(seg "$((SEGMENT-1))")"
    # the driver restarts from the previous segment's checkpoint prefix;
    # per-rank files are named <prefix>_checkpoint_r<rank>.txt (io/checkpoint.hpp).
    if ! ls "${PREV}/${PREFIX}"_checkpoint_r*.txt >/dev/null 2>&1; then
        echo "WARNING: no checkpoint found under ${PREV}/ — segment ${SEGMENT} may have nothing to restart from."
    fi
    RESTART_ARG="--restart ${PREV}/${PREFIX}"
    echo "Restarting segment ${SEGMENT} from ${PREV}/${PREFIX}"
else
    echo "Initial segment ${SEGMENT} (cold start)."
fi

echo "BASE=${BASE}  SEGMENT=${SEGMENT}  OUT=${OUT}  CKPT=${CKPT}"

ibrun ./seas_driver config/bp5_production_gmres.toml \
      --override output.output_dir="${OUT}" \
      --override output.output_prefix="${PREFIX}" \
      --checkpoint-interval "${CKPT}" \
      ${RESTART_ARG}
