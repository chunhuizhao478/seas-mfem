#!/bin/bash
#SBATCH -J bp5_iter
#SBATCH -o jobs/bp5_spatial/bp5_verify_iterative_dev_2hr_%j.out
#SBATCH -e jobs/bp5_spatial/bp5_verify_iterative_dev_2hr_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# BP5-QD iterative-solver verification gate (seas_driver — the BP5 TOML driver,
# NOT seas_spatial_dyn_driver).  Exercises the np>1 AMG path on the production km
# mesh: runs GMRES+AMG then CG+AMG, 50 steps each.  Compare both against the
# MUMPS golden to < 1e-6.
#
# For a clean, self-verifying build first (recommended after a git pull — the
# seas Makefile has no header-dependency tracking):
#   D=miniapps/seas/jobs/bp5_spatial
#   bid=$(sbatch --parsable ${D}/rebuild_seas_driver_frontera.sbatch)
#   sbatch --dependency=afterok:${bid} ${D}/bp5_verify_iterative_dev_2hr.sbatch
# (this job also does a clean in-job relink, so it is safe to submit standalone.)
#
# Acceptance per solver:
#   - every step's iterative solve CONVERGES (no "did NOT converge" abort)
#   - V_max stays physical (no NaN, traction ~25 MPa scale, never 2.99e12)
#   - station traces match the MUMPS 50-step run to < 1e-6 (regression_check.py)

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

# hypre/2.31.0 first on LD_LIBRARY_PATH (avoid old petsc libHYPRE shadowing —
# HYPRE_Init->HYPRE_Initialize rename, see frontera-hypre-initialize memory)
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

# Resolve seas-mfem root by walking UP until miniapps/seas/Makefile is found
# (works regardless of the Frontera checkout path).
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
cd "${SEAS_MFEM_ROOT}"

JOBS=8

grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas

# ensure toml11 submodule is present (SEAS_USE_TOML gating) on fresh checkouts
test -f extern/toml11/toml.hpp || test -f extern/toml11/include/toml.hpp || \
   git submodule update --init extern/toml11 || \
   { echo "ERROR: toml11 not initialized. Run: git submodule update --init miniapps/seas/extern/toml11"; exit 1; }

# Clean relink of the single driver TU (kills the stale-.o hazard after a pull;
# rebuild_seas_driver_frontera.sbatch is the fuller pre-flight build).
rm -f drivers/seas_driver.o
make seas_driver SEAS_USE_TOML=YES -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"

# GMRES first (the robust DG iterative solver, ~100 iters/solve), then CG.
# If CG converges slowly at this rank count it cannot starve the GMRES gate.
for SOLVER in gmres cg; do
  echo "============================================================"
  echo "BP5 iterative verification: solver=${SOLVER}, 50 steps, 400 ranks"
  echo "============================================================"
  # create output dir here too, so the run does not depend on the driver-side
  # mkdir fix being present on this checkout
  mkdir -p "bp5/results_verify_50step_${SOLVER}"
  ibrun ./seas_driver config/bp5_verification_50step_${SOLVER}.toml
done
