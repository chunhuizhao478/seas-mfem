#!/bin/bash
#SBATCH -J seas_qd_bp5_mumps
#SBATCH -o seas_qd_bp5_mumps_%j.out
#SBATCH -e seas_qd_bp5_mumps_%j.err
#SBATCH -p development
#SBATCH -N 4
#SBATCH -n 200
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# =========================================================================
# BP5 run on the quasi-dynamic (QD) SEAS driver with the MUMPS DIRECT solver.
# =========================================================================
# This is the MUMPS sibling of spatial_seas_bp5_50step.sbatch.  It runs
# config/safs_qd/bp5_phase5_mumps.toml = uniform friction + fault-local prestress
# on the METRES BP5 mesh (bp5_tandem_exact_meters.msh), solved with [solver]
# type="mumps" (exact LU factorization, factorized ONCE and reused per RK stage).
#
# WHY METRES + MUMPS:
#  - The QD driver consumes mesh coords in SI directly (NO km->m scaling), so the
#    BP5 mesh MUST be in metres or the SIPG penalty (~mu*p^2/h) is mis-scaled and
#    the first loaded RK step blows up at the fault edge (traction ~1e12 Pa).
#    bp5_tandem_exact_meters.msh is the km reference x1000.  Verified locally:
#    a clean physical step (V_max ~ V_init), no blow-up.
#  - MUMPS is the exact direct solver (vs cg+BoomerAMG which is slow / hard to
#    converge on DG elasticity).  Frontera's mumps/5.3 module is used here.
#
# PREREQUISITE: build the binary ONCE on Frontera first (this job does NOT build):
#   cd miniapps/seas && make seas_spatial_seas_driver
#
# Submit (from miniapps/seas):
#   sbatch jobs/safs_qd/spatial_seas_bp5_mumps.sbatch
# Override config / step count:
#   SEAS_QD_CONFIG=... SEAS_QD_MAX_STEPS=... sbatch jobs/safs_qd/spatial_seas_bp5_mumps.sbatch
#
# Expected (PASS): rank-0 stdout shows
#   [spatial_seas] SetInitialCondition done: V_max = <finite ~1e-9> m/s
#   [spatial_seas] step N  t=..  dt=..  V_max=..   (V_max finite, dt growing)
#   ... and NO "non-physical fault traction" / "Non-finite V" abort.
#
# MUMPS NOTE: 3D DG factorization is memory-heavy.  If MUMPS OOMs, raise the
# workspace via ICNTL(14) (a percentage relaxation) or reduce ranks/node.

set -euo pipefail
export LC_ALL=C
export LANG=C

# --- Frontera module stack (match the working QD sbatch; [[feedback_sbatch_modules]]) ---
module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
# Put the module HYPRE/MUMPS libs FIRST so the old petsc libHYPRE does not shadow
# HYPRE_Initialize ([[frontera-hypre-initialize-ldlibrarypath-shadowing]]).
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"

# --- Run from miniapps/seas (where the binary + config live) ---
cd "${SLURM_SUBMIT_DIR}"

BINARY=./seas_spatial_seas_driver
if [[ ! -x "${BINARY}" ]]; then
    echo "ERROR: ${BINARY} not found/executable in $(pwd)."
    echo "       Build it first:  cd miniapps/seas && make seas_spatial_seas_driver"
    exit 1
fi

CONFIG_TOML="${SEAS_QD_CONFIG:-config/safs_qd/bp5_phase5_mumps.toml}"
if [[ ! -f "${CONFIG_TOML}" ]]; then
    echo "ERROR: config ${CONFIG_TOML} not found."
    exit 1
fi

MAX_STEPS="${SEAS_QD_MAX_STEPS:-50}"

echo "[job] $(date)  ${BINARY}  config=${CONFIG_TOML}  max-steps=${MAX_STEPS}  solver=mumps"
ibrun "${BINARY}" \
    --config "${CONFIG_TOML}" \
    --max-steps "${MAX_STEPS}"
echo "[job] done $(date)"
