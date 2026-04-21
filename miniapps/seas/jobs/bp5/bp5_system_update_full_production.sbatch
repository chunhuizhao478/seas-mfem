#!/bin/bash
#SBATCH -J bp5_suf
#SBATCH -o bp5_system_update_full_production_%j.out
#SBATCH -e bp5_system_update_full_production_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# system_update branch: full production run (400 ranks, ~1800 yr).
# Final verification after all 7 refactoring phases (0-6).
# Compares against seas-mfem parallel golden (golden_parallel_400r_p1_1000m).
#
# Before submitting, on the login node:
#   cd /scratch2/10024/zhaochun/seas-project/seas-mfem
#   git fetch origin && git checkout system_update && git pull
#   cd miniapps/seas && sbatch jobs/bp5/bp5_system_update_full_production.sbatch
#
# What to check:
#   1. Recurrence interval ~240 yr (same as v62/Phase 0 golden)
#   2. Coseismic slip ~3-5 m per event
#   3. Station time series match golden within recompilation noise
#   4. No blowup, no dt→0, no bracket failures
#
# Post-run comparison:
#   python3 scripts/regression_check.py \
#       --sim-dir $RESULT_DIR --sim-prefix $OUTPUT_PREFIX \
#       --ref-dir bp5/benchmark_data/golden_parallel_400r_p1_1000m \
#       --tolerance 1e-2

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${LD_LIBRARY_PATH}"

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
make clean
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="bp5/results_system_update_full_production_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_system_update_full_production_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "system_update full production: 400 ranks, 48hr, ~1800 yr earthquake cycles"
echo "Phases 0-6 complete. Final verification run."

# Adaptive fault-surface ParaView PVD output (fault-only: no volume PVD).
# Rationale for the values below on a 1800-yr / ~8-event / ~100-s-rupture run:
#   dt_inter = 15 yr  -> ~120 interseismic frames (was 360 at 5 yr — cut 3x
#                        after observing 1.86M-inode blowup on prior run).
#   dt_nu    = 30 s   -> ~3 frames per event's nucleation phase.
#   dt_co    = 2.0 s  -> ~50 frames per ~100-s rupture (was ~200 at 0.5 s).
#   hyst     = 10     -> enter at 1e-3 / 1e-6, exit at 1e-4 / 1e-7 (no chatter).
# Expected disk footprint: ~500 fault-surface frames x ~400 ranks x ~12 fields
# ~= 200K VTU files, well under Frontera scratch inode quota (~2M).
# If post-run inode pressure still shows, tar per-cycle after the run:
#   cd "${RESULT_DIR}/ParaView/FaultSurface"
#   for p in fault_surface_c*.pvtu; do c=${p#fault_surface_c}; c=${c%.pvtu}; \
#     tar -cf fault_surface_c${c}.tar fault_surface_r*_c${c}.vtu && \
#     rm fault_surface_r*_c${c}.vtu; done
ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --ref-dir bp5/benchmark_data/golden_parallel_400r_p1_1000m \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --checkpoint-interval 5000 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --regression-tolerance 1e-2 \
      --verify \
      --paraview \
      --paraview-fault-only \
      --paraview-adaptive \
      --paraview-v-co        1e-3 \
      --paraview-v-nu        1e-6 \
      --paraview-hyst        10.0 \
      --paraview-dt-co       2.0 \
      --paraview-dt-nu       30.0 \
      --paraview-dt-inter-yr 15.0
