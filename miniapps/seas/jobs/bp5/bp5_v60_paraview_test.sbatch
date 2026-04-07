#!/bin/bash
#SBATCH -J bp5_pv
#SBATCH -o bp5_v60_paraview_%j.out
#SBATCH -e bp5_v60_paraview_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5 v60 ParaView: Full simulation to capture the ~25 yr blowup.
# Fault surface VTU output every 0.01 yr (~3.65 days) = ~2500 snapshots.
# Dense enough to see the normal stress erosion develop frame by frame
# near the failure point at ~25.27 yr.
#
# Fields: displacement (3D), slip, slip_rate, traction, state_variable,
#         normal_stress, param_a, param_Dc, fault_x2, fault_x3.
# Open FaultSurface/fault_surface.pvd in ParaView.

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

RESULT_DIR="bp5/results_v60_paraview_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_v60_pv_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo "Full run, fault surface VTU every 0.01 yr (~2500 snapshots over 25 yr)"

# 0.01 yr = 315576 s (seconds_per_year = 31557600)
ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --checkpoint-interval 5000 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --paraview-dt 315576
