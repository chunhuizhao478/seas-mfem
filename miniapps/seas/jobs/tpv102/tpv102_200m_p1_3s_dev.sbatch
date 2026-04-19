#!/bin/bash
#SBATCH -J tpv102_200m_p1
#SBATCH -o tpv102_200m_p1_%j.out
#SBATCH -e tpv102_200m_p1_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# SCEC TPV102 dynamic rupture — order-1, Frontera dev queue, 2 hr.
#
# Physics:
#   Mesh: tpv102_200m.msh (h = 200 m fault resolution, ~2.46 M tets)
#   Order: 1 (p1 L2 on tets) — matches the direct benchmark baseline:
#          Aagaard's PyLith v1.8.1a "Tet4 200m" submission on the SCEC
#          TPV102 comparison page.  Using the same element order/mesh
#          resolution lets us cross-plot with that reference directly.
#   Absorbing BCs on outer walls + bottom (attr 5), free surface top (attr 1),
#   embedded fault at Y=0 (attr 3).
#
# Run parameters:
#   tfinal = 3.0 s
#   dt     = 5 ms  →  600 explicit-RK4 steps
#            With p=1 on h=200 m tets (cp ≈ 5800 m/s) this is inside the
#            RK4 stability region.  The driver still prints the usual
#            CFL comparison line in the log; safely ignore the
#            "dt_override > dt_cfl" advisory — dt_cfl is reported at
#            cfl_factor=0.5 and the explicit RK4 envelope is larger.
#   ParaView fault-surface PVD every 0.1 s; volume PVD suppressed.
#   Station output every 0.1 s at the SCEC probe points.
#
# Parallel setup: 8 nodes × 50 ranks = 400 ranks.  Dev queue, 2 hr.
#
# Before submitting:
#   1. Build driver on Frontera:
#        cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas
#        make clean && make seas_tpv102_driver -j8
#   2. Ensure the mesh exists (gitignored, must be regenerated):
#        ls tpv102/mesh/tpv102_200m.msh
#      If missing:
#        conda activate pythonenv
#        gmsh -3 tpv102/mesh/tpv102_200m.geo -o tpv102/mesh/tpv102_200m.msh
#   3. Submit:
#        sbatch jobs/tpv102/tpv102_200m_p1_3s_dev.sbatch

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

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
make seas_tpv102_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_200m_p1_3s_dev_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "tfinal=3.0 s, dt=5 ms, order=1 (matches PyLith Tet4 200m)"
echo "Output every 0.1 s, fault-surface PVD only"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --tfinal 3.0 \
      --dt 0.005 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv
