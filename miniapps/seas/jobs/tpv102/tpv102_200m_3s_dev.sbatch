#!/bin/bash
#SBATCH -J tpv102_200m_3s
#SBATCH -o tpv102_200m_3s_%j.out
#SBATCH -e tpv102_200m_3s_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# SCEC TPV102 dynamic rupture — Frontera dev queue, 2 hr.
#
# Physics:
#   Mesh: tpv102_200m.msh (h = 200 m fault resolution, ~2.46 M tets)
#   Order: 2 (p2 L2 on tets)
#   Absorbing BCs on outer walls + bottom (attr 5), free surface top (attr 1),
#   embedded fault at Y=0 (attr 3).
#
# Run parameters:
#   tfinal = 3.0 s
#   dt     = 5 ms  →  600 steps
#   ParaView output every 0.1 s (31 snapshots incl. t=0)
#   Fault-surface PVD only; volume PVD suppressed (--no-domain-pv)
#   Station output every 0.1 s at the SCEC probe points
#
# Parallel setup: 8 nodes × 50 ranks = 400 ranks.  Same node/rank layout
# used for BP5 production.  Dev queue max = 2 hr, max = 40 nodes.
#
# Before submitting:
#   1. Build driver on Frontera:
#        cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas
#        make clean && make seas_tpv102_driver -j8
#   2. Ensure the mesh exists:
#        ls tpv102/mesh/tpv102_200m.msh
#      If missing, regenerate with Gmsh:
#        conda activate pythonenv
#        gmsh -3 tpv102/mesh/tpv102_200m.geo -o tpv102/mesh/tpv102_200m.msh
#   3. Submit:
#        sbatch jobs/tpv102/tpv102_200m_3s_dev.sbatch

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
make seas_tpv102_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_200m_3s_dev_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "tfinal=3.0 s, dt=5 ms, output every 0.1 s, fault-surface PVD only"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 2 \
      --bc-mode absorbing \
      --tfinal 3.0 \
      --dt 0.005 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv
