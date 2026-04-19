#!/bin/bash
#SBATCH -J tpv102_200m_p2
#SBATCH -o tpv102_200m_p2_%j.out
#SBATCH -e tpv102_200m_p2_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# SCEC TPV102 dynamic rupture — order-2, Frontera dev queue, 2 hr.
#
# Physics:
#   Mesh: tpv102_200m.msh (h = 200 m fault resolution, ~2.46 M tets)
#   Order: 2 (p2 L2 on tets) — higher-order variant, comparable to
#   Wenqiang Zhang's Mixed-Flux DG 200m O4 submission in SCEC TPV102.
#   Absorbing BCs on outer walls + bottom (attr 5), free surface top (attr 1),
#   embedded fault at Y=0 (attr 3).
#
# Run parameters:
#   tfinal = 3.0 s
#   dt     = CFL-limited (auto, no --dt override); the driver reports the
#            selected dt in the log.  For h=200 m, p=2, cp ~ 5800 m/s with
#            cfl_factor=0.5 this is ≈ 0.17 ms — about 18k steps for 3 s.
#   ParaView fault-surface PVD every 0.1 s; volume PVD suppressed.
#   Station output every 0.1 s at the SCEC probe points.
#
# IMPORTANT: the full 3 s run at p=2 is expected to exceed the dev-queue
# 2 hr wall-time limit on 400 ranks.  This sbatch is intended for:
#   (a) validating the CFL-stable auto-dt path and the fault-surface
#       output at early times (seeing the nucleation pulse form), or
#   (b) re-targeting the "normal" queue (48 hr) for a full 3 s production
#       run — just change -p development -> -p normal and -t 02:00:00
#       -> -t 24:00:00.
# For the p=1 variant that fits in dev-queue 2 hr, see
# tpv102_200m_p1_3s_dev.sbatch (matches PyLith Tet4 200m submission).
#
# Parallel setup: 8 nodes × 50 ranks = 400 ranks.  Same layout as BP5
# production.  Dev queue limits: 2 hr, 40 nodes, 1 job/user.
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
#        sbatch jobs/tpv102/tpv102_200m_p2_3s_dev.sbatch

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

RESULT_DIR="tpv102/results_200m_p2_3s_dev_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "tfinal=3.0 s, dt=auto (CFL-limited at cfl_factor=0.5), order=2"
echo "Output every 0.1 s, fault-surface PVD only"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 2 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 3.0 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv
