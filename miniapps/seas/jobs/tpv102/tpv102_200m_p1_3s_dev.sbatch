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
#   tfinal = 3.0 s (see WALL-TIME NOTE below — dev queue 2 hr will
#                    not complete the full run; use normal queue
#                    48 hr for production).
#   dt     = auto (CFL-limited) via the driver's ComputeMaxDt using
#            cfl_factor=0.5.  For h=200 m, p=1, cp ≈ 5800 m/s on tets
#            this is ~0.285 ms.  Do NOT override with a larger dt:
#            explicit RK4 for our Godunov-flux DG has a linear stability
#            envelope of only ~2× the CFL-1 limit on hyperbolic systems,
#            i.e. dt ≲ 0.6 ms; anything larger (e.g. the 5 ms value
#            that looked reasonable for implicit PyLith integrators)
#            causes the wave field to amplify ~10×/step once dynamics
#            engage — blow-up was observed at step ~30 in job 7664250.
#   ParaView fault-surface PVD every 0.1 s; volume PVD suppressed.
#   Station output every 0.1 s at the SCEC probe points.
#
# WALL-TIME NOTE: ~10 500 steps × ~1 s/step on 400 ranks ≈ 3 hr.  The
# 2 hr dev queue will produce results up to ~t ≈ 2 s.  For the full
# 3 s run, change `-p development` → `-p normal` and `-t 02:00:00` →
# `-t 12:00:00` (normal queue max is 48 hr).
#
# Parallel setup: 8 nodes × 50 ranks = 400 ranks.
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
echo "tfinal=3.0 s, dt=auto (CFL-limited at cfl_factor=0.5)"
echo "order=1 (matches PyLith Tet4 200m), output every 0.1 s, fault-surface PVD only"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 3.0 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv
