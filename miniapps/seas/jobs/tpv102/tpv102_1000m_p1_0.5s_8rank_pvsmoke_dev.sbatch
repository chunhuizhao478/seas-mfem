#!/bin/bash
#SBATCH -J tpv102_1000m_8r_pv
#SBATCH -o tpv102_1000m_8r_pv_%j.out
#SBATCH -e tpv102_1000m_8r_pv_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# TPV102 ParaView visual smoke test (v2 check-doc recommendation #4).
#
# Purpose: produce a fault-surface PVD at 8 ranks on the 1000 m coarse
# mesh, dense enough in time (--paraview-every 10) that if R-001 is
# incomplete the flagship visual signature — stripes of opposite sign
# in `slip_rate_strike` at partition seams — will be clearly visible
# in ParaView.
#
# Review flow after the run finishes:
#   1. scp ${RESULT_DIR}/ParaView/FaultSurface/fault_surface.pvd to a
#      local machine, open in ParaView.
#   2. Color the fault by `slip_rate_strike` with a divergent colormap
#      (blue-white-red, symmetric around 0).
#   3. At each frame, look at partition seams (use the `mpi_rank`
#      scalar to highlight them).  Seams should show continuous slip
#      rate across them — not alternating +/- stripes.
#
# If stripes ARE visible at seams: R-001 is incomplete under the current
# MFEM face-normal convention, the R-101 runtime diagnostic should have
# caught it (check the job .out for "R-101 shared-fault DOFData
# consistency FAILED" → MFEM_ABORT).  If the diagnostic did NOT catch
# it and stripes are present, the diagnostic itself has a gap — file
# that as a v3 regression.
#
# Expected wall time: ~15 min.
#
# Usage:
#   sbatch jobs/tpv102/tpv102_1000m_p1_0.5s_8rank_pvsmoke_dev.sbatch

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

RESULT_DIR="tpv102/results_1000m_p1_0.5s_8r_pv_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "tfinal=0.5 s, dt=auto, 8 ranks, paraview every 10 steps"
echo "After run: open ${RESULT_DIR}/ParaView/FaultSurface/fault_surface.pvd"
echo "and check slip_rate_strike at partition seams for opposite-sign stripes."

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_1000m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 0.5 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.05 \
      --paraview \
      --paraview-every 10 \
      --no-domain-pv \
      --debug-qnorm
