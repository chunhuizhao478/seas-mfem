#!/bin/bash
#SBATCH -J bp5_v57_pts_p1_mumps
#SBATCH -o bp5_v57_pts_p1_mumps_%j.out
#SBATCH -e bp5_v57_pts_p1_mumps_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5 v57: 1000 m, IP, p=1, PETSc TS with MUMPS BLR
#
# Purpose:
# - build MFEM/SEAS with PETSc and MFEM_USE_MUMPS enabled
# - use PETSc TS for time stepping
# - keep the BP5 elasticity solve on MUMPS BLR
#
# IMPORTANT:
# - on Frontera, the MUMPS module is coupled to PETSc 3.15
# - this job therefore builds against petsc/3.15 intentionally

export LC_ALL=C
export LANG=C

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

export PETSC_MODULE=petsc/3.15
export USE_MUMPS=1
export JOBS=8

# Reconfigure from a clean tree so we do not mix stale PETSc 3.23 objects
# with the PETSc 3.15 + MUMPS build requested by this job.
make distclean || true

bash build_frontera.sh || { echo "BUILD FAILED"; exit 1; }

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || {
    echo "ERROR: MFEM was not built with PETSc enabled"
    exit 1
}
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || {
    echo "ERROR: MFEM was not built with MUMPS enabled"
    exit 1
}

cd miniapps/seas
mkdir -p bp5/results_v57_1000m_ip_p1_petsc_ts_mumps

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v57_1000m_ip_p1_petsc_ts_mumps \
      --output-prefix bp5_v57_petsc_ts_p1_mumps \
      --checkpoint-interval 0 \
      --write-every-step \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg
