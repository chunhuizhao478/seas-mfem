#!/bin/bash
#SBATCH -J bp5_v55_nopsi
#SBATCH -o bp5_v55_nopsi_%j.out
#SBATCH -e bp5_v55_nopsi_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5 v55: IP p=1, 1000m mesh, full run with post-step psi clamp disabled.
#
# Purpose:
#   - Match Tandem's dynamic behavior more closely by leaving accepted RK states
#     untouched after each step.
#   - Test whether removing MFEM's legacy post-step psi clamp changes the
#     early nucleation decay at strk-16dp+10 and the global Vmax trend.
#
# Key dynamic change:
#   - No post-step psi clamp (Tandem-style)
#
# Configuration:
#   - Mesh: bp5_tandem.msh (1000m reference mesh)
#   - Order: p=1
#   - DG: IP
#   - Solver: MUMPS-BLR, tol=1e-12
#   - Runtime: 48h wall clock, 8 nodes / 400 MPI ranks

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v55_1000m_ip_p1_nopsi

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --no-psi-clamp \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v55_1000m_ip_p1_nopsi \
      --output-prefix bp5_v55_ip_p1_nopsi \
      --checkpoint-interval 1000 \
      --write-every-step
