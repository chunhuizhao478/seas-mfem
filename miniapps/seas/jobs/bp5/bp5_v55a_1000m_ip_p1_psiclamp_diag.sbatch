#!/bin/bash
#SBATCH -J bp5_v55a_psi
#SBATCH -o bp5_v55a_psi_%j.out
#SBATCH -e bp5_v55a_psi_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5 v55a: IP p=1, 1000m mesh, legacy psi-clamp A/B run with diagnostics.
#
# Purpose:
#   - Re-enable the legacy post-step psi clamp explicitly.
#   - Record whether and when clamp hits occur, so the unclamped v55 run can be
#     compared against the previous MFEM behavior on the same production setup.
#
# New dynamic flags:
#   --psi-clamp
#   --diag-psi-clamp

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

mkdir -p bp5/results_v55a_1000m_ip_p1_psiclamp_diag

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --psi-clamp \
      --diag-psi-clamp \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v55a_1000m_ip_p1_psiclamp_diag \
      --output-prefix bp5_v55a_ip_p1_psiclamp_diag \
      --checkpoint-interval 1000 \
      --write-every-step
