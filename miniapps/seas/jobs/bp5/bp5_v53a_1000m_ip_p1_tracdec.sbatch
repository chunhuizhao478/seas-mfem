#!/bin/bash
#SBATCH -J bp5_v53a_tracdec       # Job name
#SBATCH -o bp5_v53a_tracdec_%j.out
#SBATCH -e bp5_v53a_tracdec_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5 v53a: IP p=1, 1000m mesh, full run with station traction decomposition
#
# Purpose:
#   - Keep the same stable v53 production configuration
#   - Add station-level traction decomposition output:
#       tau_stress + tau_corr = tau_total
#   - Diagnose whether the remaining Tandem mismatch at strk-16dp+10 comes
#     mainly from the stress part, the correction part, or their balance
#
# New diagnostic flag:
#   --diag-station-traction-decomp
#
# Extra output files:
#   <prefix>_tracdec_<station>.txt
#
# Columns:
#   time(s),
#   tau_stress_strike(MPa), tau_stress_dip(MPa),
#   tau_corr_strike(MPa),   tau_corr_dip(MPa),
#   tau_total_strike(MPa),  tau_total_dip(MPa)

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

mkdir -p bp5/results_v53a_1000m_ip_p1_tracdec

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --diag-station-traction-decomp \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v53a_1000m_ip_p1_tracdec \
      --output-prefix bp5_v53a_ip_p1_tracdec \
      --checkpoint-interval 1000 \
      --write-every-step
