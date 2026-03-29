#!/bin/bash
#SBATCH -J bp5_v53b_jumpres       # Job name
#SBATCH -o bp5_v53b_jumpres_%j.out
#SBATCH -e bp5_v53b_jumpres_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5 v53b: IP p=1, 1000m mesh, full run with station traction decomposition
# plus jump-residual diagnostics
#
# Purpose:
#   - Keep the same stable v53 production configuration
#   - Retain station-level traction decomposition output
#   - Add station-level fault jump residual output:
#       R = [[u]] - delta
#   - Determine whether the oversized correction near strk-16dp+10 is driven
#     by a large jump residual itself or by penalty amplification
#
# New diagnostic flags:
#   --diag-station-traction-decomp
#   --diag-station-jump-residual
#
# Extra output files:
#   <prefix>_tracdec_<station>.txt
#   <prefix>_jumpres_<station>.txt
#
# Jump residual columns:
#   time(s),
#   jump_res_strike(m), jump_res_dip(m), jump_res_mag(m)

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

mkdir -p bp5/results_v53b_1000m_ip_p1_jumpres

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --diag-station-traction-decomp \
      --diag-station-jump-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v53b_1000m_ip_p1_jumpres \
      --output-prefix bp5_v53b_ip_p1_jumpres \
      --checkpoint-interval 1000 \
      --write-every-step
