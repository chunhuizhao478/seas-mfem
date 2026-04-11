#!/bin/bash
#SBATCH -J bp5_v55_d8_p1             # Job name
#SBATCH -o bp5_v55_d8_p1_%j.out
#SBATCH -e bp5_v55_d8_p1_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v55 D8-only: V_vec sign convention matched to Tandem
#
# D8 change only (D4 log10 solver reverted):
#   V_vec negated in SolveSlipRateVectorPsi to match Tandem's
#   -(V/|tau|)*tau convention (DieterichRuinaBase.h:174).
#   Below-Wf hardcoded rate also negated to -Vp.
#   Scalar Brent solver unchanged (linear space [0, tau/eta]).
#
# The sign chain V→S→EmbedSlip→sign*delta_u is self-consistent:
# MFEM's sign factor in slip/traction assembly handles face orientation.
# Unit tests: 105/105 cross_verify, 76/76 bp5_fault_operator, 32/32 friction,
# 23/23 state_evolution — all pass.
#
# Purpose: isolate D8 effect. The D8+D4 combined run crashed during
# coseismic (segfault in GetSharedFaceTransformations at V~200 m/s).
# This D8-only run determines whether D8 alone causes the nucleation
# change or if D4 was the trigger.
#
# Configuration identical to v53 except for the D8 sign change.

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

mkdir -p bp5/results_v55_1000m_ip_p1_d8only

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v55_1000m_ip_p1_d8only \
      --output-prefix bp5_v55_d8_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
