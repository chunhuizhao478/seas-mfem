#!/bin/bash
#SBATCH -J bp5_v40_ip_p1      # Job name
#SBATCH -o bp5_v40_ip_p1_%j.out
#SBATCH -e bp5_v40_ip_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v40: IP method — p=1, 1000m mesh, 1800 yr
#
# Fix: IP penalty correction sign in ComputeTraction.
#
# The traction formula was T = {σ·n̂_fault} - η*pj, but with
# n̂_fault = (0,-1,0), the standard DG formula becomes:
#   T = {σ·n̂_fault} + η*pj  (sign flips when normal flips)
#
# The code used correction = +penalty*sign*jump, then T = stress - correction.
# Fixed to correction = -penalty*sign*jump, so T = stress - (-η*pj) = stress + η*pj.
#
# This affected BOTH interior and shared face IP traction paths.
# BR2 was unaffected (its lifting-based correction naturally includes n̂_fault).
#
# Previous: v38a/v39 (IP) blew up at step 43, V_max=216 m/s, tau=1 GPa
#           v37 (BR2) = 295 yr, Tandem (IP) = 240 yr

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

mkdir -p bp5/results_v40_1000m_ip_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v40_1000m_ip_p1_full \
      --output-prefix bp5_v40_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
