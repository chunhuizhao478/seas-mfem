#!/bin/bash
#SBATCH -J bp5_v39_ip_p1      # Job name
#SBATCH -o bp5_v39_ip_p1_%j.out
#SBATCH -e bp5_v39_ip_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v39: IP method — p=1, 1000m mesh, 1800 yr
#
# Fixes since v38a (which blew up at step 43):
#   v38b: Shared-face Dirichlet loading (parallel bug)
#   v38c: Shared-face IP traction sign (missing sign * in correction)
#   v39:  Interior-face IP traction sign (same missing sign * bug)
#
# The interior-face fix is the critical one: ~half of all interior fault
# faces had sign=-1, making the penalty correction destabilizing (positive
# feedback instead of negative). This drove exponential blowup.
#
# All three bugs caused: correction = penalty * jump (raw DG convention)
# instead of:            correction = penalty * sign * jump (canonical)
#
# Previous: v38a (IP) blew up at 5.7e-9 yr
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

mkdir -p bp5/results_v39_1000m_ip_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v39_1000m_ip_p1_full \
      --output-prefix bp5_v39_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
