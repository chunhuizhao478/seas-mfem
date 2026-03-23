#!/bin/bash
#SBATCH -J bp5_v50a_2500m_p4      # Job name
#SBATCH -o bp5_v50a_2500m_p4_%j.out
#SBATCH -e bp5_v50a_2500m_p4_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # 200 MPI ranks
#SBATCH -t 04:00:00                  # 4 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v50a: 2500m mesh, p=4, IP method
#
# First p=4 run with order-aware CFL dt.
#   c_N_1 = 4*(4+2)/3 = 8.0
#   beta  = 4.0 * 8.0 / (8/3) = 12.0
#   dt_CFL = 2*eta*h_min / (12*mu)  ≈ smaller than p=2
#
# Tandem reference: bp5qd_tandem_p4 (4000m mesh, p=4)
# Our mesh is finer (2500m) to compensate for lower order accuracy vs Tandem.
#
# nbf_per_face = (4+1)(4+2)/2 = 15 DOFs per fault face.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# BUILD job — rebuild for v50+ order-aware CFL fix
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v50a_2500m_p4

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v50a_2500m_p4 \
      --output-prefix bp5_v50a_2500m_p4 \
      --write-every-step
