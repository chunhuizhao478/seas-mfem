#!/bin/bash
#SBATCH -J bp5_v50b_4000m_p6      # Job name
#SBATCH -o bp5_v50b_4000m_p6_%j.out
#SBATCH -e bp5_v50b_4000m_p6_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # 200 MPI ranks
#SBATCH -t 04:00:00                  # 4 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v50b: 4000m mesh, p=6, IP method
#
# Direct comparison with Tandem p6 reference (same 4000m mesh).
#   c_N_1 = 6*(6+2)/3 = 16.0
#   beta  = 4.0 * 16.0 / (8/3) = 24.0
#   dt_CFL = 2*eta*h_min / (24*mu)
#
# Tandem reference: bp5qd_tandem_p6 (4000m mesh, p=6)
# SAME mesh resolution → closest possible comparison.
#
# nbf_per_face = (6+1)(6+2)/2 = 28 DOFs per fault face.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# NO BUILD — v50a builds. Wait for binary to exist.
for i in $(seq 1 60); do
   [ -f seas_bp5_full ] && break
   echo "Waiting for v50a build... ($i/60)"
   sleep 10
done
[ ! -f seas_bp5_full ] && { echo "ERROR: seas_bp5_full not found after 10min"; exit 1; }

mkdir -p bp5/results_v50b_4000m_p6

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_4000m.msh \
      --mesh-scale 1000 \
      --order 6 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v50b_4000m_p6 \
      --output-prefix bp5_v50b_4000m_p6 \
      --write-every-step
