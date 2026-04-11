#!/bin/bash
#SBATCH -J bp5_v50d_p2_ip         # Job name
#SBATCH -o bp5_v50d_p2_ip_%j.out
#SBATCH -e bp5_v50d_p2_ip_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # 200 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v50d: 2500m mesh, p=2, IP method (SAME mesh as v50a p=4)
#
# DIAGNOSTIC: Baseline comparison on the same 2500m mesh.
# If p=2 nucleates correctly on 2500m mesh, the issue is p=4-specific.
# This isolates mesh effects from order effects.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# NO BUILD — v50c builds. Wait for binary.
for i in $(seq 1 60); do
   [ -f seas_bp5_full ] && break
   echo "Waiting for v50c build... ($i/60)"
   sleep 10
done
[ ! -f seas_bp5_full ] && { echo "ERROR: seas_bp5_full not found after 10min"; exit 1; }

mkdir -p bp5/results_v50d_p2_ip

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v50d_p2_ip \
      --output-prefix bp5_v50d_p2_ip \
      --write-every-step
