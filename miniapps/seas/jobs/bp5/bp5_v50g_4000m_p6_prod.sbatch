#!/bin/bash
#SBATCH -J bp5_v50g_p6_prod       # Job name
#SBATCH -o bp5_v50g_p6_prod_%j.out
#SBATCH -e bp5_v50g_p6_prod_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks — production scale
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v50g PRODUCTION: 4000m mesh, p=6, IP, ClosedUniform face nodes
#
# Direct comparison with Tandem p6 reference (same 4000m mesh, same p=6).
# ClosedUniform nodes fix the GaussLobatto conditioning disaster at p>=3.
# At p=6: GL cond(M)=479658 vs Equi cond(M)=215 — 2232x improvement.
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

# NO BUILD — v50g p4 prod builds. Wait for binary.
for i in $(seq 1 60); do
   [ -f seas_bp5_full ] && break
   echo "Waiting for build... ($i/60)"
   sleep 10
done
[ ! -f seas_bp5_full ] && { echo "ERROR: seas_bp5_full not found after 10min"; exit 1; }

mkdir -p bp5/results_v50g_p6_prod

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
      --output-dir bp5/results_v50g_p6_prod \
      --output-prefix bp5_v50g_p6_prod \
      --face-basis-type equi \
      --checkpoint-interval 5000 \
      --write-every-step
