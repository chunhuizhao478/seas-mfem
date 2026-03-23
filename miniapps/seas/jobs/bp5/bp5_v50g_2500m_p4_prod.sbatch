#!/bin/bash
#SBATCH -J bp5_v50g_p4_prod       # Job name
#SBATCH -o bp5_v50g_p4_prod_%j.out
#SBATCH -e bp5_v50g_p4_prod_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks — production scale
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v50g PRODUCTION: 2500m mesh, p=4, IP, ClosedUniform face nodes
#
# ROOT CAUSE FIX: GaussLobatto nodes on triangles have cond(M)=2901 at p=4.
# ClosedUniform has cond(M)=58 — 50x better conditioning.
# Phase 1 test (v50g, 116 steps) confirmed V growing at p=4.
#
# This is the full production run to capture multiple earthquake cycles.
# Compare with Tandem p4 reference (bp5qd_tandem_p4, 4000m mesh).

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# BUILD job
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v50g_p4_prod

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
      --output-dir bp5/results_v50g_p4_prod \
      --output-prefix bp5_v50g_p4_prod \
      --face-basis-type equi \
      --checkpoint-interval 5000 \
      --write-every-step
