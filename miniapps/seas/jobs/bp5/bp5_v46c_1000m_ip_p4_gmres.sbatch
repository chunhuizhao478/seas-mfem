#!/bin/bash
#SBATCH -J bp5_v46c_ip_p4_gmres # Job name
#SBATCH -o bp5_v46c_ip_p4_gmres_%j.out
#SBATCH -e bp5_v46c_ip_p4_gmres_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 64                        # 64 nodes (p=4 → 35 DOFs/tet, ~6.7M total DOFs)
#SBATCH -n 3200                      # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v46c: IP p=4, 1000m mesh — GMRES+BlockILU iterative solver
#
# MUMPS factorization fails at p=4 due to memory (35 DOFs/tet on 1000m mesh).
# Switch to GMRES with BlockILU preconditioner (element-sized blocks,
# natural for DG). Doubled nodes from 32→64 for memory headroom.
#
# v46 implements multi-DOF fault discretization (Phases 1-5):
#   face_order = 4, nbf = (4+1)(4+2)/2 = 15 DOFs per fault face
#   L2 traction projection, per-quad-point slip interpolation
#   Per-DOF friction at 15 GaussLobatto nodes per face
#
# SUCCESS CRITERIA: V_max > 1 m/s (earthquake nucleation)
# EXPECTED: ~8 events, dip slip convergence toward Tandem (p=4 default)

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

mkdir -p bp5/results_v46c_1000m_ip_p4_gmres

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method IP \
      --solver gmres \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v46c_1000m_ip_p4_gmres \
      --output-prefix bp5_v46c_ip_p4_gmres \
      --checkpoint-interval 1000 \
      --write-every-step
