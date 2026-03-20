#!/bin/bash
#SBATCH -J bp5_v45c_ip_p4_cg    # Job name
#SBATCH -o bp5_v45c_ip_p4_cg_%j.out
#SBATCH -e bp5_v45c_ip_p4_cg_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 64                        # 64 nodes (p=4 → 35 DOFs/tet, ~6.7M total DOFs)
#SBATCH -n 3200                      # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v45c: IP p=4, 1000m mesh — CG+BoomerAMG iterative solver
#
# MUMPS factorization fails at p=4 due to memory.
# CG+AMG alternative: SIPG bilinear form is SPD, so CG is valid.
# BoomerAMG uses SetElasticityOptions for vector-valued problems.
# Doubled nodes from 32→64 for memory headroom.
#
# Compare with GMRES+BlockILU variant to pick faster/more robust solver.
#
# v45 implements multi-DOF fault discretization (Phases 1-5):
#   face_order = 4, nbf = (4+1)(4+2)/2 = 15 DOFs per fault face
#
# SUCCESS CRITERIA: V_max > 1 m/s (earthquake nucleation)

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

mkdir -p bp5/results_v45c_1000m_ip_p4_cg

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method IP \
      --solver cg \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v45c_1000m_ip_p4_cg \
      --output-prefix bp5_v45c_ip_p4_cg \
      --checkpoint-interval 1000 \
      --write-every-step
