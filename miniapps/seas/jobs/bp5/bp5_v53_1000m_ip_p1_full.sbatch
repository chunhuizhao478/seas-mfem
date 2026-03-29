#!/bin/bash
#SBATCH -J bp5_v53_ip_p1          # Job name
#SBATCH -o bp5_v53_ip_p1_%j.out
#SBATCH -e bp5_v53_ip_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v53: IP p=1, 1000m mesh, 1800 yr — Tandem-style 3-node fault sampling
#
# Root fix in v53:
#   - Current MFEM p=1 IP no longer uses one face-averaged fault DOF per
#     triangular face.
#   - Instead, it now uses the same nodal fault-space structure as Tandem p=1:
#       3 fault DOFs per triangular face
#   - BP5 parameters (a, Dc/L, V_init, tau_pre) are therefore sampled at the
#     three face nodes rather than at a single centroid value.
#
# This is the first full 1800-year production run after that discretization fix.
# The purpose is to test whether p=1 near-fault 1000m now recovers Tandem-like
# nucleation and rupture evolution.
#
# Configuration:
#   - Mesh: bp5_tandem.msh (1000m reference mesh)
#   - Order: p=1
#   - DG: IP
#   - Solver: MUMPS-BLR, tol=1e-12
#   - Initialization: current Tandem-style defaults in code
#   - Runtime: 48h wall clock, 8 nodes / 400 MPI ranks
#
# Expected validation target:
#   - Compare against Tandem p=1 near-fault 1000m station outputs
#   - Focus on initial nucleation timing at strk-24dp+10, strk-16dp+10,
#     and strk+00dp+10

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

mkdir -p bp5/results_v53_1000m_ip_p1_full

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
      --output-dir bp5/results_v53_1000m_ip_p1_full \
      --output-prefix bp5_v53_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
