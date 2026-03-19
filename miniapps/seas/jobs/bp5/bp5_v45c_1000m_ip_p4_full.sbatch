#!/bin/bash
#SBATCH -J bp5_v45c_ip_p4       # Job name
#SBATCH -o bp5_v45c_ip_p4_%j.out
#SBATCH -e bp5_v45c_ip_p4_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 32                        # 32 nodes (p=4 → 35 DOFs/tet, ~6.7M total DOFs)
#SBATCH -n 1600                      # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v45c: IP p=4, 1000m mesh — HIGH-ORDER CONVERGENCE TEST
#
# v45 implements multi-DOF fault discretization (Phases 1-5, all unit tests pass):
#   - FaceQuadrature: L2 projection with nbf = (p+1)(p+2)/2 = 15 DOFs per face
#   - L2 traction projection in ComputeTraction (GalerkinProject)
#   - Per-quad-point slip interpolation in AssembleSlipContributionIP
#   - Per-DOF friction at 15 GaussLobatto nodes per face
#   - Per-DOF BP5 parameters (a, Dc, tau_pre, V_init) from FaultGeometry
#
# v42b (same config, face-averaged traction) FAILED: fault locked at V→1e-10 m/s.
# v45c should fix this — multi-DOF eliminates the p²-scaling penalty bias.
#
# SUCCESS CRITERIA: V_max > 1 m/s (earthquake nucleation)
# EXPECTED: ~8 events, strike deficit < 5% (vs ~20% at p=1)
#
# Note: p=4 on 1000m mesh has 35 DOFs/tet and 15 fault DOFs/face.
# Memory usage ~4x that of p=2. MUMPS BLR should handle it.

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

mkdir -p bp5/results_v45c_1000m_ip_p4_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v45c_1000m_ip_p4_full \
      --output-prefix bp5_v45c_ip_p4 \
      --checkpoint-interval 1000 \
      --write-every-step
