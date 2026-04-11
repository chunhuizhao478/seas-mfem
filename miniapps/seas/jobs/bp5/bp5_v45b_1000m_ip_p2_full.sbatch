#!/bin/bash
#SBATCH -J bp5_v45b_ip_p2       # Job name
#SBATCH -o bp5_v45b_ip_p2_%j.out
#SBATCH -e bp5_v45b_ip_p2_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (p=2 → 10 DOFs/tet, ~1.9M total DOFs)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v45b: IP p=2, 1000m mesh — CRITICAL MULTI-DOF TEST
#
# v45 implements multi-DOF fault discretization (Phases 1-5, all unit tests pass):
#   - FaceQuadrature: L2 projection with nbf = (p+1)(p+2)/2 = 6 DOFs per face
#   - L2 traction projection in ComputeTraction (GalerkinProject)
#   - Per-quad-point slip interpolation in AssembleSlipContributionIP
#   - Per-DOF friction at 6 GaussLobatto nodes per face
#   - Per-DOF BP5 parameters (a, Dc, tau_pre, V_init) from FaultGeometry
#
# v42a (same config, face-averaged traction) FAILED: nucleation stalled at V=4e-5 m/s.
# v44b (same multi-DOF approach, had 3x penalty bug) was reverted.
# v45b is the clean implementation: no penalty changes, only multi-DOF.
#
# The fix mechanism: multi-DOF slip develops polynomial variation matching the DG
# solution, eliminating the mismatch between constant-per-face slip and the IP
# penalty's p²-scaling that locked the fault at p≥2.
#
# SUCCESS CRITERIA: V_max > 1 m/s (earthquake nucleation and propagation)
# EXPECTED: ~8 events, ~240yr recurrence, matching Tandem (p=2 default)

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

mkdir -p bp5/results_v45b_1000m_ip_p2_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v45b_1000m_ip_p2_full \
      --output-prefix bp5_v45b_ip_p2 \
      --checkpoint-interval 1000 \
      --write-every-step
