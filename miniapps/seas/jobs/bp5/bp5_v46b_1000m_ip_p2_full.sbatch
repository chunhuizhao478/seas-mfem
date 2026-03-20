#!/bin/bash
#SBATCH -J bp5_v46b_ip_p2       # Job name
#SBATCH -o bp5_v46b_ip_p2_%j.out
#SBATCH -e bp5_v46b_ip_p2_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (p=2 → 10 DOFs/tet, ~1.9M total DOFs)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v46b: IP p=2, 1000m mesh — CRITICAL MULTI-DOF TEST (Tandem initialization)
#
# v46 changes from previous versions:
#   1. Slip indexing fix in ComputeTraction IP penalty (v45→v46):
#      Previously (v45b), penalty used slip_bc(2*fi) — face index instead of
#      DOF index — reading wrong faces' slip. Created huge spurious traction
#      forces (V_max 0.049 → 196 m/s in 56 steps, OOM). Now correctly:
#        - Reads per-DOF slip via slip_bc(2 * (fi * nbf + kk))
#        - Interpolates slip to quad points (InterpolateToQuadPoints)
#        - Uses per-quad-point slip in penalty: jump = (u1-u2) - sign*slip_q
#
#   2. Tandem initialization defaults adopted (v46 Phase 2):
#      - V_nuc = 0.01 m/s (was 0.03, Tandem bp5.lua:72)
#      - delta_tau_factor = 0.0 (was 1.0, Tandem has no delta_tau overstress)
#      - psi_init_mode = "tandem" (equilibrium psi from stress, not fixed 0.8072)
#      - No initial overstress → cleaner nucleation, matches Tandem exactly
#
# Multi-DOF fault discretization:
#   face_order = 2, nbf = (2+1)(2+2)/2 = 6 DOFs per fault face
#   L2 traction projection (GalerkinProject), per-quad-point slip interpolation
#   Per-DOF friction evaluation (automatic via expanded num_nodes)
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

mkdir -p bp5/results_v46b_1000m_ip_p2_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v46b_1000m_ip_p2_full \
      --output-prefix bp5_v46b_ip_p2 \
      --checkpoint-interval 1000 \
      --write-every-step
