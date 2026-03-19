#!/bin/bash
#SBATCH -J bp5_v46b_ip_p2       # Job name
#SBATCH -o bp5_v46b_ip_p2_%j.out
#SBATCH -e bp5_v46b_ip_p2_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (p=2 → 10 DOFs/tet, ~1.9M total DOFs)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v46b: IP p=2, 1000m mesh — CRITICAL MULTI-DOF TEST
#
# v46 fixes the ComputeTraction penalty slip indexing bug that caused v45b
# to blow up at p=2. The penalty correction now correctly:
#   1. Reads per-DOF slip via slip_bc(2 * (fi * nbf + kk))
#   2. Interpolates slip to quadrature points (InterpolateToQuadPoints)
#   3. Uses per-quad-point slip in the penalty: jump = (u1-u2) - sign*slip_q
#
# Previously (v45b), the penalty used slip_bc(2*fi) — the face index instead
# of the DOF index — reading wrong faces' slip. This created huge spurious
# traction forces that drove V_max from 0.049 to 196 m/s in 56 steps (OOM).
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
