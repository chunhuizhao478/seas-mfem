#!/bin/bash
#SBATCH -J bp5_v44b_ip_p2       # Job name
#SBATCH -o bp5_v44b_ip_p2_%j.out
#SBATCH -e bp5_v44b_ip_p2_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (p=2 → 10 DOFs/tet, ~1.9M total DOFs)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v44b: IP p=2, 1000m mesh — CRITICAL MULTI-DOF TEST
#
# v44 implements multi-DOF fault discretization (Tandem approach):
#   face_order = 2, nbf = (2+1)(2+2)/2 = 6 DOFs per fault face
#   L2 traction projection, per-quad-point slip interpolation
#   Per-DOF friction evaluation (automatic via expanded num_nodes)
#
# v42a (same config, face-averaged traction) FAILED: nucleation stalled at V=4e-5 m/s.
# v44b should fix this because multi-DOF slip develops polynomial variation matching
# the DG solution, eliminating the penalty-induced traction bias.
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

mkdir -p bp5/results_v44b_1000m_ip_p2_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v44b_1000m_ip_p2_full \
      --output-prefix bp5_v44b_ip_p2 \
      --checkpoint-interval 1000 \
      --write-every-step
