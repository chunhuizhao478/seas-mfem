#!/bin/bash
#SBATCH -J bp5_v47b_ip_p2       # Job name
#SBATCH -o bp5_v47b_ip_p2_%j.out
#SBATCH -e bp5_v47b_ip_p2_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (p=2 → 10 DOFs/tet, ~1.9M total DOFs)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v47b: IP p=2, 1000m mesh — PENALTY ×3 CORRECTION (CRITICAL TEST)
#
# v47 change: IP penalty corrected from 1/3 to full physical value.
#   MFEM's nl_q / Weight() = A/(3V) due to reference element scaling
#   (ref_tet_vol=1/6, ref_tri_area=1/2). Tandem uses physical A/V.
#   Fix: multiply by dim=3 at all 15 penalty locations.
#   See bp5_debug_v47.md for full derivation.
#
# Previous v46 changes (retained):
#   - Slip indexing fix (v46)
#   - Tandem initialization: V_nuc=0.01, delta_tau_factor=0, psi_init="tandem"
#
# This is the KEY test: v46 p=2 matched benchmark in closeup but diverged
# more than p=1 at later times. With correct penalty, p=2 should show
# BETTER long-term convergence than p=1 (higher spatial accuracy + correct
# penalty scaling with p).
#
# Multi-DOF fault discretization:
#   face_order = 2, nbf = (2+1)(2+2)/2 = 6 DOFs per fault face
#   L2 traction projection (GalerkinProject), per-quad-point slip interpolation
#   Per-DOF friction evaluation (automatic via expanded num_nodes)
#
# SUCCESS CRITERIA:
#   1. V_max > 1 m/s (earthquake nucleation)
#   2. p=2 matches Tandem BETTER than p=1 at later events
#   3. Event timing drift < v46b

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

mkdir -p bp5/results_v47b_1000m_ip_p2_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v47b_1000m_ip_p2_full \
      --output-prefix bp5_v47b_ip_p2 \
      --checkpoint-interval 1000 \
      --write-every-step
