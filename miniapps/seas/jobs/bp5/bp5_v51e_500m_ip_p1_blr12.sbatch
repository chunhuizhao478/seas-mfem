#!/bin/bash
#SBATCH -J bp5_v51e_500m_ip_p1   # Job name
#SBATCH -o bp5_v51e_500m_ip_p1_%j.out
#SBATCH -e bp5_v51e_500m_ip_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 32                        # 32 nodes (234k tets, ~4x 1000m)
#SBATCH -n 1600                      # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v51e: p=1, 500m mesh, IP, MUMPS-BLR 1e-12, ×3 penalty
#
# v51d confirmed p=1 on 1000m with ×3 penalty cannot nucleate (k_elastic > k_crit).
# This run tests whether halving h (500m) restores nucleation:
#   - DOF spacing = 500m (same as p=2 on 1000m, which nucleates)
#   - DG jump residual O(h²) = 4× smaller than 1000m
#   - k_elastic should drop below k_crit = 5.77 MPa/m
#
# Previous v41b ran 500m p=1 with 1/3 penalty (pre-v47) → nucleated OK.
# This is the first test of 500m p=1 with CORRECT ×3 penalty.
#
# Configuration:
#   - Order: p=1 (nbf=1, no multi-DOF issues)
#   - Mesh: bp5_tandem_500m.msh (233,783 tets)
#   - DG: IP with ×3 penalty (current default since v47)
#   - Solver: MUMPS-BLR with tolerance 1e-12
#   - v50 defaults: CFL-aware dt_init, V-guard ON (factor=100)
#   - v51 change: AgingLawPsi exp cap removed (matches Tandem/benchmark)
#
# Success: V_nuc grows and earthquake nucleates → p=1 viable at 500m
# Failure: V_nuc decays → p=1 needs even finer mesh or must use p≥2
#
# See bp5_debug_v51.md Section 20 (v51d result) and Section 24 (benchmark audit).

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

mkdir -p bp5/results_v51e_500m_ip_p1_blr12

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_500m.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v51e_500m_ip_p1_blr12 \
      --output-prefix bp5_v51e_500m_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
