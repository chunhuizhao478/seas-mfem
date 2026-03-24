#!/bin/bash
#SBATCH -J bp5_v51d_ip_p1_blr12   # Job name
#SBATCH -o bp5_v51d_ip_p1_blr12_%j.out
#SBATCH -e bp5_v51d_ip_p1_blr12_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v51d: p=1, 1000m, IP, MUMPS-BLR 1e-12 — MISSED TEST FROM v47
#
# v47a ran p=1 1000m with ×3 penalty + BLR 1e-10 → V_nuc decayed, no earthquake.
# The v47 investigation plan called for testing tighter BLR tolerance or exact MUMPS
# at p=1, but this was never done (only p=2 got the tight-tolerance tests).
#
# This run disambiguates:
#   - If V_nuc grows → v47a V-decay was a BLR 1e-10 artifact, NOT penalty stiffness
#   - If V_nuc decays → confirms ×3 penalty + p=1 + 1000m is genuinely too stiff
#
# Configuration:
#   - Order: p=1 (nbf=1, no multi-DOF issues)
#   - Mesh: bp5_tandem.msh (1000m)
#   - DG: IP with ×3 penalty (current default since v47)
#   - Solver: MUMPS-BLR with tolerance 1e-12 (vs v47a's 1e-10)
#   - v50 defaults: CFL-aware dt_init, V-guard ON (factor=100)
#
# See bp5_debug_v51.md Section 20 for full rationale.

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

mkdir -p bp5/results_v51d_1000m_ip_p1_blr12

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
      --output-dir bp5/results_v51d_1000m_ip_p1_blr12 \
      --output-prefix bp5_v51d_ip_p1_blr12 \
      --checkpoint-interval 1000 \
      --write-every-step
