#!/bin/bash
#SBATCH -J bp5_v46a_ip_p1       # Job name
#SBATCH -o bp5_v46a_ip_p1_%j.out
#SBATCH -e bp5_v46a_ip_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v46a: IP p=1, 1000m mesh — REGRESSION TEST (Tandem initialization)
#
# v46 changes from previous versions:
#   1. Slip indexing fix in ComputeTraction IP penalty (v45→v46):
#      At p>=2, penalty used slip_bc(2*face_index) instead of
#      slip_bc(2*(face_index*nbf + dof_index)), reading wrong faces' slip.
#      At p=1 (nbf=1): face_index = DOF_index, so no effect on p=1.
#
#   2. Tandem initialization defaults adopted (v46 Phase 2):
#      - V_nuc = 0.01 m/s (was 0.03, Tandem bp5.lua:72)
#      - delta_tau_factor = 0.0 (was 1.0, Tandem has no delta_tau overstress)
#      - psi_init_mode = "tandem" (equilibrium psi from stress, not fixed 0.8072)
#      - With delta_tau_factor=0, tau_pre is exact equilibrium:
#        tau_pre = sigma_n * a * asinh(V_init/(2*V0) * exp(psi_ss/a))
#        Initial V = V_init everywhere (no overstress-induced jump)
#      - Math proof: psi_init = f0 + b*ln(V0/Vp) = 0.8072 regardless of mode
#        (tau_pre V_nuc terms cancel in InitialStatePsi formula)
#
# EXPECTED: ~8 events, ~250yr recurrence, matching v45a/v41a reference.
# Small initial transient differences expected (V_nuc=0.01 vs 0.03, no overstress)
# but long-term dynamics should converge to same attractor.

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

mkdir -p bp5/results_v46a_1000m_ip_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v46a_1000m_ip_p1_full \
      --output-prefix bp5_v46a_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
