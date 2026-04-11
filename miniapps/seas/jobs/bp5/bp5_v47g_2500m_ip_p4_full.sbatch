#!/bin/bash
#SBATCH -J bp5_v47g_ip_p4       # Job name
#SBATCH -o bp5_v47g_ip_p4_%j.out
#SBATCH -e bp5_v47g_ip_p4_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (13.5k tets × 105 DOFs/tet → 1.4M DOFs)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v47g: IP p=4, 2500m mesh — MATCH TANDEM uphoff.2 WITH CORRECT ×3 PENALTY
#
# This matches Tandem's SCEC benchmark submission "uphoff.2" (dx=2.5km, P4).
# Tandem's penalty at this config: penalty_ip ≈ 9.4 GPa/m (LARGER than our
# p=2 h=1000m value of 7.83 GPa/m). Tandem nucleates fine because p=4 gives
# O(h^5) accuracy — the DG displacement field is accurate enough that the
# penalty correction residual is negligible.
#
# v47 changes vs v42b:
#   - IP penalty ×3 correction (ref element scaling fix, 15 locations)
#   - Slip indexing fix (v46: fi*nbf+kk instead of fi)
#   - Tandem initialization: V_nuc=0.01, delta_tau_factor=0, psi_init="tandem"
#   - dt_init fix: uses max(V_init, V_nuc) → dt ≈ 0.13s
#
# v42b had c_N_1 fix but NOT the ×3 penalty fix (penalty was 1/3 of physical).
# v47g has the CORRECT penalty. If p=4 at h=2500m matches Tandem, then:
#   1. The ×3 penalty fix is validated at the resolution where it matters
#   2. The p=1/p=2 h=1000m failures are confirmed as resolution effects
#   3. The path forward is higher p (not penalty tuning)
#
# Mesh: bp5_tandem_2500m.msh (13,503 tets)
# DOFs: ~1.4M (3 components × 35 DOFs/tet × 13.5k tets)
# c_N_1: p*(p+D-1)/D = 4*(4+3-1)/3 = 8.0
# penalty_ip ≈ (D+1)*c_N_1*(D*A/V)*(c1²/c0) / 2 ≈ 9.4 GPa/m

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

mkdir -p bp5/results_v47g_2500m_ip_p4_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v47g_2500m_ip_p4_full \
      --output-prefix bp5_v47g_ip_p4 \
      --checkpoint-interval 1000 \
      --write-every-step
