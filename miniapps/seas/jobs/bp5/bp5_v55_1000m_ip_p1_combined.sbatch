#!/bin/bash
#SBATCH -J bp5_v55_comb_p1           # Job name
#SBATCH -o bp5_v55_comb_p1_%j.out
#SBATCH -e bp5_v55_comb_p1_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v55 D8+D4+D5: Combined IP integrator matching Tandem structure
#
# Changes from baseline:
#   D8: V_vec negated + GetSlip/SetSlip negation (correct, verified)
#   D4: Log10 Brent solver (matching Tandem's search space)
#   D5: NEW - Combined DGElasticityIPCombinedIntegrator replaces the split
#       DGElasticityIntegrator + DGElasticityIPPenaltyIntegrator for K,
#       and the manual AssembleSlipContributionIP loop for b.
#       Both K and b now use the same traction operator formula, quadrature
#       (2p+1), and penalty computation, guaranteeing K-b consistency.
#       Mirrors Tandem's assembleSurface + rhsFacet kernel structure.
#
# All tests pass: 97+76+32+23+352 = 580 tests, 0 failures.

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

mkdir -p bp5/results_v55_1000m_ip_p1_combined

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
      --output-dir bp5/results_v55_1000m_ip_p1_combined \
      --output-prefix bp5_v55_comb_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
