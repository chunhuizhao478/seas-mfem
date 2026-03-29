#!/bin/bash
#SBATCH -J bp5_v53_p1_exact          # Job name
#SBATCH -o bp5_v53_p1_exact_%j.out
#SBATCH -e bp5_v53_p1_exact_%j.err
#SBATCH -p normal                    # Queue (48hr max)
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v53: p=1, 1000m, IP, Exact MUMPS — MFEM standard DG integrators
#
# Key changes from v51d:
#   1. Bilinear form uses MFEM's DGElasticityIntegrator(alpha, kappa)
#      instead of split custom penalty integrator
#   2. Boundary Dirichlet RHS uses MFEM's DGElasticityDirichletLFIntegrator
#   3. Solver: exact MUMPS (matches Tandem's Frontera config)
#
# This is the first production run with:
#   - Standard MFEM DG code (no custom penalty integrator)
#   - Same solver as Tandem (exact MUMPS, no BLR)
#   - Same mesh/order as Tandem p=1 reference
#
# Compare against:
#   - Tandem p=1 1000m (results_tandem_1000m_p1/)
#   - MFEM v51d p=1 BLR (results_v51/bp5_v51d_ip_p1_blr12_*)

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

mkdir -p bp5/results_v53_1000m_ip_p1_exact

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v53_1000m_ip_p1_exact \
      --output-prefix bp5_v53_p1_exact \
      --checkpoint-interval 1000 \
      --write-every-step \
      --diag-rhs-z
