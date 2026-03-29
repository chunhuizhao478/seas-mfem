#!/bin/bash
#SBATCH -J bp5_v55_trac_p1           # Job name
#SBATCH -o bp5_v55_trac_p1_%j.out
#SBATCH -e bp5_v55_trac_p1_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v55 traction recovery: all fixes including Tandem-style traction
#
# Changes from previous v55_combined run:
#   - D8: V_vec sign (Tandem convention + GetSlip/SetSlip negation)
#   - D4: Log10 Brent solver
#   - D5: Combined integrator (K-b consistency, 2p+1 quadrature)
#   - D5 fix: Face matrix symmetrization for CG compatibility
#   - NEW: Tandem-style traction recovery using ComputeTractionAtQuadPoints
#     + ProjectTractionToFaultDOFs with sign_flipped double negation
#     - Per-quad-point Jinv (not centroid)
#     - Per-quad-point normal
#     - Per-quad-point penalty (same formula as K)
#     - nl_q-weighted L2 projection with fault basis applied inside integral
#   - FaultBasis.sign_flipped flag matching Tandem's AdapterBase
#
# All tests pass: 587/587, 0 failures.
#
# Key question: does the Tandem-style traction recovery close the 40% V gap?

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

mkdir -p bp5/results_v55_1000m_ip_p1_traction

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
      --output-dir bp5/results_v55_1000m_ip_p1_traction \
      --output-prefix bp5_v55_trac_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
