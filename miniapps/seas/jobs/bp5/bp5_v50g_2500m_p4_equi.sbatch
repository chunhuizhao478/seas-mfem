#!/bin/bash
#SBATCH -J bp5_v50g_p4_equi       # Job name
#SBATCH -o bp5_v50g_p4_equi_%j.out
#SBATCH -e bp5_v50g_p4_equi_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # 200 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v50g: 2500m mesh, p=4, IP, ClosedUniform face nodes
#
# ROOT CAUSE TEST: GaussLobatto nodes on triangles have cond(M)=2901 at p=4
# (Lebesgue constant = 46). ClosedUniform has cond(M)=58 (Lebesgue = 3.5).
# The ill-conditioned mass matrix inversion in GalerkinProject amplifies
# traction projection errors by 46x, flipping DOF signs and killing nucleation.
#
# Expected: If node distribution is the root cause, V should GROW at p=4.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# BUILD job
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v50g_p4_equi

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v50g_p4_equi \
      --output-prefix bp5_v50g_p4_equi \
      --face-basis-type equi \
      --write-every-step
