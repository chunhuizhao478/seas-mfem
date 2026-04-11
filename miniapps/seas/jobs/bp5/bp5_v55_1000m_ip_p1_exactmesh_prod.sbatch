#!/bin/bash
#SBATCH -J bp5_v55_exct_prod         # Job name
#SBATCH -o bp5_v55_exct_prod_%j.out
#SBATCH -e bp5_v55_exct_prod_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v55: Tandem's EXACT mesh (62874 tets), 48hr production run
#
# The 2hr dev run confirmed: exact mesh produces different V trajectory
# from the old mesh (63451 tets). V is ~5% higher on the exact mesh
# at early steps, closer to Tandem. Full production run needed to
# compare earthquake timing and interseismic traction accumulation.

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

mkdir -p bp5/results_v55_1000m_ip_p1_exactmesh_prod

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v55_1000m_ip_p1_exactmesh_prod \
      --output-prefix bp5_v55_exct_prod \
      --checkpoint-interval 1000 \
      --write-every-step
