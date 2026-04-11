#!/bin/bash
#SBATCH -J bp5_v50c_p4_br2        # Job name
#SBATCH -o bp5_v50c_p4_br2_%j.out
#SBATCH -e bp5_v50c_p4_br2_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # 200 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v50c: 2500m mesh, p=4, BR2 method (NO IP penalty)
#
# DIAGNOSTIC: Tests whether p=4 nucleation failure is caused by IP penalty.
# BR2 uses lifting operator instead of penalty → no penalty over-stiffening.
#
# If BR2 nucleates correctly (V grows), IP penalty is the root cause.
# If BR2 also fails, the problem is elsewhere (multi-DOF, psi init, etc).

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# BUILD job — rebuild with penalty_factor support
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v50c_p4_br2

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method BR2 \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v50c_p4_br2 \
      --output-prefix bp5_v50c_p4_br2 \
      --write-every-step
