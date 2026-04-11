#!/bin/bash
#SBATCH -J bp5_v18_medium_gmres # Job name
#SBATCH -o bp5_v18_medium_gmres_%j.out
#SBATCH -e bp5_v18_medium_gmres_%j.err
#SBATCH -p normal               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 48:00:00             # 48 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v18: GMRES+BlockILU solver comparison
# Goal: compare with CG+AMG — BlockILU is DG-native (element-block factorization)

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

if [ ! -f bp5/mesh/bp5_1000m.msh ]; then
    gmsh -3 bp5/mesh/bp5.geo -o bp5/mesh/bp5_1000m.msh -setnumber h 1.0
fi

mkdir -p bp5/results_v18_medium_gmres

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver gmres \
      --check-residual \
      --monitor-traction 500 \
      --tfinal 9.45e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v18_medium_gmres \
      --output-prefix bp5_v18_medium_gmres \
      --checkpoint-interval 5000
