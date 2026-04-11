#!/bin/bash
#SBATCH -J bp5_v19_short_slu    # Job name
#SBATCH -o bp5_v19_short_slu_%j.out
#SBATCH -e bp5_v19_short_slu_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours (short run)
#SBATCH -A EAR20006             # TACC allocation

# BP5 v19: SuperLU_DIST direct solver
# Goal: test if SuperLU handles memory distribution better than MUMPS
# NOTE: Requires MFEM compiled with MFEM_USE_SUPERLU=YES
#       Check: module spider superlu on Frontera

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
# module load superlu-dist  # uncomment when available
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

if [ ! -f bp5/mesh/bp5_1000m.msh ]; then
    gmsh -3 bp5/mesh/bp5.geo -o bp5/mesh/bp5_1000m.msh -setnumber h 1.0
fi

mkdir -p bp5/results_v19_short_superlu

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver superlu \
      --check-residual \
      --monitor-traction 100 \
      --tfinal 3.15e6 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v19_short_superlu \
      --output-prefix bp5_v19_short_superlu \
      --checkpoint-interval 1000 \
      --write-every-step
