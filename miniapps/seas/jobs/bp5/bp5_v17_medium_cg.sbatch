#!/bin/bash
#SBATCH -J bp5_v17_medium_cg    # Job name
#SBATCH -o bp5_v17_medium_cg_%j.out
#SBATCH -e bp5_v17_medium_cg_%j.err
#SBATCH -p normal               # Queue
#SBATCH -N 4                    # nodes
#SBATCH -n 224                  # MPI tasks (56 per node)
#SBATCH -t 12:00:00             # 12 hours (medium run)
#SBATCH -A EAR20006             # TACC allocation

# BP5 v17 verification step 3: medium run (~5 yr) with CG+AMG
# Goal: stress rate ~0.001 MPa/yr, no spurious EQ#2

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

mkdir -p bp5/results_v17_medium_cg

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver cg \
      --check-residual \
      --monitor-traction 500 \
      --tfinal 1.575e8 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v17_medium_cg \
      --output-prefix bp5_v17_medium_cg \
      --checkpoint-interval 5000
