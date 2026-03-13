#!/bin/bash
#SBATCH -J bp5_v17_full_cg      # Job name
#SBATCH -o bp5_v17_full_cg_%j.out
#SBATCH -e bp5_v17_full_cg_%j.err
#SBATCH -p normal               # Queue
#SBATCH -N 4                    # nodes
#SBATCH -n 224                  # MPI tasks (56 per node)
#SBATCH -t 48:00:00             # 48 hours (full run)
#SBATCH -A EAR20006             # TACC allocation

# BP5 v17 verification step 4: full run with CG+AMG (default solver)
# Goal: correct EQ recurrence (~200 yr), no spurious events

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

make seas_bp5_full

if [ ! -f bp5/mesh/bp5_1000m.msh ]; then
    gmsh -3 bp5/mesh/bp5.geo -o bp5/mesh/bp5_1000m.msh -setnumber h 1.0
fi

mkdir -p bp5/results_v17_full_cg

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver cg \
      --check-residual \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v17_full_cg \
      --output-prefix bp5_v17_full_cg \
      --checkpoint-interval 5000
