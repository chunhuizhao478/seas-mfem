#!/bin/bash
#SBATCH -J bp5_1000m            # Job name
#SBATCH -o bp5_1000m_%j.out     # stdout (%j = job ID)
#SBATCH -e bp5_1000m_%j.err     # stderr
#SBATCH -p normal               # Queue
#SBATCH -N 2                    # nodes
#SBATCH -n 100                  # MPI tasks
#SBATCH -t 48:00:00             # wall time
#SBATCH -A EAR20006             # TACC allocation

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

# Ensure hypre/2.31.0 is found before PETSc's old bundled libHYPRE.so
export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# Rebuild with latest optimizations
make seas_bp5_full

# Generate mesh if not present
if [ ! -f bp5/mesh/bp5_1000m.msh ]; then
    gmsh -3 bp5/mesh/bp5.geo -o bp5/mesh/bp5_1000m.msh -setnumber h 1.0
fi

mkdir -p bp5/results_1000m

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_1000m \
      --checkpoint-interval 5000
