#!/bin/bash
#SBATCH -J bp1_25m            # Job name
#SBATCH -o bp1_25m_%j.out     # stdout (%j = job ID)
#SBATCH -e bp1_25m_%j.err     # stderr
#SBATCH -p normal             # Queue
#SBATCH -N 1                  # 1 node
#SBATCH -n 40                 # MPI tasks (~61k elements, 40 cores for 3000yr)
#SBATCH -t 48:00:00           # wall time
#SBATCH -A EAR20006           # TACC allocation

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# Rebuild with latest optimizations
make seas_bp1_full

# Generate mesh if not present
if [ ! -f bp1/mesh/bp1_25m.msh ]; then
    gmsh -2 bp1/mesh/bp1.geo -o bp1/mesh/bp1_25m.msh -setnumber hf 0.025
fi

mkdir -p bp1/results_25m

ibrun ./seas_bp1_full \
      --mesh bp1/mesh/bp1_25m.msh \
      --ref-dir bp1/benchmark_data \
      --output-dir bp1/results_25m \
      --checkpoint-interval 5000
