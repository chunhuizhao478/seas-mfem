#!/bin/bash
#SBATCH -J bp1_ss_50m_ip      # Job name
#SBATCH -o bp1_ss_50m_ip_%j.out  # stdout (%j = job ID)
#SBATCH -e bp1_ss_50m_ip_%j.err  # stderr
#SBATCH -p normal             # Queue
#SBATCH -N 2                  # 2 nodes
#SBATCH -n 80                 # MPI tasks (~290k elements, ~3.6k elem/core)
#SBATCH -t 48:00:00           # wall time
#SBATCH -A EAR20006           # TACC allocation

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

# Ensure hypre/2.31.0 is found before PETSc's old bundled libHYPRE.so
export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd "$SLURM_SUBMIT_DIR/../.."

# Rebuild with latest optimizations
make seas_bp1_bdrload

# Generate mesh if not present
if [ ! -f bp1/mesh/bp1_ss_50m.msh ]; then
    gmsh -2 bp1/mesh/bp1_selfsimilar.geo -o bp1/mesh/bp1_ss_50m.msh -setnumber h 0.05
fi

mkdir -p bp1/results_ss_50m_ip

ibrun ./seas_bp1_bdrload \
      --mesh bp1/mesh/bp1_ss_50m.msh \
      --ref-dir bp1/benchmark_data \
      --output-dir bp1/results_ss_50m_ip \
      --dg-method IP \
      --checkpoint-interval 5000
