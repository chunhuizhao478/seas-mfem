#!/bin/bash
#SBATCH -J bp2_50m            # Job name
#SBATCH -o bp2_50m_%j.out     # stdout (%j = job ID)
#SBATCH -e bp2_50m_%j.err     # stderr
#SBATCH -p normal             # Queue
#SBATCH -N 1                  # 1 node
#SBATCH -n 10                 # MPI tasks (31k elements, 10 cores optimal)
#SBATCH -t 48:00:00           # wall time
#SBATCH -A EAR20006           # TACC allocation

export LC_ALL=C
export LANG=C

module list

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

mkdir -p bp2/results_50m

ibrun ./seas_bp2_full \
      --mesh bp2/mesh/bp2_50m.msh \
      --ref-dir bp2/benchmark_data_200m \
      --output-dir bp2/results_50m \
      --checkpoint-interval 5000
