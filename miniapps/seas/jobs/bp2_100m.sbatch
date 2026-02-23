#!/bin/bash
#SBATCH -J bp2_100m           # Job name
#SBATCH -o bp2_100m_%j.out    # stdout (%j = job ID)
#SBATCH -e bp2_100m_%j.err    # stderr
#SBATCH -p normal             # Queue
#SBATCH -N 1                  # 1 node
#SBATCH -n 8                  # MPI tasks (~8k elements, 8 cores optimal)
#SBATCH -t 48:00:00           # wall time
#SBATCH -A EAR20006           # TACC allocation

export LC_ALL=C
export LANG=C

module list

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

mkdir -p bp2/results_100m

ibrun ./seas_bp2_full \
      --mesh bp2/mesh/bp2_100m.msh \
      --ref-dir bp2/benchmark_data_200m \
      --output-dir bp2/results_100m \
      --checkpoint-interval 5000
