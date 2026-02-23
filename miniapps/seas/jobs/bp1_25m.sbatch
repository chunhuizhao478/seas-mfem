#!/bin/bash
#SBATCH -J bp1_25m            # Job name
#SBATCH -o bp1_25m_%j.out     # stdout (%j = job ID)
#SBATCH -e bp1_25m_%j.err     # stderr
#SBATCH -p normal             # Queue
#SBATCH -N 1                  # 1 node
#SBATCH -n 20                 # MPI tasks
#SBATCH -t 96:00:00           # wall time (BP1: ~30 EQs over 3000 yr)
#SBATCH -A EAR20006           # TACC allocation

export LC_ALL=C
export LANG=C

module list

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

mkdir -p bp1/results_25m

ibrun ./seas_bp1_full \
      --mesh bp1/mesh/bp1_25m.msh \
      --ref-dir bp1/benchmark_data \
      --output-dir bp1/results_25m \
      --checkpoint-interval 5000
