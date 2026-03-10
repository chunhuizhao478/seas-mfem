#!/bin/bash
#SBATCH -J bp1_weak_np80         # Job name
#SBATCH -o bp1_weak_np80_%j.out  # stdout (%j = job ID)
#SBATCH -e bp1_weak_np80_%j.err  # stderr
#SBATCH -p normal                # Queue
#SBATCH -N 2                     # 2 nodes
#SBATCH -n 80                    # MPI tasks
#SBATCH -t 48:00:00              # wall time
#SBATCH -A EAR20006              # TACC allocation

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd "$SLURM_SUBMIT_DIR/../.."

mkdir -p bp1/scaling_analysis/weak_scaling/np80/

ibrun ./seas_bp1_bdrload \
      --mesh bp1/mesh/bp1_ss_weak_np80.msh \
      --output-dir bp1/scaling_analysis/weak_scaling/np80/ \
      --checkpoint-interval 5000
