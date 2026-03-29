#!/bin/bash
#SBATCH -J bp5_2500m_p4
#SBATCH -o bp5_2500m_p4_%j.out
#SBATCH -e bp5_2500m_p4_%j.err
#SBATCH --chdir=/scratch2/10024/zhaochun/tandem/bp5_2500m_p4
#SBATCH -p normal
#SBATCH -N 4
#SBATCH -n 224
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

module load gcc/9.1.0
module load impi/19.0.9

export OMP_NUM_THREADS=1

TANDEM=/scratch2/10024/zhaochun/tandem/build_3d_4p_staging/app/tandem
OPTIONS=/scratch2/10024/zhaochun/tandem/tandem/examples/options

echo "num_nodes: $SLURM_JOB_NUM_NODES ntasks: $SLURM_NTASKS"

ibrun $TANDEM bp5.toml --petsc \
      -options_file ${OPTIONS}/lu_mumps.cfg \
      -options_file ${OPTIONS}/rk45.cfg \
      -ts_monitor
