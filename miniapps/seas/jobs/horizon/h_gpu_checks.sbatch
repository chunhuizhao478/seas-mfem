#!/bin/bash
#SBATCH -J h_gpu_checks
#SBATCH -o h_gpu_checks_%j.out
#SBATCH -e h_gpu_checks_%j.err
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -A EAR26024

# GPU plan on Horizon: jobs/vista/g1_gpu_checks.sbatch (ex1p and ex2p on
# the GPU, a device-buffer MPI_Allreduce on 2 GPUs, ex1p with cuDSS on 1 and
# 2 GPUs) on one GB200 node, with the settings of jobs/horizon/h_env.sh.
#
#   sbatch jobs/horizon/h_gpu_checks.sbatch   (from miniapps/seas of the Horizon clone)

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/horizon/h_env.sh" || exit 1
exec bash "${SLURM_SUBMIT_DIR}/jobs/vista/g1_gpu_checks.sbatch"
