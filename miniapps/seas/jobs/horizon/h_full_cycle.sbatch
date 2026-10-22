#!/bin/bash
#SBATCH -J h_full_cycle
#SBATCH -o h_full_cycle_%j.out
#SBATCH -e h_full_cycle_%j.err
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -A ear26024

# GPU plan G5 on Horizon: jobs/vista/g5_gpu_full_cycle.sbatch (the
# 1500 yr BP5 cycle at p=1 on one GPU, deck of the CPU gold 1074099, with
# speed.txt and the stage timers) on one GB200, with the settings of
# jobs/horizon/h_env.sh.  QD_PATH picks the stage solve as there: cudss
# (default) or green (opt-in, reads GREEN_CACHE).
#
#   sbatch [--export=ALL,QD_PATH=green] jobs/horizon/h_full_cycle.sbatch
#   (from miniapps/seas of the Horizon clone, after h_gpu_gates passed)

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/horizon/h_env.sh" || exit 1
exec bash "${SLURM_SUBMIT_DIR}/jobs/vista/g5_gpu_full_cycle.sbatch"
