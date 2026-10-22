#!/bin/bash
#SBATCH -J h_green_cache
#SBATCH -o h_green_cache_%j.out
#SBATCH -e h_green_cache_%j.err
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 02:00:00
#SBATCH -A EAR26024

# GPU plan on Horizon: write the dense Green's-function cache GREEN_CACHE
# (production mesh, p=1, 37.5 GB) that the G3c gate and the dense-G full
# cycle read.  The Vista file cannot be copied here directly, and on Vista
# the QD speed job built it in about 4 minutes on one 144-rank node.
#
# Runs jobs/vista_qd_speed/qd_speed_timing.sbatch with MODES=green and
# RUN_STEPS=10 on the 144 Grace cores of one node, on the CPU (no
# --device): it builds G with the solve path, writes it with its read-back
# check, checks it against the solve path, then takes 10 steps.  If the
# file already exists, the run reads and checks it instead.
#
#   sbatch jobs/horizon/h_green_cache.sbatch   (from miniapps/seas of the Horizon clone)

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/horizon/h_env.sh" || exit 1
export MODES=green RUN_STEPS=10 ORDER=1
exec bash "${SLURM_SUBMIT_DIR}/jobs/vista_qd_speed/qd_speed_timing.sbatch"
