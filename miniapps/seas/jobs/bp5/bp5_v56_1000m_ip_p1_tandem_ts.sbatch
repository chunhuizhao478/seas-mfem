#!/bin/bash
#SBATCH -J bp5_v56_tts_p1
#SBATCH -o bp5_v56_tts_p1_%j.out
#SBATCH -e bp5_v56_tts_p1_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5 v56: 1000 m, IP, p=1, Tandem-style time stepping policy
#
# Purpose:
# - keep the current MFEM operator/fault fixes
# - switch startup/acceptance defaults to a Tandem-style policy
#   (fixed startup dt, no V-guard, no psi clamp)
# - run the full 48-hour production case for direct comparison

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v56_1000m_ip_p1_tandem_ts

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v56_1000m_ip_p1_tandem_ts \
      --output-prefix bp5_v56_tandem_ts_p1 \
      --checkpoint-interval 1000 \
      --write-every-step \
      --tandem-time-stepping \
      --tandem-dt-init 0.01
