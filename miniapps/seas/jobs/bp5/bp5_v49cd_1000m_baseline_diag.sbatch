#!/bin/bash
#SBATCH -J bp5_v49cd_diag        # Job name
#SBATCH -o bp5_v49cd_diag_%j.out
#SBATCH -e bp5_v49cd_diag_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node (192 GB on Frontera)
#SBATCH -n 1                         # 1 rank — SERIAL
#SBATCH -t 01:00:00                  # 1 hour (will crash quickly)
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v49cd: BASELINE DIAGNOSTIC ON 1000m (no fixes, will crash)
#
# Collects all diagnostic data on the 1000m mesh WITHOUT any fixes:
#   --diag-normals:          CalcOrtho vs basis.normal consistency (Candidate 4)
#   --diag-first-traction:   Per-DOF traction at first stage (seed location)
#   --diag-rk-stages:        Per-stage V_max/slip_max (cascade trace)
#
# This is the 1000m REFERENCE run. Compare with v49e (2500m reference).
# The job WILL crash (same instability as v48i). The diagnostics print
# BEFORE the crash, giving us the seed and cascade data.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# Build already done by v49a (same code). Skip if binary exists.
if [ ! -f seas_bp5_full ]; then
   make clean-build
   make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }
fi

mkdir -p bp5/results_v49cd_1000m_diag

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v49cd_1000m_diag \
      --output-prefix bp5_v49cd_diag \
      --write-every-step \
      --diag-normals \
      --diag-first-traction \
      --diag-rk-stages
