#!/bin/bash
#SBATCH -J bp5_v49i_vguard       # Job name
#SBATCH -o bp5_v49i_vguard_%j.out
#SBATCH -e bp5_v49i_vguard_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node (192 GB on Frontera)
#SBATCH -n 1                         # 1 rank — SERIAL
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v49i: V GUARD — dt_init = 0.13s + stage-level V monitoring on 1000m mesh
#
# Uses the SAME dt_init = 0.13s as the crashing runs, but with
# --v-guard 100 (reject step if any stage produces V > 100 × V_stage0).
#
# The V guard detects the cascade at stage 3 (where V amplifies 40×),
# rejects the step, and halves dt. This should eventually reach a stable dt
# without pre-computing the CFL limit.
#
# Expected: STABLE after a few rejected steps. The V guard catches the
# CFL violation and reduces dt adaptively.
# V guard factor = 100 means: if V grows by more than 100× within a
# single step's stages, something is wrong → reject.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# NO BUILD — v49f builds. Wait for binary to exist.
for i in $(seq 1 60); do
   [ -f seas_bp5_full ] && break
   echo "Waiting for v49f build... ($i/60)"
   sleep 10
done
[ ! -f seas_bp5_full ] && { echo "ERROR: seas_bp5_full not found after 10min"; exit 1; }

mkdir -p bp5/results_v49i_vguard

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
      --output-dir bp5/results_v49i_vguard \
      --output-prefix bp5_v49i_vguard \
      --write-every-step \
      --v-guard 100 \
      --diag-rk-stages
