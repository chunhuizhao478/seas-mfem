#!/bin/bash
#SBATCH -J bp5_v25_t1_uni       # Job name
#SBATCH -o bp5_v25_t1_uni_%j.out
#SBATCH -e bp5_v25_t1_uni_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v25 Test 1: Uniform fault — verify VS zone health
#
# Both fixes active:
#   H26: IP-style traction penalty (no BR2 lifting)
#   H27: Uniform fault mesh (no per-zone sizing)
#
# V_nuc = Vp, delta_tau = 0 (completely uniform, no nucleation)
# PASS: V/Vp stays in [0.5, 2.0] at z=22km for 300 years, no traction blowup

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

# Mesh bp5_1000m.msh is committed to git (uniform res_f=1, Tandem-style)
mkdir -p bp5/results_v25_test1_uniform

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --bc-mode far-field \
      --check-residual \
      --delta-tau-factor 0 \
      --V-nuc 1e-9 \
      --tfinal 9.45e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v25_test1_uniform \
      --output-prefix bp5_v25_t1_uni \
      --checkpoint-interval 1000 \
      --write-every-step
