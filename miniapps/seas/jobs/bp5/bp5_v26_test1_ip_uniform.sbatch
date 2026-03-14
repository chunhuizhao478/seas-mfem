#!/bin/bash
#SBATCH -J bp5_v26_t1_ip_uni    # Job name
#SBATCH -o bp5_v26_t1_ip_uni_%j.out
#SBATCH -e bp5_v26_t1_ip_uni_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v26 Test 1: IP method + new uniform mesh — verify VS zone health
#
# H27: Uniform fault mesh (Tandem-style, committed to git)
# H28: DG method = IP (not BR2) — avoids BR2 traction lifting issue
#
# V_nuc = Vp, delta_tau = 0 (completely uniform, no nucleation)
# PASS: V/Vp stays in [0.5, 2.0] at z=22km for 300 years, no blowup

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
mkdir -p bp5/results_v26_test1_ip_uniform

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --dg-method IP \
      --solver mumps-blr \
      --bc-mode far-field \
      --check-residual \
      --delta-tau-factor 0 \
      --V-nuc 1e-9 \
      --tfinal 9.45e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v26_test1_ip_uniform \
      --output-prefix bp5_v26_t1_ip_uni \
      --checkpoint-interval 1000 \
      --write-every-step
