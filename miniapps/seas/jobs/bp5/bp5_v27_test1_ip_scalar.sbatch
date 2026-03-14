#!/bin/bash
#SBATCH -J bp5_v27_t1_ips       # Job name
#SBATCH -o bp5_v27_t1_ips_%j.out
#SBATCH -e bp5_v27_t1_ips_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v27 Test 1: IP method + Tandem scalar traction penalty + new mesh
#
# H27: Uniform fault mesh (Tandem-style)
# H28: IP traction uses Tandem's scalar penalty formula:
#      penalty = (p0+p1)/4, p = (D+1)*c_N_1*(A/V)*(c1²/c0)
#      correction = penalty * jump (scalar × vector, isotropic)
#
# Uniform fault (V_nuc=Vp, delta_tau=0)
# PASS: V/Vp in [0.5, 2.0] at z=22km, no blowup

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

# Mesh bp5_1000m.msh committed to git (uniform res_f=1)
mkdir -p bp5/results_v27_test1_ip_scalar

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
      --output-dir bp5/results_v27_test1_ip_scalar \
      --output-prefix bp5_v27_t1_ips \
      --checkpoint-interval 1000 \
      --write-every-step
