#!/bin/bash
#SBATCH -J bp5_v26_t3_ip_scec   # Job name
#SBATCH -o bp5_v26_t3_ip_scec_%j.out
#SBATCH -e bp5_v26_t3_ip_scec_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v26 Test 3: IP method + new mesh + SCEC-correct initialization
#
# H27: Uniform fault mesh
# H28: DG method = IP
#
# Default SCEC initialization (V_nuc=0.03, delta_tau_factor=1, SCEC psi)
# Expected: immediate earthquake, then VS recovers and cycles
# Run 600 years

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
mkdir -p bp5/results_v26_test3_ip_scec

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --dg-method IP \
      --solver mumps-blr \
      --bc-mode far-field \
      --check-residual \
      --tfinal 1.89e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v26_test3_ip_scec \
      --output-prefix bp5_v26_t3_ip_scec \
      --checkpoint-interval 1000 \
      --write-every-step
