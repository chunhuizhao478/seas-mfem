#!/bin/bash
#SBATCH -J bp5_v25_t3_scec      # Job name
#SBATCH -o bp5_v25_t3_scec_%j.out
#SBATCH -e bp5_v25_t3_scec_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v25 Test 3: Full BP5 with SCEC-correct initialization
#
# Both fixes active:
#   H26: IP-style traction penalty (no BR2 lifting)
#   H27: Uniform fault mesh (no per-zone sizing)
#
# Default SCEC initialization (V_nuc=0.03, delta_tau_factor=1, SCEC psi)
# Tests whether SCEC-correct initial conditions work with both fixes
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
mkdir -p bp5/results_v25_test3_scec

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --bc-mode far-field \
      --check-residual \
      --tfinal 1.89e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v25_test3_scec \
      --output-prefix bp5_v25_t3_scec \
      --checkpoint-interval 1000 \
      --write-every-step
