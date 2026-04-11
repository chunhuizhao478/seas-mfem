#!/bin/bash
#SBATCH -J bp5_v23_t3_texact    # Job name
#SBATCH -o bp5_v23_t3_texact_%j.out
#SBATCH -e bp5_v23_t3_texact_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v23 Test 3-corrected: Match Tandem exactly
#
# V_nuc = 0.01 (Tandem's value, not MFEM's 0.03)
# psi-init-mode = tandem (absorb delta_tau into psi)
# delta_tau_factor = 1.0 (default, same as Tandem)
#
# Purpose: Match Tandem's exact initialization to check if MFEM reproduces cycling.
# Expected: no immediate earthquake, first event at ~150 yr, deep VS at ~Vp.
# Run 600 years to capture first event + verify recurrence.

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

if [ ! -f bp5/mesh/bp5_1000m.msh ]; then
    gmsh -3 bp5/mesh/bp5.geo -o bp5/mesh/bp5_1000m.msh -setnumber h 1.0
fi

mkdir -p bp5/results_v23_test3_tandem_exact

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --bc-mode far-field \
      --check-residual \
      --psi-init-mode tandem \
      --V-nuc 0.01 \
      --tfinal 1.89e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v23_test3_tandem_exact \
      --output-prefix bp5_v23_t3_texact \
      --checkpoint-interval 1000 \
      --write-every-step
