#!/bin/bash
#SBATCH -J bp5_v22_t3_tpsi      # Job name
#SBATCH -o bp5_v22_t3_tpsi_%j.out
#SBATCH -e bp5_v22_t3_tpsi_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v22 Test 3: Tandem-style psi initialization (absorb delta_tau into psi)
#
# Purpose: Check if matching Tandem's initialization produces matching behavior.
# Uses --psi-init-mode tandem to absorb delta_tau into psi via InitialStatePsi.
# System should start in near-equilibrium, first event at ~150 years.
#
# Expected: first event ~150yr, recurrence ~200yr, deep VS maintains V ≈ Vp.

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

mkdir -p bp5/results_v22_test3_tandem_psi

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --bc-mode far-field \
      --check-residual \
      --psi-init-mode tandem \
      --tfinal 9.45e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v22_test3_tandem_psi \
      --output-prefix bp5_v22_t3_tpsi \
      --checkpoint-interval 1000 \
      --write-every-step
