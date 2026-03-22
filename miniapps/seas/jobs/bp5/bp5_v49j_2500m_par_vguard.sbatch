#!/bin/bash
#SBATCH -J bp5_v49j_par_vg       # Job name
#SBATCH -o bp5_v49j_par_vg_%j.out
#SBATCH -e bp5_v49j_par_vg_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node
#SBATCH -n 48                        # 48 MPI ranks (same as v47m)
#SBATCH -t 01:00:00                  # 1 hour (v47m crashed in seconds)
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v49j: 2500m mesh, 48 ranks, V-guard ON
#
# Purpose: Determine if Issue B (2500m parallel crash from v47m) is ALSO
# a CFL/cascade issue that V-guard can fix, or a structural parallel bug.
#
# Same setup as v47m (which crashed immediately) but with:
#   --v-guard        → stage-level V rejection (factor 100×)
#   --diag-rk-stages → trace cascade through stages
#
# Expected outcomes:
#   - If V-guard fixes it: Issue B is cascade-related (shared faces may
#     have different effective stiffness), and V-guard is sufficient.
#   - If V-guard doesn't fix it: Issue B is a structural parallel bug
#     (sign, assembly, data exchange) → need K×u vs f diagnostic.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# BUILD — need Phase 2 code with V-guard
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v49j_2500m_par_vguard

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v49j_2500m_par_vguard \
      --output-prefix bp5_v49j_par_vg \
      --write-every-step \
      --v-guard \
      --diag-rk-stages
