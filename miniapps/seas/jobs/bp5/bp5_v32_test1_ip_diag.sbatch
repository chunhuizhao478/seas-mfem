#!/bin/bash
#SBATCH -J bp5_v32_t1_ip_diag   # Job name
#SBATCH -o bp5_v32_t1_ip_diag_%j.out
#SBATCH -e bp5_v32_t1_ip_diag_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 00:30:00             # 30 minutes (short diagnostic run)
#SBATCH -A EAR20006             # TACC allocation

# BP5 v32 Test 1: IP method + uniform fault + traction decomposition diagnostic
#
# PURPOSE: Diagnose IP blowup root cause by decomposing traction into
#   T_stress = {sigma}·n_hat   (average stress traction)
#   T_penalty = penalty * jump (IP penalty correction)
#   jump = [[u]] - sign*delta  (residual displacement jump)
#
# If |jump| ~ 0: instability comes from stress computation
# If |jump| >> 0: bilinear form / RHS / solver inconsistency
#
# Short run: just enough steps to see the instability seed (~100 steps)

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

mkdir -p bp5/results_v32_test1_ip_diag

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --diag-traction-decomp \
      --diag-vtk \
      --delta-tau-factor 0 \
      --V-nuc 1e-9 \
      --tfinal 1e8 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v32_test1_ip_diag \
      --output-prefix bp5_v32_t1_ip_diag \
      --checkpoint-interval 1000 \
      --write-every-step
