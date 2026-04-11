#!/bin/bash
#SBATCH -J bp5_v30_t2_ip_tan    # Job name
#SBATCH -o bp5_v30_t2_ip_tan_%j.out
#SBATCH -e bp5_v30_t2_ip_tan_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v30 Test 2: IP method + Tandem coord system + Tandem initialization
#
# All v30 fixes active (see test1 for list).
# V_nuc=0.01 (Tandem), psi-init-mode=tandem, delta_tau=0
# Expected: first event ~150yr, VS maintains V~Vp, cycling
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

mkdir -p bp5/results_v30_test2_ip_tandem

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --diag-vtk \
      --psi-init-mode tandem \
      --V-nuc 0.01 \
      --delta-tau-factor 0 \
      --tfinal 1.89e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v30_test2_ip_tandem \
      --output-prefix bp5_v30_t2_ip_tan \
      --checkpoint-interval 1000 \
      --write-every-step
