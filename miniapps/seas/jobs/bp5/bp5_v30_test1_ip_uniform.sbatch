#!/bin/bash
#SBATCH -J bp5_v30_t1_ip_uni    # Job name
#SBATCH -o bp5_v30_t1_ip_uni_%j.out
#SBATCH -e bp5_v30_t1_ip_uni_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v30 Test 1: IP method + Tandem coordinate system + uniform fault
#
# What's new in v30:
#   - Tandem coordinate system exclusively (X=strike, Y=normal, Z=depth<0)
#   - Tandem mesh directly (tags 1=Natural, 3=Fault, 5=Dirichlet)
#   - ref_normal = (0,-1,0), dip = strike x n (no negation)
#   - Fixed BR2 nor_q(0)->nor_q(1) bug
#   - Fixed tau_pre sign (negated to match Tandem)
#   - Fixed shared-face depth/coord bugs
#
# Uniform fault: V_nuc=Vp, delta_tau=0 (no nucleation, no earthquake)
# PASS criterion: V/Vp in [0.5, 2.0] for 300yr, no blowup, no VS lockup

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

mkdir -p bp5/results_v30_test1_ip_uniform

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --diag-vtk \
      --delta-tau-factor 0 \
      --V-nuc 1e-9 \
      --tfinal 9.45e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v30_test1_ip_uniform \
      --output-prefix bp5_v30_t1_ip_uni \
      --checkpoint-interval 1000 \
      --write-every-step
