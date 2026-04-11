#!/bin/bash
#SBATCH -J bp5_v38a_ip_p1      # Job name
#SBATCH -o bp5_v38a_ip_p1_%j.out
#SBATCH -e bp5_v38a_ip_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v38a: IP method — p=1, 1000m mesh, 1800 yr
#
# CRITICAL TEST: Tandem BP5 uses IP (not BR2). IP was never tested after
# v31's DG slip sign fix. All v24-v28 IP blowups occurred with the wrong
# sign (penalty correcting in wrong direction → exponential blowup).
#
# IP code path has Tandem-matching formulation:
#   - Custom penalty integrator: penalty * |nor| (not kappa * |nor|²)
#   - Scalar traction: T = {σ·n̂} - penalty * ([[u]] - δ)
#   - Tag-based fault detection: excludes z=0/Wf boundary faces
#   - Per-quad-point stress + penalty in ComputeTraction
#
# Previous (BR2): v34 = 277 yr, v37 = 295 yr, Tandem (IP) = 240 yr

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

mkdir -p bp5/results_v38a_1000m_ip_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v38a_1000m_ip_p1_full \
      --output-prefix bp5_v38a_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
