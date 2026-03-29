#!/bin/bash
#SBATCH -J bp5_v55_prod_p1           # Job name
#SBATCH -o bp5_v55_prod_p1_%j.out
#SBATCH -e bp5_v55_prod_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v55 PRODUCTION: Combined integrator, full 1800-year run
#
# All v55 fixes: D8 (V_vec sign) + D4 (log10 Brent) + D5 (combined integrator)
# Dev run (2hr) showed nucleation: V bottoms at 0.0067 at t~60s then
# recovers to 0.06 at t=143s. This production run should capture the
# full earthquake and subsequent interseismic/nucleation cycles.

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

mkdir -p bp5/results_v55_1000m_ip_p1_production

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v55_1000m_ip_p1_production \
      --output-prefix bp5_v55_prod_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
