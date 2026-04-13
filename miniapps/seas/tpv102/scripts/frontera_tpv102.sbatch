#!/bin/bash
#SBATCH -J tpv102_seas_mfem
#SBATCH -o tpv102_%j.out
#SBATCH -e tpv102_%j.err
#SBATCH -p normal
#SBATCH -N 4
#SBATCH -n 224
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# SCEC TPV102 benchmark: full 12s dynamic rupture simulation
# Fine mesh (h=500m, ~300k elements), order 2
#
# Before submitting:
#   1. Build driver:
#      cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas
#      make clean && make seas_tpv102_driver -j8
#   2. Generate mesh (if not committed):
#      conda activate pythonenv
#      gmsh -3 tpv102/mesh/tpv102_fine.geo -o tpv102/mesh/tpv102_fine.msh
#   3. Submit:
#      sbatch tpv102/scripts/frontera_tpv102.sbatch

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load fftw3/3.3.8
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

test -x ./seas_tpv102_driver || { echo "ERROR: driver not built"; exit 1; }

ibrun ./seas_tpv102_driver \
    --mesh tpv102/mesh/tpv102_fine.msh \
    --mesh-scale 1000 \
    --order 2 \
    --bc-mode absorbing \
    --tfinal 12.0 \
    --output-dir tpv102/results/${SLURM_JOB_ID} \
    --output-prefix tpv102
