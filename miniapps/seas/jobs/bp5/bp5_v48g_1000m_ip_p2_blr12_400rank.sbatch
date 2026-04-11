#!/bin/bash
#SBATCH -J bp5_v48g_1k_blr12    # Job name
#SBATCH -o bp5_v48g_1k_blr12_%j.out
#SBATCH -e bp5_v48g_1k_blr12_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 ranks
#SBATCH -t 48:00:00                  # 48 hours (production run)
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v48g: IP p=2, 1000m mesh, 400 ranks, BLR tol=1e-12 (NEW DEFAULT)
#
# ROOT CAUSE (v48): MUMPS-BLR at 1e-10 was too loose for the ×3 IP penalty.
# The BLR approximation error changed the effective fault stiffness enough
# to suppress nucleation and trigger blowup.
#
# Evidence (2500m mesh, 8 ranks):
#   v48c (exact MUMPS):  V GROWS 0.010→0.019  ✓  (job 7606315)
#   v48d (BLR 1e-12):    V GROWS 0.010→0.014  ✓  (job 7606342)
#   v48e (BLR 1e-14):    V GROWS 0.010→0.016  ✓  (job 7606340)
#   v48  (BLR 1e-10):    V DECAYS 0.010→0.004 ✗  (job 7606241)
#
# This is the PRODUCTION test: 1000m mesh (the mesh that always crashed
# with BLR 1e-10) with the new default BLR 1e-12.
#
# SUCCESS CRITERIA:
#   1. No TRACTION BLOWUP / crash
#   2. V_max grows (nucleation proceeds)
#   3. First earthquake occurs
#
# See bp5_debug_v48.md Section 13.

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

mkdir -p bp5/results_v48g_1000m_blr12_400rank

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v48g_1000m_blr12_400rank \
      --output-prefix bp5_v48g_1k_blr12 \
      --checkpoint-interval 1000 \
      --write-every-step
