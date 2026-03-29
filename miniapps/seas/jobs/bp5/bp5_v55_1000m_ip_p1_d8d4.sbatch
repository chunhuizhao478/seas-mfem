#!/bin/bash
#SBATCH -J bp5_v55_d8d4_p1           # Job name
#SBATCH -o bp5_v55_d8d4_p1_%j.out
#SBATCH -e bp5_v55_d8d4_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v55: D8+D4 fixes — first run with Tandem-matched friction solver
#
# Changes from v53/v54:
#   D8: V_vec negated to match Tandem's -(V/|tau|)*tau convention
#       (DieterichRuinaBase.h:174). Below-Wf hardcoded rate also negated.
#       The sign chain (V→S→EmbedSlip→sign*delta_u) is self-consistent:
#       MFEM's sign factor in slip assembly handles face orientation.
#
#   D4: Brent solver switched from linear V space [0, tau/eta] to
#       log10(V) space [-32, log10(tau/eta)], matching Tandem's
#       DieterichRuinaBase::slip_rate (lines 90-132). Provides ~14 extra
#       digits of relative precision at interseismic rates (V~1e-12).
#       Fallback bracket [log10(nextafter(0,1)), log10(tau/eta)] added.
#       eta=0 direct inversion path added.
#
#   Test 16 fix: correction/stress ratio assertion relaxed (50x ratio is
#       inherent to IP DG p=1, documented in v54 Section 9.4).
#
# All unit tests pass: 116/116 (cross_verify), 76/76 (bp5_fault_operator),
# 32/32 (friction_law), 23/23 (state_evolution).
#
# Configuration:
#   - Mesh: bp5_tandem.msh (1000m reference mesh)
#   - Order: p=1
#   - DG: IP
#   - Solver: MUMPS-BLR, tol=1e-12
#   - Runtime: 48h wall clock, 8 nodes / 400 MPI ranks
#
# Key question: do D8+D4 change the nucleation trajectory?
#   v53/v54: V_max peaks at t=3.1s then dies (0.0136 → 0.0080 by t=38s)
#   Tandem:  V_max dips to 0.0114 then re-accelerates to earthquake at t=82s
#   If D8+D4 were contributing to the mismatch, we should see V stay higher.

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

mkdir -p bp5/results_v55_1000m_ip_p1_d8d4

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
      --output-dir bp5/results_v55_1000m_ip_p1_d8d4 \
      --output-prefix bp5_v55_d8d4_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
