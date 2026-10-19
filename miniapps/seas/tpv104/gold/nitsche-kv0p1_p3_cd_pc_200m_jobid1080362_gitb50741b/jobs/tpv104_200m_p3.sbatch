#!/bin/bash
#SBATCH -J tpv104_200m_p3
#SBATCH -o tpv104_200m_p3_%j.out
#SBATCH -e tpv104_200m_p3_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2048
#SBATCH -t 04:30:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, tpv104 on the 200 m mesh, p = 3, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv104/configs/tpv104_nitsche.toml (Chunhui Zhao, 2026-10-19).
# Same deck, mesh and 16 x 128 ranks as the p = 2 gold (Vista 1075917); only the order differs.
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Estimate: dt = 0.591 x the p = 2 dt, 176,570 steps to 12 s; step 47 to 58 ms (the p = 2 step with
# the speedups, 11 to 14 ms, from Vista 1079640 scaled by the old TPV104 / TPV205 step ratio of the
# p = 2 golds 1075917 / 1075916 or by tets per rank, times 4.24 for p = 3): 2.3 to 2.9 h of steps,
# about 2.5 to 3.1 h = 40 to 50 node-h = 13 to 16 SU.
# Requested: 16 nodes x 4.5 h x 0.33 SU = 23.8 SU at most.
# 128 ranks per node as at p = 2 (every rank reads the serial mesh, 1.30 x the TPV205 mesh).  p = 3
# memory is unmeasured; run_nitsche_vista.sh logs the first node's peak.
PROB=tpv104
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_200m.msh
ORDER=3
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
