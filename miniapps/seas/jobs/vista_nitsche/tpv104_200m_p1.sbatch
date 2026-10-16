#!/bin/bash
#SBATCH -J tpv104_200m_p1
#SBATCH -o tpv104_200m_p1_%j.out
#SBATCH -e tpv104_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 03:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, tpv104 on the 200 m mesh, p = 1, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv104/configs/tpv104_nitsche.toml (Chunhui Zhao, 2026-10-16).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Estimate (Frontera Nitsche TPV102 200 m KV+PC timings scaled by tets and CFL dt; Grace core taken
# equal to a Frontera core, unverified): about 510 core-h = 3.6 node-h = 1.2 SU, about 0.9 h on 4 nodes.  Requested: 4 nodes x 3 h x 0.33 SU = 4.0 SU at most.
# 144 ranks per node is provisional: the 200 m load check in build_verify_vista.sbatch decides it (over about 200 GB per node -> 128 per node).
PROB=tpv104
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_200m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
