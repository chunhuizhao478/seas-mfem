#!/bin/bash
#SBATCH -J tpv104_200m_p2
#SBATCH -o tpv104_200m_p2_%j.out
#SBATCH -e tpv104_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2048
#SBATCH -t 06:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, tpv104 on the 200 m mesh, p = 2, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv104/configs/tpv104_nitsche.toml (Chunhui Zhao, 2026-10-16).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Estimate (Frontera Nitsche TPV102 200 m KV+PC timings scaled by tets and CFL dt; Grace core taken
# equal to a Frontera core, unverified): about 5,370 core-h, 2.6 h on 16 nodes = 42 node-h = 14 SU.
# Requested: 16 nodes x 6 h x 0.33 SU = 31.7 SU at most.
# 128 ranks per node: every rank loads the whole serial mesh at startup (drivers/nitsche_dyn_driver.cpp:750),
# estimated ~1.5 GB per rank for this mesh, so 144 ranks would need ~212 GB of the node's 237 GB.
PROB=tpv104
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_200m.msh
ORDER=2
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
