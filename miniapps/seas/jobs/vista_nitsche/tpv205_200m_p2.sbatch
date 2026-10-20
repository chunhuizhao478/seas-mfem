#!/bin/bash
#SBATCH -J tpv205_200m_p2
#SBATCH -o tpv205_200m_p2_%j.out
#SBATCH -e tpv205_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 05:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, tpv205 on the 200 m mesh, p = 2, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv205/configs/tpv205_nitsche.toml (Chunhui Zhao, 2026-10-16).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Estimate (Frontera Nitsche TPV102 200 m KV+PC timings scaled by tets and CFL dt; Grace core taken
# equal to a Frontera core, unverified): about 4,370 core-h, 1.9 h on 16 nodes = 30 node-h = 10 SU.
# Requested: 16 nodes x 5 h x 0.33 SU = 26.4 SU at most.
# 144 ranks per node: the 200 m load check in build_verify_vista.sbatch confirms it fits.
PROB=tpv205
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=2
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
