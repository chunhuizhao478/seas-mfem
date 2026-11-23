#!/bin/bash
#SBATCH -J tpv205_200m_p1
#SBATCH -o tpv205_200m_p1_%j.out
#SBATCH -e tpv205_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 03:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, tpv205 on the 200 m mesh, p = 1, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv205/configs/tpv205_nitsche.toml (Chunhui Zhao, 2026-10-16).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Estimate (Frontera Nitsche TPV102 200 m KV+PC timings scaled by tets and CFL dt; Grace core taken
# equal to a Frontera core, unverified): about 420 core-h, 0.7 h on 4 nodes = 2.9 node-h = 1.0 SU.
# Requested: 4 nodes x 3 h x 0.33 SU = 4.0 SU at most.
# 144 ranks per node: the 200 m load check in build_verify_vista.sbatch confirms it fits.
PROB=tpv205
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
