#!/bin/bash
#SBATCH -J tpv31_200m_p3
#SBATCH -o tpv31_200m_p3_%j.out
#SBATCH -e tpv31_200m_p3_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 04:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, tpv31 on the 200 m mesh, p = 3, 0 -> 15.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv31/configs/tpv31_nitsche.toml (Chunhui Zhao, TPV6/TPV7/TPV31 plan
# rev 1, 2026-10-19).  Rank layout of the TPV205 p = 3 gold (16 x 144).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Estimate: 202,590 steps to 15.0 s (dt from the min 6V/A / c_p of the cloud mesh against TPV205's);
# step cost from the TPV205 p = 3 gold (Vista 1080361, 31.1 ms of driver wall per step) scaled by the
# tet count, plus 2 min of setup; about 1.9 h = 10.0 SU.
# The built-in kernel mesh (cloud: 2,031,963 tets) has dt 0.883 x TPV205's against 0.926 for the
# first (OpenCASCADE) mesh this estimate used: about 5 % more steps, inside the time limit.
# Requested: 16 nodes x 4 h x 0.33 SU = 21.1 SU at most.
PROB=tpv31
DECK=tpv31/configs/tpv31_nitsche.toml
MESH=tpv31/mesh/tpv31_nitsche_200m.msh
ORDER=3
TFINAL=15.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
