#!/bin/bash
#SBATCH -J tpv31_200m_p1
#SBATCH -o tpv31_200m_p1_%j.out
#SBATCH -e tpv31_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 01:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, tpv31 on the 200 m mesh, p = 1, 0 -> 15.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv31/configs/tpv31_nitsche.toml (Chunhui Zhao, TPV6/TPV7/TPV31 plan
# rev 1, 2026-10-19).  Rank layout of the TPV205 p = 1 gold (4 x 144).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Estimate: 56,180 steps to 15.0 s (dt from the min 6V/A / c_p of the cloud mesh against TPV205's);
# step cost from the TPV205 p = 1 gold (Vista 1076378, 1102 s for 41,614 steps on 4 x 144, before the
# speedups) scaled by the tet count; about 0.45 h = 0.6 SU.
# The built-in kernel mesh (cloud: 2,031,963 tets) has dt 0.883 x TPV205's against 0.926 for the
# first (OpenCASCADE) mesh this estimate used: about 5 % more steps, inside the time limit.
# Requested: 4 nodes x 1 h x 0.33 SU = 1.3 SU at most.
PROB=tpv31
DECK=tpv31/configs/tpv31_nitsche.toml
MESH=tpv31/mesh/tpv31_nitsche_200m.msh
ORDER=1
TFINAL=15.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
