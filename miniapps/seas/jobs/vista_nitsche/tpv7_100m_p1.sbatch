#!/bin/bash
#SBATCH -J tpv7_100m_p1
#SBATCH -o tpv7_100m_p1_%j.out
#SBATCH -e tpv7_100m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 8
#SBATCH -n 384
#SBATCH -t 03:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, tpv7 on the 100 m mesh, p = 1, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv7/configs/tpv7_nitsche.toml (Chunhui Zhao, TPV6/TPV7/TPV31
# plan rev 1.1, 2026-10-19: 100 m, the node spacing the description asks for).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Mesh: tpv6/mesh/tpv6_nitsche_100m.msh, made on Vista by build_verify_vista.sbatch (TPV6 and TPV7
# share it).  The cloud gmsh 4.15.2 mesh of the same .geo: 6,037,351 tets, min 6V/A 12.5 m on the
# near (6000 m/s) side, so dt is 0.411 x the TPV205 200 m dt at every p.
# 48 ranks per node: every rank reads the whole serial mesh at startup (drivers/nitsche_dyn_driver.cpp:829),
# and the first-node peaks of the p = 3 golds (145.1 GB at 144 ranks per node for TPV205's 1.89 M tets,
# Vista 1080361; 177.7 GB at 128 per node for TPV104's 2.47 M tets, 1080362) come to 0.53-0.56 GB per
# rank per million tets, so 3.2-3.4 GB per rank here: 154-163 GB per node at 48 per node, of 237.6 GB
# (205-218 GB at 64).
# Estimate: 101,288 steps to 12.0 s.
# Step cost from the TPV205 100 m p = 1 gold (Vista 1076419, 4318 s of stepping for 110,004 steps on
# 8 x 64, 5.44 M tets, before the speedups) scaled by the tet count: about 44 ms per step if the
# step cost goes with tets per node, 58 ms if it goes with tets per rank; 1.2 to 1.6 h, plus about
# 0.3 h of fault output and 5 min of setup = 1.6 to 2.0 h, 8 nodes x 0.33 SU = 4.2 to 5.3 SU.
# Requested: 8 nodes x 3 h x 0.33 SU = 7.9 SU at most.
PROB=tpv7
RES=100m
DECK=tpv7/configs/tpv7_nitsche.toml
MESH=tpv6/mesh/tpv6_nitsche_100m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
