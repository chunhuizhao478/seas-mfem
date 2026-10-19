#!/bin/bash
#SBATCH -J tpv6_100m_p2
#SBATCH -o tpv6_100m_p2_%j.out
#SBATCH -e tpv6_100m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 32
#SBATCH -n 1536
#SBATCH -t 04:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, tpv6 on the 100 m mesh, p = 2, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv6/configs/tpv6_nitsche.toml (Chunhui Zhao, TPV6/TPV7/TPV31
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
# Estimate: 215,955 steps to 12.0 s.
# Step cost from the TPV205 200 m p = 2 rerun with the speedups (Vista 1079640, about 10.7 ms of
# driver wall per step on 16 x 144, 118 k tets per node): 17 ms per step if it goes with tets per
# node (189 k here), 51 ms if it goes with tets per rank (3,930 here, 820 there); 1.0 to 3.1 h,
# plus 10 min of setup and output, 32 nodes x 0.33 SU = 12 to 35 SU.
# Requested: 32 nodes x 4 h x 0.33 SU = 42 SU at most.
PROB=tpv6
RES=100m
DECK=tpv6/configs/tpv6_nitsche.toml
MESH=tpv6/mesh/tpv6_nitsche_100m.msh
ORDER=2
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
