#!/bin/bash
#SBATCH -J tpv104_100m_p1
#SBATCH -o tpv104_100m_p1_%j.out
#SBATCH -e tpv104_100m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 8
#SBATCH -n 512
#SBATCH -t 06:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, tpv104 on the 100 m mesh, p = 1, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv104/configs/tpv104_nitsche.toml (Chunhui Zhao, 2026-10-16).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Mesh: tpv102/mesh/tpv102_100m.msh, the mesh of the condo 100 m jobs (Mac gmsh 4.15.0 output of
# tpv102_100m.geo; TPV104 has the TPV102 geometry): 7,578,859 tets, 122,120 steps to 12 s.  Vista does
# not generate it: copy it into tpv102/mesh/ first; the run checks it against meshes.sha256.
# 64 ranks per node: every rank loads the whole serial mesh (drivers/nitsche_dyn_driver.cpp:750).
# Vista 1075919 (TPV104 200 m, 2.47 M tets) used about 0.72 GB per rank and condo job 6459052 peaked at
# 2.9 GB per rank on this mesh, so 64 per node needs about 140-190 GB of the node's 237 GB.
# Estimate from Vista 1075919 (6.0e-6 core-s per tet-step, Kelvin-Voigt + Prakash-Clifton): 1,550 core-h,
# 3.0 h on 512 ranks + about 0.5 h of fault output = 3.5 h, 8 nodes x 3.5 h x 0.33 = 9 SU.
# Requested: 8 nodes x 6 h x 0.33 SU = 15.8 SU at most.
PROB=tpv104
RES=100m
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv102/mesh/tpv102_100m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
