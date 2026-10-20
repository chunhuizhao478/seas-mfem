#!/bin/bash
#SBATCH -J tpv205_100m_p1
#SBATCH -o tpv205_100m_p1_%j.out
#SBATCH -e tpv205_100m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 8
#SBATCH -n 512
#SBATCH -t 04:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, tpv205 on the 100 m mesh, p = 1, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv205/configs/tpv205_nitsche.toml (Chunhui Zhao, 2026-10-16).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Mesh: tpv205/mesh/tpv2053d_100m.msh, the mesh of condo job 6490984 (Mac gmsh 4.15.0 output of
# tpv2053d_100m.geo): 5,438,053 tets, 110,004 steps to 12 s.  Vista does not generate it: copy it into
# tpv205/mesh/ first; the run checks it against meshes.sha256.
# 64 ranks per node: every rank loads the whole serial mesh (drivers/nitsche_dyn_driver.cpp:750).
# Vista 1076378 (TPV205 200 m, 1.89 M tets) used about 0.81 GB per rank (MaxRSS 117 GB at 144 per node),
# so about 2.3 GB per rank here and about 150 GB per node at 64 per node, of 237 GB.
# Estimate from Vista 1076378 (5.9e-6 core-s per tet-step, 0.48 s per fault frame): 980 core-h,
# 1.9 h on 512 ranks + about 0.3 h of fault output = 2.3 h, 8 nodes x 2.3 h x 0.33 = 6 SU.
# Requested: 8 nodes x 4 h x 0.33 SU = 10.6 SU at most.
PROB=tpv205
RES=100m
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_100m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
