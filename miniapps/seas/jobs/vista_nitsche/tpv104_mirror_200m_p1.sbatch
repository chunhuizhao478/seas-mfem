#!/bin/bash
#SBATCH -J tpv104_mirror_200m_p1
#SBATCH -o tpv104_mirror_200m_p1_%j.out
#SBATCH -e tpv104_mirror_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 512
#SBATCH -t 02:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, tpv104 on the mirror-symmetric 200 m mesh with the 22 km fault, p = 1, 0 -> 12.0 s,
# Kelvin-Voigt eta_bar 0.1 + Prakash-Clifton as committed in tpv104/configs/tpv104_nitsche.toml
# (Chunhui Zhao, 2026-10-17).  The fault extends 4 km past the 18 km transition edge into the
# velocity-strengthening border the deck's boxcar_taper already gives (a = 0.02, V_w = 1 m/s).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch and
# make_mirror_meshes_vista.sbatch passed at this commit.
# Mesh (cloud gmsh 4.15.2; Vista's differs slightly): 3,453,310 tets, 56,506 fault faces, min 6V/A
# 24.24 m; the driver's dt on it is 2.2826e-4 s, 52,600 steps to 12 s (old 200 m mesh: 2.47 M tets,
# 26.81 m, 2.45096e-4 s).  Estimate from Vista 1075919 (TPV104 200 m p1: 7.1e-6 core-s per tet-step,
# whole job): 360 core-h, 0.7 h on 512 ranks, 4 nodes x 0.7 h x 0.33 = 1 SU.
# Requested: 4 nodes x 2 h x 0.33 SU = 2.6 SU at most.
# 128 ranks per node: every rank loads the whole serial mesh (drivers/nitsche_dyn_driver.cpp:750); at
# the largest measured 0.43 GB per million tets (Vista 1076378) that is 1.5 GB per rank, 190 GB of 237 GB.
PROB=tpv104
RES=mirror200m
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_mirror_ext22_200m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
