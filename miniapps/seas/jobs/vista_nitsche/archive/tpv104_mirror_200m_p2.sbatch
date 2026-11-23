#!/bin/bash
#SBATCH -J tpv104_mirror_200m_p2
#SBATCH -o tpv104_mirror_200m_p2_%j.out
#SBATCH -e tpv104_mirror_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2048
#SBATCH -t 03:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, tpv104 on the mirror-symmetric 200 m mesh with the 22 km fault, p = 2, 0 -> 12.0 s,
# Kelvin-Voigt eta_bar 0.1 + Prakash-Clifton as committed in tpv104/configs/tpv104_nitsche.toml
# (Chunhui Zhao, 2026-10-17).  The fault extends 4 km past the 18 km transition edge into the
# velocity-strengthening border the deck's boxcar_taper already gives (a = 0.02, V_w = 1 m/s).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch and
# make_mirror_meshes_vista.sbatch passed at this commit.
# Mesh as tpv104_mirror_200m_p1.sbatch: 3,453,310 tets, min 6V/A 24.24 m, so dt about 1.04e-4 s and
# 115,500 steps to 12 s (Vista 1075917, TPV104 200 m p2 on the old mesh: 1.1496e-4 s at 26.81 m).
# Estimate from 1075917 (2.47 M tets, 104,388 steps, 3,809 s on 2,048 ranks: 3.0e-5 core-s per
# tet-step, whole job): 3,400 core-h, 1.6 h on 2,048 ranks, 16 nodes x 1.6 h x 0.33 = 8.7 SU.
# Requested: 16 nodes x 3 h x 0.33 SU = 15.8 SU at most.
# 128 ranks per node, as 1075917, which peaked at 91 GB per node; scaled by tets about 128 GB here.
PROB=tpv104
RES=mirror200m
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_mirror_ext22_200m.msh
ORDER=2
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
