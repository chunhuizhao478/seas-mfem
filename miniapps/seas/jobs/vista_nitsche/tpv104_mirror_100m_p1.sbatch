#!/bin/bash
#SBATCH -J tpv104_mirror_100m_p1
#SBATCH -o tpv104_mirror_100m_p1_%j.out
#SBATCH -e tpv104_mirror_100m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 640
#SBATCH -t 04:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, tpv104 on the mirror-symmetric 100 m mesh with the 22 km fault, p = 1, 0 -> 12.0 s,
# Kelvin-Voigt eta_bar 0.1 + Prakash-Clifton as committed in tpv104/configs/tpv104_nitsche.toml
# (Chunhui Zhao, 2026-10-17).  The fault extends 4 km past the 18 km transition edge into the
# velocity-strengthening border the deck's boxcar_taper already gives (a = 0.02, V_w = 1 m/s).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch and
# make_mirror_meshes_vista.sbatch passed at this commit.
# Mesh (cloud gmsh 4.15.2; Vista's differs slightly): 11,473,766 tets, 225,318 fault faces, min 6V/A
# 14.22 m, so dt about 1.34e-4 s and 89,600 steps to 12 s (the driver's dt is 9.42e-6 s per metre of
# min 6V/A on both 200 m mirror meshes; old 100 m mesh tpv102_100m.msh: 7.58 M tets, 122,120 steps).
# Estimate from Vista 1076418 (TPV104 100 m p1, 5.0e-6 core-s per tet-step, whole job): 1,440 core-h,
# 2.3 h on 640 ranks, 16 nodes x 2.3 h x 0.33 = 12 SU.
# Requested: 16 nodes x 4 h x 0.33 SU = 21.1 SU at most.
# 40 ranks per node: every rank loads the whole serial mesh (drivers/nitsche_dyn_driver.cpp:750);
# 1076418 peaked at 194 GB per node with 64 ranks on 7.58 M tets (0.40 GB per rank per million tets),
# so about 4.6 GB per rank here and 184 GB per node at 40 per node, of 237 GB.
PROB=tpv104
RES=mirror100m
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_mirror_ext22_100m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
