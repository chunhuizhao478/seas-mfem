#!/bin/bash
#SBATCH -J tpv7_mirror_100m_p3
#SBATCH -o tpv7_mirror_100m_p3_%j.out
#SBATCH -e tpv7_mirror_100m_p3_%j.err
#SBATCH -p gg
#SBATCH -N 32
#SBATCH -n 4608
#SBATCH -t 10:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, tpv7 on the 100 m MIRROR mesh, p = 3, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv7/configs/tpv7_nitsche.toml: the SCEC submission run (Chunhui
# Zhao, 2026-11-22: "yes prepare mirror mesh").  tpv7_100m_p3.sbatch with the mesh
# tpv6/mesh/tpv6_nitsche_mirror_100m.msh (make_tpv6_mirror_mesh_vista.sbatch) instead of the gmsh mesh,
# on which 1159765 aborted at t = 5.99 s (sigma_n_eff -55 kPa at (14264, 0, -49) m, fault end at the
# free surface).  Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at
# this commit and make_tpv6_mirror_mesh_vista.sbatch made the mesh.
# Layout: 32 x 144 ranks (the gg maximum), --dense-fill direct, as tpv7_100m_p3.sbatch.
# Time: 1159765 reached t = 5.99 s in 2:15 on 32 nodes, so about 4.5 h for 12 s (if the mirror mesh has
# about the same tets) = about 48 SU.  Memory: 1159765 peaked well inside the node limit.
# Stall guard as in TPV34 p = 4.
# Requested: 32 nodes x 10 h x 0.33 SU = 105.6 SU at most.
PROB=tpv7
RES=mirror100m
DECK=tpv7/configs/tpv7_nitsche.toml
MESH=tpv6/mesh/tpv6_nitsche_mirror_100m.msh
ORDER=3
TFINAL=12.0
DRIVER_ARGS="--dense-fill direct"
STALL_SAMPLE_AFTER=900
STALL_ABORT=1
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
