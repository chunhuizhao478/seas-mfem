#!/bin/bash
#SBATCH -J tpv104_mirror_200m_p4
#SBATCH -o tpv104_mirror_200m_p4_%j.out
#SBATCH -e tpv104_mirror_200m_p4_%j.err
#SBATCH -p gg
#SBATCH -N 32
#SBATCH -n 4608
#SBATCH -t 10:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, SCEC TPV104 at p = 4 on the mirror-symmetric 200 m mesh with the 22 km fault
# (tpv104/mesh/tpv104_mirror_ext22_200m.msh, make_mirror_meshes_vista.sbatch), 0 -> 12.0 s, the SCEC
# submission run (Chunhui Zhao, 2026-11-21: "Mirror mesh for TPV104 p4"; "for benchmark case, mirror
# mesh is fine and standard").  tpv104_200m_p4.sbatch (1156772) on the gmsh mesh aborted at t = 8.185 s
# with sigma_n_eff = -1.1e4 Pa at the fault end at the free surface (18000, 0, -200): that mesh is not
# mirror-symmetric about the fault, so sigma_n swings although the exact change is zero; on the mirror
# mesh it stays at 120 MPa (1000 m cloud runs, 2026-10-17).  The fault extends 4 km past the 18 km
# transition edge into the velocity-strengthening border the deck's boxcar_taper already gives.
# Deck as committed in tpv104/configs/tpv104_nitsche.toml (off-fault and on-fault stations, rupture times).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch and
# make_mirror_meshes_vista.sbatch passed at this commit.
# Layout: 32 x 144 ranks (the gg maximum); 3,453,310 tets: 107.9 k per node.
# Memory (estimate): 1156772 peaked at 70.2 GiB committed per node with 77.1 k tets per node; scaled by
# tets per node about 98 GiB, under the 118.8 GiB limit.  --dense-fill direct as 1156772.
# Time (estimate, not measured): 1156772 ran 181,860 steps in 10,852 s driver wall (59.7 ms per step,
# setup included), so about 266,600 steps to 12 s on the gmsh mesh; the mirror mesh has min 6V/A
# 24.24 m (gmsh 26.81 m), so x 1.106 steps, and x 1.40 tets: about 6.8 h = 72 SU.
# Stall guard as in TPV34 p = 4.
# Requested: 32 nodes x 10 h x 0.33 SU = 105.6 SU at most.
PROB=tpv104
RES=mirror200m
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_mirror_ext22_200m.msh
ORDER=4
TFINAL=12.0
DRIVER_ARGS="--dense-fill direct"
STALL_SAMPLE_AFTER=900
STALL_ABORT=1
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
