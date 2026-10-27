#!/bin/bash
#SBATCH -J tpv205_mirror_200m_p1
#SBATCH -o tpv205_mirror_200m_p1_%j.out
#SBATCH -e tpv205_mirror_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 02:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, tpv205 on the mirror-symmetric 200 m mesh, p = 1, 0 -> 12.0 s, Kelvin-Voigt
# eta_bar 0.1 + Prakash-Clifton as committed in tpv205/configs/tpv205_nitsche.toml (Chunhui Zhao,
# 2026-10-17).  Same fault and stations as tpv205_200m_p1.sbatch; only the mesh is mirror-symmetric.
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch and
# make_mirror_meshes_vista.sbatch passed at this commit.
# Mesh (cloud gmsh 4.15.2; Vista's differs slightly): 1,976,500 tets, 26,688 fault faces, min 6V/A
# 30.60 m; the driver's dt on it is 2.8817e-4 s, 41,600 steps to 12 s (old 200 m mesh: 1.89 M tets,
# 41,614 steps).  Estimate from Vista 1076378 (TPV205 200 m p1: 8.3e-6 core-s per tet-step, whole
# job): 190 core-h, 0.35 h on 576 ranks, 4 nodes x 0.35 h x 0.33 = 0.5 SU.
# Requested: 4 nodes x 2 h x 0.33 SU = 2.6 SU at most.
# 144 ranks per node: 1076378 used 0.81 GB per rank on 1.89 M tets, so about 0.85 GB per rank, 122 GB per node.
PROB=tpv205
RES=mirror200m
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_mirror_200m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
