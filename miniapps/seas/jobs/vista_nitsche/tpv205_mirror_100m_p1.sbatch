#!/bin/bash
#SBATCH -J tpv205_mirror_100m_p1
#SBATCH -o tpv205_mirror_100m_p1_%j.out
#SBATCH -e tpv205_mirror_100m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 8
#SBATCH -n 512
#SBATCH -t 03:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, tpv205 on the mirror-symmetric 100 m mesh, p = 1, 0 -> 12.0 s, Kelvin-Voigt
# eta_bar 0.1 + Prakash-Clifton as committed in tpv205/configs/tpv205_nitsche.toml (Chunhui Zhao,
# 2026-10-17).  Same fault and stations as tpv205_100m_p1.sbatch; only the mesh is mirror-symmetric.
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch and
# make_mirror_meshes_vista.sbatch passed at this commit.
# Mesh (cloud gmsh 4.15.2; Vista's differs slightly): 6,036,560 tets, 105,898 fault faces, min 6V/A
# 13.58 m, so dt about 1.28e-4 s and 93,800 steps to 12 s (the driver's dt is 9.42e-6 s per metre of
# min 6V/A on both 200 m mirror meshes; old 100 m mesh: 5.44 M tets, 110,004 steps).
# Estimate from Vista 1076419 (TPV205 100 m p1, 4.5e-6 core-s per tet-step, whole job): 710 core-h,
# 1.4 h on 512 ranks, 8 nodes x 1.4 h x 0.33 = 3.7 SU.
# Requested: 8 nodes x 3 h x 0.33 SU = 7.9 SU at most.
# 64 ranks per node, as tpv205_100m_p1.sbatch: at 0.43 GB per rank per million tets (Vista 1076378)
# about 2.6 GB per rank, 166 GB per node, of 237 GB.
PROB=tpv205
RES=mirror100m
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_mirror_100m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
