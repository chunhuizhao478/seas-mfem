#!/bin/bash
#SBATCH -J tpv205_mirror_200m_p2
#SBATCH -o tpv205_mirror_200m_p2_%j.out
#SBATCH -e tpv205_mirror_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 02:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, tpv205 on the mirror-symmetric 200 m mesh, p = 2, 0 -> 12.0 s, Kelvin-Voigt
# eta_bar 0.1 + Prakash-Clifton as committed in tpv205/configs/tpv205_nitsche.toml (Chunhui Zhao,
# 2026-10-17).  Same fault and stations as tpv205_200m_p1.sbatch; only the mesh is mirror-symmetric.
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch and
# make_mirror_meshes_vista.sbatch passed at this commit.
# Mesh as tpv205_mirror_200m_p1.sbatch: 1,976,500 tets, min 6V/A 30.60 m, about the old mesh's, so dt
# about 1.35e-4 s and 88,700 steps to 12 s as in Vista 1075916 (TPV205 200 m p2 on the old mesh).
# Estimate from 1075916 (1.89 M tets, 88,725 steps, 2,680 s on 2,304 ranks: 3.7e-5 core-s per
# tet-step, whole job): 1,800 core-h, 0.8 h on 2,304 ranks, 16 nodes x 0.8 h x 0.33 = 4.1 SU.
# Requested: 16 nodes x 2 h x 0.33 SU = 10.6 SU at most.
# 144 ranks per node, as 1075916, which peaked at 117 GB per node; scaled by tets about 122 GB here.
PROB=tpv205
RES=mirror200m
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_mirror_200m.msh
ORDER=2
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
