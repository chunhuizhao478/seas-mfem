#!/bin/bash
#SBATCH -J tpv26_mirror_200m_p2
#SBATCH -o tpv26_mirror_200m_p2_%j.out
#SBATCH -e tpv26_mirror_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 24
#SBATCH -n 3072
#SBATCH -t 02:00:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# B2 of rheology plan rev 1 (Chunhui Zhao, 2026-10-20): SCEC TPV26 (forced rupture, depth-dependent
# effective stress, linear elastic) through the Nitsche driver on the mirror-symmetric 200 m mesh, p = 2,
# 0 -> 13.0 s, deck tpv26/configs/tpv26_nitsche_kv0p1.toml (Kelvin-Voigt eta_bar 0.1 + Prakash-Clifton, as
# the TPV205/TPV104 golds).  Compared with DRDG3D TPV26 over 0-12 s by the plan's criteria; a gold if it
# passes.  Submit from miniapps/seas/jobs/vista_rheology/ AFTER build_tpv2627_vista.sbatch passed at this
# commit (it builds the driver and writes the mesh and its sha256 to MESH_SHA256_TPV2627_VISTA).
# Mesh (cloud, gmsh 4.15.2, built-in kernel; Vista's can differ slightly): 2,996,482 tets, 46,334 fault
# faces, min 6V/A 29.94 m, so dt about 1.32e-4 s and 98,400 steps to 13 s (TPV205 mirror 200 m: 30.60 m,
# 1.35e-4 s).
# Estimate from the trunk's 1 s run in LTS job 1082384 (TPV205 200 m p2, 1.89 M tets, 7,394 steps,
# 199 s on 2,304 ranks: 3.3e-5 core-s per tet-step with setup): 2,700 core-h, 0.9 h on 3,072 ranks,
# 24 nodes x 0.9 h x 0.33 = 7 SU; the larger fault adds rank-0 frame output.
# Requested: 24 nodes x 2 h x 0.33 SU = 15.8 SU at most.
# 128 ranks per node: 125,000 tets per node, about the TPV205 gold's 118,000 (first node 138 GB of 238 GB
# with the trunk code), while every rank still reads the whole serial mesh at startup (1.6x the TPV205 mesh);
# the load check of build_tpv2627_vista.sbatch gives the startup memory to confirm this before submitting.
PROB=tpv26
RES=mirror200m
DECK=tpv26/configs/tpv26_nitsche_kv0p1.toml
MESH=tpv26/mesh/tpv2627_mirror_200m.msh
ORDER=2
TFINAL=13.0
source "$SLURM_SUBMIT_DIR/../vista_nitsche/run_nitsche_vista.sh"
