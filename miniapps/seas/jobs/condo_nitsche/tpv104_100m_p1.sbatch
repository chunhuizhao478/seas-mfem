#!/bin/bash
#SBATCH --job-name=tpv104_100m_p1
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=19 --ntasks-per-node=10 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=30:00:00
#SBATCH --output=tpv104_100m_p1_%j.out
#SBATCH --error=tpv104_100m_p1_%j.err
# Nitsche dynamic rupture, tpv104 on the 100 m mesh, p = 1, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv104/configs/tpv104_nitsche.toml (Chunhui Zhao, 2026-10-16).
# Submit from miniapps/seas/jobs/condo_nitsche/ AFTER build_verify_condo.sbatch passed at this commit.
# Mesh: tpv102/mesh/tpv102_100m.msh (7,578,859 tets, h_min 10.44 m at two slivers), the mesh of condo
# job 6459052; TPV104 has the TPV102 geometry (tpv104_200m.geo is tpv102_200m.geo apart from comments).
# 10 ranks per node: every rank loads the whole serial mesh and runs METIS on it
# (drivers/nitsche_dyn_driver.cpp:750; MFEM ParMesh partitions on every rank).  At 14 per node with
# --mem=56G this setup failed twice (condo jobs 6490983, 6491016: METIS_PartGraphKway errors and
# SIGBUS during partitioning), although 6459052 (older build b0c0228, TPV102 deck) fit at 14 per node.
# Estimate: TPV205 100 m (6490984, same build, Kelvin-Voigt + Prakash-Clifton, 266 ranks) steps at
# 0.214 s per step on 5.44 M tets, ~1.4x the undamped 6459052 rate per tet.  Scaled to 7.58 M
# tets, 122,120 steps and 190 ranks: about 16 h.  Limit 30 h.
PROB=tpv104
RES=100m
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv102/mesh/tpv102_100m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_condo.sh"
