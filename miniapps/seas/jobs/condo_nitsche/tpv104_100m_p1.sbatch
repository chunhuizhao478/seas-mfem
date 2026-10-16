#!/bin/bash
#SBATCH --job-name=tpv104_100m_p1
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=19 --ntasks-per-node=14 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=30:00:00
#SBATCH --output=tpv104_100m_p1_%j.out
#SBATCH --error=tpv104_100m_p1_%j.err
# Nitsche dynamic rupture, tpv104 on the 100 m mesh, p = 1, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv104/configs/tpv104_nitsche.toml (Chunhui Zhao, 2026-10-16).
# Submit from miniapps/seas/jobs/condo_nitsche/ AFTER build_verify_condo.sbatch passed at this commit.
# Mesh: tpv102/mesh/tpv102_100m.msh (7,578,859 tets, h_min 10.44 m at two slivers), the mesh of condo
# job 6459052; TPV104 has the TPV102 geometry (tpv104_200m.geo is tpv102_200m.geo apart from comments).
# Layout as 6459052: 14 ranks per node because every rank loads the whole serial mesh
# (drivers/nitsche_dyn_driver.cpp:750), 2.44 GB per rank for this mesh; 20 per node was OOM-killed
# (6459047); 6459052 peaked at 41.1 GB per node.
# Estimate: 6459052 (same mesh and layout, Prakash-Clifton only) stepped 122,120 steps in 26,006 s plus
# 948 s for 119 fault frames.  Kelvin-Voigt costs ~2x per step (Frontera 200 m damping runs) and this deck
# writes 600 fault frames: about 52,000 + 4,800 s = 16 h.  Limit 30 h.
PROB=tpv104
RES=100m
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv102/mesh/tpv102_100m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_condo.sh"
