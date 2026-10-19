#!/bin/bash
#SBATCH -J tpv102_200m_p3
#SBATCH -o tpv102_200m_p3_%j.out
#SBATCH -e tpv102_200m_p3_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2048
#SBATCH -t 04:30:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, tpv102 on the 200 m mesh, p = 3, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton (Chunhui Zhao, 2026-10-19).  Deck tpv102/configs/tpv102_nitsche_pc_p2_kv0p1.toml,
# the deck of the p = 2 gold (Frontera 7995488, which stopped at 9 s); the order comes from ORDER.
# Mesh: tpv104/mesh/tpv104_200m.msh, made by build_verify_vista.sbatch.  tpv104_200m.geo equals
# tpv102_200m.geo apart from comments, so this is the Vista mesh of tpv102_200m.geo; the TPV102 golds
# used the Frontera mesh of the same .geo (2,466,377 tets, sha256 eea25a27...).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Estimate: as tpv104_200m_p3.sbatch (same mesh and ranks), 176,570 steps to 12 s, about 2.5 to 3.1 h
# = 13 to 16 SU.  Requested: 16 nodes x 4.5 h x 0.33 SU = 23.8 SU at most.
PROB=tpv102
DECK=tpv102/configs/tpv102_nitsche_pc_p2_kv0p1.toml
MESH=tpv104/mesh/tpv104_200m.msh
ORDER=3
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
