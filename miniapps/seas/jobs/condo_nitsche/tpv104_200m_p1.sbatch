#!/bin/bash
#SBATCH --job-name=tpv104_200m_p1
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=19 --ntasks-per-node=14 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=06:00:00
#SBATCH --output=tpv104_200m_p1_%j.out
#SBATCH --error=tpv104_200m_p1_%j.err
# Debug run for the condo Bus errors (Chunhui Zhao, 2026-10-16): TPV104 on the 200 m mesh, p = 1,
# 0 -> 12.0 s, deck tpv104/configs/tpv104_nitsche.toml, on the layout of the failed 100 m jobs
# (19 nodes x 14 ranks, --mem=56G).  TPV104 100 m died with SIGBUS at setup (6490983, 6491016; 6491334
# at 10 per node), TPV205 100 m at t = 4.15 s (6490984).  Only the mesh changes: 2.47 M instead of
# 7.58 M tets, so about a third of the memory per rank.  Same mesh as Vista job 1075919
# (tpv104_200m.msh made by build_verify_vista.sbatch, staged in condo-meshes-20261016), so the stations
# can be compared with that run directly.
# Estimate: 0.26 s per step at 5.44 M tets on 266 ranks (6490984) scaled to 2.47 M tets and 48,961
# steps: about 1.6 h.  Limit 6 h.
PROB=tpv104
RES=200m
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_200m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_condo.sh"
