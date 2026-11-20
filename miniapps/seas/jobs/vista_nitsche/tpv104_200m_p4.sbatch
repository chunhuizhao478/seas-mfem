#!/bin/bash
#SBATCH -J tpv104_200m_p4
#SBATCH -o tpv104_200m_p4_%j.out
#SBATCH -e tpv104_200m_p4_%j.err
#SBATCH -p gg
#SBATCH -N 32
#SBATCH -n 4608
#SBATCH -t 07:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, SCEC TPV104 on the 200 m mesh, p = 4 (about 50 m node spacing, the spacing
# the description asks for), 0 -> 12.0 s, the SCEC submission run (Chunhui Zhao, 2026-11-20: p = 4 for
# TPV104 and TPV31 only, p = 3 for the other benchmarks; "yes I agree, let's prepare them").  Deck as
# committed in tpv104/configs/tpv104_nitsche.toml, with six free-surface off-fault stations and nine on-fault; the
# rupture-time files feed the SCEC contour plot (scripts/export_scec.py).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# The account can be changed at submission (sbatch -A EAR26024 ...) after the balances are read.
# Layout: 32 x 144 ranks, the most nodes one gg job may have (QOS qgg MaxTRES node=32); 2,466,549 tets: 77.1 k per node.
# Memory: --dense-fill direct as in TPV34 p = 4 (1152065: 63.4 k tets per node, peak 61.1 GiB
# committed per node, limit 118.8 GiB); scaled by tets per node this run should stay under the limit
# (estimate).
# Estimate (not measured): the p = 3 gold (16 nodes, 176,571 steps, 6839 s driver wall, 10.1 SU) times
# 4.2 in node-hours (TPV34 p = 3 -> p = 4, measured) on 32 nodes: about 4.0 h = 42 SU.
# Stall guard as in TPV34 p = 4.
# Requested: 32 nodes x 7 h x 0.33 SU = 73.9 SU at most.
PROB=tpv104
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_200m.msh
ORDER=4
TFINAL=12.0
DRIVER_ARGS="--dense-fill direct"
STALL_SAMPLE_AFTER=900
STALL_ABORT=1
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
