#!/bin/bash
#SBATCH -J tpv31_200m_p4
#SBATCH -o tpv31_200m_p4_%j.out
#SBATCH -e tpv31_200m_p4_%j.err
#SBATCH -p gg
#SBATCH -N 32
#SBATCH -n 4608
#SBATCH -t 07:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, SCEC TPV31 on the 200 m mesh, p = 4 (about 50 m node spacing, the spacing
# the description asks for), 0 -> 15.0 s, the SCEC submission run (Chunhui Zhao, 2026-11-20: p = 4 for
# TPV104 and TPV31 only, p = 3 for the other benchmarks; "yes I agree, let's prepare them").  Deck as
# committed in tpv31/configs/tpv31_nitsche.toml, with 18 off-fault and 30 on-fault stations; the
# rupture-time files feed the SCEC contour plot (scripts/export_scec.py).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# The account can be changed at submission (sbatch -A EAR26024 ...) after the balances are read.
# Layout: 32 x 144 ranks, the most nodes one gg job may have (QOS qgg MaxTRES node=32); 2,029,790 tets: 63.4 k per node.
# Memory: --dense-fill direct as in TPV34 p = 4 (1152065: 63.4 k tets per node, peak 61.1 GiB
# committed per node, limit 118.8 GiB); scaled by tets per node this run should stay under the limit
# (estimate).
# Estimate (not measured): the p = 3 gold (16 nodes, 211,157 steps, 6585 s driver wall, 9.7 SU) times
# 4.2 in node-hours (TPV34 p = 3 -> p = 4, measured) on 32 nodes: about 3.8 h = 41 SU.
# Stall guard as in TPV34 p = 4.
# Requested: 32 nodes x 7 h x 0.33 SU = 73.9 SU at most.
PROB=tpv31
DECK=tpv31/configs/tpv31_nitsche.toml
MESH=tpv31/mesh/tpv31_nitsche_200m.msh
ORDER=4
TFINAL=15.0
DRIVER_ARGS="--dense-fill direct"
STALL_SAMPLE_AFTER=900
STALL_ABORT=1
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
