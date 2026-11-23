#!/bin/bash
#SBATCH -J tpv34_200m_p3
#SBATCH -o tpv34_200m_p3_%j.out
#SBATCH -e tpv34_200m_p3_%j.err
#SBATCH -p gg
#SBATCH -N 24
#SBATCH -n 3456
#SBATCH -t 05:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, SCEC TPV34 on the 200 m mesh, p = 3, 0 -> 20.0 s, Kelvin-Voigt eta_bar 0.1
# + Prakash-Clifton as committed in tpv34/configs/tpv34_nitsche.toml; 3-D CVM-H material table, 35
# on-fault and 56 off-fault stations (Chunhui Zhao, TPV34 plan rev 2, step S6, 2026-11-12).  Rank
# layout: 24 x 144 ranks (84.6k tets per node).  The first try, job 1140834 on the TPV31 p = 3 gold
# layout (16 x 144, 126.9k tets per node), reached its 5 h limit with no step line; the cause is not
# known.  Its first node peaked at 98.3 of 118.8 GiB committed (the other 15 nodes were not sampled),
# and this layout cuts the per-node memory by a third (Chunhui Zhao, 2026-11-16: "it might just be
# memory, you need to request more nodes").
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_tpv34_vista.sbatch passed at this
# commit.
# Estimate (not measured: the TPV34 step is set by the table's fastest element): about 290,000 steps;
# the TPV31 p = 3 gold 1096106 did 211,157 steps in 1:50:08 on 16 nodes, so at most about 2.6 h here.
# Requested: 24 nodes x 5 h x 0.33 SU = 39.6 SU at most.
PROB=tpv34
DECK=tpv34/configs/tpv34_nitsche.toml
MESH=tpv34/mesh/tpv34_nitsche_200m.msh
TABLE=$WORK/tpv34_material/tpv34_cvmh_table.h5
ORDER=3
TFINAL=20.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
