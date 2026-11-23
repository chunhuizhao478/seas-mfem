#!/bin/bash
#SBATCH -J tpv34_200m_p1
#SBATCH -o tpv34_200m_p1_%j.out
#SBATCH -e tpv34_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 01:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, SCEC TPV34 on the 200 m mesh, p = 1, 0 -> 20.0 s, Kelvin-Voigt eta_bar 0.1
# + Prakash-Clifton as committed in tpv34/configs/tpv34_nitsche.toml; 3-D CVM-H material table, 35
# on-fault and 56 off-fault stations (Chunhui Zhao, TPV34 plan rev 2, step S6, 2026-11-12).  Rank
# layout of the TPV31 p = 1 gold (4 x 144).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_tpv34_vista.sbatch passed at this
# commit.
# Estimate (not measured: the TPV34 step is set by the table's fastest element): about 80,000 steps (TPV31 p = 1 gold 1085893: 58,551 steps to 15 s in 9:41 of job wall on this layout; 20 s and a similar step give about 1.4 x): about 15 min = 0.3 SU.
# Requested: 4 nodes x 1 h x 0.33 SU = 1.3 SU at most.
PROB=tpv34
DECK=tpv34/configs/tpv34_nitsche.toml
MESH=tpv34/mesh/tpv34_nitsche_200m.msh
TABLE=$WORK/tpv34_material/tpv34_cvmh_table.h5
ORDER=1
TFINAL=20.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
