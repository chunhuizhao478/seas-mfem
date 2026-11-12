#!/bin/bash
#SBATCH -J tpv34_200m_p2
#SBATCH -o tpv34_200m_p2_%j.out
#SBATCH -e tpv34_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 01:30:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, SCEC TPV34 on the 200 m mesh, p = 2, 0 -> 20.0 s, Kelvin-Voigt eta_bar 0.1
# + Prakash-Clifton as committed in tpv34/configs/tpv34_nitsche.toml; 3-D CVM-H material table, 35
# on-fault and 56 off-fault stations (Chunhui Zhao, TPV34 plan rev 2, step S6, 2026-11-12).  Rank
# layout of the TPV31 p = 2 gold (16 x 144).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_tpv34_vista.sbatch passed at this
# commit.
# Estimate (not measured: the TPV34 step is set by the table's fastest element): about 170,000 steps (TPV31 p = 2 gold 1096105: 124,835 steps to 15 s in 24:15 of job wall on this layout; about 1.4 x): about 35 min = 3.1 SU.
# Requested: 16 nodes x 1.5 h x 0.33 SU = 7.9 SU at most.
PROB=tpv34
DECK=tpv34/configs/tpv34_nitsche.toml
MESH=tpv34/mesh/tpv34_nitsche_200m.msh
TABLE=$WORK/tpv34_material/tpv34_cvmh_table.h5
ORDER=2
TFINAL=20.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
