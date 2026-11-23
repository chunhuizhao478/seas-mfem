#!/bin/bash
#SBATCH -J tpv34_200m_p3_check
#SBATCH -o tpv34_200m_p3_check_%j.out
#SBATCH -e tpv34_200m_p3_check_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 00:30:00
#SBATCH -A EAR24034
# Short check of the TPV34 200 m p = 3 run (Chunhui Zhao, TPV34 plan rev 2, step S6, 2026-11-15).
# Job 1140834 (tpv34_200m_p3.sbatch, this layout) ran into its 5 h limit with no step line in
# driver.log and empty station files; the cause is not known.  This job repeats exactly that run
# (same deck, mesh, table, order, 16 x 144 ranks) but only to t = 0.02 s, logs every 10 steps, and
# if no step line appears within 8 min of the start it records what the ranks on the first node are
# doing (run_nitsche_vista.sh stall sampler: CPU, state and call stacks, three samples 60 s apart).
# Not a result run: output is diagnostics only.
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_tpv34_vista.sbatch passed at this
# commit.
# Requested: 16 nodes x 0.5 h x 0.33 SU = 2.6 SU at most.
PROB=tpv34
DECK=tpv34/configs/tpv34_nitsche.toml
MESH=tpv34/mesh/tpv34_nitsche_200m.msh
TABLE=$WORK/tpv34_material/tpv34_cvmh_table.h5
ORDER=3
TFINAL=0.02
RES=200m_p3check
LOG_EVERY=10
STALL_SAMPLE_AFTER=480
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
