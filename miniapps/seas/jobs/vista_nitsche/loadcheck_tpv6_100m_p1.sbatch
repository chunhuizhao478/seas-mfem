#!/bin/bash
#SBATCH -J loadcheck_tpv6_100m_p1
#SBATCH -o loadcheck_tpv6_100m_p1_%j.out
#SBATCH -e loadcheck_tpv6_100m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 8
#SBATCH -n 384
#SBATCH -t 00:30:00
#SBATCH -A EAR24034
# Memory check on the run layout of tpv6_100m_p1.sbatch and tpv7_100m_p1.sbatch (Chunhui Zhao,
# 2026-10-20): the same nodes, ranks, deck, mesh and order, to t = 0.005 s (about 40 steps) instead
# of 12 s.  TPV7 has the same mesh and layout as TPV6, so this check covers both.
# Why: build_verify_vista.sbatch's one-node load check of this mesh (job 1083417) ran out of memory
# after the driver's setup had ended, with the whole problem on one node; the runs spread it over 8.
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Pass: driver exit 0, its "=== done" line, and a first-node peak (rank 0's node, sampled as in
# every run) of at most 214 GB, 90 % of the node's 237.6 GB.  Output in
# $SCRATCH/nitsche_benchmarks/tpv6_100m_loadcheck_p1_<job>.
# Requested: 8 nodes x 0.5 h x 0.33 SU = 1.3 SU at most.
PROB=tpv6
RES=100m_loadcheck
DECK=tpv6/configs/tpv6_nitsche.toml
MESH=tpv6/mesh/tpv6_nitsche_100m.msh
ORDER=1
TFINAL=0.005
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
