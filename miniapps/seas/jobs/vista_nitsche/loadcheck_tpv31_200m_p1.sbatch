#!/bin/bash
#SBATCH -J loadcheck_tpv31_200m_p1
#SBATCH -o loadcheck_tpv31_200m_p1_%j.out
#SBATCH -e loadcheck_tpv31_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 00:30:00
#SBATCH -A EAR24034
# Memory check on the run layout of tpv31_200m_p1.sbatch (Chunhui Zhao, 2026-10-20): the same
# nodes, ranks, deck, mesh and order, to t = 0.005 s (about 10 steps) instead of 15 s.
# Why: build_verify_vista.sbatch's one-node load check of this mesh (p = 2, job 1083417) ran out of
# memory after the driver's setup had ended, with the whole problem on one node; the run spreads
# it over 4.
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Pass: driver exit 0, its "=== done" line, and a first-node peak (rank 0's node, sampled as in
# every run) of at most 214 GB, 90 % of the node's 237.6 GB.  Output in
# $SCRATCH/nitsche_benchmarks/tpv31_200m_loadcheck_p1_<job>.
# Requested: 4 nodes x 0.5 h x 0.33 SU = 0.7 SU at most.
PROB=tpv31
RES=200m_loadcheck
DECK=tpv31/configs/tpv31_nitsche.toml
MESH=tpv31/mesh/tpv31_nitsche_200m.msh
ORDER=1
TFINAL=0.005
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
