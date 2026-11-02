#!/bin/bash
#SBATCH -J sst_gts_tpv104_p1
#SBATCH -o sst_gts_tpv104_p1_%j.out
#SBATCH -e sst_gts_tpv104_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH --ntasks-per-node=144
#SBATCH -t 01:00:00
#SBATCH -A EAR24034
#SBATCH --exclude=i617-131,i617-101
# SeisSol and MFEM speed test on Vista, equal global step (Chunhui Zhao 2026-11-02: damping stays in MFEM,
# "use global time stepping for both cases", equal time steps, no station output): tpv104 on the 200 m mesh,
# 0 to 12.0 s, 4 nodes, one allocation: MFEM p = 1 (Nitsche, central difference, --lts off, Kelvin-Voigt
# and Prakash-Clifton on as the gold deck, 144 ranks per node), then SeisSol o2 (ClusteredLTS = 1,
# FixTimestep = FIX_DT, SS_LAYOUTS ranks per node).  Both take the step FIX_DT; the job fails if either
# log shows another.  No station output in either code.
# Submit from miniapps/seas/jobs/seissol_speed_vista/ AFTER build_speed_test_vista.sbatch passed at this
# commit.  MFEM_RPN and SS_LAYOUTS are placeholders until probe_gts.sbatch has run (a follow-up commit).
# Estimate (inferred from 1114562: MFEM stepping x 1/0.445 without LTS, SeisSol x the
# update ratio): MFEM about 460 s, SeisSol about 1050 s, 2.0 node-h = 0.7 SU; requested 4 x 1 h x 0.33 = 1.3 SU at most.
TAG=gts
PROB=tpv104
MESH=tpv104/mesh/tpv104_200m.msh
PUML=tpv104_200m
ORDER=1
TFINAL=12.0
MFEM_RUN=1
MFEM_RPN=144
SS_LAYOUTS="8"        # probes 1112916, 1114522: 8 ranks per node, 17 + 1 threads each
GTS=1                 # global time stepping in both codes
FIX_DT=2.4509e-4      # both codes' step [s]: MFEM's own CFL step 2.45096e-4 s (1114562), rounded down
NO_STATIONS=1         # no station output in either code
source "$SLURM_SUBMIT_DIR/run_pair_vista.sh"
