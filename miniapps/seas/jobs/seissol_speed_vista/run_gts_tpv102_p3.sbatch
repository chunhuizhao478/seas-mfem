#!/bin/bash
#SBATCH -J sst_gts_tpv102_p3
#SBATCH -o sst_gts_tpv102_p3_%j.out
#SBATCH -e sst_gts_tpv102_p3_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH --ntasks-per-node=144
#SBATCH -t 05:00:00
#SBATCH -A EAR24034
#SBATCH --exclude=i617-131,i617-101
# SeisSol and MFEM speed test on Vista, equal global step (Chunhui Zhao 2026-11-02: damping stays in MFEM,
# "use global time stepping for both cases", equal time steps, no station output): tpv102 on the 200 m mesh,
# 0 to 12.0 s, 16 nodes, one allocation: MFEM p = 3 (Nitsche, central difference, --lts off, Kelvin-Voigt
# and Prakash-Clifton on as the gold deck, 128 ranks per node), then SeisSol o4 (ClusteredLTS = 1,
# FixTimestep = FIX_DT, SS_LAYOUTS ranks per node).  Both take the step FIX_DT; the job fails if either
# log shows another.  No station output in either code.
# Submit from miniapps/seas/jobs/seissol_speed_vista/ AFTER build_speed_test_vista.sbatch passed at this
# commit.  MFEM_RPN and SS_LAYOUTS are placeholders until probe_gts.sbatch has run (a follow-up commit).
# Estimate (inferred from 1114567: MFEM stepping x 1/0.445 without LTS, SeisSol x the
# update ratio): MFEM about 6320 s, SeisSol about 2890 s, 42.3 node-h = 13.9 SU; requested 16 x 5 h x 0.33 = 26.4 SU at most.
TAG=gts
PROB=tpv102
MESH=tpv104/mesh/tpv104_200m.msh
PUML=tpv104_200m
ORDER=3
TFINAL=12.0
MFEM_RUN=1
MFEM_RPN=128
SS_LAYOUTS="8"        # probes 1112916, 1114522: 8 ranks per node, 17 + 1 threads each
GTS=1                 # global time stepping in both codes
FIX_DT=6.7961e-5      # both codes' step [s]: MFEM's own CFL step 6.79616e-5 s (1114567), rounded down
NO_STATIONS=1         # no station output in either code
source "$SLURM_SUBMIT_DIR/run_pair_vista.sh"
