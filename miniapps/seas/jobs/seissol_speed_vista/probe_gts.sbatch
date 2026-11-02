#!/bin/bash
#SBATCH -J sst_probe_gts
#SBATCH -o sst_probe_gts_%j.out
#SBATCH -e sst_probe_gts_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH --ntasks-per-node=144
#SBATCH -t 01:00:00
#SBATCH -A EAR24034
#SBATCH --exclude=i617-131,i617-101
# SeisSol and MFEM speed test on Vista, equal global step (Chunhui Zhao 2026-11-02): probe before the nine
# run_gts_* pairs.  TPV205 p = 2 / o3 on the 200 m mesh, 16 nodes (as run_gts_tpv205_p2), 0 to 2.0 s,
# global time stepping in both codes at the one step FIX_DT (MFEM's own CFL step 1.35251e-4 s in
# 1114560, rounded down), no station output.  It shows
#   1. that both codes take FIX_DT (the job fails if either log shows another step) and that both accept
#      the edited deck copies (SeisSol FixTimestep, ClusteredLTS = 1, no pickpoints or receivers; MFEM
#      dt_max = FIX_DT, no [output.stations]) and write no station files;
#   2. MFEM's best layout: 144, 96 and 72 ranks per node (it has never been probed; 72 per node computed
#      1.48x faster at p = 1 in S1 1102100);
#   3. SeisSol's best layout under global stepping: 8 and 16 ranks per node (probed so far under LTS only).
# The fastest of each picks MFEM_RPN and SS_LAYOUTS for the run_gts_* jobs (a follow-up commit).
# Submit from miniapps/seas/jobs/seissol_speed_vista/ AFTER build_speed_test_vista.sbatch passed at this
# commit.
# Estimate (inferred from 1114560 on 16 nodes, scaled to 2 s): MFEM about 125 s per layout (3 runs),
# SeisSol about 140 s per layout (2 runs), plus about 90 s setup each: about 20 min = 5.6 node-h = 1.9 SU;
# requested 16 x 1 h x 0.33 = 5.3 SU at most.
TAG=probe_gts
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
PUML=tpv2053d_200m
ORDER=2
TFINAL=2.0
MFEM_RUN=1
MFEM_RPN=144
MFEM_LAYOUTS="144 96 72"
SS_LAYOUTS="8 16"
GTS=1
FIX_DT=1.3525e-4
NO_STATIONS=1
source "$SLURM_SUBMIT_DIR/run_pair_vista.sh"
