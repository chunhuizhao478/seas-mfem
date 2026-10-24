#!/bin/bash
#SBATCH -J v2_lts_tpv205_200m_p2
#SBATCH -o v2_lts_tpv205_200m_p2_%j.out
#SBATCH -e v2_lts_tpv205_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 00:30:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131,i617-101
# Fault-writer note rev1 check V2 (Chunhui Zhao, 2026-10-24): as v2_lts_tpv205_200m_p1.sbatch at p = 2,
# the rank layout of LTS gold 1082384 and V1 p2 1095225 (16 nodes x 144 ranks), 0 to 12.0 s, the
# spread variant only (the defaults, W0 + W1; the p = 1 job compares it with rank0), stations compared
# with V1 o2's (1095225, kept on scratch) over the 12 s (meant to be bitwise equal; V1 differs from gold
# 1082384 by the 8.3e-11 rank-0 mesh setup round-off); wall time against V1 o2's 537.5 s.  Then the
# spread file's values against V1 o2's kept file (other chunk lengths, so values only).
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_fault_writer_vista.sbatch passed at this
# commit.  Estimate: about 7 min for the run and 3 min for the comparison of 15 GB of values, with
# setup about 12 min = 3.2 node-h = 1.1 SU.  Requested: 16 nodes x 0.5 h x 0.33 SU = 2.6 SU at most.
TAG=v2
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=2
TFINAL=12.0
# the station reference: V1 o2 (run_m0_vista.sh compares the variants marked gold with it)
GOLD=$SCRATCH/nitsche_profiling/v1_tpv205_200m_p2_1095225/o2
VARIANTS=(
   "spread|tpv205/configs/tpv205_nitsche.toml||gold"
)
FAULT_COMPARE=(
   "@BASE@/spread/fault.vtkhdf|$SCRATCH/nitsche_profiling/v1_tpv205_200m_p2_1095225/o2/fault.vtkhdf|--values-only"
)
source "$SLURM_SUBMIT_DIR/run_m0_vista.sh"
