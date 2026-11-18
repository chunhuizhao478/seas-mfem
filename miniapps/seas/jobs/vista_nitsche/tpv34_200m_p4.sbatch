#!/bin/bash
#SBATCH -J tpv34_200m_p4
#SBATCH -o tpv34_200m_p4_%j.out
#SBATCH -e tpv34_200m_p4_%j.err
#SBATCH -p gg
#SBATCH -N 32
#SBATCH -n 4608
#SBATCH -t 08:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, SCEC TPV34 on the 200 m mesh, p = 4, 0 -> 20.0 s, Kelvin-Voigt eta_bar 0.1
# + Prakash-Clifton as committed in tpv34/configs/tpv34_nitsche.toml; 3-D CVM-H material table, 35
# on-fault and 56 off-fault stations (Chunhui Zhao, 2026-11-18: "prepare the TPV34 p4 run").
# Layout: 32 x 144 ranks, the most nodes one gg job may have (QOS qgg MaxTRES node=32, measured
# 2026-11-18); 63.4k tets per node.
# Memory: with the default fill the first time step assembles K and K_eta next to the dense blocks.
# The p = 3 run 1147574 (24 nodes, 84.6k tets per node) peaked at 75.7 GiB committed on its first
# node and settled at 43.5 GiB, 21.8 GiB of it there before the first step; scaled by the element
# DOFs squared ((35/20)^2) and the tets per node, p = 4 here would peak near 145 GiB against the
# node limit of 118.8 GiB (estimate).  So this job runs with --dense-fill direct, which fills the
# blocks from the element and face matrices and never assembles K: the same block entries (bitwise,
# make test-dense-fill-direct and step 5b of build_verify_tpv34_vista.sbatch), and a peak near the
# settled value, about 72 GiB per node (estimate, not measured at p = 4).
# Stall guard: when driver.log has no step line 15 min after the start, run_nitsche_vista.sh
# records what the first node's ranks are doing three times and then cancels the job (the 16-node
# p = 3 job 1140834 sat at its first steps for 5 h; cause not known).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_tpv34_vista.sbatch passed at this
# commit.
# Estimate (not measured): dt = 4.67e-5 s (p = 3 dt x C(4)/C(3) x 9/16), about 428,000 steps; the
# p = 2 and p = 3 runs cost 172 and 534 node-ms per step, the same growth gives about 1,340 at
# p = 4 = 42 ms per step on 32 nodes: about 5.0 h = 53 SU.
# Requested: 32 nodes x 8 h x 0.33 SU = 84.5 SU at most.
PROB=tpv34
DECK=tpv34/configs/tpv34_nitsche.toml
MESH=tpv34/mesh/tpv34_nitsche_200m.msh
TABLE=$WORK/tpv34_material/tpv34_cvmh_table.h5
ORDER=4
TFINAL=20.0
DRIVER_ARGS="--dense-fill direct"
STALL_SAMPLE_AFTER=900
STALL_ABORT=1
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
