#!/bin/bash
#SBATCH -J tpv7_100m_p3
#SBATCH -o tpv7_100m_p3_%j.out
#SBATCH -e tpv7_100m_p3_%j.err
#SBATCH -p gg
#SBATCH -N 32
#SBATCH -n 4608
#SBATCH -t 10:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, tpv7 on the 100 m mesh, p = 3, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv7/configs/tpv7_nitsche.toml: the SCEC submission run (Chunhui
# Zhao, 2026-11-21: "Prepare a TPV7 100 m p3 rerun", after the p = 2 run 1098968 missed the small-error
# bar against DRDG3D: median arrival 0.228 s, peak slip rate 9.0 %, x2_12_x3_7.5 +30.8 %).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit
# (it makes tpv6/mesh/tpv6_nitsche_100m.msh, which TPV6 and TPV7 share: 6,199,161 tets on Vista).
# Layout: 32 x 144 ranks (the gg maximum, QOS qgg MaxTRES node=32), 193.7 k tets per node, with
# --dense-fill direct (fills the dense blocks from the element and face matrices, no global CSR; the
# fix that let TPV34 p = 4 run, 1152065).  The old 48-node layout of this script was over the gg limit.
# Memory (estimate, not measured): TPV34 p = 4 with direct fill peaked at 61.1 GiB committed per node
# with 63.4 k tets per node; per tet p = 3 has (60/105)^2 = 0.33 of the p = 4 block size, so about
# 61 GiB per node here, under the 118.8 GiB limit (the p = 3 memory note gave the same ~61 GiB).
# Time (estimate, not measured): the p = 2 run 1098968 took 1:07 on 32 nodes; p = 2 -> p = 3 costs
# 4.6 x in node-hours (TPV31, measured), so about 5.1 h on 32 nodes = 54 SU.
# Stall guard as in TPV34 p = 4.
# Requested: 32 nodes x 10 h x 0.33 SU = 105.6 SU at most.
PROB=tpv7
RES=100m
DECK=tpv7/configs/tpv7_nitsche.toml
MESH=tpv6/mesh/tpv6_nitsche_100m.msh
ORDER=3
TFINAL=12.0
DRIVER_ARGS="--dense-fill direct"
STALL_SAMPLE_AFTER=900
STALL_ABORT=1
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
