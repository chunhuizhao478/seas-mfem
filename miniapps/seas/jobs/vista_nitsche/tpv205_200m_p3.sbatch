#!/bin/bash
#SBATCH -J tpv205_200m_p3
#SBATCH -o tpv205_200m_p3_%j.out
#SBATCH -e tpv205_200m_p3_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 03:30:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, tpv205 on the 200 m mesh, p = 3, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv205/configs/tpv205_nitsche.toml (Chunhui Zhao, 2026-10-19).
# Same deck, mesh and 16 x 144 ranks as the p = 2 gold (Vista 1075916); only the order differs.
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Estimate: dt = 0.591 x the p = 2 dt (driver C(p) 0.282 / 9 against 0.212 / 4), 150,076 steps to 12 s;
# step cost 4.24 x p = 2 (cloud, TPV205 5000 m, p = 3 against p = 2, both with dense blocks) applied to
# the 9.36 ms per step of Vista job 1079640: 39.7 ms, 1.7 h of steps, about 1.9 h = 30 node-h = 10 SU.
# Requested: 16 nodes x 3.5 h x 0.33 SU = 18.5 SU at most.
# Memory: p = 3 setup takes about 4 x the p = 2 setup in the cloud; the first node peaked at 138 GB of
# 238 GB at p = 2, mostly from every rank reading the serial mesh.  Unmeasured at p = 3:
# run_nitsche_vista.sh logs the first node's peak.
PROB=tpv205
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=3
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
