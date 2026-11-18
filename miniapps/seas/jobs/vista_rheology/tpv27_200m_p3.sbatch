#!/bin/bash
#SBATCH -J tpv27_mirror_200m_p3
#SBATCH -o tpv27_mirror_200m_p3_%j.out
#SBATCH -e tpv27_mirror_200m_p3_%j.err
#SBATCH -p gg
#SBATCH -N 36
#SBATCH -n 4608
#SBATCH -t 04:00:00
#SBATCH -A EAR24034
#SBATCH --exclude=i617-131
# TPV27 at p = 3 (Chunhui Zhao, 2026-11-02): the convergence check after TPV27 p = 2 (1114571) and p = 1
# (1114572), whose sigma_n after rupture sits 0.05 to 1.5 MPa above DRDG3D at p = 2 and changes sign at
# several stations from p = 1 to p = 2.  Same deck (tpv27/configs/tpv27_nitsche_kv0p1.toml: the plastic
# state at the 10 nodes of degree 2 in each tet), mesh and output as tpv27_200m_p2.sbatch, 0 -> 13.0 s.
# Submit from miniapps/seas/jobs/vista_rheology/ AFTER build_tpv2627_vista.sbatch passed at this commit.
# Layout: 36 nodes x 128 ranks (83.3 k tets per node).  Memory: the step-1 operator assembly peaks at about
# 0.968 MiB per tet at p = 3 (cloud VmData; TPV31 p = 3 on Vista at 127 k tets per node peaked at 115.2 GiB
# committed), so about 79 GiB plus the rank-0 mesh setup's ~18 GiB, about 97 GiB per node, under the 118.8 GiB
# limit (24 nodes would need about 137 GiB).  The plastic state (7 values at 10 nodes per tet) adds 1.7 GB in all.
# Estimate: dt scales by C(3)/9 / (C(2)/4) = 0.591, so about 203,000 steps; p = 3 costs 2.77 x p = 2 per
# step at the same layout (TPV31 1096106 vs 1096105) and TPV27 p = 2 ran 17.1 ms/step on 24 nodes, so about
# 32 ms/step on 36 nodes and 1.8 h, 36 nodes x 2 h x 0.33 = 24 SU expected.  Requested: 36 nodes x 4 h x
# 0.33 SU = 47.5 SU at most, on EAR24034 (the p = 3 account).
PROB=tpv27
RES=mirror200m
DECK=tpv27/configs/tpv27_nitsche_kv0p1.toml
MESH=tpv26/mesh/tpv2627_mirror_200m.msh
ORDER=3
TFINAL=13.0
# first-node memory, 2 s samples (kB): MemTotal MemAvailable CommitLimit Committed_AS
MEMLOG=$SCRATCH/nitsche_benchmarks/${PROB}_${RES}_p${ORDER}_${SLURM_JOB_ID}_meminfo_kB.log
mkdir -p $SCRATCH/nitsche_benchmarks
( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable|^CommitLimit|^Committed_AS/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $MEMLOG &
MON=$!
( source "$SLURM_SUBMIT_DIR/../vista_nitsche/run_nitsche_vista.sh" )
rc=$?
kill $MON 2>/dev/null; wait $MON 2>/dev/null
awk '{u=($2-$3)/1048576; if(u>m)m=u; if($5/1048576>c)c=$5/1048576; t=$2/1048576; l=$4/1048576} END{printf "first node memory: peak in use %.1f GiB of %.1f GiB, peak committed %.1f GiB of the %.1f GiB limit (2 s samples, %s)\n", m, t, c, l, FILENAME}' $MEMLOG
exit $rc
