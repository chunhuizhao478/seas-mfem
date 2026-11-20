#!/bin/bash
#SBATCH -J tpv26_mirror_200m_p3
#SBATCH -o tpv26_mirror_200m_p3_%j.out
#SBATCH -e tpv26_mirror_200m_p3_%j.err
#SBATCH -p gg
#SBATCH -N 32
#SBATCH -n 4096
#SBATCH -t 04:00:00
#SBATCH -A EAR24034
#SBATCH --exclude=i617-131
# TPV26 at p = 3, the SCEC submission run (Chunhui Zhao, 2026-11-20: p = 3 for every benchmark but TPV104
# and TPV31; "yes I agree, let's prepare them").  The gold is p = 2 (1099061).  Same deck
# (tpv26/configs/tpv26_nitsche_kv0p1.toml: 12 on-fault and the 6 Part 9 off-fault surface stations), mesh
# and output as tpv26_200m_p2.sbatch, 0 -> 13.0 s; the layout of tpv27_200m_p3.sbatch (same mesh).
# Submit from miniapps/seas/jobs/vista_rheology/ AFTER build_tpv2627_vista.sbatch passed at this commit.
# The account can be changed at submission (sbatch -A EAR26024 ...) after the balances are read.
# Layout: 32 nodes x 128 ranks (93.7 k tets per node), the most the gg QOS allows per job.  Memory as
# estimated in tpv27_200m_p3.sbatch: about 74 GiB per node, upper bound 107 GiB, limit 118.8 GiB (TPV26 is
# elastic, so it has no plastic state on top).
# Estimate (not measured): TPV26 p = 2 ran 26 min on 24 nodes (about 3.4 SU); p2 -> p3 costs 4.6x in
# node-hours (TPV31 1096105 -> 1096106), so about 1.5 h on 32 nodes = 16 SU.
# Requested: 32 nodes x 4 h x 0.33 SU = 42.2 SU at most.
PROB=tpv26
RES=mirror200m
DECK=tpv26/configs/tpv26_nitsche_kv0p1.toml
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
