#!/bin/bash
#SBATCH -J tpv27_mirror_200m_p1
#SBATCH -o tpv27_mirror_200m_p1_%j.out
#SBATCH -e tpv27_mirror_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 8
#SBATCH -n 1024
#SBATCH -t 01:00:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# TPV27 at p = 1 (Chunhui Zhao, 2026-10-25): tpv27_200m_p2.sbatch at order 1, to see the Drucker-Prager
# formulation with the plastic state at one node per tet (degree 0) next to p = 2 (four nodes, degree 1).
# The same deck, mesh and comparison with DRDG3D TPV27 over 0-12 s; there is no TPV26 p = 1 run on this
# mesh, so TPV27 minus TPV26 is measured at p = 2 only.  Submit from miniapps/seas/jobs/vista_rheology/
# AFTER build_tpv2627_vista.sbatch passed at this commit.
# Layout: 8 nodes x 128 ranks (1,024), rank-0 mesh setup (the default), first-node committed memory
# logged.  Memory: the first time step peaks at about 0.04 MiB per tet at p = 1 (TPV31 runs), 14 GiB per
# node here, plus about 18 GiB of setup.
# Estimate (inferred): dt about 2.7e-4 s, 2.1 times the p = 2 dt (the driver's CFL walk at p = 1 and
# p = 2 in the cloud, smallest 6V/A 29 m), so about 47,500 steps to 13 s.  Per tet, a p = 1 step costs
# 0.26 of a p = 2 step in the cloud (TPV26 1000 m p = 1 against 4 km p = 2, 3 ranks), so 3.5e-6 core-s
# per tet-step from TPV26 1094998's 1.36e-5: about 8 min of stepping on 1,024 ranks, plus setup and the
# rheology.vtkhdf writes.  Requested: 8 nodes x 1 h x 0.33 SU = 2.6 SU at most.
PROB=tpv27
RES=mirror200m
DECK=tpv27/configs/tpv27_nitsche_kv0p1.toml
MESH=tpv26/mesh/tpv2627_mirror_200m.msh
ORDER=1
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
