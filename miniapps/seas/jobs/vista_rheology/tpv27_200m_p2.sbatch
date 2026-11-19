#!/bin/bash
#SBATCH -J tpv27_mirror_200m_p2
#SBATCH -o tpv27_mirror_200m_p2_%j.out
#SBATCH -e tpv27_mirror_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 24
#SBATCH -n 3072
#SBATCH -t 02:00:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# B4 of rheology plan rev 1 (Chunhui Zhao, 2026-10-25): SCEC TPV27 (TPV26 with a non-associative
# Drucker-Prager viscoplastic bulk, Duvaut-Lions Tv = 0.03 s) through the Nitsche driver on the
# mirror-symmetric 200 m mesh, p = 2, 0 -> 13.0 s, deck tpv27/configs/tpv27_nitsche_kv0p1.toml (the TPV26
# deck with [rheology.drucker_prager] enabled: the plastic state at the 4 nodes of degree 1 in each tet;
# plastic strain and eta every 0.5 s to <output>/rheology.vtkhdf).  Compared with DRDG3D TPV27 over 0-12 s
# and with TPV26 p = 2 on the same mesh (tpv26_200m_p2.sbatch); a gold if it passes.  Submit from
# miniapps/seas/jobs/vista_rheology/ AFTER build_tpv2627_vista.sbatch passed at this commit (it builds the
# driver and writes the mesh and its sha256 to MESH_SHA256_TPV2627_VISTA).
# Layout as tpv26_200m_p2.sbatch: 24 nodes x 128 ranks, rank-0 mesh setup (the default), the first node's
# committed memory logged (Vista nodes commit at most 118.8 GiB each; TPV26 1094998 peaked at 56.6 GiB).
# The plastic state adds 7 values at 4 nodes per tet, 0.7 GB in all.
# Estimate: TPV26 1094998 ran 97,755 steps in 22 min (12.8 ms/step).  On cloud p = 2 runs (4 km mesh,
# 3 ranks) TPV27 cost the same per step as TPV26 before anything yielded; the 27 rheology.vtkhdf
# snapshots (each about 0.2-0.6 GB) add their write time.  So about 25-40 min, 24 nodes x 0.7 h x 0.33
# = 5 SU at most expected.  Requested: 24 nodes x 2 h x 0.33 SU = 15.8 SU at most.
PROB=tpv27
RES=mirror200m
DECK=tpv27/configs/tpv27_nitsche_kv0p1.toml
MESH=tpv26/mesh/tpv2627_mirror_200m.msh
ORDER=2
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
