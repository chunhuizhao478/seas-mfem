#!/bin/bash
#SBATCH -J tpv26_mirror_200m_p2
#SBATCH -o tpv26_mirror_200m_p2_%j.out
#SBATCH -e tpv26_mirror_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 24
#SBATCH -n 3072
#SBATCH -t 02:00:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# B2 of rheology plan rev 1 (Chunhui Zhao, 2026-10-20): SCEC TPV26 (forced rupture, depth-dependent
# effective stress, linear elastic) through the Nitsche driver on the mirror-symmetric 200 m mesh, p = 2,
# 0 -> 13.0 s, deck tpv26/configs/tpv26_nitsche_kv0p1.toml (Kelvin-Voigt eta_bar 0.1 + Prakash-Clifton, as
# the TPV205/TPV104 golds).  Compared with DRDG3D TPV26 over 0-12 s by the plan's criteria; a gold if it
# passes.  Submit from miniapps/seas/jobs/vista_rheology/ AFTER build_tpv2627_vista.sbatch passed at this
# commit (it builds the driver and writes the mesh and its sha256 to MESH_SHA256_TPV2627_VISTA).
# Mesh (Vista gmsh 4.15.0-git, built-in kernel, build job 1086883): 2,995,976 tets, 46,328 fault faces,
# dt 1.330e-4 s, 97,755 steps to 13 s.
# Estimate from the trunk's 1 s run in LTS job 1082384 (TPV205 200 m p2, 1.89 M tets, 7,394 steps,
# 199 s on 2,304 ranks: 3.3e-5 core-s per tet-step with setup): 2,700 core-h, 0.9 h on 3,072 ranks,
# 24 nodes x 0.9 h x 0.33 = 7 SU; the larger fault adds rank-0 frame output.
# Requested: 24 nodes x 2 h x 0.33 SU = 15.8 SU at most.
# Memory: Vista compute nodes use strict overcommit (mode 2, ratio 50, no swap, read on i617-101 by the
# benchmarks thread, 2026-10-21): a node can commit at most 118.8 GiB, and past that every allocation
# fails.  Job 1087808 (2f4ede1, every rank read the whole serial mesh and partitioned it) stopped at setup
# with std::bad_alloc on node i617-101.  The driver now defaults to --mesh-setup rank0 (from the Nitsche
# profiling branch, bc91ac9): rank 0 reads and partitions the mesh and sends each rank its part.  On TPV205
# 200 m p2, 16 x 144, that took the first node's peak committed memory from 114.1 GiB to 37.6 GiB (profiling
# job 1090566).  Here: 125,000 tets per node, about that run's 118,000, and a serial mesh 1.6x larger on
# rank 0.  The first node's committed memory is logged below.
PROB=tpv26
RES=mirror200m
DECK=tpv26/configs/tpv26_nitsche_kv0p1.toml
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
