#!/bin/bash
#SBATCH --job-name=tpv102_100m_p1_pc
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=19 --ntasks-per-node=14 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch1/zhaochun/nitsche-100m-20261015/logs/tpv102_100m_p1_pc_%j.out
# SCEC TPV102 (RS aging), Nitsche driver (Chunhui Zhao): p = 1 on the 100 m on-fault mesh,
# Prakash-Clifton (v_star 1 m/s, L 0.01 m), 0 -> 12 s, on the SCEC condo nodes
# (Xeon E5-2640 v4, 20 cores, 63 GB).
#   code   seas-mfem b0c0228 built on Discovery (/scratch1/zhaochun/nitsche-p3-20261015/seas-mfem);
#          the pre-flight below checks that it runs on these nodes
#   mesh   tpv102_100m.msh (gmsh default optimizer; the Netgen pass crashed), h_min 10.44 m from 2 slivers (next 15.0 m), sha256 cb6719b1d0fa6fbb5a7dcf1d9f830081aafeeb2480a7c7cb2452220ac5c9f4d6
#   deck   tpv102_nitsche_pc_100m.toml = the 7994410 deck unchanged (order via --order 1)
#   dt     driver-computed, ~9.83e-5 s = 0.5 x 0.113 x 10.44 m / 6000 m/s (~122,000 steps); dumps every 203 steps (~0.02 s); stations every step
#   sizing ESTIMATE: 200 m p=1 (7994409) 20 CLX core-s/step x 3.07 (tets) = 61 core-s/step;
#          Broadwell ~0.6x CLX per core (inferred) -> ~0.39 s/step on 266 ranks -> ~13 h of the 24 h limit
#   layout 14 ranks/node, not 20: 6459047 (20/node) was OOM-killed at mesh load (every rank reads the
#          serial mesh + METIS: 48.9 GB MaxRSS for 20 ranks = 2.44 GB/rank) -> ~34 GB/node at load; ppr:7:package crashes prterun silently here (maptest 6459051), ppr:14:node maps
set -o pipefail
module purge; module load gcc/13.3.0 openmpi/5.0.5
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # one thread per MPI rank
export HDF5_PLUGIN_PATH=/scratch1/zhaochun/nitsche-discovery-20261014/extern/h5z-zfp/install/plugin
P=/scratch1/zhaochun/nitsche-100m-20261015
ROOT=/scratch1/zhaochun/nitsche-p3-20261015/seas-mfem
MESH=$P/mesh/tpv102_100m.msh
DECK=$P/tpv102_nitsche_pc_100m.toml
OUT=$P/runs/tpv102_nitsche_100m_p1_pc_$SLURM_JOB_ID
mkdir -p $OUT
cd $ROOT/miniapps/seas
echo "host $(hostname)  cpu $(grep -m1 "model name" /proc/cpuinfo)"
echo "commit $(git rev-parse HEAD)  dirty: $(git status --short --untracked-files=no | wc -l) files"
[ "$(sha256sum $MESH | cut -d" " -f1)" = cb6719b1d0fa6fbb5a7dcf1d9f830081aafeeb2480a7c7cb2452220ac5c9f4d6 ] || { echo "ERROR: mesh checksum mismatch"; exit 2; }
[ -x ./seas_nitsche_dyn_driver ] && [ -x ./seas_test_nitsche_pc ] || { echo "ERROR: binaries not reachable"; exit 2; }
cp $DECK $OUT/
# pre-flight: the binary runs on these nodes and the Prakash-Clifton unit test passes
mpirun -np 4 --bind-to core ./seas_test_nitsche_pc > $OUT/preflight_pc.log 2>&1
rc=$?; grep -E "Summary" $OUT/preflight_pc.log
[ $rc -eq 0 ] || { echo "ERROR: pre-flight seas_test_nitsche_pc failed (rc $rc); see $OUT/preflight_pc.log"; exit 3; }
T0=$(date +%s)
mpirun -np $SLURM_NTASKS --map-by ppr:14:node --bind-to core \
       ./seas_nitsche_dyn_driver --config $DECK --mesh $MESH \
       --order 1 --tfinal 12.0 --output-dir $OUT --station-every 1 --dump-every 203 --log-every 1000 \
       > $OUT/driver.log 2>&1
rc=$?
echo "driver exit $rc after $(( $(date +%s) - T0 )) s"; tail -3 $OUT/driver.log
du -sh $OUT
exit $rc
