#!/bin/bash
#SBATCH -J tpv102_nitsche_100m_p2_pc
#SBATCH -o tpv102_nitsche_100m_p2_pc_%j.out
#SBATCH -e tpv102_nitsche_100m_p2_pc_%j.err
#SBATCH -p normal
#SBATCH -N 48
#SBATCH -n 2304
#SBATCH -t 16:00:00
#SBATCH -A EAR20006
# SCEC TPV102 (RS aging) through the Nitsche dynamic-rupture driver (Chunhui Zhao): p = 2 on the
# 100 m on-fault mesh, Prakash-Clifton (v_star 1 m/s, L 0.01 m), 0 -> 12 s, normal queue.
#   code   $WORK/seas-mfem-nitsche-b0c = fresh checkout of b0c0228 (PC + station copy averaging),
#          built by build_nitsche_b0c.sbatch (7995286)
#   mesh   $WORK/seas-mfem-nitsche-meshes/tpv102_100m.msh
#          sha256 cb6719b1d0fa6fbb5a7dcf1d9f830081aafeeb2480a7c7cb2452220ac5c9f4d6
#          7,578,859 tets, 151,564 fault faces, h_min 10.44 m (2 slivers; next 15.0 m)
#   deck   tpv102_nitsche_pc_100m.toml = the Discovery 100 m deck = tpv102_nitsche_pc_p2.toml (identical)
#   dt     ~4.61e-5 s (0.5 x 0.212 x 10.44 m / (4 x 6000 m/s)), ~260,000 steps; dumps every 434 steps (~0.02 s)
#   sizing 48 ranks/node (mesh load 2.44 GB/rank -> ~117 GB of 192 GB); measured CLX rate (7993430)
#          p=2 2.57e-5 core-s per tet-step -> ~14,100 core-h / 2304 ranks = ~6.1 h (+ output); ~300-350 SU
set -o pipefail
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # one thread per MPI rank
ROOT=$WORK/seas-mfem-nitsche-b0c
MESH=$WORK/seas-mfem-nitsche-meshes/tpv102_100m.msh
DECK=$WORK/seas-mfem-nitsche-jobs/tpv102_nitsche_pc_100m.toml
module load intel/19.1.1 impi/19.0.9 hypre/2.31.0 mumps/5.3 parmetis petsc/3.15 fftw3/3.3.8
export LD_LIBRARY_PATH="$ROOT/extern/hdf5/install/lib:$ROOT/extern/h5z-zfp/install/lib:${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"
cd $ROOT/miniapps/seas
OUT=$SCRATCH/nitsche_tpv102/tpv102_nitsche_100m_p2_$SLURM_JOB_ID
mkdir -p $OUT
echo "commit $(git rev-parse HEAD)  dirty: $(git status --short --untracked-files=no | wc -l) files"
echo "mesh sha256 $(sha256sum $MESH | cut -d" " -f1)"
[ "$(sha256sum $MESH | cut -d" " -f1)" = cb6719b1d0fa6fbb5a7dcf1d9f830081aafeeb2480a7c7cb2452220ac5c9f4d6 ] \
   || { echo "ERROR: mesh checksum mismatch"; exit 3; }
[ -x ./seas_nitsche_dyn_driver ] && [ -x ./seas_test_nitsche_pc ] || { echo "ERROR: binaries not built"; exit 3; }
cp $DECK $OUT/
# pre-flight: the Prakash-Clifton unit test at 4 ranks
ibrun -n 4 -o 0 ./seas_test_nitsche_pc > $OUT/preflight_pc.log 2>&1
rc=$?; grep -E "Summary" $OUT/preflight_pc.log
[ $rc -eq 0 ] && grep -q "0 failed" $OUT/preflight_pc.log || { echo "ERROR: pre-flight seas_test_nitsche_pc failed (rc $rc)"; exit 3; }
T0=$(date +%s)
ibrun ./seas_nitsche_dyn_driver --config $DECK --mesh $MESH \
      --order 2 --tfinal 12.0 --output-dir $OUT --station-every 1 --dump-every 434 --log-every 1000 \
      > $OUT/driver.log 2>&1
rc=$?
echo "driver exit $rc after $(( $(date +%s) - T0 )) s"; tail -3 $OUT/driver.log
du -sh $OUT
exit $rc
