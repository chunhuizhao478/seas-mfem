#!/bin/bash
#SBATCH -J tpv102_nitsche_200m_p2_pc_kv0p1
#SBATCH -o tpv102_nitsche_200m_p2_pc_kv0p1_%j.out
#SBATCH -e tpv102_nitsche_200m_p2_pc_kv0p1_%j.err
#SBATCH -p development
#SBATCH -N 40
#SBATCH -n 2240
#SBATCH -t 02:00:00
#SBATCH -A EAR20006
# SCEC TPV102 (RS aging) through the Nitsche dynamic-rupture driver, plan N-07 / P7
# (Chunhui Zhao).  p = 2 on the exact 200 m Gmsh mesh, 0 -> 9 s (Chunhui Zhao 2026-10-15 11:05: "run to 9s is enough").
#   mesh   $WORK/seas-mfem-nitsche-meshes/tpv102_200m.msh
#          sha256 eea25a27813dafd60a5cb5bc3bf1c8de4674e9b07858c3f9b7b8d159fa5be219
#          2,466,377 tets, 37,930 fault faces, h_min 26.81 m
#   dt     0.11842 ms (0.5 x C(p) x h_min / (p^2 c_p)), 101337 steps to 12 s
#   dumps  every 169 steps (~0.02 s of model time); stations every step
#   sizing development queue (40 nodes/user, 2 h max); measured CLX rate, smoke job 7993430: p=2 2.57e-5 core-s per tet-step -> 1786 core-h / 2016 = 0.89 h
#   damped  Kelvin-Voigt eta_bar 0.1, node scaling, h_e = 6V/A (kv_eta/kv_scale/kv_length, the ONLY deck
#           change vs the undamped PC reference 7994410); cost x1.5/step (Mac gate f).  Rate 3.5e-5
#           core-s/tet-step (7994410 steps) x 2,466,377 tets x 101337 steps x 1.5 = 1.31e7 core-s / 2240
#           ranks = 5850 s = 1.63 h; from the measured 4717 s driver of 7994410 (2016 ranks): x1.5 x
#           2016/2240 = 6370 s = 1.77 h.  40 nodes (dev max) instead of 36 to fit in 2 h.
#           tfinal 9 s by Chunhui Zhao's choice: 12 s at the measured x1.9 damped cost exceeds the 2 h dev limit;
#           tfinal 9 s: 4717 s x 1.9 x 2016/2240 x 9/12 = ~6050 s = 1.68 h.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # one thread per MPI rank
ROOT=$WORK/seas-mfem-nitsche-damp   # b0c0228 + kv_node.patch (uncommitted), Prakash-Clifton + Kelvin-Voigt
MESH=$WORK/seas-mfem-nitsche-meshes/tpv102_200m.msh
module load intel/19.1.1 impi/19.0.9 hypre/2.31.0 mumps/5.3 parmetis petsc/3.15 fftw3/3.3.8
export LD_LIBRARY_PATH="$ROOT/extern/hdf5/install/lib:$ROOT/extern/h5z-zfp/install/lib:${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"
cd $ROOT/miniapps/seas
OUT=$SCRATCH/nitsche_tpv102/tpv102_nitsche_200m_p2_pc_kv0p1_$SLURM_JOB_ID
mkdir -p $OUT
echo "commit $(git rev-parse HEAD)  dirty: $(git status --short --untracked-files=no | wc -l) files"
echo "mesh sha256 $(sha256sum $MESH | cut -d' ' -f1)"
[ "$(sha256sum $MESH | cut -d' ' -f1)" = eea25a27813dafd60a5cb5bc3bf1c8de4674e9b07858c3f9b7b8d159fa5be219 ] \
   || { echo "ERROR: mesh checksum mismatch"; exit 2; }
[ -x ./seas_nitsche_dyn_driver ] || { echo "ERROR: seas_nitsche_dyn_driver not built"; exit 2; }
grep -q 'kv_length", "inscribed"' drivers/nitsche_dyn_driver.cpp || { echo "ERROR: not the KV code"; exit 2; }
DECK=$ROOT/damp_jobs/tpv102_nitsche_pc_p2_kv0p1.toml   # reference deck of 7994409/7994410 + [numerics] kv keys
cp $DECK $OUT/
T0=$(date +%s)
# --station-copies first: the undamped references 7994409/7994410 ran ec092fa+pc.patch, which writes single-copy stations; b0c0228 defaults to averaged copies, so match the reference.
ibrun ./seas_nitsche_dyn_driver --config $DECK --mesh $MESH \
      --order 2 --tfinal 9.0 --output-dir $OUT --station-every 1 --dump-every 169 --log-every 1000 --station-copies first \
      > $OUT/driver.log 2>&1
rc=$?
echo "driver exit $rc after $(( $(date +%s) - T0 )) s"; tail -3 $OUT/driver.log
# post-processing is run BY HAND after the job
du -sh $OUT
exit $rc
