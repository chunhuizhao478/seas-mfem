#!/bin/bash
#SBATCH -J tpv102_nitsche_200m_p1_pc_kv0p1
#SBATCH -o tpv102_nitsche_200m_p1_pc_kv0p1_%j.out
#SBATCH -e tpv102_nitsche_200m_p1_pc_kv0p1_%j.err
#SBATCH -p development
#SBATCH -N 4
#SBATCH -n 224
#SBATCH -t 02:00:00
#SBATCH -A EAR20006
# SCEC TPV102 (RS aging) through the Nitsche dynamic-rupture driver, plan N-07 / P7
# (Chunhui Zhao).  p = 1 on the exact 200 m Gmsh mesh, 0 -> 12 s.
#   mesh   $WORK/seas-mfem-nitsche-meshes/tpv102_200m.msh
#          sha256 eea25a27813dafd60a5cb5bc3bf1c8de4674e9b07858c3f9b7b8d159fa5be219
#          2,466,377 tets, 37,930 fault faces, h_min 26.81 m
#   dt     0.25248 ms (0.5 x C(p) x h_min / (p^2 c_p)), 47530 steps to 12 s
#   dumps  every 79 steps (~0.02 s of model time); stations every step
#   sizing development queue (40 nodes/user, 2 h max); measured CLX rate, smoke job 7993430: p=1 5.54e-6 core-s per tet-step -> 180 core-h / 224 = 0.81 h
#   damped  Kelvin-Voigt eta_bar 0.1, node scaling, h_e = 6V/A (kv_eta/kv_scale/kv_length, the ONLY change
#           vs the undamped PC reference 7994409 besides the code); cost x1.2/step (Mac gate f) on the
#           measured 4244 s driver of 7994409 -> ~5100 s = 1.42 h; 2 h walltime
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # one thread per MPI rank
ROOT=$WORK/seas-mfem-nitsche-damp   # b0c0228 + kv_node.patch (uncommitted), Prakash-Clifton + Kelvin-Voigt
MESH=$WORK/seas-mfem-nitsche-meshes/tpv102_200m.msh
module load intel/19.1.1 impi/19.0.9 hypre/2.31.0 mumps/5.3 parmetis petsc/3.15 fftw3/3.3.8
export LD_LIBRARY_PATH="$ROOT/extern/hdf5/install/lib:$ROOT/extern/h5z-zfp/install/lib:${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"
cd $ROOT/miniapps/seas
OUT=$SCRATCH/nitsche_tpv102/tpv102_nitsche_200m_p1_pc_kv0p1_$SLURM_JOB_ID
mkdir -p $OUT
echo "commit $(git rev-parse HEAD)  dirty: $(git status --short --untracked-files=no | wc -l) files"
echo "mesh sha256 $(sha256sum $MESH | cut -d' ' -f1)"
[ "$(sha256sum $MESH | cut -d' ' -f1)" = eea25a27813dafd60a5cb5bc3bf1c8de4674e9b07858c3f9b7b8d159fa5be219 ] \
   || { echo "ERROR: mesh checksum mismatch"; exit 2; }
[ -x ./seas_nitsche_dyn_driver ] || { echo "ERROR: seas_nitsche_dyn_driver not built"; exit 2; }
grep -q 'kv_length", "inscribed"' drivers/nitsche_dyn_driver.cpp || { echo "ERROR: not the KV code"; exit 2; }
DECK=$ROOT/damp_jobs/tpv102_nitsche_pc_p1_kv0p1.toml   # reference deck of 7994409/7994410 + [numerics] kv keys
cp $DECK $OUT/
T0=$(date +%s)
# --station-copies first: the undamped references 7994409/7994410 ran ec092fa+pc.patch, which writes single-copy stations; b0c0228 defaults to averaged copies, so match the reference.
ibrun ./seas_nitsche_dyn_driver --config $DECK --mesh $MESH \
      --order 1 --tfinal 12.0 --output-dir $OUT --station-every 1 --dump-every 79 --log-every 1000 --station-copies first \
      > $OUT/driver.log 2>&1
rc=$?
echo "driver exit $rc after $(( $(date +%s) - T0 )) s"; tail -3 $OUT/driver.log
# post-processing is run BY HAND after the job
du -sh $OUT
exit $rc
