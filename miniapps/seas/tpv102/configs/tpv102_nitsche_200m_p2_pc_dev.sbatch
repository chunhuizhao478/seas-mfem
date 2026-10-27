#!/bin/bash
#SBATCH -J tpv102_nitsche_200m_p2_pc
#SBATCH -o tpv102_nitsche_200m_p2_pc_%j.out
#SBATCH -e tpv102_nitsche_200m_p2_pc_%j.err
#SBATCH -p development
#SBATCH -N 36
#SBATCH -n 2016
#SBATCH -t 01:30:00
#SBATCH -A EAR20006
# SCEC TPV102 (RS aging) through the Nitsche dynamic-rupture driver, plan N-07 / P7
# (Chunhui Zhao).  p = 2 on the exact 200 m Gmsh mesh, 0 -> 12 s.
#   mesh   $WORK/seas-mfem-nitsche-meshes/tpv102_200m.msh
#          sha256 eea25a27813dafd60a5cb5bc3bf1c8de4674e9b07858c3f9b7b8d159fa5be219
#          2,466,377 tets, 37,930 fault faces, h_min 26.81 m
#   dt     0.11842 ms (0.5 x C(p) x h_min / (p^2 c_p)), 101337 steps to 12 s
#   dumps  every 169 steps (~0.02 s of model time); stations every step
#   sizing development queue (40 nodes/user, 2 h max); measured CLX rate, smoke job 7993430: p=2 2.57e-5 core-s per tet-step -> 1786 core-h / 2016 = 0.89 h
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # one thread per MPI rank
ROOT=$WORK/seas-mfem-nitsche-pc   # ec092fa + pc.patch (sha256 2891bdd9...317a), Prakash-Clifton
MESH=$WORK/seas-mfem-nitsche-meshes/tpv102_200m.msh
module load intel/19.1.1 impi/19.0.9 hypre/2.31.0 mumps/5.3 parmetis petsc/3.15 fftw3/3.3.8
export LD_LIBRARY_PATH="$ROOT/extern/hdf5/install/lib:$ROOT/extern/h5z-zfp/install/lib:${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"
cd $ROOT/miniapps/seas
OUT=$SCRATCH/nitsche_tpv102/tpv102_nitsche_200m_p2_pc_$SLURM_JOB_ID
mkdir -p $OUT
echo "commit $(git rev-parse HEAD)  dirty: $(git status --short --untracked-files=no | wc -l) files"
echo "mesh sha256 $(sha256sum $MESH | cut -d' ' -f1)"
[ "$(sha256sum $MESH | cut -d' ' -f1)" = eea25a27813dafd60a5cb5bc3bf1c8de4674e9b07858c3f9b7b8d159fa5be219 ] \
   || { echo "ERROR: mesh checksum mismatch"; exit 2; }
[ -x ./seas_nitsche_dyn_driver ] || { echo "ERROR: seas_nitsche_dyn_driver not built"; exit 2; }
DECK=$WORK/seas-mfem-nitsche-jobs/tpv102_nitsche_pc_p2.toml   # 7993802 deck + [friction.prakash_clifton] v_star 1, length 0.01
cp $DECK $OUT/
T0=$(date +%s)
ibrun ./seas_nitsche_dyn_driver --config $DECK --mesh $MESH \
      --order 2 --tfinal 12.0 --output-dir $OUT --station-every 1 --dump-every 169 --log-every 1000 \
      > $OUT/driver.log 2>&1
rc=$?
echo "driver exit $rc after $(( $(date +%s) - T0 )) s"; tail -3 $OUT/driver.log
# post-processing is run BY HAND after the job
du -sh $OUT
exit $rc
