#!/bin/bash
#SBATCH --job-name=nitsche_build_verify
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=1 --ntasks-per-node=20 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=02:00:00
#SBATCH --output=nitsche_build_verify_%j.out
#SBATCH --error=nitsche_build_verify_%j.err
# Discovery build + verification on an SCEC condo node for the Nitsche TPV205 / TPV104 100 m
# benchmarks (Chunhui Zhao, 2026-10-16).  Network work on the Discovery login node:
#   cd /scratch1/$USER && git clone --recurse-submodules -b nitsche-benchmarks \
#       https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-benchmarks
#   cd seas-mfem-benchmarks/miniapps/seas/jobs/condo_nitsche && sbatch build_verify_condo.sbatch
# Discovery has no gmsh, so the meshes are made elsewhere from the committed .geo files and
# staged in $MESH_STAGE; step 2 copies them in (*.msh is git-ignored) and checks meshes.sha256.
# Steps:
#   0. the repo CLAUDE.md guard (GitHub first), plus no untracked sources
#   1. build_discovery.sh: MFEM + seas_nitsche_dyn_driver + six gate tests (modules only)
#   2. meshes from $MESH_STAGE, sha256 against meshes.sha256:
#      tpv205/mesh/tpv2053d_1000m.msh  the Mac mesh behind the smoke reference below
#      tpv205/mesh/tpv2053d_100m.msh   Mac gmsh 4.15.0 output of tpv2053d_100m.geo
#      tpv102/mesh/tpv102_100m.msh     Mac gmsh 4.15.0 output of tpv102_100m.geo, the mesh of condo
#                                      job 6459052; TPV104 runs on it (tpv104_200m.geo is
#                                      tpv102_200m.geo apart from comments; there is no tpv104_100m.geo)
#      tpv104/mesh/tpv104_200m.msh     Vista gmsh 4.15.0 mesh of Vista job 1075919, for the condo debug
#                                      run tpv104_200m_p1.sbatch
#   3. gate tests, as on the Mac at c906a8a and on Vista (job 1075901): gates 74/0, pf 22/0,
#      hook-spy at np 1, 2, 4 + abort, lsw 13/13, srw 9/9, lsw-inline at np 1, 2, 4
#   4. driver smoke: TPV205 1000 m p=1 to 1.2 s on 8 ranks, on the same mesh file as the Mac
#      reference "816 steps, t = 1.20137 s; max|V| 1.52293 m/s; max|V|dt/Dc 0.0056054", so the step
#      count and t must match and max|V| should agree to a few digits (x86 gcc vs Apple clang).
# Cost: 1 node x <= 2 h on the SCEC condo (scec_608).
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MESH_STAGE=/scratch1/$USER/condo-meshes-20261016
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd $ROOT
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "commit $(git rev-parse HEAD)"
UNTRACKED=$(git ls-files --others --exclude-standard -- '*.cpp' '*.hpp' '*.inl' Makefile '*.toml' '*.sbatch')
[ -z "$UNTRACKED" ] || { echo "ERROR: untracked sources in $ROOT:"; echo "$UNTRACKED"; exit 2; }
echo "host $(hostname)  cpu $(grep -m1 'model name' /proc/cpuinfo)"
T0=$(date +%s)
JOBS=20 bash build_discovery.sh > build_discovery_$SLURM_JOB_ID.log 2>&1
rc=$?; echo "build_discovery.sh exit $rc after $(( $(date +%s) - T0 )) s"; tail -20 build_discovery_$SLURM_JOB_ID.log
[ $rc -eq 0 ] || exit $rc
# same module environment for the runs below
module purge
module load gcc/13.3.0 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0 hdf5/1.14.3
cd miniapps/seas
echo "=== meshes from $MESH_STAGE ==="
for m in tpv205/mesh/tpv2053d_1000m.msh tpv205/mesh/tpv2053d_100m.msh tpv102/mesh/tpv102_100m.msh \
         tpv104/mesh/tpv104_200m.msh; do
   [ -f $m ] || cp -v $MESH_STAGE/$(basename $m) $m || { echo "ERROR: cannot copy $(basename $m) from $MESH_STAGE"; exit 2; }
done
sha256sum -c jobs/condo_nitsche/meshes.sha256 || { echo "ERROR: mesh checksum mismatch"; exit 2; }
echo "=== gate tests ==="
for t in test-nitsche-gates test-nitsche-pf test-nitsche-hook-spy test-nitsche-lsw test-nitsche-srw test-nitsche-lsw-inline; do
   s=$(date +%s); make MFEM_DIR=../.. HDF5_DIR="$HDF5_ROOT" SEAS_FS_LIB= $t > ${t}_$SLURM_JOB_ID.log 2>&1; rc=$?
   echo "$t exit $rc in $(( $(date +%s) - s )) s: PASSED $(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c PASSED), FAILED $(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c FAILED)"
done
echo "=== driver smoke: TPV205 1000 m p=1 1.2 s, 8 ranks ==="
OUT=/scratch1/$USER/nitsche_benchmarks/verify_tpv205_1000m_p1_1p2s_$SLURM_JOB_ID; mkdir -p $OUT
s=$(date +%s)
mpirun -np 8 --bind-to core ./seas_nitsche_dyn_driver --config tpv205/configs/tpv205_nitsche.toml \
   --mesh tpv205/mesh/tpv2053d_1000m.msh --order 1 --tfinal 1.2 --output-dir $OUT --log-every 1000 \
   > $OUT/driver.log 2>&1
echo "driver exit $? in $(( $(date +%s) - s )) s"; grep "=== done" $OUT/driver.log
echo "Mac reference: 816 steps, t = 1.20137 s; max|V| 1.52293 m/s; max|V|dt/Dc 0.0056054"
echo "total $(( $(date +%s) - T0 )) s"
