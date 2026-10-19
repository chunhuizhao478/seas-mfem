#!/bin/bash
#SBATCH -J nitsche_profiling_build
#SBATCH -o nitsche_profiling_build_%j.out
#SBATCH -e nitsche_profiling_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 01:30:00
#SBATCH -A EAR26024
# Vista build + checks for the Nitsche profiling speedups (plan items 1, 3, 4, 5, 6;
# Chunhui Zhao, 2026-10-17).  Network work on the login node, compiling here:
#   cd $WORK && git clone --recurse-submodules -b claude/nitsche-profiling-w336yq \
#       https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-profiling
#   cd seas-mfem-profiling/miniapps/seas/jobs/vista_profiling && sbatch build_profiling_vista.sbatch
# Steps:
#   0. the repo CLAUDE.md guard (GitHub first), plus no untracked sources
#   1. build_vista.sh (MFEM + seas_nitsche_dyn_driver + six gate tests), then seas_bench_nitsche_kernels
#   2. the six Nitsche gate tests, as build_verify_vista.sbatch runs them
#   3. meshes: tpv2053d_200m.msh is the mesh of the TPV205 200 m golds (made on Vista by job 1075901 in
#      $WORK/seas-mfem-benchmarks); it is copied from there and its sha256 checked against
#      jobs/vista_profiling/meshes.sha256, so the runs can be compared with gold 1075916.
#      tpv2053d_1000m.msh (bench only) is made here with TACC gmsh, msh 2.2.
#   4. bench_nitsche_kernels: the checks of each new call against the call it replaces on 16 ranks
#      (1000 m, p = 1), and the fault writers at the 200 m fault size on 1 rank ([H]).
# Cost: 1 node x <= 1.5 h x 0.33 SU = <= 0.5 SU.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
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
T0=$(date +%s)
JOBS=64 bash build_vista.sh > build_vista_$SLURM_JOB_ID.log 2>&1
rc=$?; echo "build_vista.sh exit $rc after $(( $(date +%s) - T0 )) s"; tail -5 build_vista_$SLURM_JOB_ID.log
[ $rc -eq 0 ] || exit $rc
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6 gmsh/4.15.0
cd miniapps/seas
make -j64 MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= seas_bench_nitsche_kernels > bench_build_$SLURM_JOB_ID.log 2>&1
rc=$?; echo "bench build exit $rc"; [ $rc -eq 0 ] || { tail -20 bench_build_$SLURM_JOB_ID.log; exit $rc; }
echo "=== gate tests ==="
for t in test-nitsche-gates test-nitsche-pf test-nitsche-hook-spy test-nitsche-lsw test-nitsche-srw test-nitsche-lsw-inline; do
   s=$(date +%s); make MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= $t > ${t}_$SLURM_JOB_ID.log 2>&1; rc=$?
   echo "$t exit $rc in $(( $(date +%s) - s )) s: PASSED $(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c PASSED), FAILED $(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c FAILED)"
done
echo "=== meshes ==="
SRC=$WORK/seas-mfem-benchmarks/miniapps/seas/tpv205/mesh/tpv2053d_200m.msh
[ -f tpv205/mesh/tpv2053d_200m.msh ] || cp "$SRC" tpv205/mesh/tpv2053d_200m.msh
grep -q "^$(sha256sum tpv205/mesh/tpv2053d_200m.msh | cut -d' ' -f1)  tpv205/mesh/tpv2053d_200m.msh\$" jobs/vista_profiling/meshes.sha256 \
   && echo "tpv2053d_200m.msh: sha256 matches gold 1075916" \
   || { echo "ERROR: tpv2053d_200m.msh sha256 $(sha256sum tpv205/mesh/tpv2053d_200m.msh | cut -d' ' -f1) is not the gold mesh"; exit 2; }
gmsh -3 tpv205/mesh/tpv2053d_1000m.geo -format msh22 -o tpv205/mesh/tpv2053d_1000m.msh > gmsh_1000m_$SLURM_JOB_ID.log 2>&1
echo "tpv2053d_1000m.msh: gmsh exit $?, sha256 $(sha256sum tpv205/mesh/tpv2053d_1000m.msh | cut -d' ' -f1)"
echo "=== bench: checks on 16 ranks, TPV205 1000 m, p = 1 ==="
ibrun -n 16 -o 0 ./seas_bench_nitsche_kernels --mesh tpv205/mesh/tpv2053d_1000m.msh --order 1 --reps 20 \
   > bench_np16_$SLURM_JOB_ID.log 2>&1
echo "bench exit $?"; grep -a "checks passed\|FAILED\|of triad\|speedup\|x vs\|union pattern" bench_np16_$SLURM_JOB_ID.log
echo "=== bench [H]: fault writers at the 200 m fault size, 1 rank ==="
HDF=$SCRATCH/nitsche_profiling/hdf_bench_$SLURM_JOB_ID; mkdir -p $HDF
ibrun -n 1 -o 0 ./seas_bench_nitsche_kernels --hdf-dir $HDF --hdf-cells 105000 --hdf-frames 20 \
   > bench_h_$SLURM_JOB_ID.log 2>&1
echo "bench [H] exit $?"; cat bench_h_$SLURM_JOB_ID.log; rm -rf $HDF
echo "total $(( $(date +%s) - T0 )) s"
