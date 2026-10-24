#!/bin/bash
#SBATCH -J fault_writer_build
#SBATCH -o fault_writer_build_%j.out
#SBATCH -e fault_writer_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 01:00:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131,i617-101
# Vista build + checks for the fault writer (fault-writer note rev1, W0/W1, and the per-rank balance
# report L1; Chunhui Zhao, 2026-10-24), before the v2_*.sbatch runs; build_mesh_parts_vista.sbatch with
# the fault-frame spread test in place of the gold-mesh parts runs.  Network work on the login node,
# compiling here:
#   cd $WORK/seas-mfem-profiling && git fetch origin && git checkout <the pushed commit>
#   cd miniapps/seas/jobs/vista_profiling && sbatch build_fault_writer_vista.sbatch
# Steps:
#   0. the repo CLAUDE.md guard (GitHub first), plus no untracked sources
#   1. build_vista.sh (MFEM, whose mesh/vtkhdf.cpp gained the chunk writes, + seas_nitsche_dyn_driver
#      + six gate tests; writes BUILD_COMMIT_VISTA), then seas_test_nitsche_lts,
#      seas_test_station_output, seas_test_mesh_parts and seas_test_fault_frame_spread
#   2. the six Nitsche gate tests, the LTS suite (L-1 to L-10, about 13 min on one core), the station
#      output test, the mesh parts test and the fault-frame spread test (FS-1 to FS-4), each on 1, 2, 3
#      and 4 ranks
#   3. the fault-frame spread test on 16 ranks (15 helpers; the TPV205 200 m frame length)
#   4. meshes: the TPV205 200 m gold mesh is copied from $WORK/seas-mfem-benchmarks if absent and its
#      sha256 checked against jobs/vista_profiling/meshes.sha256
# Cost: 1 node x <= 1 h x 0.33 SU = <= 0.33 SU.
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
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
cd miniapps/seas
TESTS="seas_test_nitsche_lts seas_test_station_output seas_test_mesh_parts seas_test_fault_frame_spread"
make -j64 MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= $TESTS > lts_test_build_$SLURM_JOB_ID.log 2>&1
rc=$?; echo "$TESTS build exit $rc"; [ $rc -eq 0 ] || { tail -20 lts_test_build_$SLURM_JOB_ID.log; exit $rc; }
echo "=== gate tests ==="
fail=0
for t in test-nitsche-gates test-nitsche-pf test-nitsche-hook-spy test-nitsche-lsw test-nitsche-srw test-nitsche-lsw-inline test-nitsche-lts test-station-output test-mesh-parts test-fault-frame-spread; do
   s=$(date +%s); make MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= $t > ${t}_$SLURM_JOB_ID.log 2>&1; rc=$?
   nf=$(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c FAILED)
   echo "$t exit $rc in $(( $(date +%s) - s )) s: PASSED $(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c PASSED), FAILED $nf"
   [ $rc -eq 0 ] && [ $nf -eq 0 ] || fail=1
done
grep -a "Summary" test-nitsche-lts_$SLURM_JOB_ID.log test-station-output_$SLURM_JOB_ID.log
grep -a "passed," test-mesh-parts_$SLURM_JOB_ID.log test-fault-frame-spread_$SLURM_JOB_ID.log
echo "=== fault-frame spread test on 16 ranks ==="
t=fault_frame_spread_16
s=$(date +%s); ibrun -n 16 -o 0 ./seas_test_fault_frame_spread > ${t}_$SLURM_JOB_ID.log 2>&1; rc=$?
nf=$(grep -c FAILED ${t}_$SLURM_JOB_ID.log)
echo "16 ranks exit $rc in $(( $(date +%s) - s )) s: PASSED $(grep -c PASSED ${t}_$SLURM_JOB_ID.log), FAILED $nf"
grep -a "helpers\|passed," ${t}_$SLURM_JOB_ID.log
[ $rc -eq 0 ] && [ $nf -eq 0 ] || fail=1
rm -rf test_fault_frame_spread_out
echo "=== meshes ==="
m=tpv205/mesh/tpv2053d_200m.msh
[ -f $m ] || cp "$WORK/seas-mfem-benchmarks/miniapps/seas/$m" $m
grep -q "^$(sha256sum $m | cut -d' ' -f1)  $m\$" jobs/vista_profiling/meshes.sha256 \
   && echo "$m: sha256 matches the gold mesh" \
   || { echo "ERROR: $m sha256 $(sha256sum $m | cut -d' ' -f1) is not the gold mesh"; fail=1; }
echo "total $(( $(date +%s) - T0 )) s"
[ $fail -eq 0 ] || { echo "ERROR: a check failed; do not submit the v2_*.sbatch jobs"; exit 1; }
echo "all checks passed"
