#!/bin/bash
#SBATCH -J mesh_parts_build
#SBATCH -o mesh_parts_build_%j.out
#SBATCH -e mesh_parts_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 01:30:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131,i617-101
# Vista build + checks for the mesh setup change (partition-memory report rev1 options D and E,
# speed-next S1; Chunhui Zhao, 2026-10-22), before the mp_*.sbatch runs; build_speed_next_vista.sbatch
# with the mesh parts test added.  Network work on the login node, compiling here:
#   cd $WORK/seas-mfem-profiling && git fetch origin && git checkout <the pushed commit>
#   cd miniapps/seas/jobs/vista_profiling && sbatch build_mesh_parts_vista.sbatch
# Steps:
#   0. the repo CLAUDE.md guard (GitHub first), plus no untracked sources
#   1. build_vista.sh (MFEM + seas_nitsche_dyn_driver + six gate tests; writes BUILD_COMMIT_VISTA),
#      then seas_test_nitsche_lts, seas_test_station_output and seas_test_mesh_parts
#   2. the six Nitsche gate tests, the LTS suite (L-1 to L-10, about 13 min on one core), the
#      station output test (O2-1 to O2-4 on 1, 2, 3 and 4 ranks) and the mesh parts test (MP-1 to
#      MP-8 on 1, 2, 3 and 4 ranks, on the built-in mesh)
#   3. meshes: the TPV205 and TPV104 200 m gold meshes are copied from $WORK/seas-mfem-benchmarks and
#      their sha256 checked against jobs/vista_profiling/meshes.sha256
#   4. the mesh parts test on both gold meshes, 16 ranks each (METIS partition; every rank also
#      reads the serial mesh for the comparison: 3.4 GB peak per rank on 4 ranks, less on 16).
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
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
cd miniapps/seas
make -j64 MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= seas_test_nitsche_lts seas_test_station_output seas_test_mesh_parts > lts_test_build_$SLURM_JOB_ID.log 2>&1
rc=$?; echo "seas_test_nitsche_lts, seas_test_station_output and seas_test_mesh_parts build exit $rc"; [ $rc -eq 0 ] || { tail -20 lts_test_build_$SLURM_JOB_ID.log; exit $rc; }
echo "=== gate tests ==="
fail=0
for t in test-nitsche-gates test-nitsche-pf test-nitsche-hook-spy test-nitsche-lsw test-nitsche-srw test-nitsche-lsw-inline test-nitsche-lts test-station-output test-mesh-parts; do
   s=$(date +%s); make MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= $t > ${t}_$SLURM_JOB_ID.log 2>&1; rc=$?
   nf=$(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c FAILED)
   echo "$t exit $rc in $(( $(date +%s) - s )) s: PASSED $(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c PASSED), FAILED $nf"
   [ $rc -eq 0 ] && [ $nf -eq 0 ] || fail=1
done
grep -a "Summary" test-nitsche-lts_$SLURM_JOB_ID.log test-station-output_$SLURM_JOB_ID.log
grep -a "passed," test-mesh-parts_$SLURM_JOB_ID.log
echo "=== meshes ==="
for m in tpv205/mesh/tpv2053d_200m.msh tpv104/mesh/tpv104_200m.msh; do
   [ -f $m ] || cp "$WORK/seas-mfem-benchmarks/miniapps/seas/$m" $m
   grep -q "^$(sha256sum $m | cut -d' ' -f1)  $m\$" jobs/vista_profiling/meshes.sha256 \
      && echo "$m: sha256 matches the gold mesh" \
      || { echo "ERROR: $m sha256 $(sha256sum $m | cut -d' ' -f1) is not the gold mesh"; fail=1; }
done
echo "=== mesh parts test on the gold meshes (16 ranks) ==="
for spec in tpv205/mesh/tpv2053d_200m.msh:103 tpv104/mesh/tpv104_200m.msh:3; do
   m=${spec%%:*}; a=${spec##*:}; t=mesh_parts_$(basename $m .msh)
   s=$(date +%s); ibrun -n 16 -o 0 ./seas_test_mesh_parts --mesh $m --fault-attr $a > ${t}_$SLURM_JOB_ID.log 2>&1; rc=$?
   nf=$(grep -c FAILED ${t}_$SLURM_JOB_ID.log)
   echo "$m (fault attribute $a) exit $rc in $(( $(date +%s) - s )) s: PASSED $(grep -c PASSED ${t}_$SLURM_JOB_ID.log), FAILED $nf"
   grep -a "DG operator\|passed," ${t}_$SLURM_JOB_ID.log
   [ $rc -eq 0 ] && [ $nf -eq 0 ] || fail=1
done
rm -rf test_mesh_parts_out
echo "total $(( $(date +%s) - T0 )) s"
[ $fail -eq 0 ] || { echo "ERROR: a check failed; do not submit the mp_*.sbatch jobs"; exit 1; }
echo "all checks passed"
