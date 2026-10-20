#!/bin/bash
#SBATCH -J nitsche_map_build
#SBATCH -o nitsche_map_build_%j.out
#SBATCH -e nitsche_map_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 00:30:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# MAP build (speed-next, Chunhui Zhao 2026-10-20): seas_nitsche_dyn_driver_g, the driver compiled
# with -g for Linaro MAP and perf, AFTER build_lts_vista.sbatch passed at the same commit (it reuses that
# build's libmfem.a and objects).  -g does not change the generated code: the job strips the debug
# sections from drivers/nitsche_dyn_driver_g.o and compares it with drivers/nitsche_dyn_driver.o byte
# for byte, and prints the result.  Writes miniapps/seas/BUILD_COMMIT_VISTA_G.
# Submit from miniapps/seas/jobs/vista_profiling/.  Cost: 1 node x about 10 min = 0.06 SU (0.17 SU at most).
export LC_ALL=C LANG=C
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
HEAD=$(git rev-parse HEAD)
cd miniapps/seas
[ "$(cat BUILD_COMMIT_VISTA 2>/dev/null)" = "$HEAD" ] \
   || { echo "ERROR: build_lts_vista.sbatch built '$(cat BUILD_COMMIT_VISTA 2>/dev/null)', checkout is $HEAD"; exit 2; }
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
rm -f drivers/nitsche_dyn_driver_g.o seas_nitsche_dyn_driver_g BUILD_COMMIT_VISTA_G
s=$(date +%s)
make -j16 MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= seas_nitsche_dyn_driver_g > map_build_$SLURM_JOB_ID.log 2>&1
rc=$?; echo "seas_nitsche_dyn_driver_g build exit $rc in $(( $(date +%s) - s )) s"
[ $rc -eq 0 ] && [ -x seas_nitsche_dyn_driver_g ] || { tail -30 map_build_$SLURM_JOB_ID.log; exit 2; }
T=$(mktemp -d)
objcopy --strip-debug drivers/nitsche_dyn_driver.o $T/a.o
objcopy --strip-debug drivers/nitsche_dyn_driver_g.o $T/b.o
if cmp -s $T/a.o $T/b.o; then
   echo "code check: drivers/nitsche_dyn_driver_g.o without its debug sections is byte-identical to drivers/nitsche_dyn_driver.o"
else
   echo "WARNING: the stripped -g object differs from drivers/nitsche_dyn_driver.o; section sizes:"
   size -A $T/a.o | sort > $T/a.txt; size -A $T/b.o | sort > $T/b.txt; diff $T/a.txt $T/b.txt | head -20
fi
rm -rf $T
readelf -S seas_nitsche_dyn_driver_g | grep -c "\.debug_" | awk '{print "debug sections in seas_nitsche_dyn_driver_g: " $1}'
echo "$HEAD" > BUILD_COMMIT_VISTA_G
echo "MAP build OK: commit $HEAD"
