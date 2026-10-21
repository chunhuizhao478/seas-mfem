#!/bin/bash
#SBATCH -J complex_build_verify
#SBATCH -o complex_build_verify_%j.out
#SBATCH -e complex_build_verify_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 03:00:00
#SBATCH -A EAR26024
# Vista build + verification for the complex fault geometries TPV12, TPV22 and TPV15 (complex-geometries
# plan rev 2, phase P5; Chunhui Zhao, 2026-10-21).  Network work on the login node, compiling here:
#   cd $WORK && git clone --recurse-submodules -b nitsche-dynamic-rupture-complex-geometries \
#       https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-complex-geometries
#   cd seas-mfem-complex-geometries/miniapps/seas/jobs/vista_complex && sbatch build_verify_complex_vista.sbatch
# Steps:
#   0. the repo CLAUDE.md guard (GitHub first), plus no untracked sources; the base commit 64e2b92 of step 4
#      is in the clone and has the same MFEM sources and toml11 commit as HEAD
#   1. build_vista.sh (MFEM + seas_nitsche_dyn_driver + six gate tests; writes BUILD_COMMIT_VISTA), then
#      seas_test_nitsche_options
#   2. meshes from the committed .geo with TACC gmsh, msh 2.2, sha256 in MESH_SHA256_COMPLEX_VISTA:
#      the 1000 m and 200 m meshes of TPV12, TPV22 and TPV15, and tpv205/mesh/tpv2053d_1000m.msh for step 4
#      (the 100 m TPV12 mesh has its own job and checksum file, make_tpv12_100m_mesh_vista.sbatch)
#   3. tests: the six gate suites of build_verify_vista.sbatch, and test-nitsche-options (52 checks, 1 and 2 ranks)
#   4. G-off on Vista: the driver of the base commit 64e2b92 (its miniapps/seas sources by git archive into
#      $SCRATCH, this build's MFEM by symlink) and this driver run TPV205 1000 m, p = 1, 1.2 s, 8 ranks,
#      under global stepping and --lts rate2, with no option table; both runs must finish, and every output
#      file must be identical byte for byte (the cloud gate, 18 of 18 cells, repeated on this hardware),
#      fault.vtkhdf by h5diff (HDF5 stores object times; the cloud build had no HDF5, so this file is new
#      to the gate).  The base driver is a comparison binary only, never a result.
#   5. option smoke: TPV12, TPV22 and TPV15 at 1000 m, p = 1, 0.5 s, 8 ranks (TPV15 with --lts rate2), with
#      the committed decks: exit 0 and "=== done", the option and tip-lock lines, the station values at t = 0
#   6. 200 m setup check: each problem at p = 1 on one node with the ranks per node its production jobs use
#      (TPV12 144, TPV15 96, TPV22 72), a few steps: tets, dt, node memory peak.  Every rank reads the whole
#      serial mesh (drivers/nitsche_dyn_driver.cpp:927), about 0.81 GB per rank for the 1.89 M tets of
#      TPV205 200 m (Vista 1076378), so about 1.08 GB for TPV15 (2.53 M tets) and 1.55 GB for TPV22
#      (3.63 M), hence the fewer ranks per node (inferred, which this step measures).  That part of a node's
#      memory is the same at every order; the rest of the problem is here all on one node, more than a
#      production node holds.  p = 2 does not fit on one node at 144 ranks (the TPV205 and TPV104 checks
#      of build_verify_vista.sbatch ran out of memory in every build), so it is not tried.  A run that
#      dies before its dt line counts as a failure; a later death and the ms per step (the whole problem
#      on one node, not the production speed) are only reported.
# Exit 0 and VERIFIED_COMMIT_COMPLEX_VISTA = HEAD only when every check passed; run_complex_vista.sh
# refuses to run otherwise.  Exit 1 when a check failed (all are still run and printed), 2 for the guard
# or the build, 3 when a mesh was not made.
# Cost: 1 node x <= 3 h x 0.33 SU = <= 1.0 SU.
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
UNTRACKED=$(git ls-files --others --exclude-standard -- '*.cpp' '*.hpp' '*.inl' Makefile '*.toml' '*.sbatch' '*.geo' '*.csv')
[ -z "$UNTRACKED" ] || { echo "ERROR: untracked sources in $ROOT:"; echo "$UNTRACKED"; exit 2; }
HEAD=$(git rev-parse HEAD)
rm -f miniapps/seas/VERIFIED_COMMIT_COMPLEX_VISTA
BASE_COMMIT=64e2b921
git cat-file -e $BASE_COMMIT^{commit} || { echo "ERROR: base commit $BASE_COMMIT not in this clone"; exit 2; }
# step 4 links this build's MFEM into the base tree: the MFEM library sources and the toml11 commit must match
[ -z "$(git diff --name-only $BASE_COMMIT HEAD -- config fem general linalg mesh makefile ':(glob)*.hpp')" ] \
   || { echo "ERROR: HEAD's MFEM sources differ from $BASE_COMMIT's"; exit 2; }
[ "$(git rev-parse $BASE_COMMIT:miniapps/seas/extern/toml11)" = "$(git rev-parse HEAD:miniapps/seas/extern/toml11)" ] \
   || { echo "ERROR: the toml11 submodule commit differs from $BASE_COMMIT's"; exit 2; }
T0=$(date +%s)
JOBS=64 bash build_vista.sh > build_vista_$SLURM_JOB_ID.log 2>&1
rc=$?; echo "build_vista.sh exit $rc after $(( $(date +%s) - T0 )) s"; tail -20 build_vista_$SLURM_JOB_ID.log
[ $rc -eq 0 ] || exit 2
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6 gmsh/4.15.0
command -v h5diff >/dev/null || { echo "ERROR: h5diff not on PATH (phdf5 module)"; exit 2; }
cd miniapps/seas
make -j64 MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= seas_test_nitsche_options > build_options_$SLURM_JOB_ID.log 2>&1 \
   || { echo "ERROR: seas_test_nitsche_options build failed"; tail -20 build_options_$SLURM_JOB_ID.log; exit 2; }

echo "=== meshes (gmsh $(gmsh --version 2>&1)) ==="
MESHES="tpv12/mesh/tpv12_1000m tpv12/mesh/tpv12_200m tpv22/mesh/tpv22_1000m tpv22/mesh/tpv22_200m
        tpv15/mesh/tpv15_1000m tpv15/mesh/tpv15_200m tpv205/mesh/tpv2053d_1000m"
for g in $MESHES; do
   t=$(date +%s); gmsh -3 $g.geo -format msh22 -o $g.msh > $(basename $g)_gmsh_$SLURM_JOB_ID.log 2>&1; rc=$?
   ne=$(grep -a -m1 -A1 '^\$Elements' $g.msh 2>/dev/null | tail -1)
   echo "$g.msh: gmsh exit $rc in $(( $(date +%s) - t )) s, ${ne:-0} elements (all types)"
   [ $rc -eq 0 ] && [ "${ne:-0}" -gt 0 ] \
      || { echo "ERROR: $g.msh not made; $(basename $g)_gmsh_$SLURM_JOB_ID.log:"; grep -a "Error" $(basename $g)_gmsh_$SLURM_JOB_ID.log | head -5; exit 3; }
done
sha256sum $(for g in $MESHES; do echo $g.msh; done) | tee MESH_SHA256_COMPLEX_VISTA
FAIL=0

echo "=== tests ==="
for t in test-nitsche-gates test-nitsche-pf test-nitsche-hook-spy test-nitsche-lsw test-nitsche-srw test-nitsche-lsw-inline \
         test-nitsche-options; do
   s=$(date +%s); make MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= $t > ${t}_$SLURM_JOB_ID.log 2>&1; rc=$?
   nf=$(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c FAILED)
   echo "$t exit $rc in $(( $(date +%s) - s )) s: PASSED $(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c PASSED), FAILED $nf"
   [ $rc -eq 0 ] && [ $nf -eq 0 ] || FAIL=$((FAIL + 1))
done
grep -a "passed," test-nitsche-options_$SLURM_JOB_ID.log || { echo "test-nitsche-options printed no result line"; FAIL=$((FAIL + 1)); }

echo "=== G-off on Vista: base $BASE_COMMIT against this commit, TPV205 1000 m p = 1, 1.2 s, 8 ranks, no option table ==="
# The base tree: its miniapps/seas sources (git archive, without the benchmark data folders, which hold no
# source), and this checkout's MFEM level by symlink (the same sources, checked in step 0, and this
# build's libmfem.a, config.mk and _config.hpp) and toml11 (the same commit).  Tried in the cloud: the
# driver built this way is byte for byte the 64e2b92 driver built the usual way.
X=$SCRATCH/complex_geometries/base_$SLURM_JOB_ID
mkdir -p $X/miniapps
for e in $ROOT/*; do [ "$(basename $e)" = miniapps ] || ln -s $e $X/; done
SEAS_PATHS=$(git -C $ROOT ls-tree --name-only $BASE_COMMIT miniapps/seas/ \
   | grep -v -E '^miniapps/seas/(tpv[0-9]+|bp[0-9]+|safs|document|debug_document|archive)$')
git -C $ROOT archive $BASE_COMMIT $SEAS_PATHS | tar -x -C $X
rmdir $X/miniapps/seas/extern/toml11 2>/dev/null; ln -s $ROOT/miniapps/seas/extern/toml11 $X/miniapps/seas/extern/toml11
if ( cd $X/miniapps/seas && make -j64 MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= seas_nitsche_dyn_driver ) \
      > build_base_$SLURM_JOB_ID.log 2>&1; then
   for v in gts lts; do
      FLAGS=""; [ $v = lts ] && FLAGS="--lts rate2"
      ok=1
      for d in base new; do
         DRV=./seas_nitsche_dyn_driver; [ $d = base ] && DRV=$X/miniapps/seas/seas_nitsche_dyn_driver
         OUT=$SCRATCH/complex_geometries/goff_${v}_${d}_$SLURM_JOB_ID; mkdir -p $OUT
         s=$(date +%s)
         ibrun -n 8 -o 0 $DRV --config tpv205/configs/tpv205_nitsche.toml --mesh tpv205/mesh/tpv2053d_1000m.msh \
            --order 1 --tfinal 1.2 --output-dir $OUT --dump-every 20 --log-every 1000 $FLAGS > $OUT.log 2>&1
         rc=$?
         echo "$v $d: driver exit $rc in $(( $(date +%s) - s )) s; $(grep -a '=== done' $OUT.log)"
         [ $rc -eq 0 ] && grep -q "=== done" $OUT.log && [ -f $OUT/fault.vtkhdf ] \
            && ls $OUT/tpv205_station_*.dat >/dev/null 2>&1 || { ok=0; tail -10 $OUT.log; }
      done
      A=$SCRATCH/complex_geometries/goff_${v}_base_$SLURM_JOB_ID; B=$SCRATCH/complex_geometries/goff_${v}_new_$SLURM_JOB_ID
      n=0; nd=0
      for f in $(ls $A); do
         n=$((n + 1))
         case $f in
            *.vtkhdf) h5diff -q $A/$f $B/$f || { nd=$((nd + 1)); echo "   DIFF $f (h5diff)"; } ;;
            *) cmp -s $A/$f $B/$f || { nd=$((nd + 1)); echo "   DIFF $f"; } ;;
         esac
      done
      [ "$(ls $B | wc -l)" = "$n" ] || { echo "   file lists differ: base $n, new $(ls $B | wc -l)"; nd=$((nd + 1)); }
      echo "G-off $v: $n files compared, $nd differ"
      [ $ok = 1 ] && [ $n -gt 0 ] && [ $nd = 0 ] && echo "G-off $v PASSED" || { echo "G-off $v FAILED"; FAIL=$((FAIL + 1)); }
   done
else
   echo "ERROR: base driver build failed"; tail -20 build_base_$SLURM_JOB_ID.log; FAIL=$((FAIL + 1))
fi

echo "=== option smoke: 1000 m, p = 1, 0.5 s, 8 ranks, committed decks ==="
for case in "tpv12 tpv12/configs/tpv12_nitsche.toml tpv12/mesh/tpv12_1000m.msh" \
            "tpv22 tpv22/configs/tpv22_nitsche.toml tpv22/mesh/tpv22_1000m.msh" \
            "tpv15 tpv15/configs/tpv15_nitsche.toml tpv15/mesh/tpv15_1000m.msh --lts rate2"; do
   set -- $case; PROB=$1; DECK=$2; MESH=$3; shift 3
   OUT=$SCRATCH/complex_geometries/smoke_${PROB}_1000m_$SLURM_JOB_ID; mkdir -p $OUT
   s=$(date +%s)
   ibrun -n 8 -o 0 ./seas_nitsche_dyn_driver --config $DECK --mesh $MESH --order 1 --tfinal 0.5 \
      --output-dir $OUT --log-every 1000 "$@" > $OUT/driver.log 2>&1
   rc=$?
   echo "$PROB [$*]: driver exit $rc in $(( $(date +%s) - s )) s"
   grep -a "prestress: \[\|fault options\|fault tips:\|forced rupture:\|tension floor:\|LTS (plan" $OUT/driver.log
   if [ $rc -eq 0 ] && grep -a "=== done" $OUT/driver.log; then
      for f in $OUT/${PROB}_station_*[0-9].dat; do
         echo "   $(basename $f .dat): t = 0 tau1 tau2 sigma_n = $(awk '!/^#/{print $6, $7, $8; exit}' $f)"
      done
   else
      tail -20 $OUT/driver.log; FAIL=$((FAIL + 1))
   fi
done

echo "=== 200 m setup check: p = 1, one node, the production ranks per node ==="
# Node memory in use (MemTotal - MemAvailable, /proc/meminfo) every 2 s; ms per step from the driver's
# wall column between the log lines after the first, which carries the setup.
for case in "tpv12 144 tpv12/configs/tpv12_nitsche.toml tpv12/mesh/tpv12_200m.msh 0.012" \
            "tpv22 72 tpv22/configs/tpv22_nitsche.toml tpv22/mesh/tpv22_200m.msh 0.01" \
            "tpv15 96 tpv15/configs/tpv15_nitsche.toml tpv15/mesh/tpv15_200m.msh 0.0025 --lts rate2"; do
   set -- $case; PROB=$1; NP=$2; DECK=$3; MESH=$4; TF=$5; shift 5
   echo "--- $PROB 200 m p = 1 on $NP ranks [$*], tfinal $TF s"
   OUT=$SCRATCH/complex_geometries/setupcheck_${PROB}_200m_p1_$SLURM_JOB_ID; mkdir -p $OUT
   ( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $OUT/meminfo_kB.log &
   MON=$!
   s=$(date +%s)
   ibrun -n $NP -o 0 ./seas_nitsche_dyn_driver --config $DECK --mesh $MESH --order 1 --tfinal $TF \
      --output-dir $OUT --log-every 5 "$@" > $OUT/driver.log 2>&1
   rc=$?; kill $MON 2>/dev/null; wait $MON 2>/dev/null
   echo "driver exit $rc in $(( $(date +%s) - s )) s"
   grep -a "tets;\|  dt = \|LTS (plan" $OUT/driver.log | head -3
   grep -q -a "  dt = " $OUT/driver.log || { echo "FAILED: no dt line (the driver died in setup)"; FAIL=$((FAIL + 1)); }
   [ $rc -eq 0 ] && grep -a "^  step " $OUT/driver.log | awk -v np=$NP '{for(i=1;i<=NF;i++){if($i=="step")st=$(i+1); if($i=="wall")w=$(i+1)} n++; if(n==1){s0=st; w0=w} s1=st; w1=w}
      END{if(n>1 && s1>s0) printf "%.2f ms per step over steps %d-%d on %d ranks (whole problem on one node; not the production speed)\n", 1000*(w1-w0)/(s1-s0), s0, s1, np}'
   grep -a "=== done" $OUT/driver.log || tail -8 $OUT/driver.log
   awk -v np=$NP '{u=($2-$3)/1048576; if(u>m)m=u; t=$2/1048576} END{printf "node memory: peak in use %.1f GB of %.1f GB, %.2f GB per rank\n", m, t, m/np}' $OUT/meminfo_kB.log
done
echo "total $(( $(date +%s) - T0 )) s"
echo "failed checks: $FAIL"
[ $FAIL -eq 0 ] || exit 1
echo $HEAD > VERIFIED_COMMIT_COMPLEX_VISTA
echo "VERIFIED_COMMIT_COMPLEX_VISTA = $HEAD"
