#!/bin/bash
#SBATCH -J tpv34_build_verify
#SBATCH -o tpv34_build_verify_%j.out
#SBATCH -e tpv34_build_verify_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 02:00:00
#SBATCH -A EAR24034
# Vista build + verification for the Nitsche TPV34 runs (Chunhui Zhao, TPV34 plan rev 2, step S6,
# 2026-11-12).  The TPV34 part of build_verify_vista.sbatch: it does not remake the other
# benchmarks' meshes.  Submit from miniapps/seas/jobs/vista_nitsche/ of the clone
# $WORK/seas-mfem-prepare-rh-aftershock after `git pull` on the login node.
# Steps:
#   0. the repo CLAUDE.md guard (GitHub first), plus no untracked sources
#   1. build_vista.sh: MFEM + seas_nitsche_dyn_driver (modules only); writes BUILD_COMMIT_VISTA
#   2. meshes from the committed .geo with TACC gmsh, msh 2.2: tpv34/mesh/tpv34_nitsche_1000m.msh and
#      tpv34_nitsche_200m.msh; sha256 in MESH_SHA256_TPV34_VISTA; the job stops if a mesh fails
#   3. the material table $WORK/tpv34_material/tpv34_cvmh_table.h5: sha256 against
#      jobs/vista_nitsche/tables.sha256; the job stops if it differs
#   4. unit tests: test-nitsche-tpv-setup (30/30), test-nitsche-table-material (13 and 10),
#      test-body-stations (11 on 1, 2, 3 ranks)
#   5. driver smoke: TPV34 1000 m p = 1 to 20 s on 8 ranks.  Mac reference (d0b287e7, its own 1000 m
#      mesh, 234,727 tets): "14662 steps, t = 20.0009 s; max|V| 1.22367 m/s; max|V|dt/Dc 0.00927357".
#      The Vista mesh can differ (gmsh output depends on the platform): compare numbers, not bytes.
#   6. 200 m load check: p = 2 on 144 ranks, tfinal 0.003 s: tets, dt, steps to 20 s for p = 2.  On
#      one node this checks only the driver's setup (the TPV31 one ran out of memory after it, job
#      1083417); it counts as failed only when there is no dt line.
# Exit status: 3 when a mesh or the table fails, 1 when a test or the smoke run fails, else 0.
# Cost: 1 node x <= 2 h x 0.33 SU = <= 0.7 SU.
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
rc=$?; echo "build_vista.sh exit $rc after $(( $(date +%s) - T0 )) s"; tail -20 build_vista_$SLURM_JOB_ID.log
[ $rc -eq 0 ] || exit $rc
# same module environment for the runs below
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6 gmsh/4.15.0
cd miniapps/seas
echo "=== meshes (gmsh $(gmsh --version 2>&1)) ==="
for g in tpv34/mesh/tpv34_nitsche_1000m tpv34/mesh/tpv34_nitsche_200m; do
   t=$(date +%s); gmsh -3 $g.geo -format msh22 -o $g.msh > $(basename $g)_gmsh_$SLURM_JOB_ID.log 2>&1
   rc=$?; ne=$(grep -a -A1 '^\$Elements' $g.msh 2>/dev/null | tail -1)
   echo "$g.msh: gmsh exit $rc in $(( $(date +%s) - t )) s, ${ne:-0} elements"
   [ $rc -eq 0 ] && [ "${ne:-0}" -gt 0 ] \
      || { echo "ERROR: $g.msh not made; $(basename $g)_gmsh_$SLURM_JOB_ID.log:"; grep -a "Error" $(basename $g)_gmsh_$SLURM_JOB_ID.log | head -5; exit 3; }
done
sha256sum tpv34/mesh/tpv34_nitsche_1000m.msh tpv34/mesh/tpv34_nitsche_200m.msh | tee MESH_SHA256_TPV34_VISTA
echo "=== material table ==="
TABLE=$WORK/tpv34_material/tpv34_cvmh_table.h5
[ -f "$TABLE" ] || { echo "ERROR: material table $TABLE missing"; exit 3; }
TSHA=$(sha256sum "$TABLE" | cut -d' ' -f1)
grep -q "^$TSHA  $(basename "$TABLE")\$" jobs/vista_nitsche/tables.sha256 \
   || { echo "ERROR: $TABLE sha256 $TSHA not in jobs/vista_nitsche/tables.sha256"; exit 3; }
echo "material table $TABLE sha256 $TSHA"
FAIL=0
echo "=== unit tests ==="
for t in test-nitsche-tpv-setup test-nitsche-table-material test-body-stations; do
   s=$(date +%s); make MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= $t > ${t}_$SLURM_JOB_ID.log 2>&1; rc=$?
   nf=$(grep -a -v echo ${t}_$SLURM_JOB_ID.log | grep -a -c FAILED)
   echo "$t exit $rc in $(( $(date +%s) - s )) s: PASSED $(grep -a -v echo ${t}_$SLURM_JOB_ID.log | grep -a -c PASSED), FAILED $nf"
   [ $rc -eq 0 ] && [ $nf -eq 0 ] || FAIL=$((FAIL + 1))
done
echo "=== driver smoke: TPV34 1000 m p=1 20 s, 8 ranks ==="
OUT=$SCRATCH/nitsche_benchmarks/verify_tpv34_1000m_p1_$SLURM_JOB_ID; mkdir -p $OUT
s=$(date +%s)
ibrun -n 8 -o 0 ./seas_nitsche_dyn_driver --config tpv34/configs/tpv34_nitsche.toml \
   --mesh tpv34/mesh/tpv34_nitsche_1000m.msh --order 1 --tfinal 20.0 --material-table "$TABLE" \
   --output-dir $OUT --log-every 1000 > $OUT/driver.log 2>&1
echo "driver exit $? in $(( $(date +%s) - s )) s"; grep -a "=== done" $OUT/driver.log || FAIL=$((FAIL + 1))
grep -a "material:\|table box\|off-fault stations" $OUT/driver.log
echo "off-fault station files: $(ls $OUT/*_offfault_*.dat 2>/dev/null | wc -l) (56 expected)"
echo "Mac reference (d0b287e7): 14662 steps, t = 20.0009 s; max|V| 1.22367 m/s; max|V|dt/Dc 0.00927357"
echo "=== load check: TPV34 200 m p=2, 144 ranks, tfinal 0.003 s ==="
OUT=$SCRATCH/nitsche_benchmarks/loadcheck_tpv34_p2_$SLURM_JOB_ID; mkdir -p $OUT
( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $OUT/meminfo_kB.log &
MON=$!
s=$(date +%s)
ibrun -n 144 -o 0 ./seas_nitsche_dyn_driver --config tpv34/configs/tpv34_nitsche.toml \
   --mesh tpv34/mesh/tpv34_nitsche_200m.msh --order 2 --tfinal 0.003 --material-table "$TABLE" \
   --output-dir $OUT --log-every 1 > $OUT/driver.log 2>&1
rc=$?; kill $MON 2>/dev/null; wait $MON 2>/dev/null
echo "driver exit $rc in $(( $(date +%s) - s )) s"
[ $rc -eq 0 ] || grep -q -a "  dt = " $OUT/driver.log || FAIL=$((FAIL + 1))
grep -a "tets;" $OUT/driver.log | head -1
grep -a "  dt = " $OUT/driver.log | head -1
grep -a "=== done" $OUT/driver.log || tail -15 $OUT/driver.log
awk '{u=($2-$3)/1048576; if(u>m)m=u; t=$2/1048576} END{printf "node memory: peak in use %.1f GB of %.1f GB (MemTotal - MemAvailable, 2 s samples)\n", m, t}' $OUT/meminfo_kB.log
echo "total $(( $(date +%s) - T0 )) s"
echo "failed checks: $FAIL"
[ $FAIL -eq 0 ] || exit 1
