#!/bin/bash
#SBATCH -J tpv2627_build
#SBATCH -o tpv2627_build_%j.out
#SBATCH -e tpv2627_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 02:00:00
#SBATCH -A EAR26024
# Vista build + checks for TPV26/27, rheology plan rev 1, phases R1 and R3 (Chunhui Zhao, 2026-10-20),
# before the B0 jobs (b0_*.sbatch) and the TPV26 run (tpv26_200m_p2.sbatch).  Network work on the login node,
# compiling here:
#   cd $WORK && git clone --recurse-submodules -b nitsche-dynamic-rupture-nonlinearity \
#       https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-nonlinearity
#   (or, in that clone: git fetch origin && git checkout <the pushed commit>)
#   cd seas-mfem-nonlinearity/miniapps/seas/jobs/vista_rheology && sbatch build_tpv2627_vista.sbatch
# Steps:
#   0. the repo CLAUDE.md guard (GitHub first), plus no untracked sources
#   1. build_vista.sh (MFEM + seas_nitsche_dyn_driver + six gate tests; writes BUILD_COMMIT_VISTA), then
#      seas_test_nitsche_lts, the rheology unit tests and the operator tests O1-O7
#   2. tests: rheology U1-U11 (make test-rheology: 24, 13, 36 and 25 checks in the cloud), the six
#      Nitsche gate test sets as build_verify_vista.sbatch runs them, the LTS suite (about 13 min;
#      the fault closure it calls changed in R1) and O1-O7 (make test-nitsche-rheology: 36 checks per
#      rank count, 1, 2 and 4 ranks; the operator changed in R3)
#   3. meshes: TPV26/27 1000 m and 200 m from tpv26/mesh/tpv2627_half.geo with TACC gmsh, sha256 to
#      miniapps/seas/MESH_SHA256_TPV2627_VISTA (run_nitsche_vista.sh accepts it); the TPV205 and TPV104
#      200 m gold meshes for the B0 jobs, copied from $WORK/seas-mfem-benchmarks and checked against
#      jobs/vista_profiling/meshes.sha256, as build_lts_vista.sbatch does
#   4. B1 on Vista: TPV26 1000 m p=1, 0 -> 2.5 s, 8 ranks.  Pass: no abort, and fault nodes farther than
#      r_crit = 4 km from the hypocenter rupture (scripts/make_cplot.py), so the rupture leaves the forced
#      zone.  Cloud reference (its own built-in-kernel 1000 m mesh, 3 ranks, 0 -> 2.5 s): 1772 steps, dt
#      1.4112e-03 s, no abort; all 59 nodes inside r_crit ruptured, and 90 outside it by 2.5 s, the
#      farthest 7.1 km away (on the first, OpenCASCADE mesh: 1 node by 1.5 s, the plan's B1 time, too
#      close to call; 88-89 by 2.5 s).  The off-fault station files (tpv26_offfault_*.dat) are listed.
#   4b. R3 output check: TPV27 deck as committed, 1000 m p=1, 0 -> 1.0 s, 8 ranks.  Pass: no abort and
#      <output_dir>/rheology.vtkhdf holds 3 snapshots (output_dt_s = 0.5: t = 0, 0.5, 1.0 s).  The
#      writer needs parallel HDF5, which the cloud does not have, so this is its first run.
#   5. 200 m load check: TPV26 200 m, p = 2, all 144 ranks, tfinal 0.003 s; node memory peak, tets,
#      implied min 6V/A and dt, steps to 13 s, as build_verify_vista.sbatch does for TPV205/TPV104.  It
#      decides the ranks per node of tpv26_200m_p2.sbatch.
# Cost: 1 node x <= 2 h x 0.33 SU = <= 0.66 SU.
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
MK="make MFEM_DIR=../.. HDF5_DIR=$TACC_PHDF5_DIR SEAS_FS_LIB="
$MK -j64 seas_test_nitsche_lts seas_test_rheology_pointwise seas_test_rheology_initial_stress \
   seas_test_rheology_config seas_test_rheology_forced_rupture seas_test_nitsche_rheology \
   > tests_build_$SLURM_JOB_ID.log 2>&1
rc=$?; echo "test build exit $rc"; [ $rc -eq 0 ] || { tail -20 tests_build_$SLURM_JOB_ID.log; exit $rc; }
fail=0
echo "=== rheology unit tests U1-U11 ==="
$MK test-rheology > test-rheology_$SLURM_JOB_ID.log 2>&1; rc=$?
grep -a "passed, " test-rheology_$SLURM_JOB_ID.log
nf=$(grep -a "passed, " test-rheology_$SLURM_JOB_ID.log | awk '{s += $5} END {print s + 0}')
echo "test-rheology exit $rc, failed $nf"
[ $rc -eq 0 ] && [ $nf -eq 0 ] || fail=1
echo "=== Nitsche gate tests and LTS suite ==="
for t in test-nitsche-gates test-nitsche-pf test-nitsche-hook-spy test-nitsche-lsw test-nitsche-srw test-nitsche-lsw-inline test-nitsche-lts test-nitsche-rheology; do
   s=$(date +%s); $MK $t > ${t}_$SLURM_JOB_ID.log 2>&1; rc=$?
   nf=$(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c FAILED)
   echo "$t exit $rc in $(( $(date +%s) - s )) s: PASSED $(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c PASSED), FAILED $nf"
   [ $rc -eq 0 ] && [ $nf -eq 0 ] || fail=1
done
grep -a "Summary" test-nitsche-lts_$SLURM_JOB_ID.log
echo "=== meshes (gmsh $(gmsh --version 2>&1)) ==="
s=$(date +%s)
SHA_FILE=MESH_SHA256_TPV2627_VISTA bash tpv26/mesh/make_tpv2627_meshes.sh 1000 200; rc=$?
echo "make_tpv2627_meshes.sh exit $rc in $(( $(date +%s) - s )) s"
[ $rc -eq 0 ] || { echo "ERROR: TPV26/27 meshes failed"; exit 1; }
for m in tpv205/mesh/tpv2053d_200m.msh tpv104/mesh/tpv104_200m.msh; do
   [ -f $m ] || cp "$WORK/seas-mfem-benchmarks/miniapps/seas/$m" $m
   grep -q "^$(sha256sum $m | cut -d' ' -f1)  $m\$" jobs/vista_profiling/meshes.sha256 \
      && echo "$m: sha256 matches the gold mesh" \
      || { echo "ERROR: $m sha256 $(sha256sum $m | cut -d' ' -f1) is not the gold mesh"; fail=1; }
done
echo "=== B1: TPV26 1000 m p=1, 0 -> 2.5 s, 8 ranks ==="
OUT=$SCRATCH/rheology/b1_tpv26_1000m_p1_$SLURM_JOB_ID; mkdir -p $OUT
s=$(date +%s)
ibrun -n 8 -o 0 ./seas_nitsche_dyn_driver --config tpv26/configs/tpv26_nitsche_kv0p1.toml \
   --order 1 --tfinal 2.5 --output-dir $OUT --log-every 200 > $OUT/driver.log 2>&1
rc=$?; echo "driver exit $rc in $(( $(date +%s) - s )) s"; grep -a "=== done" $OUT/driver.log || { tail -15 $OUT/driver.log; fail=1; }
python3 -I scripts/make_cplot.py --run-dir $OUT --prefix tpv26 --problem TPV26 \
   --code-version "$(git rev-parse HEAD)" --element-size "1000 m, p = 1 (smoke)" --out $OUT/cplot.dat
python3 -I - $OUT/cplot.dat <<'EOF' || fail=1
import math, sys
rows = [tuple(map(float, l.split())) for l in open(sys.argv[1])
        if len(l.split()) == 3 and l[0] in '-0123456789']
r = lambda j, k: math.hypot(j + 5000.0, k - 10000.0)      # hypocenter j = -5 km, k = 10 km
inside = [t for j, k, t in rows if r(j, k) < 4000.0]
outside = [(r(j, k), t) for j, k, t in rows if r(j, k) >= 4000.0]
hit = [x for x in outside if x[1] < 1e9]
print(f"B1: {sum(t < 1e9 for t in inside)} of {len(inside)} nodes inside r_crit and {len(hit)} of "
      f"{len(outside)} outside ruptured by 2.5 s" +
      (f"; farthest {max(hit)[0] / 1e3:.1f} km at {max(hit)[1]:.3f} s" if hit else ""))
sys.exit(0 if hit else 1)
EOF
ls $OUT/tpv26_offfault_*.dat 2>/dev/null | sed 's/^/  /'
echo "=== R3 output check: TPV27 1000 m p=1, 0 -> 1.0 s, 8 ranks ==="
OUT=$SCRATCH/rheology/r3out_tpv27_1000m_p1_$SLURM_JOB_ID; mkdir -p $OUT
s=$(date +%s)
ibrun -n 8 -o 0 ./seas_nitsche_dyn_driver --config tpv27/configs/tpv27_nitsche_kv0p1.toml \
   --order 1 --tfinal 1.0 --output-dir $OUT --log-every 200 > $OUT/driver.log 2>&1
rc=$?; echo "driver exit $rc in $(( $(date +%s) - s )) s"; grep -a "=== done" $OUT/driver.log || { tail -15 $OUT/driver.log; fail=1; }
grep -a "rheology" $OUT/driver.log | head -3
command -v h5dump >/dev/null || echo "ERROR: h5dump (phdf5 module) not on PATH"
ns=$(h5dump -a /VTKHDF/Steps/NSteps $OUT/rheology.vtkhdf 2>/dev/null | awk '/\(0\)/{print $2}')
echo "rheology.vtkhdf: NSteps ${ns:-missing}"
[ "$ns" = "3" ] || { echo "ERROR: rheology.vtkhdf should hold 3 snapshots"; fail=1; }
echo "=== 200 m load check: TPV26 p=2, 144 ranks, tfinal 0.003 s ==="
OUT=$SCRATCH/rheology/loadcheck_tpv26_200m_p2_$SLURM_JOB_ID; mkdir -p $OUT
( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $OUT/meminfo_kB.log &
MON=$!
s=$(date +%s)
ibrun ./seas_nitsche_dyn_driver --config tpv26/configs/tpv26_nitsche_kv0p1.toml \
   --mesh tpv26/mesh/tpv2627_mirror_200m.msh --order 2 --tfinal 0.003 --output-dir $OUT --log-every 1 \
   > $OUT/driver.log 2>&1
rc=$?; kill $MON 2>/dev/null; wait $MON 2>/dev/null
echo "driver exit $rc in $(( $(date +%s) - s )) s"
grep -a "tets;" $OUT/driver.log | head -1
grep -a "  dt = " $OUT/driver.log | head -1
grep -a "  dt = " $OUT/driver.log | head -1 | awk '{for(i=1;i<=NF;i++){if($i=="dt"&&$(i+1)=="=")dt=$(i+2); if($i=="C(p)")C=$(i+1)}} END{h=dt*4*6000/(0.5*C); printf "implied h_min (min 6V/A) %.2f m; steps to 13 s %d\n", h, int(13/dt+0.999)}'
grep -a "=== done" $OUT/driver.log || tail -15 $OUT/driver.log
awk '{u=($2-$3)/1048576; if(u>m)m=u; t=$2/1048576} END{printf "node memory: peak in use %.1f GB of %.1f GB (MemTotal - MemAvailable, 2 s samples)\n", m, t}' $OUT/meminfo_kB.log
echo "total $(( $(date +%s) - T0 )) s"
[ $fail -eq 0 ] || { echo "ERROR: a check failed; do not submit the b0_*.sbatch or tpv26_*.sbatch jobs"; exit 1; }
echo "all checks passed"
