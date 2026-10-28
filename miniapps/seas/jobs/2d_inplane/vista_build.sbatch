#!/bin/bash
#SBATCH -J bp3fd_build
#SBATCH -o bp3fd_build_%j.out
#SBATCH -e bp3fd_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 04:00:00
#SBATCH -A EAR20006
# 2D in-plane BP3-FD on Vista: build seas_fd_cycle_driver on one gg node in
# a checkout of branch 2d-inplane and check it on the BP3-FD path (Chunhui
# Zhao, 2026-10-27: "cancel this run, add snapshots, submit it on cluster").
# On the login node first (it has internet; the compute nodes may not):
#   cd $WORK
#   git clone --recurse-submodules -b 2d-inplane \
#       https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-2d-inplane
#   cd seas-mfem-2d-inplane && STAGE=fetch bash miniapps/seas/jobs/fd_cycle/build_vista.sh
#   bash miniapps/seas/jobs/2d_inplane/submit_vista_bp3_fd.sh
# Steps (logs in <checkout>-runs/build_<jobid>, next to the checkout):
#   0. python3 with numpy, h5py and tomllib (CYCLE_PY, as
#      jobs/fd_cycle/vista_build.sbatch): step 4 reads the snapshots back
#   1. jobs/fd_cycle/build_vista.sh, stages deps and mfem: PETSc 3.24.6 with
#      MUMPS and HDF5 into <checkout>-deps, MFEM in place (the GitHub-first
#      guard runs inside, before anything is built)
#   2. seas_fd_cycle_driver from fresh objects, then BUILD_COMMIT_VISTA.
#      build_vista.sh's seas stage is not used: its cycle-all ends with
#      fully-dynamic-cycle's rule M1 check (cycle-owned-check), which refuses
#      this branch's 2D files by design
#   3. the 2D tests of the BP3-FD path (tests/2d_inplane/2d_inplane.mk):
#      test-2d-deck, test-prestress-bp3 (B7: the five BP3-FD decks on their
#      1 km meshes, 1 and 2 ranks), test-2d-absorbing and test-absorbing-all
#      (the dynamic regime's walls, 1 and 2 ranks)
#   4. snapshots on a 2D mesh: config/fd_cycle/bp3_fd_d60_thrust.toml on its
#      1 km mesh to the end of the first dynamic phase on 2 ranks (A), and
#      again stopped at QD step 1200 and restarted from that checkpoint (B);
#      tests/cycle/check_snapshots.py reads both back (layout, the station
#      rows, the schedule) and compares B with A (MUMPS on 2 ranks: 3e-9 in
#      the cloud, tolerance 1e-8).  No python vtk here: --no-vtk
# Exits non-zero when any step fails.  Cost: 1 node x <= 4 h on gg
# = <= 1.3 SU (EAR20006).
# Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
DEPS_ROOT=${DEPS_ROOT:-$ROOT-deps}
PETSC_PREFIX=$DEPS_ROOT/petsc-${PETSC_VERSION:-3.24.6}
RUNS=$ROOT-runs/build_$SLURM_JOB_ID
mkdir -p "$RUNS" || exit 2
echo "host $(hostname)  cpu $(lscpu | grep -m1 'Model name' | xargs)  logs $RUNS"
T0=$(date +%s)

# 0. The tests' python (jobs/fd_cycle/vista_build.sbatch, step 0)
CYCLE_PY=${CYCLE_PY:-${WORK:?}/seissol/venv-post/bin/python3}
CYCLE_PY_LIB=${CYCLE_PY_LIB:-/opt/apps/gcc14/cuda12/python3/3.11.8/lib}
mkdir -p "$RUNS/pybin" || exit 2
printf '#!/bin/bash\nexport LD_LIBRARY_PATH="%s${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"\nexec "%s" "$@"\n' \
   "$CYCLE_PY_LIB" "$CYCLE_PY" > "$RUNS/pybin/python3" && chmod +x "$RUNS/pybin/python3" || exit 2
"$RUNS/pybin/python3" -c "import numpy, h5py, tomllib; print('python', __import__('sys').version.split()[0], 'numpy', numpy.__version__, 'h5py', h5py.__version__)" \
   || { echo "ERROR: $CYCLE_PY (lib $CYCLE_PY_LIB) cannot import numpy, h5py and tomllib; set CYCLE_PY"; exit 2; }

# 1. PETSc and MFEM (the GitHub-first guard runs inside)
JOBS=64 STAGE="deps mfem" DEPS_ROOT="$DEPS_ROOT" \
   bash "$ROOT/miniapps/seas/jobs/fd_cycle/build_vista.sh" > "$RUNS/build_deps_mfem.log" 2>&1
rc=$?; echo "build_vista.sh deps mfem exit $rc after $(( $(date +%s) - T0 )) s"; tail -5 "$RUNS/build_deps_mfem.log"
[ $rc -eq 0 ] || exit $rc

# 2. The driver, from fresh objects (the shared Makefile has no header
# dependencies), with the toolchain of step 1
for m in $MODULES; do module load $m || exit 2; done
export PATH="$RUNS/pybin:$PATH"
cd "$ROOT/miniapps/seas" || exit 2
case "$(git -C "$ROOT" submodule status miniapps/seas/extern/toml11)" in
   " "*) ;;
   *) echo "ERROR: toml11 submodule missing or not at its recorded commit"; exit 2 ;;
esac
find . -name '*.o' -not -path './extern/*' -delete
MK="make -f tests/2d_inplane/2d_inplane.mk SEAS_FS_LIB= HDF5_DIR=$PETSC_PREFIX"
s=$(date +%s)
$MK -j64 seas_fd_cycle_driver > "$RUNS/build_seas.log" 2>&1; rc=$?
echo "seas_fd_cycle_driver build exit $rc in $(( $(date +%s) - s )) s, $(grep -c warning "$RUNS/build_seas.log") warning lines"
[ $rc -eq 0 ] || { tail -30 "$RUNS/build_seas.log"; exit $rc; }
git -C "$ROOT" rev-parse HEAD > BUILD_COMMIT_VISTA
echo "commit $(cat BUILD_COMMIT_VISTA)"
fails=0

# 3. The 2D tests of the BP3-FD path
for t in test-2d-deck test-prestress-bp3 test-2d-absorbing test-absorbing-all; do
   s=$(date +%s)
   $MK $t > "$RUNS/$t.log" 2>&1; rc=$?
   echo "$t exit $rc in $(( $(date +%s) - s )) s"
   grep -E "passed ===|[0-9]+ / [0-9]+ passed|: PASSED \(|: FAILED|^  FAIL|  FAILED" "$RUNS/$t.log" | tail -8
   [ $rc -eq 0 ] || fails=$((fails + 1))
done

# 4. Snapshots on a 2D mesh: a run to the end of event 1, and the same run
# stopped and restarted
S4=$RUNS/snapshots
mkdir -p "$S4/a" "$S4/b" || exit 2
sed -e 's|^path  = "bp3/mesh/bp3_d60_25m.msh"$|path  = "bp3/mesh/bp3_d60_1km.msh"|' \
    -e 's/^max_dynamic_phases = 0$/max_dynamic_phases = 1/' \
    config/fd_cycle/bp3_fd_d60_thrust.toml > "$S4/d60_1km.toml"
diff config/fd_cycle/bp3_fd_d60_thrust.toml "$S4/d60_1km.toml"
grep -q '^path  = "bp3/mesh/bp3_d60_1km.msh"$' "$S4/d60_1km.toml" \
   && grep -q '^max_dynamic_phases = 1$' "$S4/d60_1km.toml" \
   || { echo "ERROR: the 1 km deck lacks its change (the sed patterns no longer match)"; exit 2; }
P=bp3_fd_d60_thrust
s=$(date +%s)
mpirun -n 2 ./seas_fd_cycle_driver --deck "$S4/d60_1km.toml" --output-dir "$S4/a" > "$S4/a/run.log" 2>&1
ra=$?
mpirun -n 2 ./seas_fd_cycle_driver --deck "$S4/d60_1km.toml" --output-dir "$S4/b" --max-steps 1200 \
   > "$S4/b/run1.log" 2>&1
rb1=$?
mpirun -n 2 ./seas_fd_cycle_driver --deck "$S4/d60_1km.toml" --output-dir "$S4/b" --restart "$S4/b/$P" \
   > "$S4/b/run2.log" 2>&1
rb2=$?
echo "snapshot runs: A exit $ra, B exit $rb1 then $rb2 (restart), $(( $(date +%s) - s )) s"
grep -E "snapshots:|absorbing boundary:|QD -> dynamic|end: " "$S4/a/run.log"
grep -E "restarted from" "$S4/b/run2.log"
DT=$(grep -m1 -oE "steps of [0-9.e+-]+ s" "$S4/a/run.log" | awk '{print $3}')
rc=1
if [ $ra -eq 0 ] && [ $rb1 -eq 0 ] && [ $rb2 -eq 0 ] && [ -n "$DT" ]; then
   python3 tests/cycle/check_snapshots.py --dir "$S4/a" --prefix $P --deck "$S4/d60_1km.toml" \
      --dt-dyn-step "$DT" --no-vtk > "$S4/check_a.log" 2>&1 \
   && python3 tests/cycle/check_snapshots.py --dir "$S4/b" --prefix $P --deck "$S4/d60_1km.toml" \
      --dt-dyn-step "$DT" --no-vtk --against "$S4/a" --tol 1e-8 > "$S4/check_b.log" 2>&1
   rc=$?
fi
echo "snapshot check exit $rc"
grep -hE "layout|rows matched| vs |ALL OK|FAIL" "$S4/check_a.log" "$S4/check_b.log" 2>/dev/null
[ $rc -eq 0 ] || fails=$((fails + 1))

echo "total $(( $(date +%s) - T0 )) s, failing steps: $fails"
[ $fails -eq 0 ] || exit 1
