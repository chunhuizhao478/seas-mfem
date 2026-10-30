#!/bin/bash
#SBATCH -J l2_build
#SBATCH -o l2_build_%j.out
#SBATCH -e l2_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 02:00:00
#SBATCH -A EAR20006
# The 2.5D note's Level 2 test on Vista, part 1 (Chunhui Zhao, 2026-10-30;
# /mnt/project-files/2d-inplane/nonplanar-literature/
# inplane-2p5d-derivation-rev8.pdf, Section 7.2): build seas_nitsche_dyn_driver
# in the checkout of branch 2d-inplane ($WORK/seas-mfem-2d-inplane), mesh the
# 3-D reference with TACC gmsh, and check both before the run
# (elongated_3d_vista.sbatch).  Submit from miniapps/seas/elongated2p5d/jobs/.
# Steps (logs in <checkout>-runs/l2_build_<jobid>):
#   0. the repo CLAUDE.md guard (GitHub first; no untracked sources)
#   1. jobs/fd_cycle/build_vista.sh, stages deps and mfem (PETSc and MFEM as
#      jobs/2d_inplane/vista_build.sbatch; MFEM without HDF5, so the decks
#      keep paraview_fault = "off")
#   2. seas_nitsche_dyn_driver from fresh objects with
#      tests/2d_inplane/2d_inplane.mk; the commit goes to
#      BUILD_COMMIT_VISTA_NITSCHE (BUILD_COMMIT_VISTA stays the cycle
#      driver's)
#   3. make test-2d5-dyn-driver: the [finite_width] refusals and accepted
#      tables on the TPV205-2D deck (the 3-D refusal is skipped: no MESHES)
#   4. meshes from elongated2p5d/mesh/elongated_3d.geo (built-in kernel; TACC
#      gmsh has no OpenCASCADE) at lc_fault 1000 m and 400 m, msh 2.2; sha256
#      to elongated2p5d/MESH_SHA256_VISTA, which the run checks.  Cloud gmsh
#      4.15.2 gave 536 k tets (1000 m) and 2,409,373 tets, 25,660 fault faces,
#      min 6V/A 58.4 m (400 m, 2 min 48 s); Vista's output can differ slightly
#   5. smoke: elongated_3d_buried.toml on the 1000 m mesh, p = 2, 0.5 s, 144
#      ranks (the cloud printed dt 6.90e-4 s for this mesh)
#   6. 400 m load check: the same deck on the 400 m mesh, p = 2, 0.003 s, 144
#      ranks on this one node; node memory every 2 s (Committed_AS against
#      CommitLimit, jobs/fd_cycle/vista_memlog.sh).  The run uses 16 nodes,
#      so one node holding it leaves a margin of 16.  Prints dt and the steps
#      to 40 s
# Exits non-zero when any step fails.  Cost: 1 node x <= 2 h on gg
# = <= 0.66 SU (EAR20006).
# Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
DEPS_ROOT=${DEPS_ROOT:-$ROOT-deps}
PETSC_PREFIX=$DEPS_ROOT/petsc-${PETSC_VERSION:-3.24.6}
RUNS=$ROOT-runs/l2_build_$SLURM_JOB_ID
mkdir -p "$RUNS" || exit 2
cd "$ROOT" || exit 2

# 0. GitHub-first guard (CLAUDE.md); no fetch on a compute node
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
UNTRACKED=$(git ls-files --others --exclude-standard -- '*.cpp' '*.hpp' '*.inl' Makefile '*.mk' '*.toml' '*.sbatch')
[ -z "$UNTRACKED" ] || { echo "ERROR: untracked sources in $ROOT:"; echo "$UNTRACKED"; exit 2; }
HEAD=$(git rev-parse HEAD)
echo "commit $HEAD"
echo "host $(hostname)  logs $RUNS"
T0=$(date +%s)
fails=0

# 1. PETSc and MFEM (the guard runs inside build_vista.sh too)
JOBS=64 STAGE="deps mfem" DEPS_ROOT="$DEPS_ROOT" \
   bash "$ROOT/miniapps/seas/jobs/fd_cycle/build_vista.sh" > "$RUNS/build_deps_mfem.log" 2>&1
rc=$?; echo "build_vista.sh deps mfem exit $rc after $(( $(date +%s) - T0 )) s"; tail -5 "$RUNS/build_deps_mfem.log"
[ $rc -eq 0 ] || exit $rc

# 2. The Nitsche driver, from fresh objects
module purge
for m in $MODULES; do module load $m || exit 2; done
cd "$ROOT/miniapps/seas" || exit 2
case "$(git -C "$ROOT" submodule status miniapps/seas/extern/toml11)" in
   " "*) ;;
   *) echo "ERROR: toml11 submodule missing or not at its recorded commit"; exit 2 ;;
esac
find . -name '*.o' -not -path './extern/*' -delete
MK="make -f tests/2d_inplane/2d_inplane.mk SEAS_FS_LIB= HDF5_DIR=$PETSC_PREFIX"
s=$(date +%s)
$MK -j64 seas_nitsche_dyn_driver > "$RUNS/build_seas.log" 2>&1; rc=$?
echo "seas_nitsche_dyn_driver build exit $rc in $(( $(date +%s) - s )) s, $(grep -c warning "$RUNS/build_seas.log") warning lines"
[ $rc -eq 0 ] || { tail -30 "$RUNS/build_seas.log"; exit $rc; }
echo "$HEAD" > BUILD_COMMIT_VISTA_NITSCHE

# 3. The [finite_width] tests of the driver
s=$(date +%s)
TMPDIR=$RUNS $MK test-2d5-dyn-driver > "$RUNS/test-2d5-dyn-driver.log" 2>&1; rc=$?
echo "test-2d5-dyn-driver exit $rc in $(( $(date +%s) - s )) s"
grep -E "PASSED|FAILED|SKIPPED" "$RUNS/test-2d5-dyn-driver.log" | cut -c1-140
[ $rc -eq 0 ] || fails=$((fails + 1))

# 4. The meshes, with TACC gmsh in its own module set
for lc in 1000 400; do
   m=elongated2p5d/mesh/elongated_3d_${lc}m.msh
   s=$(date +%s)
   ( module purge; module load TACC
     module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
     module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6 gmsh/4.15.0
     gmsh --version 2>&1 | head -1
     gmsh -3 elongated2p5d/mesh/elongated_3d.geo -format msh22 -setnumber lc_fault $lc -o "$m" ) \
      > "$RUNS/gmsh_${lc}m.log" 2>&1
   rc=$?
   echo "$m: gmsh exit $rc in $(( $(date +%s) - s )) s; $(grep -m1 -E '^[0-9]+\.[0-9]+' "$RUNS/gmsh_${lc}m.log"); $(grep -E 'Info *: [0-9]+ nodes [0-9]+ elements' "$RUNS/gmsh_${lc}m.log" | tail -1)"
   [ $rc -eq 0 ] || { tail -10 "$RUNS/gmsh_${lc}m.log"; exit $rc; }
done
sha256sum elongated2p5d/mesh/elongated_3d_1000m.msh elongated2p5d/mesh/elongated_3d_400m.msh \
   | tee elongated2p5d/MESH_SHA256_VISTA

# 5. Smoke: 1000 m, p = 2, 0.5 s, 144 ranks
DECK=elongated2p5d/configs/elongated_3d_buried.toml
OUT=$RUNS/smoke_1000m
mkdir -p "$OUT" || exit 2
s=$(date +%s)
ibrun ./seas_nitsche_dyn_driver --config "$DECK" --mesh elongated2p5d/mesh/elongated_3d_1000m.msh \
   --order 2 --tfinal 0.5 --output-dir "$OUT" --log-every 200 > "$OUT/driver.log" 2>&1
rc=$?
echo "smoke 1000 m exit $rc in $(( $(date +%s) - s )) s"
grep -m1 -E "^  dt = " "$OUT/driver.log"; grep "=== done" "$OUT/driver.log" || tail -15 "$OUT/driver.log"
[ $rc -eq 0 ] || fails=$((fails + 1))

# 6. 400 m load check: p = 2, 0.003 s, 144 ranks on this node
OUT=$RUNS/load_400m
mkdir -p "$OUT" || exit 2
"$ROOT/miniapps/seas/jobs/fd_cycle/vista_memlog.sh" 2 > "$OUT/memlog.txt" 2>&1 &
MEMLOG=$!
s=$(date +%s)
ibrun ./seas_nitsche_dyn_driver --config "$DECK" --mesh elongated2p5d/mesh/elongated_3d_400m.msh \
   --order 2 --tfinal 0.003 --output-dir "$OUT" --log-every 1 > "$OUT/driver.log" 2>&1
rc=$?
kill $MEMLOG 2>/dev/null
echo "load check 400 m exit $rc in $(( $(date +%s) - s )) s"
grep "=== done" "$OUT/driver.log" || tail -15 "$OUT/driver.log"
awk '{ if ($6 > p) { p = $6; lim = $4 } } END { printf "node peak committed %.1f GiB of %.1f GiB\n", p, lim }' "$OUT/memlog.txt"
DT=$(grep -m1 -E "^  dt = " "$OUT/driver.log" | awk '{print $3}')
[ -n "$DT" ] && awk -v dt="$DT" 'BEGIN { printf "dt %.4g s: %d steps to 40 s\n", dt, 40 / dt + 0.5 }'
[ $rc -eq 0 ] || fails=$((fails + 1))

echo "total $(( $(date +%s) - T0 )) s, failing steps: $fails"
[ $fails -eq 0 ] || exit 1
