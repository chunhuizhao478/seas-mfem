#!/bin/bash
#SBATCH -J d2_mf_bench
#SBATCH -o d2_mf_bench_%j.out
#SBATCH -e d2_mf_bench_%j.err
#SBATCH -p gh
#SBATCH -N 1
#SBATCH -n 72
#SBATCH -t 04:00:00
#SBATCH -A EAR20006

# GPU dynamic-rupture plan D2 (plan rev2 §4.3), first step: the matrix-free
# SIPG product (domain/nitsche_sipg_operator.hpp) on one GH200 and on its
# Grace CPU, with tests/unit/test_nitsche_sipg_mf.cpp on TPV205 1000 m.
#   checks  y = K u, K_eta v and K u + K_eta v against the assembled matrices
#           (hypre on the host) at 1e-13 (plan rev2 §5), p = 1 and 2, both
#           kernel layouts, heterogeneous moduli and the absorbing boundary
#           made Dirichlet (every side kind): 1 rank on the GPU; 4 ranks
#           sharing it (kernel B and the halo), only when the GPU's compute
#           mode is Default; 72 ranks on Grace.
#   bench   the volume pass, the face pass and both, median of 50 calls, as
#           the decks run (no Dirichlet faces): the GPU with both layouts and
#           team batches 2, 4, 8; Grace, 72 ranks, the element layout and the
#           dense blocks (item 8a) on the same matrices.
# The decision number is the GPU's GFLOP/s of K u + K_eta v at p = 2 (model
# count, NitscheSipgOperator::Flops): below about 1500 the plan reports and
# asks (plan rev2 §4.3).  SUMMARY.txt in $OUT collects the tables.
#
# It builds in the clone it is submitted from ($ROOT) at its HEAD with
# build_horizon.sh MODE=cuda (MFEM configured in place; the CUDA PETSc of
# $WORK/seas-gpu-deps and the HDF5 of $WORK/seas-gpu-dyn-deps, both read
# only), and stops if another job of yours is queued or running there.  Do
# not pull into a clone that queued jobs still use: the D1 GPU checks in
# $WORK/seas-mfem-gpu-dyn-cuda fail when its HEAD moves.  Login node, in a
# clone of its own:
#   git clone -b gpu-dyn-nitsche https://github.com/chunhuizhao478/seas-mfem.git \
#      $WORK/seas-mfem-gpu-dyn-d2
#   cd $WORK/seas-mfem-gpu-dyn-d2/miniapps/seas && sbatch jobs/vista/d2_mf_bench.sbatch
# or, after the D1 GPU checks have finished, in the D1 CUDA clone:
#   cd $WORK/seas-mfem-gpu-dyn-cuda && git pull --ff-only
#   cd miniapps/seas && sbatch jobs/vista/d2_mf_bench.sbatch
# The meshes are D0a's, copied from $WORK/seas-mfem-gpu-dyn-nit when absent
# (sha256 checked).  A failing check or build stops nothing else: every case
# runs and the job exits 1 when any failed.
# Cost: 1 gh node x <= 4 h (1 SU per node-hour); the build is most of it.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ROOT=$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)              # submitted from miniapps/seas
NIT=${NIT:-$WORK/seas-mfem-gpu-dyn-nit}
# no other job of this user may use the clone this job builds in (the D1
# checks compare BUILD_COMMIT_CUDA with HEAD and run the clone's binaries)
Q=$(squeue -h -u "$USER" -o "%i %Z") || { echo "ERROR: squeue failed"; exit 2; }
OTHERS=$(awk -v me="$SLURM_JOB_ID" -v r="$ROOT" '$1 != me && (($2 == r) || index($2, r "/") == 1) {print $1}' <<< "$Q" | tr '\n' ' ')
[ -z "$OTHERS" ] || { echo "ERROR: jobs $OTHERS are queued or running in $ROOT; use another clone or wait"; exit 2; }
# shellcheck source=d1_env.sh
source "$ROOT/miniapps/seas/jobs/vista/d1_env.sh"
T0=$(date +%s)

d1_guard "$ROOT"
[ "$(git -C "$ROOT" rev-parse HEAD)" = "$(git -C "$ROOT" rev-parse '@{u}' 2>/dev/null)" ] \
   || { echo "ERROR: $ROOT HEAD is not its upstream; pull first"; exit 2; }
SHA=$(git -C "$ROOT" rev-parse HEAD)
OUT=$SCRATCH/gpu_dyn/d2_mf_bench_$SLURM_JOB_ID
mkdir -p "$OUT"
echo "job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); checkout $ROOT at $SHA"

cd "$ROOT" || exit 2
if ! MODE=cuda GPU_ARCH=90 STAGE="hdf5 mfem seas" WITH_HDF5=1 \
      DEPS_ROOT="$WORK/seas-gpu-deps" HDF5_ROOT="$WORK/seas-gpu-dyn-deps" \
      PROVENANCE_FILE="$OUT/build_provenance_cuda.txt" \
      SEAS_TARGETS="seas_test_nitsche_sipg_mf" JOBS=72 bash build_horizon.sh > "$OUT/build.log" 2>&1; then
   echo "BUILD: FAIL ($OUT/build.log)"; tail -40 "$OUT/build.log"; exit 1
fi
grep -q '^MFEM_USE_CUDA *= YES' config/config.mk || { echo "ERROR: $ROOT is not a CUDA build"; exit 2; }
echo "BUILD: done at $SHA in $(( $(date +%s) - T0 )) s"

cd "$ROOT/miniapps/seas" || exit 2
for m in tpv205/mesh/tpv2053d_1000m.msh tpv104/mesh/tpv104_1000m.msh; do
   [ -f "$m" ] || cp -v "$NIT/miniapps/seas/$m" "$m"
done
d1_check_meshes "$ROOT/miniapps/seas" || exit 2
d1_gh_modules
CM=$(nvidia-smi --query-gpu=compute_mode --format=csv,noheader | head -1)
{ echo "job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ), commit $SHA ($ROOT)"
  echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
  nvidia-smi --query-gpu=name,compute_mode,driver_version,clocks.max.sm --format=csv,noheader
  lscpu | grep -E 'Model name|^CPU\(s\)'; } | tee "$OUT/RUN_INFO.txt"

MESH=tpv205/mesh/tpv2053d_1000m.msh
fail=0
# run NAME NRANKS ARGS...: one test run, its log in $OUT/NAME.log
run () {
   local name=$1 n=$2 t rc
   shift 2
   t=$(date +%s)
   ibrun -n "$n" -o 0 ./seas_test_nitsche_sipg_mf --mesh "$MESH" "$@" > "$OUT/$name.log" 2>&1
   rc=$?
   echo "$name ($n rank(s)): exit $rc in $(( $(date +%s) - t )) s; $(grep -a 'PASSED, ' "$OUT/$name.log" | tail -1)"
   [ "$rc" = 0 ] || fail=1
}

CHECK=(--hetero 0.2 --dirichlet 105 --layout both)
for P in 1 2; do
   run "gpu_check_p$P" 1 --order "$P" --device cuda "${CHECK[@]}" --batch 2,4,8
   if [ "$CM" = "Default" ]; then
      run "gpu_share_check_p$P" 4 --order "$P" --device cuda "${CHECK[@]}" --batch 4
   else
      echo "gpu_share_check_p$P: NOT RUN (compute mode $CM: ranks cannot share the GPU); kernel B and the halo on the GPU wait for the multi-GPU runs"
   fi
   run "grace_check_p$P" 72 --order "$P" --device cpu "${CHECK[@]}" --batch 4
   run "gpu_bench_p$P" 1 --order "$P" --device cuda --no-check --reps 50 --layout both --batch 2,4,8
   run "grace_bench_p$P" 72 --order "$P" --device cpu --reps 50 --layout element --dense
done

{ echo "D2 matrix-free microbenchmark, job $SLURM_JOB_ID, commit $SHA, $MESH"
  for f in gpu_bench_p1 gpu_bench_p2 grace_bench_p1 grace_bench_p2; do
     echo; echo "== $f"
     grep -a -E '^matrix-free setup|^microbenchmark|layout: product|^  (element|team[0-9]): |dense blocks' "$OUT/$f.log"
  done
  echo; echo "== checks"
  grep -a -H -E 'PASSED:|FAILED:|skipped|PASSED, ' "$OUT"/*check*.log | sed "s|$OUT/||"; } > "$OUT/SUMMARY.txt"
cat "$OUT/SUMMARY.txt"
echo "D2 bench: $([ $fail = 0 ] && echo PASS || echo FAIL); commit $SHA; logs in $OUT; total $(( $(date +%s) - T0 )) s"
exit $fail
