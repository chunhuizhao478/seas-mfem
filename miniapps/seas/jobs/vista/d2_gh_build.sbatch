#!/bin/bash
#SBATCH -J d2_gh_build
#SBATCH -o d2_gh_build_%j.out
#SBATCH -e d2_gh_build_%j.err
#SBATCH -p gh
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -A EAR20006

# GPU dynamic-rupture plan D2 (plan rev2 §4.3), the device step: the CUDA
# build of seas_nitsche_dyn_driver and seas_test_nitsche_sipg_mf in the clone
# it is submitted from ($ROOT, at its HEAD), then the unit checks on the GPU:
#   [K] [E] [KE]  the matrix-free product against the assembled matrices
#                 (1e-13), both layouts, team batches 2 and 4
#   [M] [MC]      the device mass inverse, plain and with damped blocks,
#                 against the host NitscheMass (1e-14, plan rev2 §4.3)
#   [Gd] [T]      the device slip load against AddSlipLoadFast, and the
#                 fault traction from the values gathered off the GPU
# on TPV205 1000 m, p = 1 and 2, heterogeneous moduli, with the absorbing
# boundary made Dirichlet (every face kind of the product): 1 rank, and 4
# ranks sharing the GPU (the halo and shared fault faces) when its compute
# mode is Default, or each on its own GPU on a node of 4 (Horizon's GB200).
# BUILD_COMMIT_CUDA records the commit for d2_gh_gate.sbatch.
# Then the transfer trace (plan rev2 §4.3 acceptance, "an Nsight trace shows
# no volume-sized copy inside a step, only the planned fault gather"):
# TPV205 1000 m on 1 GPU, p = 1 and 2, 200 steps of the device step under
# Nsight Systems, checked by jobs/vista/d2_trace_check.py --gate.  Nsight
# Systems is NSYS if set, else nsys on PATH, else the CUDA module's.
# GPU dynamic-rupture plan D3 (plan rev2 §4.4), the fault on the GPU: the
# unit checks also run [DM] [DT] (the device fault maps and traction against
# ExpandOwnedToLocalFault + ToKernel and ComputeTractionFast, 1e-13), at 2
# ranks too; seas_test_device_closures runs the device Nitsche closures
# (slip weakening, aging, slip law, strong rate weakening) on the GPU against
# the host closures, 10^6 random points per law and an edge set (1e-12 in
# V); and the trace is repeated with --device-fault 1, where
# d2_trace_check.py --gate fails on any copy in a step of 8 bytes per fault
# DOF or more (plan rev2 §6, "no copy of 8 bytes per fault DOF or more").
#
# MFEM is configured in place with build_horizon.sh MODE=cuda (the CUDA
# PETSc of DEPS_ROOT, read only, and the HDF5 of HDF5_ROOT, built there when
# absent; on Vista $WORK/seas-gpu-deps and $WORK/seas-gpu-dyn-deps).  It stops
# if another job of yours is queued or running in the clone, unless that job
# waits on this one (a gate submitted after it).  Login node, in the D2 clone:
#   cd $WORK/seas-mfem-gpu-dyn-d2 && git pull --ff-only
#   cd miniapps/seas && sbatch jobs/vista/d2_gh_build.sbatch
# The meshes are D0a's, copied from MESH_SRC (default the miniapps/seas of
# $WORK/seas-mfem-gpu-dyn-nit) when absent (sha256 checked).  On Horizon,
# jobs/horizon/h_d2_build.sbatch runs this script with the settings of
# jobs/horizon/h_dyn_env.sh (GPU_ARCH, MFEM_CUDA_ARCH, DEPS_ROOT, HDF5_ROOT,
# MESH_SRC, the modules and NSYS).  Every check runs; the job exits 1 when
# any failed.
# Cost: 1 gh node x <= 4 h (1 SU per node-hour); the build is most of it.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ROOT=$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)              # submitted from miniapps/seas
NIT=${NIT:-$WORK/seas-mfem-gpu-dyn-nit}
MESH_SRC=${MESH_SRC:-$NIT/miniapps/seas}
# Other jobs in this clone, except those still waiting on this one (their
# dependency names this job: "afterok:<id>(unfulfilled)"), which start on
# the binary this job builds.
Q=$(squeue -h -u "$USER" -o "%i %E %Z") || { echo "ERROR: squeue failed"; exit 2; }
OTHERS=$(awk -v me="$SLURM_JOB_ID" -v r="$ROOT" '$1 != me && index($2, ":" me "(") == 0 && (($3 == r) || index($3, r "/") == 1) {print $1}' <<< "$Q" | tr '\n' ' ')
[ -z "$OTHERS" ] || { echo "ERROR: jobs $OTHERS are queued or running in $ROOT; use another clone or wait"; exit 2; }
# shellcheck source=d1_env.sh
source "$ROOT/miniapps/seas/jobs/vista/d1_env.sh"
T0=$(date +%s)

d1_guard "$ROOT"
[ "$(git -C "$ROOT" rev-parse HEAD)" = "$(git -C "$ROOT" rev-parse '@{u}' 2>/dev/null)" ] \
   || { echo "ERROR: $ROOT HEAD is not its upstream; pull first"; exit 2; }
SHA=$(git -C "$ROOT" rev-parse HEAD)
OUT=$SCRATCH/gpu_dyn/d2_gh_build_$SLURM_JOB_ID
mkdir -p "$OUT"
echo "job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); checkout $ROOT at $SHA"

cd "$ROOT" || exit 2
rm -f miniapps/seas/BUILD_COMMIT_CUDA
if ! MODE=cuda GPU_ARCH="${GPU_ARCH:-90}" STAGE="hdf5 mfem seas" WITH_HDF5=1 \
      DEPS_ROOT="${DEPS_ROOT:-$WORK/seas-gpu-deps}" HDF5_ROOT="${HDF5_ROOT:-$WORK/seas-gpu-dyn-deps}" \
      PROVENANCE_FILE="$OUT/build_provenance_cuda.txt" \
      SEAS_TARGETS="seas_nitsche_dyn_driver seas_test_nitsche_sipg_mf seas_test_device_closures" JOBS=72 \
      bash build_horizon.sh > "$OUT/build.log" 2>&1; then
   echo "BUILD: FAIL ($OUT/build.log)"; tail -40 "$OUT/build.log"; exit 1
fi
grep -q '^MFEM_USE_CUDA *= YES' config/config.mk || { echo "ERROR: $ROOT is not a CUDA build"; exit 2; }
git rev-parse HEAD > miniapps/seas/BUILD_COMMIT_CUDA
echo "BUILD: done at $SHA in $(( $(date +%s) - T0 )) s"

cd "$ROOT/miniapps/seas" || exit 2
for m in tpv205/mesh/tpv2053d_1000m.msh tpv104/mesh/tpv104_1000m.msh; do
   [ -f "$m" ] || cp -v "$MESH_SRC/$m" "$m"
done
d1_check_meshes "$ROOT/miniapps/seas" || exit 2
d1_gh_modules
CM=$(nvidia-smi --query-gpu=compute_mode --format=csv,noheader | head -1)
NGPU=$(nvidia-smi -L | grep -c '^GPU ')
{ echo "job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ), commit $SHA ($ROOT)"
  echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
  nvidia-smi --query-gpu=name,compute_mode,driver_version --format=csv,noheader; } | tee "$OUT/RUN_INFO.txt"

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
   grep -a -E 'PASSED:|FAILED:|skipped' "$OUT/$name.log" | sed 's/^/   /'
   [ "$rc" = 0 ] || fail=1
}

CHECK=(--hetero 0.2 --dirichlet 105 --layout both --batch 2,4)
for P in 1 2; do
   run "gpu_check_p$P" 1 --order "$P" --device cuda "${CHECK[@]}"
   if [ "$NGPU" -ge 4 ] || [ "$CM" = "Default" ]; then
      run "gpu_share2_check_p$P" 2 --order "$P" --device cuda "${CHECK[@]}"
      run "gpu_share_check_p$P" 4 --order "$P" --device cuda "${CHECK[@]}"
   else
      echo "gpu_share_check_p$P: NOT RUN ($NGPU GPU(s), compute mode $CM: 2 and 4 ranks cannot share them); the halo on the GPU waits for d2_gh_gate.sbatch at 2 and 4 GPUs"
   fi
done

# D3: the device closures against the host closures, on the GPU
t=$(date +%s)
ibrun -n 1 -o 0 ./seas_test_device_closures --device cuda --n 1000000 > "$OUT/device_closures.log" 2>&1
rc=$?
echo "device_closures (1 rank): exit $rc in $(( $(date +%s) - t )) s; $(tail -n 1 "$OUT/device_closures.log")"
grep -a -E 'points;|FAIL' "$OUT/device_closures.log" | sed 's/^/   /'
[ "$rc" = 0 ] || fail=1

# The transfer trace: one rank, one GPU, the device step after step 1
NSYS_BIN="${NSYS:-$(command -v nsys 2>/dev/null || true)}"
if [ -z "$NSYS_BIN" ]; then
   for r in "${CUDA_HOME:-}" "${TACC_CUDA_DIR:-}" "${CUDA_PATH:-}"; do
      [ -n "$r" ] || continue
      if [ -x "$r/bin/nsys" ]; then NSYS_BIN="$r/bin/nsys"; break; fi
      c=$(ls -d "$r"/nsight-systems-*/bin/nsys 2>/dev/null | sort -V | tail -n 1)
      if [ -n "$c" ]; then NSYS_BIN=$c; break; fi
   done
fi
if [ -z "$NSYS_BIN" ] || [ ! -x "$NSYS_BIN" ]; then
   echo "TRACE: FAIL (no Nsight Systems: set NSYS=/path/to/nsys or load a module that has it)"; fail=1
else
   echo "Nsight Systems: $NSYS_BIN ($("$NSYS_BIN" --version 2>&1 | head -1))"
   # the D2 device step (the fault on the host), then the D3 device fault
   for DF in 0 1; do
   for P in 1 2; do
      O="$OUT/trace_p$P"; [ "$DF" = 0 ] || O="$OUT/trace_df1_p$P"
      mkdir -p "$O/run"
      t=$(date +%s)
      ibrun -n 1 -o 0 "$NSYS_BIN" profile --trace=cuda,nvtx --sample=none --cpuctxsw=none \
            --force-overwrite=true -o "$O/trace" \
            ./seas_nitsche_dyn_driver --config tpv205/configs/tpv205_nitsche.toml --mesh "$MESH" \
            --order "$P" --output-dir "$O/run" --log-every 50 --device cuda --matrix-free 1 \
            --stage-timers 50 --max-steps 200 --device-fault "$DF" > "$O/driver.log" 2>&1
      rc=$?
      echo "$(basename "$O"): exit $rc in $(( $(date +%s) - t )) s; $(grep -a "=== done" "$O/driver.log" || echo 'no done line')"
      if [ "$rc" != 0 ]; then echo "TRACE $(basename "$O"): FAIL (run exit $rc, $O/driver.log)"; fail=1; continue; fi
      if [ "$DF" = 1 ] && ! grep -a -q "fault on the device: " "$O/driver.log"; then
         echo "TRACE $(basename "$O"): FAIL (the summary does not say the fault ran on the device)"; fail=1; continue
      fi
      if ! "$NSYS_BIN" export --type=sqlite --force-overwrite=true -o "$O/trace.sqlite" \
            "$O/trace.nsys-rep" > "$O/export.log" 2>&1; then
         echo "TRACE $(basename "$O"): FAIL (nsys export, $O/export.log)"; fail=1; continue
      fi
      python3 jobs/vista/d2_trace_check.py "$O/trace.sqlite" --log "$O/driver.log" --gate > "$O/check.txt" 2>&1
      crc=$?
      sed 's/^/   /' "$O/check.txt"
      echo "TRACE $(basename "$O"): $([ $crc = 0 ] && echo PASS || echo "FAIL (check exit $crc)")"
      [ "$crc" = 0 ] || fail=1
   done
   done
fi

echo "D2 and D3 gh build and checks: $([ $fail = 0 ] && echo PASS || echo FAIL); commit $SHA; logs in $OUT; total $(( $(date +%s) - T0 )) s"
exit $fail
