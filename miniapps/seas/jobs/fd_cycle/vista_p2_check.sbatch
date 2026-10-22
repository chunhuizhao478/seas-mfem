#!/bin/bash
#SBATCH -J fdc_p2_check
#SBATCH -o fdc_p2_check_%j.out
#SBATCH -e fdc_p2_check_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 05:00:00
#SBATCH -A EAR20006
# Fully dynamic BP5 cycle on Vista: the checks of the p = 2 run's speed-ups
# at its own size (4 nodes, 576 ranks), before vista_p2_1800yr.sbatch.
# Each step runs config/fd_cycle/bp5_fd_p2_1800yr.toml or a variant of it
# written into the log folder (the sed lines below are the whole change):
#   1. green   QD only (V_on = 100 m/s), 300 QD steps: stage operators,
#              the Green's function as an H-matrix, the exact QD velocity.
#              Builds and writes FDC_HCACHE when it is missing (about
#              80 min), else reads it.
#   2. loops   the same 300 steps on the face loops (no Green's function),
#              576 ranks: step 1 against it, compare_interp.py --tol 1e-5
#              on the station, off-fault and source rows (--scale quantity
#              in every comparison: each column against the size of its
#              quantity over the stations of its kind, as the dip
#              components stay at round-off)
#   3. loops144  step 2 on one node (144 ranks, its own partition): the
#              station read-out against another rank count (the QD
#              thread's open item: station -16dp+00 off 2.4e-3 m at 576
#              ranks in 1084134), same comparison against step 1
#   4. lts     the deck with [cycle.lts] mode = "rate2" (the first event
#              dynamic from t = 0), 800 sync steps (--max-dyn-steps),
#              reading the H-matrix
#   5. gts     the deck with mode = "off", the deck's own setting (global
#              stepping), 6200 steps of 0.81 ms (5 s): both regimes and
#              the H-matrix in memory at once, and the dynamic step's cost.
#              Step 4 against it in their common time range after t = 0:
#              on-fault stations and source rows --tol 1e-3, off-fault
#              stations --tol 5e-2 (test_cycle_bench_run.sh's measures).
#              That comparison counts only when the deck runs LTS; with
#              global stepping it is printed for the record.
# The runs start with TACC's ibrun; every node's committed memory goes to
# memlog.txt every 2 s (vista_memlog.sh), and each run prints its peak.
# Logs and runs in <checkout>-runs/p2check_<jobid>; the summary at the end.
# Exits non-zero when a run or a comparison fails.
# Cost: 4 nodes x <= 5 h on gg = <= 6.7 SU (EAR20006).
# Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
DECK=config/fd_cycle/bp5_fd_p2_1800yr.toml
HCACHE=${FDC_HCACHE:-$SCRATCH/seas-qd-green/bp5_fd_p2_h1e-10.bin}
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd "$ROOT" || exit 2

# ---- GitHub-first guard, and the binaries match this code ----
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
BUILT=$(cat miniapps/seas/BUILD_COMMIT_VISTA 2>/dev/null)
git cat-file -e "$BUILT^{commit}" 2>/dev/null \
   || { echo "ERROR: binaries built at '$BUILT', not a commit here; run vista_build.sbatch"; exit 2; }
CODE=$(git diff --name-only "$BUILT" HEAD | grep -v -E '^miniapps/seas/(config|jobs)/fd_cycle/')
[ -z "$CODE" ] || { echo "ERROR: code changed since the build at $BUILT:"; echo "$CODE"; exit 2; }
echo "commit $(git rev-parse HEAD), binaries built at $BUILT"

RUNS=$ROOT-runs/p2check_$SLURM_JOB_ID
mkdir -p "$RUNS" "$(dirname "$HCACHE")" || exit 2
module purge
for m in $MODULES; do module load $m || exit 2; done
cd miniapps/seas
P=bp5_fd_p2_1800yr

# The variants (prefix kept, so the files pair by name).
sed -e 's/^V_on               = 1.0e-2$/V_on               = 1.0e2/' "$DECK" > "$RUNS/qd_green.toml"
sed -e 's/^operators     = "matrix"$/operators     = "loops"/' \
    -e 's/^\[cycle.qd.green\]$/[cycle.qd.green]\nenabled          = false/' \
    "$RUNS/qd_green.toml" > "$RUNS/qd_loops.toml"
sed -e 's/^mode         = "\(off\|rate2\)"$/mode         = "rate2"/' "$DECK" > "$RUNS/lts.toml"
sed -e 's/^mode         = "\(off\|rate2\)"$/mode         = "off"/' "$DECK" > "$RUNS/gts.toml"
for f in qd_green qd_loops lts gts; do
   echo "--- $f.toml against $DECK"; diff "$DECK" "$RUNS/$f.toml"
done
grep -q '^V_on               = 1.0e2$' "$RUNS/qd_green.toml" && grep -q '^enabled          = false$' \
   "$RUNS/qd_loops.toml" && grep -q '^mode         = "rate2"$' "$RUNS/lts.toml" \
   && grep -q '^mode         = "off"$' "$RUNS/gts.toml" \
   || { echo "ERROR: a variant deck lacks its change (the sed patterns no longer match the deck)"; exit 2; }
GATE_LTS=0
grep -q '^mode         = "rate2"$' "$DECK" && GATE_LTS=1

# Every node's committed memory every 2 s, in the background (vista_memlog.sh);
# each run prints its peak per node.
srun --overlap -N "$SLURM_NNODES" -n "$SLURM_NNODES" --ntasks-per-node=1 \
   "$ROOT/miniapps/seas/jobs/fd_cycle/vista_memlog.sh" 2 > "$RUNS/memlog.txt" 2>&1 &
MEMLOG=$!
trap 'kill $MEMLOG 2>/dev/null' EXIT

fails=0
run() {   # name, ranks, deck, driver options...
   local n=$1 np=$2 deck=$3; shift 3
   mkdir -p "$RUNS/$n"
   local s=$(date +%s)
   # TACC's ibrun, as the QD thread's 576-rank p = 2 runs (jobs/vista_qd_speed):
   # with a bare mpirun -np 576 --map-by ppr:144:node --bind-to core (1089963)
   # every 576-rank run died in its first MUMPS factorization ("Error
   # allocating RECV buffers for matrix distribution"), while 144 ranks on one
   # node ran.  ibrun -n N -o 0: N ranks from the first node on, 144 per node.
   ibrun -n "$np" -o 0 ./seas_fd_cycle_driver \
      --deck "$deck" --output-dir "$RUNS/$n" --log-every 100 "$@" \
      > "$RUNS/$n/run.log" 2>&1
   local rc=$?
   local e=$(date +%s)
   echo "$n: exit $rc in $(( e - s )) s on $np ranks"
   awk -v s="$s" -v e="$e" '$2 >= s && $2 <= e + 2 && $3 == "limit" {
         if ($6 > m[$1]) m[$1] = $6; L = $4 }
      END { for (h in m) printf "    %s: peak committed %.1f of %.1f GiB\n", h, m[h], L }' \
      "$RUNS/memlog.txt" | sort
   grep -E "QD stages|\[QD-GREEN\]|LTS:|benchmark output|profile points|=== fd_cycle summary" \
      -A3 "$RUNS/$n/run.log" | head -40
   [ $rc -eq 0 ] || fails=$((fails + 1))
}
cmp() {   # label, gate (1: a failure fails the job), run A, run B, options
   local l=$1 g=$2 a=$3 b=$4; shift 4
   echo "--- $l: $a against $b"
   python3 tests/cycle/compare_interp.py "$RUNS/$a" $P "$RUNS/$b" $P --scale quantity "$@" \
      && return
   if [ "$g" = 1 ]; then fails=$((fails + 1)); else echo "    (for the record only: the deck runs global stepping)"; fi
}

H=(--qd-green-hcache "$HCACHE")
run green 576 "$RUNS/qd_green.toml" "${H[@]}" --max-steps 300
run loops 576 "$RUNS/qd_loops.toml" --max-steps 300
run loops144 144 "$RUNS/qd_loops.toml" --max-steps 300
cmp "Green's function" 1 green loops --tol 1e-5
cmp "576 against 144 ranks" 1 green loops144 --tol 1e-5
run lts 576 "$RUNS/lts.toml" "${H[@]}" --max-dyn-steps 800
run gts 576 "$RUNS/gts.toml" "${H[@]}" --max-dyn-steps 6200
cmp "LTS against global stepping, on the fault" $GATE_LTS lts gts --t-min 1e-9 \
    --files "{fltst,source}" --tol 1e-3
cmp "LTS against global stepping, off the fault" $GATE_LTS lts gts --t-min 1e-9 \
    --files blkst --tol 5e-2
for n in lts gts; do
   echo "--- $n"; grep -E "end: |steps of|ms/step|s/step" "$RUNS/$n/run.log" | tail -6
done
echo "p2 check: $fails failed (runs in $RUNS)"
[ $fails -eq 0 ]
