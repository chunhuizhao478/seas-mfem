#!/bin/bash
#SBATCH -J fdc_absorb_e1
#SBATCH -o fdc_absorb_e1_%j.out
#SBATCH -e fdc_absorb_e1_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 03:00:00
#SBATCH -A EAR20006
# Fully dynamic BP5 cycle on Vista, P5: event 1 at p = 1 with reflecting and
# with absorbing outer faces ([cycle.absorbing] faces = "off" and "outer"),
# each through its whole first dynamic phase ([cycle.switch]
# max_dynamic_phases = 1: the run stops when that phase ends).  Both decks
# are config/fd_cycle/bp5_fd_p2_1800yr.toml at p = 1 without the snapshots
# and the Green's function, written into the log folder (the sed lines
# below are the whole change; "outer" adds the [cycle.absorbing] section).
# 144 ranks each, one after the other, on binaries vista_build.sbatch made
# at this commit.  The summary prints each phase's end and the final slip
# at the on-fault stations; jobs/fd_cycle/plot_absorbing_ab.py draws the
# station pages from the two run folders (numpy and matplotlib; on the Mac
# or in the cloud).
# Logs and runs in <checkout>-runs/absorb_e1_<jobid>.  Exits non-zero when a
# run fails.  Cost: 1 node x <= 3 h on gg = <= 1 SU (EAR20006).
# Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
DECK=config/fd_cycle/bp5_fd_p2_1800yr.toml
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

RUNS=$ROOT-runs/absorb_e1_$SLURM_JOB_ID
mkdir -p "$RUNS" || exit 2
module purge
for m in $MODULES; do module load $m || exit 2; done
cd miniapps/seas

sed -e '/^# Whole-field snapshots/,$d' \
    -e 's/^order = 2$/order = 1/' \
    -e 's/^operators     = "matrix"$/operators     = "loops"/' \
    -e 's/^\[cycle.qd.green\]$/[cycle.qd.green]\nenabled          = false/' \
    -e 's/^max_dynamic_phases = 0$/max_dynamic_phases = 1/' \
    -e 's/^prefix       = "bp5_fd_p2_1800yr"$/prefix       = "ab"/' \
    "$DECK" > "$RUNS/off.toml"
{ cat "$RUNS/off.toml"; printf '\n[cycle.absorbing]\nfaces = "outer"\n'; } > "$RUNS/outer.toml"
for f in off outer; do
   echo "--- $f.toml against $DECK"; diff "$DECK" "$RUNS/$f.toml"
done
grep -q '^order = 1$' "$RUNS/off.toml" && grep -q '^operators     = "loops"$' "$RUNS/off.toml" \
   && grep -q '^enabled          = false$' "$RUNS/off.toml" \
   && grep -q '^max_dynamic_phases = 1$' "$RUNS/off.toml" \
   && grep -q '^prefix       = "ab"$' "$RUNS/off.toml" \
   && ! grep -q 'snapshots' "$RUNS/off.toml" \
   || { echo "ERROR: a deck lacks its change (the sed patterns no longer match the deck)"; exit 2; }

fails=0
for n in off outer; do
   mkdir -p "$RUNS/$n"
   s=$(date +%s)
   ibrun -n 144 -o 0 ./seas_fd_cycle_driver --deck "$RUNS/$n.toml" --output-dir "$RUNS/$n" \
      --log-every 1000 > "$RUNS/$n/run.log" 2>&1
   rc=$?
   echo "$n: exit $rc in $(( $(date +%s) - s )) s on 144 ranks"
   grep -E "absorbing boundary:|time per dynamic step|end: " "$RUNS/$n/run.log"
   [ $rc -eq 0 ] || fails=$((fails + 1))
done

echo "--- final slip at the on-fault stations (m): off, outer, outer / off - 1"
for f in "$RUNS"/off/ab_fltst_*.txt; do
   g="$RUNS/outer/$(basename "$f")"
   [ -f "$g" ] || { echo "  $(basename "$f"): no outer file"; continue; }
   a=$(awk '!/^#/ && NF { s = $2 } END { print s }' "$f")
   b=$(awk '!/^#/ && NF { s = $2 } END { print s }' "$g")
   awk -v n="$(basename "$f" .txt)" -v a="$a" -v b="$b" \
      'BEGIN { printf "  %-26s %12.6f %12.6f %+11.3e\n", n, a, b, (a != 0 ? b / a - 1 : 0) }'
done
echo "plots: python3 jobs/fd_cycle/plot_absorbing_ab.py $RUNS/off $RUNS/outer ab ab.pdf \"BP5-FD event 1, p = 1\""
echo "failing runs: $fails"
[ $fails -eq 0 ] || exit 1
