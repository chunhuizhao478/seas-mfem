#!/bin/bash
# =============================================================================
# build_seissol_o2_o3_expanse.sbatch
#
# TPV104 200 m GTS ORDER BASELINE -- SeisSol build leg.
#
# SeisSol's ConvergenceOrder is COMPILE-TIME
# (SeisSol/src/Initializer/TimeStepping/GlobalTimestep.cpp:46 uses
#  `seissol::ConvergenceOrder`), so each order needs its own binary.  The
# baseline pairs MFEM p_k with SeisSol o_(k+1) (R2, accuracy-matched):
#
#     MFEM p1 (ader_order 2)  <->  SeisSol ORDER=2   <-- MISSING, built here
#     MFEM p2 (ader_order 3)  <->  SeisSol ORDER=3   <-- MISSING, built here
#     MFEM p3 (ader_order 4)  <->  SeisSol ORDER=4   <-- ALREADY EXISTS, NOT rebuilt
#
# CLONED FROM (read in full before authoring):
#   $SEISSOL_ROOT/build_o6_pgv.sbatch          -- the most recent working per-order
#                                                 build job (job 52490896, exit 0)
#   $SEISSOL_ROOT/build_seissol_expanse.sbatch -- the older generic build job
#   $SEISSOL_ROOT/build_seissol_expanse.sh     -- the shared, ORDER-parametric recipe
#
# WHY ONE JOB FOR BOTH ORDERS:
#   The recipe is ORDER-PARAMETRIC ONLY.  Verified by diffing the non-INTERNAL
#   CMakeCache entries of build-rome-o4-elastic vs build-rome-o5-elastic: the
#   entire diff is `ORDER:STRING=4` vs `ORDER:STRING=5`.  Confirmed on-machine
#   2026-07-24, o4 cache reads:
#     CMAKE_BUILD_TYPE=Release  HOST_ARCH=rome  EQUATIONS=elastic
#     PRECISION=double  NUMBER_OF_MECHANISMS=0  ASAGI=ON  NETCDF=ON
#     GEMM_TOOLS_LIST=LIBXSMM,PSpaMM  GRAPH_PARTITIONING_LIBS=parmetis
#     NUMA_AWARE_PINNING=ON
#   The two builds MUST be sequential in ONE job, not parallel:
#     * both invocations run `git submodule update --init --recursive` on the
#       SAME checkout (build_seissol_expanse.sh sec.12) -- concurrent runs race;
#     * both pip-install into the SAME conda env (sec.2) -- concurrent runs race;
#     * one job guarantees o2 and o3 see an IDENTICAL toolchain state.
#
# NO `module load` HERE ON PURPOSE:
#   build_seissol_expanse.sh sec.1 does `module purge` then loads the full stack
#   itself (slurm + cpu/0.17.3b + gcc/10.2.0 + cmake/3.21.4 + openmpi/4.1.3,
#   fallback 4.1.1) and sets CC/CXX/FC/CFLAGS/CXXFLAGS.  Anything loaded in this
#   sbatch would be purged on the recipe's first line.  build_o6_pgv.sbatch loads
#   nothing either.  `#SBATCH --export=ALL` is what carries Lmod's `module`
#   shell function into the batch shell -- do not drop it.
#
# COST / WALLTIME:
#   All nine dependencies (hdf5, netcdf, eigen, yamlcpp, lua, asagi, easi,
#   libxsmm, parmetis) are already built and stamp-guarded in
#   $SEISSOL_ROOT/deps/.stamps (all present, dated 2026-07-13 20:04-20:15), so
#   this job compiles SeisSol ONLY.  Measured fresh single-order SeisSol build
#   (o5, NPROC=16): ~8 min (cmake configure 20:15:19 -> last object 20:23:11).
#   Do NOT schedule off build_o6_pgv's 1:57 elapsed -- that job was INCREMENTAL
#   (its GeneratedCode/*.cpp kept their 2026-07-13 mtimes).  Lower orders emit
#   SMALLER yateto kernels than o5/o6, so o2+o3 should land well inside 30 min.
#   2 h is pure headroom.  `shared` bills max(cpu-term, mem-term); with --mem=32G
#   the cpu term wins, so the charge is ~16 core-hours per wall-hour: ~8 SU for a
#   30-minute build, ~32 SU if it ran the full 2 h.  (An earlier draft requested
#   64G, which would have billed as 32 cores -- double, for nothing.)
#   Disk: 110 MB per build dir (measured on build-rome-o4-elastic); 23 T free.
#
# SUBMIT (from anywhere -- the script cd's to SEISSOL_ROOT itself, so the .out
# lands beside the other build logs only if you submit from there; recommended):
#   scp <this file> czhao1@login.expanse.sdsc.edu:/expanse/projects/qstore/usc143/qwxdev/apps/expanse/rocky8.8/seisol_build/
#   ssh czhao1@login.expanse.sdsc.edu
#   cd /expanse/projects/qstore/usc143/qwxdev/apps/expanse/rocky8.8/seisol_build
#   sbatch build_seissol_o2_o3_expanse.sbatch
# WATCH:
#   squeue -u $USER ; tail -f build-o2o3-elastic.<jobid>.out
#
# EXPECTED PRODUCTS
#   $SEISSOL_ROOT/SeisSol/build-rome-o2-elastic/seissol-elastic-o2-f64
#   $SEISSOL_ROOT/SeisSol/build-rome-o3-elastic/seissol-elastic-o3-f64
#   (+ a sibling proxyseissol-elastic-oN-f64 next to each)
#   Naming is deterministic: SeisSol/CMakeLists.txt NEW_BINARY_NAMING builds
#   NAME_SUFFIX = "${EQUATIONS}-o${ORDER}-f64"; the build dir name comes from
#   build_seissol_expanse.sh: "$SEISSOL_SRC/build-${HOST_ARCH}-o${ORDER}-${EQUATIONS}".
#
# EXIT CODES
#   0  both binaries present, executable, freshly built, caches order-clean
#   2  cannot cd to SEISSOL_ROOT / recipe or source tree missing
#   3  source-commit fairness gate: the working tree HEAD moved off the pin, OR a
#      resulting BINARY (o2, o3 or the pre-existing o4) embeds a different commit
#   4  a required executable is missing or not executable
#   5  a produced binary is STALE (pre-dates this job) -> the build silently failed
#   6  codegen python interpreter missing
#   7  a build dir's CMakeCache ORDER != the order requested
#   8  a build dir's CMakeCache differs from o4's in something OTHER than ORDER
#   9  codegen toolchain (numpy / PSpaMM) is not at the pinned version, or drifted
#      during this job
#   NOTE: exit codes are FIRST-writer-wins, so the most fairness-relevant cause
#         survives to the exit status even when several checks fire.
# =============================================================================
#SBATCH --job-name=build-o2o3-elastic
#SBATCH --account=lbl107          # VALID (verified 2026-07-24: 199000 alloc / 194262 spent,
                                  # ~4.7k SU left, expires 2026-08-27). This build costs ~7 SU.
                                  # NOTE: -A EAR20006 and -p development are FRONTERA names;
                                  # they are INVALID on Expanse and produce a QOS failure.
#SBATCH --partition=shared        # billed only for the cores requested; `compute` would bill
                                  # all 128 cores/node for a 16-core compile.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16        # == NPROC below; yateto codegen is memory-heavy, keep -j modest
#SBATCH --mem=32G                 # BILLING: `shared` charges max(cpu-term, mem-term) with
                                  # TRESBillingWeights cpu=3600, mem=1800/G — so 64G would
                                  # have billed as 32 cores (64*1800 > 16*3600), DOUBLING the
                                  # cost for no benefit.  32G bills at parity with the 16 cores
                                  # actually used and is still ~3.8x the only measured peak
                                  # (8.4 GB MaxRSS, job 52490896 — though that was an
                                  # INCREMENTAL build, so treat it as a soft floor).
#SBATCH -t 02:00:00               # generous; measured need is <30 min for both orders
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.out
#SBATCH --export=ALL              # REQUIRED: carries Lmod's `module` function into the job

# pipefail only -- NOT `set -u`: the recipe sources conda's GFORTRAN hook, which
# references unbound vars (documented at the top of build_seissol_expanse.sh).
set -o pipefail

JOB_START_EPOCH=$(date +%s)

# Exit status is FIRST-writer-wins, not last.  Every failure message is printed
# either way, but the exit CODE must retain the most fairness-relevant cause: a
# plain `rc_total=N` chain let a late, milder condition (e.g. 4 = missing binary)
# silently overwrite 9 = codegen-toolchain drift, which the documented exit-code
# table tells the operator means "rebuild o4, do NOT proceed".
rc_total=0
fail () { echo "$2"; [ "$rc_total" -eq 0 ] && rc_total="$1"; return 0; }

# --- ROOT + SUB-PATH PINS (do NOT drop any of these exports) ------------------
# build_seissol_expanse.sh sec.0 defaults:
#     SEISSOL_ROOT=${SEISSOL_ROOT:-$HOME/seissol}
#     PREFIX=${SEISSOL_PREFIX:-$SEISSOL_ROOT/deps}
#     SEISSOL_SRC=${SEISSOL_SRC:-$SEISSOL_ROOT/SeisSol}
#     CONDA_DIR=${CONDA_DIR:-$SEISSOL_ROOT/conda}
# Without the SEISSOL_ROOT export it would bootstrap a WHOLE SECOND Miniforge +
# dependency tree in $HOME (hours of wall time, and a binary linked against a
# DIFFERENT dependency stack than the existing o4 -- silently voiding R2).
# The other three are pinned explicitly because `--export=ALL` propagates the
# submitting shell's environment: a stale SEISSOL_SRC / SEISSOL_PREFIX /
# CONDA_DIR left over from another project would be honoured by the recipe.
export SEISSOL_ROOT=/expanse/projects/qstore/usc143/qwxdev/apps/expanse/rocky8.8/seisol_build
export SEISSOL_SRC="$SEISSOL_ROOT/SeisSol"
export SEISSOL_PREFIX="$SEISSOL_ROOT/deps"
export CONDA_DIR="$SEISSOL_ROOT/conda"
RECIPE="$SEISSOL_ROOT/build_seissol_expanse.sh"

cd "$SEISSOL_ROOT" || { echo "FATAL: cannot cd to SEISSOL_ROOT=$SEISSOL_ROOT"; exit 2; }
[ -f "$RECIPE" ]            || { echo "FATAL: no recipe at $RECIPE"; exit 2; }
[ -d "$SEISSOL_SRC/.git" ]  || { echo "FATAL: $SEISSOL_SRC is not a git checkout"; exit 2; }
[ -d "$SEISSOL_PREFIX/.stamps" ] || { echo "FATAL: no dep stamp dir at $SEISSOL_PREFIX/.stamps"; exit 2; }

echo "############ SeisSol o2+o3 elastic build for the TPV104 GTS order baseline"
echo "=== job  : ${SLURM_JOB_ID:-<none>}  node=$(hostname)  cpus=${SLURM_CPUS_PER_TASK:-?}"
echo "=== root : $SEISSOL_ROOT"
echo "=== src  : $SEISSOL_SRC"
echo "=== start: $(date)"
df -h "$SEISSOL_ROOT" | tail -1

# --- SOURCE PIN (fairness gate, R2) ------------------------------------------
# The existing o4 binary that serves the p3 leg was compiled 2026-07-14 10:31
# from commit dc6db65133a999eb38c578c14cc443df1000dfd0
# (v1.3.1-2135-gdc6db6513, branch feature/peak-ground-velocity); `strings` on
# that binary contains both the SHA and the version string.
# o2/o3 MUST come from the SAME commit, or the "order study" is a source-version
# study.  Do NOT `git pull` / `git checkout` in this tree.
EXPECT_COMMIT=dc6db65133a999eb38c578c14cc443df1000dfd0
HAVE_COMMIT=$(git -C "$SEISSOL_SRC" rev-parse HEAD 2>/dev/null)
HAVE_DESCRIBE=$(git -C "$SEISSOL_SRC" describe --tags 2>/dev/null)
echo "=== source commit: $HAVE_COMMIT  ($HAVE_DESCRIBE)"
if [ "$HAVE_COMMIT" != "$EXPECT_COMMIT" ]; then
  echo "FATAL: SeisSol source HEAD moved."
  echo "       expected $EXPECT_COMMIT  (the commit the existing o4 binary was built from)"
  echo "       found    $HAVE_COMMIT"
  echo "       Either restore the pin, or rebuild ALL THREE orders (o2, o3 AND o4"
  echo "       into a fresh build-rome-o4-elastic) in one job so they match."
  exit 3
fi
echo "=== dep stamps present:"; ls -1 "$SEISSOL_PREFIX/.stamps" | tr '\n' ' '; echo

# --- CODEGEN-TOOLCHAIN PIN (fairness gate) -----------------------------------
# build_seissol_expanse.sh sec.2 unconditionally runs, with NO stamp guard:
#     "$PY" -m pip install --quiet --upgrade pip setuptools numpy
#     "$PY" -m pip install --quiet "git+https://github.com/SeisSol/PSpaMM.git"   (UNPINNED master)
# PSpaMM emits the GEMM microkernels.  If it drifts between the 2026-07-14 o4
# build and this job, p1/p2 would run DIFFERENT microkernels from p3 -- a silent
# speed bias in exactly the quantity this baseline measures.
# Empirically a no-op as of 2026-07-24 (numpy 2.5.1 / PSpaMM 0.3.1, both
# dist-info dirs still carrying their original 2026-07-13 19:46 mtimes, and the
# 2026-07-24 o6 rebuild left them untouched).  Recorded before, asserted after.
PYBIN="$CONDA_DIR/bin/python"
[ -x "$PYBIN" ] || { echo "FATAL: no codegen python at $PYBIN"; exit 6; }
ver_of(){ "$PYBIN" -m pip show "$1" 2>/dev/null | awk '/^Version:/{print $2}'; }
PRE_NUMPY=$(ver_of numpy); PRE_PSPAMM=$(ver_of pspamm)
echo "=== codegen toolchain BEFORE: numpy=$PRE_NUMPY pspamm=$PRE_PSPAMM   (expect 2.5.1 / 0.3.1)"
# ABSOLUTE pin, not just a before/after delta.  A delta check only catches drift
# that happens DURING this job; if some intervening job already upgraded PSpaMM
# after the 2026-07-14 o4 build, PRE == POST here and the check would be green
# while o2/o3 carry DIFFERENT GEMM microkernels from o4 -- a silent speed bias in
# exactly the quantity this baseline measures, in the direction hardest to detect
# after the fact.
EXPECT_NUMPY=2.5.1
EXPECT_PSPAMM=0.3.1
if [ "$PRE_NUMPY" != "$EXPECT_NUMPY" ] || [ "$PRE_PSPAMM" != "$EXPECT_PSPAMM" ]; then
  echo "!!! FAIRNESS ALERT: the codegen toolchain ALREADY differs from the state the"
  echo "!!!   pre-existing o4 binary was generated with (2026-07-14)."
  echo "!!!   numpy  found $PRE_NUMPY  expected $EXPECT_NUMPY"
  echo "!!!   pspamm found $PRE_PSPAMM expected $EXPECT_PSPAMM"
  echo "!!! o2/o3 built now would carry different GEMM microkernels from o4."
  echo "!!! Correct response: rebuild ALL THREE orders (including a fresh o4) in one"
  echo "!!! job so they share a toolchain state.  Do NOT proceed with a mixed set."
  [ $rc_total -eq 0 ] && rc_total=9
fi

# --- BUILDS ------------------------------------------------------------------
# HOST_ARCH=rome and PRECISION=double are hard/defaulted inside the recipe.
# EQUATIONS / NUMBER_OF_MECHANISMS / ORDER / NPROC are passed EXPLICITLY: the
# recipe's own ORDER default is 5 (sec.0), so an omitted ORDER would silently
# rebuild o5 instead of the order we need, and the failure would only surface at
# run time as a wrong `Minimum timestep:`.
# (rc_total is initialised at the top of the job -- NOT reset here, or the
#  codegen-toolchain alarm above would be discarded.)
BUILT=""            # orders this job actually compiled (see the staleness gate)
for ORD in 2 3; do
  BUILD_DIR="$SEISSOL_SRC/build-rome-o${ORD}-elastic"
  EXE="$BUILD_DIR/seissol-elastic-o${ORD}-f64"

  echo "================= ORDER=$ORD  start=$(date) ================="
  if [ -x "$EXE" ]; then
    # Already complete (e.g. this job is a resubmit after the other order failed).
    # NOT recorded in $BUILT, so the mtime staleness gate below skips it — a
    # binary produced by an EARLIER job necessarily pre-dates JOB_START_EPOCH,
    # and flagging that as "the build did not run/complete" would make the one
    # documented recovery workflow (resubmit after a single-order failure) always
    # exit 5 and invite an rm -rf of a perfectly good build dir.  Skipped orders
    # are instead validated by the embedded-commit assertion, which is the check
    # that actually matters for R2.
    echo ">>> SKIP: $EXE already exists and is executable (built by an earlier job)."
    echo ">>>       Staleness check waived for this order; commit parity still asserted."
  else
    if [ -d "$BUILD_DIR" ]; then
      # A build dir with no executable == a truncated / failed previous attempt.
      # The dep stamps protect the dependencies, but NOTHING protects a partial
      # SeisSol build dir, and cmake would happily reuse a half-written cache.
      echo ">>> WARNING: $BUILD_DIR exists but has no executable -> treating as a"
      echo ">>>          truncated previous build and removing it for a clean configure."
      rm -rf "$BUILD_DIR" || { echo "FATAL: cannot remove $BUILD_DIR"; exit 2; }
    fi
    t0=$SECONDS
    # PRECISION is pinned EXPLICITLY.  build_seissol_expanse.sh hard-codes
    # HOST_ARCH=rome but reads PRECISION=${PRECISION:-double} from the
    # ENVIRONMENT and passes it straight through as -DPRECISION=$PRECISION.
    # With #SBATCH --export=ALL a stale PRECISION=single in the submitting shell
    # would silently build o2/o3 at single precision against a double-precision
    # o4 — a precision study mislabelled as an order study, caught only after
    # ~16 wasted minutes by the CMakeCache diff below.
    EQUATIONS=elastic NUMBER_OF_MECHANISMS=0 PRECISION=double ORDER=$ORD NPROC=16 \
      bash "$RECIPE"
    rc=$?
    BUILT="$BUILT $ORD"     # this job actually compiled this order
    echo ">>> recipe exit=$rc  elapsed=$((SECONDS - t0))s"
    # NOTE: build_seissol_expanse.sh ENDS IN `echo`, so its exit status is 0 even
    # when `cmake --build` failed.  rc is logged for the record but is NOT the
    # success criterion -- the artifact checks below are.
  fi
  echo "================= ORDER=$ORD  end=$(date) ================="
done

# --- VERIFY: codegen toolchain did not drift ---------------------------------
POST_NUMPY=$(ver_of numpy); POST_PSPAMM=$(ver_of pspamm)
echo "=== codegen toolchain AFTER : numpy=$POST_NUMPY pspamm=$POST_PSPAMM"
if [ "$PRE_PSPAMM" != "$POST_PSPAMM" ] || [ "$PRE_NUMPY" != "$POST_NUMPY" ]; then
  echo "!!! FAIRNESS ALERT: the codegen toolchain changed during this job."
  echo "!!!   numpy  $PRE_NUMPY -> $POST_NUMPY"
  echo "!!!   pspamm $PRE_PSPAMM -> $POST_PSPAMM"
  echo "!!! o2/o3 GEMM microkernels may now differ from the o4 kernels built 2026-07-14."
  echo "!!! Correct response: rebuild o4 from scratch in the same state, do NOT proceed."
  [ $rc_total -eq 0 ] && rc_total=9
fi

# --- VERIFY: CMakeCache differs from o4 ONLY in ORDER (fairness artifact) -----
# Normalise the -oN- path substrings and the ORDER value, drop INTERNAL entries
# and comments, then diff.  Any surviving difference means the order study has
# become a build study.
REF_CACHE="$SEISSOL_SRC/build-rome-o4-elastic/CMakeCache.txt"
norm_cache(){ grep -v ':INTERNAL=' "$1" | grep -vE '^(//|#|$)' \
              | sed -E 's/build-rome-o[0-9]+-/build-rome-oN-/g; s/^ORDER:STRING=[0-9]+/ORDER:STRING=N/' \
              | sort; }
if [ -f "$REF_CACHE" ]; then
  for ORD in 2 3; do
    C="$SEISSOL_SRC/build-rome-o${ORD}-elastic/CMakeCache.txt"
    [ -f "$C" ] || continue
    GOT_ORDER=$(awk -F= '/^ORDER:STRING=/{print $2}' "$C")
    if [ "$GOT_ORDER" != "$ORD" ]; then
      echo "!!! FATAL: build-rome-o${ORD}-elastic/CMakeCache.txt has ORDER=$GOT_ORDER (expected $ORD)"
      [ $rc_total -eq 0 ] && rc_total=7
    fi
    D=$(diff <(norm_cache "$REF_CACHE") <(norm_cache "$C"))
    if [ -n "$D" ]; then
      echo "!!! FAIRNESS ALERT: o${ORD} CMakeCache differs from o4 in more than ORDER:"
      echo "$D" | head -40
      [ $rc_total -eq 0 ] && rc_total=8
    else
      echo "=== CMakeCache o${ORD} vs o4: ORDER-only difference  (OK)"
    fi
  done
else
  echo "WARN: reference cache $REF_CACHE missing -- CMakeCache fairness diff SKIPPED"
fi

# --- VERIFY: the executables, and print their paths --------------------------
# Checks per binary: exists, is executable, and (for the two we built) has an
# mtime at or after this job's start -- a stale o2/o3 binary would mean the
# recipe short-circuited and we would be about to time the WRONG build.
echo "=== resulting executables ==="
for ORD in 2 3 4; do
  B="$SEISSOL_SRC/build-rome-o${ORD}-elastic/seissol-elastic-o${ORD}-f64"
  if [ -x "$B" ]; then
    printf "  OK   %s\n" "$B"
    printf "       size   : %s bytes\n" "$(stat -c %s "$B" 2>/dev/null)"
    printf "       mtime  : %s\n"       "$(date -d @"$(stat -c %Y "$B")" 2>/dev/null)"
    EMB_VER="$(strings "$B" | grep -m1 -E '^v1\.3\.1-')"
    EMB_SHA="$(strings "$B" | grep -m1 -E '^[0-9a-f]{40}$')"
    printf "       version: %s\n"       "${EMB_VER:-<none>}"
    printf "       commit : %s\n"       "${EMB_SHA:-<none>}"
    # R2 SOURCE-PROVENANCE GATE, asserted on the ARTIFACT not on the working tree.
    # A clean `git rev-parse HEAD` proves nothing about a binary compiled weeks
    # ago or by another job -- and o4, which serves the entire p3 leg (the
    # headline of the baseline), is not built here at all.  If o4 came from a
    # different revision the "order study" is silently a source-version study
    # with every other gate green.  Checked for ALL THREE orders.
    if [ "$EMB_SHA" != "$EXPECT_COMMIT" ]; then
      echo "  !!!  COMMIT MISMATCH: o${ORD} binary embeds '${EMB_SHA:-<none>}'"
      echo "  !!!    expected $EXPECT_COMMIT"
      echo "  !!!    The three order legs would compare SOURCE REVISIONS, not orders."
      echo "  !!!    Correct response: rebuild ALL THREE orders (o2, o3 AND a fresh o4)"
      echo "  !!!    from one commit in one job."
      [ $rc_total -eq 0 ] && rc_total=3
    fi
    # Staleness applies ONLY to orders this job actually compiled -- see the SKIP
    # branch above.  o4 is never built here; a skipped o2/o3 is validated by the
    # commit assertion instead.
    case " $BUILT " in
      *" $ORD "*)
        MT=$(stat -c %Y "$B")
        if [ "$MT" -lt "$JOB_START_EPOCH" ]; then
          echo "  !!!  STALE: o${ORD} was compiled by this job but its binary pre-dates the"
          echo "  !!!         job start -- the recipe short-circuited and did not relink."
          [ $rc_total -eq 0 ] && rc_total=5
        fi ;;
      *) echo "       (not built by this job -- staleness check waived, commit asserted above)" ;;
    esac
  elif [ -e "$B" ]; then
    printf "  NOEXEC %s  (exists but is not executable)\n" "$B"
    [ $rc_total -eq 0 ] && rc_total=4
  else
    printf "  MISS %s\n" "$B"
    [ $rc_total -eq 0 ] && rc_total=4
  fi
done

# --- The three paths the run legs consume ------------------------------------
# The run sbatch derives these itself from $SEISSOL_ROOT and the leg order
# (run_seissol_gts_orders.sbatch, exe_for()); it does NOT read $SEISSOL_EXE.
# Leg selection is via BENCH_SS_LEGS.
echo
echo "=== SeisSol executables for the TPV104 GTS order baseline ==="
echo "  p1 leg (MFEM order=1, ader_order=2)  ->  $SEISSOL_SRC/build-rome-o2-elastic/seissol-elastic-o2-f64"
echo "  p2 leg (MFEM order=2, ader_order=3)  ->  $SEISSOL_SRC/build-rome-o3-elastic/seissol-elastic-o3-f64"
echo "  p3 leg (MFEM order=3, ader_order=4)  ->  $SEISSOL_SRC/build-rome-o4-elastic/seissol-elastic-o4-f64   [pre-existing]"
echo
echo "  Run the legs with (nothing is submitted by this job) — from the staged"
echo "  Lustre run directory, NOT from here:"
echo "    cd /expanse/lustre/projects/lbl107/czhao1/safs/tpv104_gts_order"
echo "    sbatch run_seissol_gts_orders.sbatch                              # all three, o2->o3->o4"
echo "    sbatch --export=ALL,BENCH_SS_LEGS=\"2\" run_seissol_gts_orders.sbatch   # a single leg"
echo "  Each leg picks its own deck (parameters_gts_o<N>.par) and its own binary;"
echo "  there is no SEISSOL_EXE override and no parameters.par rename step."
echo
echo "  PRE-REGISTERED dt gate (GlobalTimestep.cpp:45-46, dt ~ 1/(2*ORDER-1), CFL 0.5,"
echo "  2,464,689-tet tpv104 mesh; anchored on the verified o4 value 3.66128e-4 s)."
echo "  Each leg's first 'Minimum timestep:' line MUST read (to a RELATIVE tolerance"
echo "  of 1e-5 — the o2/o3 values are DERIVED from a printed o4 value that carries"
echo "  only 4 decimals, so the last digit is uncertain by +-1):"
echo "    o2 -> 854.299 us      o3 -> 512.579 us      o4 -> 366.128 us"
echo "  Any deviation beyond that means the wrong binary or a drifted mesh/CFL -> VOID."

echo
echo "=== finished: $(date)  total_elapsed=$((SECONDS))s  exit=$rc_total"
exit $rc_total
