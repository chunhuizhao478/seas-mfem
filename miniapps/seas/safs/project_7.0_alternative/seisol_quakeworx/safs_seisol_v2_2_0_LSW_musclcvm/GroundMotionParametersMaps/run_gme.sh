#!/usr/bin/env bash
#
# Compute ground-motion parameter maps (PGA, PGV, PGD, SA(T)) from SeisSol
# free-surface output using SeisSol's official postprocessing tool
# (ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py).
#
# Output: writes  safs-GME-surface.xdmf  (+ binaries) NEXT TO each input
#         surface.xdmf. Open in ParaView. Cell fields: PGA, PGV, PGD, SA<T>s.
#         Values are GMRotD50 (median rotated geometric-mean of the 2 horizontals).
#
# CAVEAT: the surface output here is sampled at dt = 0.5 s (Nyquist = 1 Hz).
#         => PGA and SA at periods < ~2 s are aliased/under-resolved.
#            PGV and PGD are far less affected.
#            For trustworthy PGA, re-run SeisSol with a finer free-surface
#            output interval, or use the receiver .dat files (dt = 0.005 s).
#
# Usage:
#   conda activate pythonenv
#   ./run_gme.sh                 # setup (if needed) + run on both result dirs
#   ./run_gme.sh /path/to/another_output_dir [...]   # run on custom dir(s)
#
set -o pipefail   # NOTE: not `set -u` (conda activate trips on unbound vars)

# macOS: allow numpy/Accelerate to run in forked worker processes
# (the GME tool uses multiprocessing with the 'fork' start method).
export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES

# ---- locations -------------------------------------------------------------
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$HERE/ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py"
SMTK_DIR="$HERE/gmpe-smtk"           # the tool appends <script_dir>/gmpe-smtk to sys.path
SMTK_COMMIT="4f008173e89f6e4ba4450fb43e95ffdf51a7c2ba^"   # deprecated, openquake-free

# Default result directories (override by passing dirs as arguments)
DEFAULT_DIRS=(
  "/Users/chunhuizhao/Downloads/output_safs_v2.2.0_cvm_mat"
  "/Users/chunhuizhao/Downloads/output_safs_v2.2.0_constant_mat_onfaultpoints"
)

# Runtime knobs
MP="${MP:-4}"                        # worker processes (must be <= core count)
PERIODS="${PERIODS:-}"               # e.g. PERIODS="2 5" to limit SA periods
EXTRA="${EXTRA:-}"                   # e.g. EXTRA="--CAV" or "--lowpass 1.0"

# ---- one-time setup --------------------------------------------------------
echo "==> python dependencies (installed into the SAME interpreter that runs the tool)"
PY="${PY:-python3}"
"$PY" -m pip install seissolxdmf seissolxdmfwriter scipy h5py lxml numpy || {
  echo "pip install failed"; exit 1; }

if [ ! -d "$SMTK_DIR" ]; then
  echo "==> cloning gmpe-smtk (mandatory dependency)"
  git clone https://github.com/GEMScienceTools/gmpe-smtk "$SMTK_DIR" || {
    echo "git clone failed"; exit 1; }
  ( cd "$SMTK_DIR" && git checkout "$SMTK_COMMIT" ) || {
    echo "git checkout failed"; exit 1; }
else
  echo "==> gmpe-smtk already present ($SMTK_DIR)"
fi

# ---- runs ------------------------------------------------------------------
DIRS=("$@")
if [ "${#DIRS[@]}" -eq 0 ]; then
  DIRS=("${DEFAULT_DIRS[@]}")
fi

PERIOD_ARG=""
[ -n "$PERIODS" ] && PERIOD_ARG="--periods $PERIODS"

for d in "${DIRS[@]}"; do
  xdmf="$d/safs-surface.xdmf"
  if [ ! -f "$xdmf" ]; then
    echo "!! skip: $xdmf not found"
    continue
  fi
  echo ""
  echo "============================================================"
  echo "==> $xdmf"
  echo "============================================================"
  # The SeisSol tool writes its output (fixed name safs-GME-surface.{xdmf,h5})
  # to the CURRENT directory, ignoring the input file's location. So cd into the
  # target output dir (subshell) before running -> result lands next to its input
  # and runs on different dirs don't overwrite each other.
  ( cd "$d" && "$PY" "$SCRIPT" --noMPI --MP "$MP" $PERIOD_ARG $EXTRA safs-surface.xdmf ) || {
    echo "!! GME computation failed for $d"; continue; }
  echo "--> wrote $d/safs-GME-surface.xdmf  (+ .h5)"
done

echo ""
echo "Done. Open the *-GME-surface.xdmf files in ParaView."
