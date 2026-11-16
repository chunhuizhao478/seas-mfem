#!/bin/bash
#SBATCH -J cycmovie
#SBATCH -o cycmovie_%j.out
#SBATCH -e cycmovie_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:00:00
# The standard videos and the magnitude catalog of a finished (or running)
# fully dynamic cycle run, on one Lonestar6 node (2 x 64 EPYC 7763 cores)
# (Chunhui Zhao, 2026-11-16: "we need to find a way to speed this up";
# "let's develop this, later this could be a standardlized video generation
# process, and I want the branch to save this as well").
# postprocess/cycle_movies/README.md describes the scripts this job runs.
#
# Submit from this folder of a checkout on Lonestar6, after `git fetch` and
# a checkout of the pushed commit on the login node:
#   MOVIE_RUN=/scratch/.../<run folder> MOVIE_VENV=/scratch/.../venv-post-ls6 \
#      sbatch -A <allocation> ls6_cycle_movies.sbatch
# Environment:
#   MOVIE_RUN (required)   run folder of seas_fd_cycle_driver, absolute.  It is
#                          only read (files opened read-only, HDF5 locking off).
#   MOVIE_VENV (required)  Python environment with numpy, h5py, matplotlib and
#                          imageio-ffmpeg (ffmpeg with libx264): its bin/python3
#                          runs the scripts.
#   MOVIE_OUT              output folder, absolute, under $SCRATCH (results
#                          never go to $HOME or $WORK; default
#                          $SCRATCH/cycle_movies/<run folder name>_<job id>).
#                          Must not exist: a job never writes over videos.
#   MOVIE_WORKERS          render processes (default 120 of the node's 128
#                          cores; each also drives one 1-thread ffmpeg).
#   MOVIE_MAP_DATA         folder with the Natural Earth files
#                          (postprocess/cycle_movies/fetch_natural_earth.sh,
#                          run once on a login node).  When set, the videos
#                          are also written on the line map (*_map.mp4).
#   MOVIE_REGION           region file of the line map (default
#                          postprocess/cycle_movies/regions/socal_ucerf2d.json).
#   MOVIE_ARGS             more options of cycle_movies.py for every video,
#                          e.g. "--max-mb 15 --last-phase 40".
#   MAG_MU, MAG_WIDTH_M    shear modulus at the fault (Pa) and rupture width
#                          (m).  When both are set, magnitude_catalog.py also
#                          runs; MAG_FAULT_NAMES (comma-separated) names the
#                          faults.
# The job writes RENDER_INFO.txt (commit, run folder, options, node) beside
# the videos and DONE at exit 0.
# No account line: sbatch -A <allocation>.  Cost: 1 node x <= 1 h on
# development = <= 1 SU.
# Submit only on Chunhui Zhao's typed go for these jobs.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE MPLBACKEND=Agg
RUN=${MOVIE_RUN:?set MOVIE_RUN to the run folder (absolute)}
VENV=${MOVIE_VENV:?set MOVIE_VENV to the Python environment}
WORKERS=${MOVIE_WORKERS:-120}
read -r -a EXTRA <<< "${MOVIE_ARGS:-}"
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd "$ROOT" || exit 2

# ---- GitHub-first guard (CLAUDE.md) ----
# No fetch: the compute nodes may have no internet (jobs/2d_inplane/
# ls6_build.sbatch); fetch on the login node before sbatch.
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
HEAD=$(git rev-parse HEAD)
echo "commit $HEAD"

TOOL=$ROOT/miniapps/seas/postprocess/cycle_movies
PY=$VENV/bin/python3
[ -x "$PY" ] || { echo "ERROR: no $PY"; exit 2; }
case "$RUN" in
   /*) ;;
   *) echo "ERROR: MOVIE_RUN '$RUN' is not an absolute path"; exit 2 ;;
esac
[ -d "$RUN" ] || { echo "ERROR: no run folder $RUN"; exit 2; }
[ -n "$SCRATCH" ] || { echo "ERROR: SCRATCH is not set"; exit 2; }
OUT=${MOVIE_OUT:-$SCRATCH/cycle_movies/$(basename "$RUN")_$SLURM_JOB_ID}
case "$OUT/" in
   "$SCRATCH"/*) ;;
   *) echo "ERROR: MOVIE_OUT '$OUT' is not under \$SCRATCH ($SCRATCH)"; exit 2 ;;
esac
case "$OUT/" in
   "$RUN"/*) echo "ERROR: MOVIE_OUT '$OUT' is inside the run folder"; exit 2 ;;
esac
[ ! -e "$OUT" ] || { echo "ERROR: $OUT exists"; exit 2; }
MAP=()
if [ -n "${MOVIE_MAP_DATA:-}" ]; then
   REGION=${MOVIE_REGION:-$TOOL/regions/socal_ucerf2d.json}
   [ -f "$REGION" ] || { echo "ERROR: no region file $REGION"; exit 2; }
   [ -d "$MOVIE_MAP_DATA" ] || { echo "ERROR: no map data folder $MOVIE_MAP_DATA"; exit 2; }
   MAP=(--region "$REGION" --map-data "$MOVIE_MAP_DATA")
fi
mkdir -p "$OUT" || exit 2
{ echo "commit $HEAD"; echo "run $RUN"; echo "run finished $([ -e "$RUN/DONE" ] && echo yes || echo no)"
  echo "job $SLURM_JOB_ID node $SLURM_JOB_NODELIST workers $WORKERS"
  echo "python $("$PY" -c 'import sys, numpy, h5py, matplotlib; print(sys.version.split()[0], "numpy", numpy.__version__, "h5py", h5py.__version__, "matplotlib", matplotlib.__version__)')"
  echo "movie args ${MOVIE_ARGS:-none}"; echo "map ${MAP[*]:-none}"
  echo "magnitudes mu ${MAG_MU:-none} W ${MAG_WIDTH_M:-none} names ${MAG_FAULT_NAMES:-none}"; } > "$OUT/RENDER_INFO.txt"
cat "$OUT/RENDER_INFO.txt"

s=$(date +%s)
"$PY" "$TOOL/cycle_movies.py" --run "$RUN" --out "$OUT" --workers "$WORKERS" "${EXTRA[@]}" \
   || { echo "ERROR: cycle_movies.py failed"; exit 3; }
echo "videos: $(( $(date +%s) - s )) s"
if [ ${#MAP[@]} -gt 0 ]; then
   s=$(date +%s)
   "$PY" "$TOOL/cycle_movies.py" --run "$RUN" --out "$OUT" --workers "$WORKERS" "${MAP[@]}" "${EXTRA[@]}" \
      || { echo "ERROR: cycle_movies.py (line map) failed"; exit 3; }
   echo "videos on the line map: $(( $(date +%s) - s )) s"
fi
if [ -n "${MAG_MU:-}" ] && [ -n "${MAG_WIDTH_M:-}" ]; then
   NAMES=()
   [ -z "${MAG_FAULT_NAMES:-}" ] || NAMES=(--fault-names "$MAG_FAULT_NAMES")
   "$PY" "$TOOL/magnitude_catalog.py" --run "$RUN" --out "$OUT" --mu "$MAG_MU" --width-m "$MAG_WIDTH_M" \
      "${NAMES[@]}" || { echo "ERROR: magnitude_catalog.py failed"; exit 3; }
fi
ls -la "$OUT"
touch "$OUT/DONE"
