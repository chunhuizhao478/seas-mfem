#!/bin/bash
# =============================================================================
# Pull TPV31 SPATIAL dynamic-rupture result files from Frontera into a local
# archive tree for ParaView / SCEC-station inspection.
#
# Run this on your LOCAL machine (Mac), not on Frontera.
#
# The spatial driver (seas_spatial_dyn_driver) writes one output dir per run,
#   <REMOTE_ROOT>/out_<tag>_<jobid>/
# For TPV31 the jobs name it  out_p2_<jobid>  (p2/ADER-O3),  out_j1_<jobid>
# (p1/ADER-O2 job 1), or  out_<jobid>  (the bare dev smoke).  This script
# auto-detects every `out_*` dir and rsync's the THREE artefact families the
# on-fault analysis needs into a same-named local subdir:
#
#   (1) On-fault VTU   FaultSurface/fault_surface_c<cycle>.vtu   (--paraview-fault-vtu)
#   (2) Fault PVD      FaultSurface/fault_surface.pvd            (+ any *.pvd)
#   (3) Station traces tpv31_*.dat                               (30 SCEC on-fault stations)
#
# TPV31 is the ONLY spatial benchmark that emits SCEC station .dat traces (the
# TPV31StationWriter, one .dat per station) — compare these against
# tpv31/benchmark_data/scec_{eqdyna,seisol}/ via tpv31/visualize_results.py.
#
# The LARGE volume field (volume.vtkhdf / volume/*.vtu from --paraview) is NOT
# transferred by default — pass --volume to add it.
#
# Default destination:
#   $HOME/Downloads/seas-mfem/tpv31_spatial
# Override with --dest <path> or the LOCAL_DEST environment variable.
#
# Usage:
#   ./sync_results.sh                                # all out_* dirs (fault vtu + pvd + .dat)
#   ./sync_results.sh 7757861                        # only matching job(s)
#   ./sync_results.sh out_p2                          # only the p2/O3 run(s)
#   ./sync_results.sh --volume                       # also pull the volume field (large)
#   ./sync_results.sh --dry-run                      # preview, transfer nothing
#   ./sync_results.sh --dest /Volumes/SSD/seas/tpv31_spatial
#   ./sync_results.sh --host zhaochun@frontera.tacc.utexas.edu
#   ./sync_results.sh --remote-root /scratch2/.../miniapps/seas/tpv31    # non-default checkout
#
# If --host is not given, $FRONTERA_HOST is used, falling back to the ssh
# alias `frontera` (configure in ~/.ssh/config).
#
# Whitelist (default): fault_surface*.vtu + *.pvd + *.dat
# Whitelist (--volume): + *.vtu + *.pvtu + *.vtkhdf  (adds the volume field)
# Never transferred:   cp/ checkpoints, .msh meshes, *.sbatch, SLURM .out/.err, source.
# =============================================================================

set -u

# Silence the "perl: warning: Setting locale failed ... LANG = C.UTF-8" spam:
# macOS ssh forwards LANG/LC_*, but Frontera's Lmod login shell has no
# C.UTF-8 locale.  Force a locale the remote DOES have (plain C) for this
# script's children.  ASCII-only artefact names, so C is safe.
export LC_ALL=C LANG=C

BENCH="tpv31"
REMOTE_ROOT="/scratch2/10024/zhaochun/seas-project/seas-mfem-spatial-dyn-driver/miniapps/seas/${BENCH}"
DEFAULT_DEST="$HOME/Downloads/seas-mfem/${BENCH}_spatial"
REMOTE_GLOB="out_*"

# --- parse args ---
REMOTE=""
DRY=""
DEST_ARG=""
ROOT_ARG=""
WANT_VOLUME=""
JOB_FILTERS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)        DRY="--dry-run"; shift ;;
        --volume|--all)   WANT_VOLUME="1"; shift ;;
        --dest)           DEST_ARG="$2"; shift 2 ;;
        --dest=*)         DEST_ARG="${1#--dest=}"; shift ;;
        --host)           REMOTE="$2"; shift 2 ;;
        --host=*)         REMOTE="${1#--host=}"; shift ;;
        --remote-root)    ROOT_ARG="$2"; shift 2 ;;
        --remote-root=*)  ROOT_ARG="${1#--remote-root=}"; shift ;;
        -h|--help)        awk 'NR>1 && /^#/{print} NR>1 && !/^#/{exit}' "$0"; exit 0 ;;
        *)                JOB_FILTERS+=("$1"); shift ;;
    esac
done
REMOTE="${REMOTE:-${FRONTERA_HOST:-frontera}}"
LOCAL_DEST="${DEST_ARG:-${LOCAL_DEST:-$DEFAULT_DEST}}"
REMOTE_ROOT="${ROOT_ARG:-$REMOTE_ROOT}"

# Refuse to write under /Volumes if the drive is not mounted.
if [[ "$LOCAL_DEST" == /Volumes/* ]]; then
    vol_name="${LOCAL_DEST#/Volumes/}"
    vol_name="${vol_name%%/*}"
    if [[ ! -d "/Volumes/$vol_name" ]]; then
        echo "ERROR: /Volumes/$vol_name is not mounted." >&2
        exit 1
    fi
fi

mkdir -p "$LOCAL_DEST" || {
    echo "ERROR: cannot create $LOCAL_DEST." >&2
    exit 1
}

# Enumerate result dirs on the cluster (bare names, no path).
echo "Querying cluster for ${REMOTE_GLOB} directories under ${REMOTE_ROOT}..."
ALL_SUBDIRS=()
while IFS= read -r _line; do
    [[ -n "$_line" ]] && ALL_SUBDIRS+=("$_line")
done < <(
    ssh "$REMOTE" "cd '$REMOTE_ROOT' 2>/dev/null && ls -d ${REMOTE_GLOB}/ 2>/dev/null | sed 's|/\$||'"
)
if [[ ${#ALL_SUBDIRS[@]} -eq 0 ]]; then
    echo "No ${REMOTE_GLOB} dirs found under $REMOTE:$REMOTE_ROOT" >&2
    echo "(wrong checkout? pass --remote-root <path-to>/miniapps/seas/${BENCH})" >&2
    exit 1
fi

# Filter to matching jobs if positional args were given.
SUBDIRS=()
if [[ ${#JOB_FILTERS[@]} -gt 0 ]]; then
    for d in "${ALL_SUBDIRS[@]}"; do
        for f in "${JOB_FILTERS[@]}"; do
            if [[ "$d" == *"$f"* ]]; then
                SUBDIRS+=("$d")
                break
            fi
        done
    done
    if [[ ${#SUBDIRS[@]} -eq 0 ]]; then
        echo "No result dirs matched filters: ${JOB_FILTERS[*]}" >&2
        echo "Available: ${ALL_SUBDIRS[*]}" >&2
        exit 1
    fi
else
    SUBDIRS=("${ALL_SUBDIRS[@]}")
fi

# Build the rsync include set.  Fault first (always); volume only with --volume.
FILE_INCLUDES=(
    --include='fault_surface*.vtu'   # on-fault VTU pieces (binary + legacy ASCII)
    --include='fault_surface*.pvtu'  # legacy ASCII per-cycle index (if used)
    --include='*.pvd'                # collection files (fault_surface.pvd; tiny)
    --include='*.dat'                # SCEC on-fault station traces (tpv31_*.dat)
)
WHAT="on-fault VTU + PVD + station .dat"
if [[ -n "$WANT_VOLUME" ]]; then
    FILE_INCLUDES+=(
        --include='*.vtu'            # volume pieces too (fault already matched above)
        --include='*.pvtu'
        --include='*.vtkhdf'
    )
    WHAT="$WHAT + volume field (--volume)"
fi

echo "Remote:  $REMOTE:$REMOTE_ROOT"
echo "Local:   $LOCAL_DEST"
echo "Files:   $WHAT"
echo "Mode:    ${DRY:-live}"
echo "Subdirs: ${#SUBDIRS[@]}"
for s in "${SUBDIRS[@]}"; do echo "  - $s"; done
echo

for sub in "${SUBDIRS[@]}"; do
    echo "==> $sub"
    mkdir -p "$LOCAL_DEST/$sub"
    rsync -avz --progress $DRY \
        --prune-empty-dirs \
        --include='*/' \
        "${FILE_INCLUDES[@]}" \
        --exclude='*' \
        "$REMOTE:$REMOTE_ROOT/$sub/" \
        "$LOCAL_DEST/$sub/"
    echo
done

echo "Done."
