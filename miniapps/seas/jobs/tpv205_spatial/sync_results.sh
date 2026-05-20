#!/bin/bash
# =============================================================================
# Pull TPV205 spatial-dyn STATION data files from Frontera into a local
# archive tree.
#
# Run this on your LOCAL machine (Mac), not on Frontera.
#
# Auto-detects every `results_spatial_dyn_*` directory under
#   <REMOTE_ROOT>/miniapps/seas/tpv205/
# and rsync's each one to a same-named local subdir.  Reproduces the cluster's
# directory layout but transfers ONLY the station *.dat traces (no ParaView
# fault/volume .vtkhdf — those are large and intentionally skipped here).
#
# Default destination:
#   $HOME/Downloads/seas-mfem/tpv205_spatial
# Override with --dest <path> or the LOCAL_DEST environment variable.
#
# Usage:
#   ./sync_results.sh                                # all result dirs
#   ./sync_results.sh 7736493                        # only matching jobs
#   ./sync_results.sh 7736493 7740000                # multiple
#   ./sync_results.sh --dry-run                      # preview
#   ./sync_results.sh --dest /Volumes/SSD/seas/tpv205_spatial
#   ./sync_results.sh --host zhaochun@frontera.tacc.utexas.edu
#
# If --host is not given, $FRONTERA_HOST is used, falling back to the
# ssh alias `frontera` (configure in ~/.ssh/config).
#
# Whitelist (transferred):
#   *.dat     (station outputs — e.g. tpv205_x2_0_x3_7.5.dat)
#
# Not transferred:
#   *.vtkhdf ParaView output, SLURM stdout/stderr, .msh meshes, RESULT.txt,
#   source code, anything else.
# =============================================================================

set -u

REMOTE_ROOT="/scratch2/10024/zhaochun/seas-project/seas-mfem-heterogeneous-riemann-solver/miniapps/seas/tpv205"
DEFAULT_DEST="$HOME/Downloads/seas-mfem/tpv205_spatial"

# --- parse args ---
REMOTE=""
DRY=""
DEST_ARG=""
JOB_FILTERS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)   DRY="--dry-run"; shift ;;
        --dest)      DEST_ARG="$2"; shift 2 ;;
        --dest=*)    DEST_ARG="${1#--dest=}"; shift ;;
        --host)      REMOTE="$2"; shift 2 ;;
        --host=*)    REMOTE="${1#--host=}"; shift ;;
        -h|--help)   sed -n '2,34p' "$0"; exit 0 ;;
        *)           JOB_FILTERS+=("$1"); shift ;;
    esac
done
REMOTE="${REMOTE:-${FRONTERA_HOST:-frontera}}"
LOCAL_DEST="${DEST_ARG:-${LOCAL_DEST:-$DEFAULT_DEST}}"

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

# Enumerate result dirs on the cluster.  Each entry of SUBDIRS is a bare
# directory name (no path) like "results_spatial_dyn_p1_O2_dev_job7736493".
echo "Querying cluster for results_spatial_dyn_* directories..."
ALL_SUBDIRS=()
while IFS= read -r _line; do
    [[ -n "$_line" ]] && ALL_SUBDIRS+=("$_line")
done < <(
    ssh "$REMOTE" "cd '$REMOTE_ROOT' 2>/dev/null && ls -d results_spatial_dyn_*/ 2>/dev/null | sed 's|/$||'"
)
if [[ ${#ALL_SUBDIRS[@]} -eq 0 ]]; then
    echo "No results_spatial_dyn_* dirs found under $REMOTE:$REMOTE_ROOT" >&2
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

echo "Remote:  $REMOTE:$REMOTE_ROOT"
echo "Local:   $LOCAL_DEST"
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
        --include='*.dat' \
        --exclude='*' \
        "$REMOTE:$REMOTE_ROOT/$sub/" \
        "$LOCAL_DEST/$sub/"
    echo
done

echo "Done."
