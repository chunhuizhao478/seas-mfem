#!/bin/bash
# =============================================================================
# Pull TPV104 Phase 6 result files from Frontera into a local archive tree.
#
# Run this on your LOCAL machine (Mac), not on Frontera.
#
# Auto-detects every `results_phase6_*` directory under
#   <REMOTE_ROOT>/miniapps/seas/tpv104/
# and rsync's each one to a same-named local subdir.  Reproduces the cluster's
# directory layout (ParaView/, ParaView_bulk/, top-level fault_surface.vtkhdf
# and station *.dat files).
#
# Default destination:
#   $HOME/Downloads/seas-mfem/tpv104
# Override with --dest <path> or the LOCAL_DEST environment variable.
#
# Usage:
#   ./sync_results.sh                                # all result dirs
#   ./sync_results.sh 7708789                        # only matching jobs
#   ./sync_results.sh 7708789 7710000                # multiple
#   ./sync_results.sh --dry-run                      # preview
#   ./sync_results.sh --dest /Volumes/SSD/seas/tpv104
#   ./sync_results.sh --host zhaochun@frontera.tacc.utexas.edu
#
# If --host is not given, $FRONTERA_HOST is used, falling back to the
# ssh alias `frontera` (configure in ~/.ssh/config).
#
# Whitelist (transferred):
#   *.vtkhdf  (fault_surface, ParaView/volume, ParaView_bulk/wave_bulk)
#   *.dat     (station outputs)
#
# Not transferred:
#   SLURM stdout/stderr, .msh meshes, source code, anything else.
# =============================================================================

set -u

REMOTE_ROOT="/scratch2/10024/zhaochun/seas-project/seas-mfem-paraview-compaction/miniapps/seas/tpv104"
DEFAULT_DEST="$HOME/Downloads/seas-mfem/tpv104"

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
        -h|--help)   sed -n '2,32p' "$0"; exit 0 ;;
        *)           JOB_FILTERS+=("$1"); shift ;;
    esac
done
REMOTE="${REMOTE:-${FRONTERA_HOST:-frontera}}"
LOCAL_DEST="${DEST_ARG:-${LOCAL_DEST:-$DEFAULT_DEST}}"

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

echo "Querying cluster for results_phase6_* directories..."
ALL_SUBDIRS=()
while IFS= read -r _line; do
    [[ -n "$_line" ]] && ALL_SUBDIRS+=("$_line")
done < <(
    ssh "$REMOTE" "cd '$REMOTE_ROOT' 2>/dev/null && ls -d results_phase6_*/ 2>/dev/null | sed 's|/$||'"
)
if [[ ${#ALL_SUBDIRS[@]} -eq 0 ]]; then
    echo "No results_phase6_* dirs found under $REMOTE:$REMOTE_ROOT" >&2
    exit 1
fi

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
        --include='*.vtkhdf' \
        --include='*.dat' \
        --exclude='*' \
        "$REMOTE:$REMOTE_ROOT/$sub/" \
        "$LOCAL_DEST/$sub/"
    echo
done

echo "Done."
