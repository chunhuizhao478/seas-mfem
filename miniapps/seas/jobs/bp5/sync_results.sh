#!/bin/bash
# =============================================================================
# Pull BP5 Phase 6 result files from Frontera into a local archive tree.
#
# Run this on your LOCAL machine (Mac), not on Frontera.
#
# Auto-detects every `results_phase6_*` directory under
#   <REMOTE_ROOT>/miniapps/seas/bp5/
# and rsync's each one to a same-named local subdir.  BP5 has no bulk
# collection (no ParaView_bulk/), but otherwise mirrors the same layout
# as the TPV runs: ParaView/, top-level fault_surface.vtkhdf, and station
# *.dat probe outputs.
#
# Default destination:
#   $HOME/Downloads/seas-mfem/bp5
# Override with --dest <path> or the LOCAL_DEST environment variable.
#
# Usage:
#   ./sync_results.sh                                # all result dirs
#   ./sync_results.sh 7708789                        # only matching jobs
#   ./sync_results.sh 7708789 7710000                # multiple
#   ./sync_results.sh --dry-run                      # preview
#   ./sync_results.sh --dest /Volumes/SSD/seas/bp5
#   ./sync_results.sh --host zhaochun@frontera.tacc.utexas.edu
#
# If --host is not given, $FRONTERA_HOST is used, falling back to the
# ssh alias `frontera` (configure in ~/.ssh/config).
#
# Whitelist (transferred):
#   *.vtkhdf  (fault_surface, ParaView/volume)
#   *.dat     (station outputs)
#
# Not transferred:
#   SLURM stdout/stderr, .msh meshes, source code, anything else.
#   BP5 also emits benchmark CSVs (*.csv) — not included in the default
#   whitelist; add `--include='*.csv'` above `--exclude='*'` below if
#   you also want those.
# =============================================================================

set -u

REMOTE_ROOT="/scratch2/10024/zhaochun/seas-project/seas-mfem-paraview-compaction/miniapps/seas/bp5"
DEFAULT_DEST="$HOME/Downloads/seas-mfem/bp5"

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
        -h|--help)   sed -n '2,36p' "$0"; exit 0 ;;
        *)
            # Bug-fix: catch the common mistake of passing an SSH
            # hostname as the first positional argument (which would
            # otherwise be silently swallowed into JOB_FILTERS and the
            # script would try ssh to the default `frontera` alias).
            # If `$1` looks like a hostname (contains `@` or a dot
            # plus no leading digit), bail with a hint.
            if [[ "$1" == *"@"* ]] || [[ "$1" == *.*.* && ! "$1" =~ ^[0-9] ]]; then
                echo "ERROR: '$1' looks like an SSH hostname, not a SLURM job ID." >&2
                echo "       Use --host '$1' or set FRONTERA_HOST='$1'." >&2
                echo "       SLURM job IDs are bare integers like 7725553." >&2
                exit 2
            fi
            JOB_FILTERS+=("$1"); shift ;;
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

echo "Querying cluster for results_phase6_* directories on $REMOTE..."

# Bug-fix: do the ssh call in a separate step so we can distinguish
# "ssh failed (network/auth/DNS)" from "ssh succeeded but no dirs
# found".  Previously both produced the same misleading "No
# results_phase6_* dirs found" message, hiding the actual ssh error
# (e.g., "Could not resolve hostname frontera") in the noise above.
SSH_OUTPUT=$(ssh -o BatchMode=yes -o ConnectTimeout=10 "$REMOTE" \
    "cd '$REMOTE_ROOT' 2>/dev/null && ls -d results_phase6_*/ 2>/dev/null | sed 's|/$||'" 2>&1)
SSH_RC=$?
if [[ $SSH_RC -ne 0 ]]; then
    echo "ERROR: ssh to '$REMOTE' failed (exit $SSH_RC):" >&2
    echo "  $SSH_OUTPUT" >&2
    echo >&2
    if [[ "$REMOTE" == "frontera" ]]; then
        echo "Hint: the default host alias 'frontera' requires an entry in" >&2
        echo "  ~/.ssh/config like:" >&2
        echo "    Host frontera" >&2
        echo "        HostName frontera.tacc.utexas.edu" >&2
        echo "        User <your-tacc-username>" >&2
        echo "  OR pass --host directly:" >&2
        echo "    ./sync_results.sh --host <user>@frontera.tacc.utexas.edu $*" >&2
        echo "  OR export FRONTERA_HOST=<user>@frontera.tacc.utexas.edu" >&2
    fi
    exit 1
fi

ALL_SUBDIRS=()
while IFS= read -r _line; do
    [[ -n "$_line" ]] && ALL_SUBDIRS+=("$_line")
done <<< "$SSH_OUTPUT"

if [[ ${#ALL_SUBDIRS[@]} -eq 0 ]]; then
    echo "ssh to '$REMOTE' succeeded but found no results_phase6_* dirs under" >&2
    echo "  $REMOTE_ROOT" >&2
    echo "(check that the job has actually produced output; the dir is created" >&2
    echo " by the sbatch's `mkdir -p \"\${RESULT_DIR}\"` call near the top.)" >&2
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
