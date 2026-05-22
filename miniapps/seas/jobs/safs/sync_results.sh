#!/bin/bash
# =============================================================================
# Pull SAFS spatial-dynamic-rupture result files from Frontera into a local
# archive tree for ParaView inspection.
#
# Run this on your LOCAL machine (Mac), not on Frontera.
#
# Auto-detects every `safs_dyn_*` directory under
#   <REMOTE_ROOT>   (default: the smoke output tree, see below)
# and rsync's each one to a same-named local subdir.  This covers both the
# smoke runs (safs_dyn_smoke_<jobid>) and the mixed-flux A/B runs
# (safs_dyn_{none,adjacent,all_continuous}_<jobid>).
#
# By DEFAULT only `fault.vtkhdf` is transferred (the fault-surface artefact
# you inspect).  Use --all to also grab volume.vtkhdf + ParaView_bulk/*.vtkhdf
# (much larger — ~5-80 GB).
#
# Default destination:
#   $HOME/Downloads/seas-mfem/safs
# Override with --dest <path> or the LOCAL_DEST environment variable.
#
# Usage:
#   ./sync_results.sh                                # fault.vtkhdf, all smoke dirs
#   ./sync_results.sh 7738686                        # only matching job(s)
#   ./sync_results.sh 7738686 7738999                # multiple
#   ./sync_results.sh none                           # mixed-flux A/B: --mixed-flux none run(s)
#   ./sync_results.sh adjacent                       # mixed-flux A/B: adjacent run(s)
#   ./sync_results.sh --all                          # fault + volume + bulk
#   ./sync_results.sh --dry-run                      # preview, transfer nothing
#   ./sync_results.sh --dest /Volumes/SSD/seas/safs
#   ./sync_results.sh --host zhaochun@frontera.tacc.utexas.edu
#   ./sync_results.sh --remote-root /scratch2/.../<JOBNAME>   # production run on $SCRATCH
#
# If --host is not given, $FRONTERA_HOST is used, falling back to the ssh
# alias `frontera` (configure in ~/.ssh/config).
#
# Whitelist (default): fault.vtkhdf
# Whitelist (--all):   *.vtkhdf  (fault, volume, ParaView_bulk/stress)
# Never transferred:   cp_* checkpoints, SLURM out/err, .msh meshes, source.
# =============================================================================

set -u

# Output tree.  Both sbatches write under
#   ${SEAS_MFEM_ROOT}/miniapps/seas/safs/ :
#   - spatial_dyn_smoke_*_safs.sbatch         -> safs_dyn_smoke_${SLURM_JOB_ID}
#   - spatial_dyn_smoke_nomixedflux_*.sbatch  -> safs_dyn_${MIXED_FLUX}_${SLURM_JOB_ID}
#       (MIXED_FLUX = none | adjacent | all_continuous)
# The glob `safs_dyn_*` matches both families.  For a PRODUCTION run
# (output on $SCRATCH/<JOBNAME>), pass --remote-root.
REMOTE_ROOT="/scratch2/10024/zhaochun/seas-project/seas-mfem-safs/miniapps/seas/safs"
REMOTE_GLOB="safs_dyn_*"
DEFAULT_DEST="$HOME/Downloads/seas-mfem/safs"

# --- parse args ---
REMOTE=""
DRY=""
DEST_ARG=""
ROOT_ARG=""
GRAB_ALL=""
JOB_FILTERS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)        DRY="--dry-run"; shift ;;
        --all)            GRAB_ALL="1"; shift ;;
        --dest)           DEST_ARG="$2"; shift 2 ;;
        --dest=*)         DEST_ARG="${1#--dest=}"; shift ;;
        --host)           REMOTE="$2"; shift 2 ;;
        --host=*)         REMOTE="${1#--host=}"; shift ;;
        --remote-root)    ROOT_ARG="$2"; shift 2 ;;
        --remote-root=*)  ROOT_ARG="${1#--remote-root=}"; shift ;;
        -h|--help)        sed -n '2,40p' "$0"; exit 0 ;;
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

# Build the rsync include set.
if [[ -n "$GRAB_ALL" ]]; then
    FILE_INCLUDES=( --include='*.vtkhdf' )
    WHAT="all *.vtkhdf (fault + volume + bulk)"
else
    FILE_INCLUDES=( --include='fault.vtkhdf' )
    WHAT="fault.vtkhdf only (use --all for volume + bulk)"
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
