#!/bin/bash
# =============================================================================
# Pull SAFS spatial-dynamic-rupture result files from Frontera into a local
# archive tree for ParaView inspection.
#
# Run this on your LOCAL machine (Mac), not on Frontera.
#
# Two things are transferred:
#
#   (1) ParaView artefacts (*.vtkhdf): auto-detects every `safs_dyn_*`
#       directory under <REMOTE_ROOT> and rsync's each to a same-named local
#       subdir.  This covers the smoke runs (safs_dyn_smoke_<jobid>), the
#       mixed-flux A/B runs (safs_dyn_{none,adjacent,all_continuous}_<jobid>),
#       the resolution sweep (safs_dyn_resDc2_<jobid>), and the dip-prestress
#       diagnostic (safs_dyn_{control,zerodip}_<jobid>).  By DEFAULT only
#       `fault.vtkhdf` is transferred; use --all for volume + bulk too.
#
#   (2) Diagnostic logs (default ON): the driver log
#       `spatial_dyn_<tag>_<jobid>.log` plus the SLURM `<name>_<jobid>.{out,err}`
#       from <REMOTE_LOG_DIR> (the cluster's miniapps/seas/jobs/safs), into a
#       local `logs/` subdir.  The scalar diagnostics (V_max + [DIAG] +
#       [R-101 nonfatal] lines) live in the .log; the fault.vtkhdf carries the
#       slip-rate field for ParaView.  Use --no-logs to skip the logs, or
#       --logs-only to grab just the .log/.out/.err and skip the (possibly
#       large) fault.vtkhdf.
#
# Default destination:
#   $HOME/Downloads/seas-mfem/safs   ( *.vtkhdf in <dest>/<dir>/, logs in <dest>/logs/ )
# Override with --dest <path> or the LOCAL_DEST environment variable.
#
# Usage:
#   ./sync_results.sh                                # fault.vtkhdf (all dirs) + logs
#   ./sync_results.sh 7738686                        # only matching job(s)
#   ./sync_results.sh 7738686 7738999                # multiple
#   ./sync_results.sh none                           # mixed-flux A/B: --mixed-flux none run(s)
#   ./sync_results.sh adjacent                       # mixed-flux A/B: adjacent run(s)
#   ./sync_results.sh zerodip                         # dip diagnostic: zerodip run(s)
#   ./sync_results.sh control                         # dip diagnostic: control run(s)
#   ./sync_results.sh --logs-only zerodip            # just the .log/.out/.err (skip fault.vtkhdf)
#   ./sync_results.sh --no-logs                       # vtkhdf only, skip logs
#   ./sync_results.sh --all                          # fault + volume + bulk (+ logs)
#   ./sync_results.sh --dry-run                      # preview, transfer nothing
#   ./sync_results.sh --dest /Volumes/SSD/seas/safs
#   ./sync_results.sh --host zhaochun@frontera.tacc.utexas.edu
#   ./sync_results.sh --remote-root /scratch2/.../<JOBNAME>      # production run on $SCRATCH
#   ./sync_results.sh --remote-log-dir /path/to/jobs/safs        # logs from a non-default checkout
#
# Note on filters: a JOB ID is the most reliable filter — it appears in every
# artefact name (safs_dyn_*_<jobid>, spatial_dyn_*_<jobid>.log, *_<jobid>.out).
# A tag like `control` matches the .log but NOT the SLURM .out/.err for the
# zerodip job (whose -o is fixed at safs_zerodip_<jobid> regardless of tag).
#
# If --host is not given, $FRONTERA_HOST is used, falling back to the ssh
# alias `frontera` (configure in ~/.ssh/config).
#
# Whitelist (default): fault.vtkhdf  + spatial_dyn_*.log + SLURM *.out/*.err
# Whitelist (--all):   *.vtkhdf      + the same logs  (fault, volume, bulk/stress)
# Never transferred:   cp_* checkpoints, .msh meshes, *.sbatch, this script, source.
# =============================================================================

set -u

# Silence the "perl: warning: Setting locale failed ... LANG = C.UTF-8" spam:
# macOS ssh forwards the local LANG/LC_* (SendEnv LANG LC_* in the default
# ssh_config), but Frontera's Perl-based login shell (Lmod) has no C.UTF-8
# locale, so every ssh/rsync connection prints the warning.  Forcing a locale
# the remote DOES have (plain C) for this script's children makes the forwarded
# value valid -> no warning.  ASCII-only artefact names, so C is safe.
export LC_ALL=C LANG=C

# Output tree.  All dev-queue sbatches write under
#   ${SEAS_MFEM_ROOT}/miniapps/seas/safs/ :
#   - spatial_dyn_smoke_*_safs.sbatch          -> safs_dyn_smoke_${SLURM_JOB_ID}
#   - spatial_dyn_smoke_nomixedflux_*.sbatch   -> safs_dyn_${MIXED_FLUX}_${SLURM_JOB_ID}
#       (MIXED_FLUX = none | adjacent | all_continuous)
#   - spatial_dyn_resolution_Dc2_*.sbatch      -> safs_dyn_resDc2_${SLURM_JOB_ID}
#   - spatial_dyn_zerodip_*.sbatch             -> safs_dyn_${TAG}_${SLURM_JOB_ID}
#       (TAG = control | zerodip; writes fault.vtkhdf like the smoke job)
# and ALL of them write the driver log + SLURM out/err to
#   ${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs/ .
# The glob `safs_dyn_*` matches every vtkhdf family.  For a PRODUCTION run
# (output on $SCRATCH/<JOBNAME>), pass --remote-root.
REMOTE_ROOT="/scratch2/10024/zhaochun/seas-project/seas-mfem-safs/miniapps/seas/safs"
REMOTE_GLOB="safs_dyn_*"
# Logs always live in the repo's jobs/safs (independent of where vtkhdf go);
# default is the sibling of the default REMOTE_ROOT.  Override --remote-log-dir.
REMOTE_LOG_DIR_DEFAULT="${REMOTE_ROOT%/safs}/jobs/safs"
LOG_GLOB="spatial_dyn_*.log"
DEFAULT_DEST="$HOME/Downloads/seas-mfem/safs"

# --- parse args ---
REMOTE=""
DRY=""
DEST_ARG=""
ROOT_ARG=""
LOGDIR_ARG=""
GRAB_ALL=""
SKIP_LOGS=""
LOGS_ONLY=""
JOB_FILTERS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)          DRY="--dry-run"; shift ;;
        --all)              GRAB_ALL="1"; shift ;;
        --no-logs)          SKIP_LOGS="1"; shift ;;
        --logs-only)        LOGS_ONLY="1"; shift ;;
        --dest)             DEST_ARG="$2"; shift 2 ;;
        --dest=*)           DEST_ARG="${1#--dest=}"; shift ;;
        --host)             REMOTE="$2"; shift 2 ;;
        --host=*)           REMOTE="${1#--host=}"; shift ;;
        --remote-root)      ROOT_ARG="$2"; shift 2 ;;
        --remote-root=*)    ROOT_ARG="${1#--remote-root=}"; shift ;;
        --remote-log-dir)   LOGDIR_ARG="$2"; shift 2 ;;
        --remote-log-dir=*) LOGDIR_ARG="${1#--remote-log-dir=}"; shift ;;
        -h|--help)          awk 'NR>1 && /^#/{print} NR>1 && !/^#/{exit}' "$0"; exit 0 ;;
        *)                  JOB_FILTERS+=("$1"); shift ;;
    esac
done
REMOTE="${REMOTE:-${FRONTERA_HOST:-frontera}}"
LOCAL_DEST="${DEST_ARG:-${LOCAL_DEST:-$DEFAULT_DEST}}"
REMOTE_ROOT="${ROOT_ARG:-$REMOTE_ROOT}"
REMOTE_LOG_DIR="${LOGDIR_ARG:-$REMOTE_LOG_DIR_DEFAULT}"

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

# =============================================================================
# (1) ParaView vtkhdf artefacts (skipped entirely with --logs-only).
# =============================================================================
SUBDIRS=()
if [[ -z "$LOGS_ONLY" ]]; then
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
        echo "(no vtkhdf yet -- job may have aborted before the first snapshot; try --logs-only)" >&2
        [[ -n "$SKIP_LOGS" ]] && exit 1
    else
        # Filter to matching jobs if positional args were given.
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
            fi
        else
            SUBDIRS=("${ALL_SUBDIRS[@]}")
        fi
    fi
fi

# Build the rsync include set for vtkhdf.
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
echo "Logs:    $([[ -n "$SKIP_LOGS" ]] && echo 'skipped (--no-logs)' || echo "$REMOTE_LOG_DIR -> $LOCAL_DEST/logs/")"
echo "Mode:    ${DRY:-live}$([[ -n "$LOGS_ONLY" ]] && echo '  (logs-only)')"
if [[ -z "$LOGS_ONLY" ]]; then
    echo "Subdirs: ${#SUBDIRS[@]}"
    for s in "${SUBDIRS[@]}"; do echo "  - $s"; done
fi
echo

if [[ -z "$LOGS_ONLY" ]]; then
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
fi

# =============================================================================
# (2) Diagnostic logs (driver .log + SLURM .out/.err), default ON.
# =============================================================================
if [[ -z "$SKIP_LOGS" ]]; then
    echo "==> logs (${LOG_GLOB} + SLURM .out/.err)"
    mkdir -p "$LOCAL_DEST/logs"
    if [[ ${#JOB_FILTERS[@]} -gt 0 ]]; then
        LOG_INCLUDES=()
        for f in "${JOB_FILTERS[@]}"; do
            LOG_INCLUDES+=( --include="*${f}*.log" \
                            --include="*${f}*.out" \
                            --include="*${f}*.err" )
        done
    else
        LOG_INCLUDES=( --include="${LOG_GLOB}" \
                       --include='*.out' \
                       --include='*.err' )
    fi
    rsync -avz --progress $DRY \
        "${LOG_INCLUDES[@]}" \
        --exclude='*' \
        "$REMOTE:$REMOTE_LOG_DIR/" \
        "$LOCAL_DEST/logs/"
    echo
fi

echo "Done."
