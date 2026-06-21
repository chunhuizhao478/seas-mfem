#!/bin/bash
# =============================================================================
# Pull SAFS SPATIAL dynamic-rupture result files from EXPANSE (SDSC) into a
# local archive tree for ParaView inspection.  (Follows the tpv6/tpv7_spatial
# sync_results.sh pattern, adapted for Expanse + the safs_expanse VTU jobs.)
#
# Run this on your LOCAL machine (Mac), not on Expanse.
#
# Two things are transferred:
#
#   (1) ParaView artefacts (VTU): the safs_expanse run jobs
#       (jobs/safs/safs_expanse/*.sbatch) use --paraview-fault-vtu +
#       --paraview-free-surface vtu, writing ONE output dir per run UNDER the
#       repo on qstore:
#         <REMOTE_ROOT>/safs_dyn_<case>_<jobid>/
#       e.g.  safs_dyn_sw_alt_exp_51210086     (slipweakening ALT, mixed flux)
#             safs_dyn_rs_alt_exp_51210073     (ratestate   ALT, mixed flux)
#             safs_dyn_sw_alt_uw_exp_51211002  (slipweakening ALT, upwind+ADER)
#       This auto-detects every `safs_dyn_*` dir and rsync's each to a same-named
#       local subdir.  By DEFAULT the on-fault VTU (FaultSurface/fault_surface_
#       c<cycle>.vtu) + free-surface VTU + the *.pvd collections are transferred;
#       use --all to also grab any volume *.vtu / *.vtkhdf (SAFS runs
#       paraview_volume="off", so there usually is none).  SAFS has no
#       [problem].tag, so there are NO per-side on-fault station .dat traces
#       (unlike TPV6/7) — the on-fault VTU IS the rupture record.
#
#   (2) Diagnostic logs (default ON): the driver log
#       `spatial_dyn_<tag>_<jobid>.log` (V_max / [DIAG] / [R-101] lines) plus the
#       SLURM `*_<jobid>.{out,err}` from <REMOTE_LOG_DIR> (the cluster's
#       jobs/safs/safs_expanse), into a local `logs/` subdir.  Use --no-logs to
#       skip them, or --logs-only to grab just the .log/.out/.err and skip the
#       (possibly large) VTU.
#
# Default destination:
#   $HOME/Downloads/seas-mfem/safs_expanse  ( VTU in <dest>/<dir>/, logs in <dest>/logs/ )
# Override with --dest <path> or the LOCAL_DEST environment variable.
#
# Usage:
#   ./sync_results.sh                                # fault+free-surface VTU (all dirs) + logs
#   ./sync_results.sh 51210086                       # only matching job id(s)
#   ./sync_results.sh sw_alt                          # only the slipweakening ALT run(s)
#   ./sync_results.sh _uw_                            # only the upwind+ADER run(s)
#   ./sync_results.sh --logs-only sw_alt             # just the .log/.out/.err (skip VTU)
#   ./sync_results.sh --no-logs                       # VTU only, skip logs
#   ./sync_results.sh --all                          # + volume/vtkhdf (if any) (+ logs)
#   ./sync_results.sh --dry-run                      # preview, transfer nothing
#   ./sync_results.sh --dest /Volumes/SSD/seas/safs_expanse
#   ./sync_results.sh --host czhao1@login.expanse.sdsc.edu
#   ./sync_results.sh --remote-root /expanse/.../seas-mfem/miniapps/seas/safs
#   ./sync_results.sh --remote-log-dir /path/to/jobs/safs/safs_expanse
#
# A JOB ID is the most reliable filter — it appears in every artefact name
# (safs_dyn_*_<jobid>, spatial_dyn_*_<jobid>.log, *_<jobid>.out).
#
# If --host is not given, $EXPANSE_HOST is used, falling back to the explicit
# czhao1@login.expanse.sdsc.edu (set up an ssh alias `expanse` in ~/.ssh/config
# and pass --host expanse to avoid typing the full login).
#
# Whitelist (default): fault_surface*.vtu + free_surface*.vtu + *.pvd
#                      + spatial_dyn_*.log + SLURM *.out/*.err
# Whitelist (--all):   + *.vtu + *.vtkhdf  (adds the volume field, if any)
# Never transferred:   cp_checkpoint_r*.txt, .msh meshes, *.sbatch, this script, source.
# =============================================================================

set -u

# Silence the "perl: warning: Setting locale failed ... LANG = C.UTF-8" spam:
# macOS ssh forwards LANG/LC_*, but Expanse's Lmod login shell has no C.UTF-8
# locale, so every ssh/rsync prints the warning.  Force plain C (ASCII-only
# artefact names, so it is safe) for this script's children.
export LC_ALL=C LANG=C

# Expanse checkout (a <checkout>/miniapps/seas/safs dir).  The safs_expanse run
# jobs write OUT here (OUT_BASE = <root>/miniapps/seas/safs).  --remote-root
# overrides the whole path.
REMOTE_ROOT="/expanse/projects/qstore/usc143/qwxdev/apps/expanse/rocky8.8/mfem_seas_versions/seas-mfem/miniapps/seas/safs"
REMOTE_GLOB="safs_dyn_*"
# Driver .log + SLURM .out/.err live in the safs_expanse job dir (LOG_DIR there).
# Default = <checkout>/miniapps/seas/jobs/safs/safs_expanse.  Override
# --remote-log-dir.
REMOTE_LOG_DIR_DEFAULT="${REMOTE_ROOT%/safs}/jobs/safs/safs_expanse"
LOG_GLOB="spatial_dyn_*.log"
DEFAULT_DEST="$HOME/Downloads/seas-mfem/safs_expanse"

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
        --all|--volume)     GRAB_ALL="1"; shift ;;
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
REMOTE="${REMOTE:-${EXPANSE_HOST:-czhao1@login.expanse.sdsc.edu}}"
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
# (1) ParaView VTU artefacts (skipped entirely with --logs-only).
# =============================================================================
SUBDIRS=()
if [[ -z "$LOGS_ONLY" ]]; then
    echo "Querying Expanse for ${REMOTE_GLOB} directories under ${REMOTE_ROOT}..."
    ALL_SUBDIRS=()
    while IFS= read -r _line; do
        [[ -n "$_line" ]] && ALL_SUBDIRS+=("$_line")
    done < <(
        ssh "$REMOTE" "cd '$REMOTE_ROOT' 2>/dev/null && ls -d ${REMOTE_GLOB}/ 2>/dev/null | sed 's|/\$||'"
    )
    if [[ ${#ALL_SUBDIRS[@]} -eq 0 ]]; then
        echo "No ${REMOTE_GLOB} dirs found under $REMOTE:$REMOTE_ROOT" >&2
        echo "(no VTU yet -- job may not have reached the first snapshot; try --logs-only." >&2
        echo " wrong checkout? --remote-root <path>/miniapps/seas/safs ; wrong host? --host ...)" >&2
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

# Build the rsync include set for VTU.  Fault + free surface always; the large
# volume field only with --all.
FILE_INCLUDES=(
    --include='fault_surface*.vtu'   # on-fault VTU, one fault_surface_c<cycle>.vtu per cycle
    --include='free_surface*.vtu'    # free-surface VTU, one free_surface_c<cycle>.vtu per cycle
                                     #   (both are SINGLE consolidated files, rank-0 gathered — no
                                     #   per-rank proc*.vtu anymore).
    --include='*.pvtu'               # (legacy per-rank index, if an old run is still around)
    --include='*.pvd'                # collection files (fault_surface.pvd, free_surface.pvd; tiny)
)
WHAT="on-fault VTU + free-surface VTU + PVD (use --all for volume)"
if [[ -n "$GRAB_ALL" ]]; then
    FILE_INCLUDES+=(
        --include='*.vtu'            # any remaining volume pieces
        --include='*.vtkhdf'
    )
    WHAT="all VTU/VTKHDF (fault + free surface + volume)"
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
