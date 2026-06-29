#!/bin/bash
# =============================================================================
# Pull BP5-SPATIAL quasi-dynamic (seas_spatial_seas_driver) result files from
# EXPANSE (SDSC) into a local archive tree for inspection.  Sibling of
# jobs/safs/safs_expanse/sync_results.sh, adapted for the bp5_spatial/expanse
# Caliper jobs (seas_spatial_seas_driver, CG-AMG, --max-steps 100).
#
# Run this on your LOCAL machine (Mac), not on Expanse.
#
# Two things are transferred:
#
#   (1) RUN ARTEFACTS: the run job writes ONE output dir per run on qstore:
#         <REMOTE_ROOT>/bp5_seas_cgamg_cali_<jobid>/
#       Because [output].bp5_stations = true, the SCEC station traces
#       (spatial_seas_*.txt + the global V_max probe) land there — these are the
#       science record.  By DEFAULT the station .txt/.dat traces + any *.pvd are
#       transferred; use --all to also grab fault/volume *.vtu / *.vtkhdf (the
#       Caliper job runs ParaView "off" by default, so there usually is none).
#
#   (2) CALIPER REPORTS + DIAGNOSTIC LOGS (default ON): from the run job's
#       LOG_DIR (the cluster's miniapps/seas/bp5_spatial/expanse):
#         bp5_spatial_seas_cgamg_cali_<jobid>.log         (driver stdout/stderr)
#         bp5_spatial_seas_cgamg_cali_<jobid>.cali-region-report.txt   (sparse)
#         bp5_spatial_seas_cgamg_cali_<jobid>.cali-sample-report.txt   (PRIMARY)
#         bp5_spatial_seas_cgamg_cali_<jobid>.cali                     (Spot/Hatchet)
#         bp5seas_cg_cali_<jobid>.{out,err}               (SLURM)
#       into a local `logs/` subdir.  Use --no-logs to skip them, or --logs-only
#       to grab just the reports/logs and skip the (possibly large) run dir.
#
# Default destination:
#   $HOME/Downloads/seas-mfem/bp5_spatial_expanse
#     ( run dirs in <dest>/<dir>/, reports+logs in <dest>/logs/ )
# Override with --dest <path> or the LOCAL_DEST environment variable.
#
# Usage:
#   ./sync_results.sh                       # station traces (all dirs) + cali reports + logs
#   ./sync_results.sh 51443560              # only matching job id(s)
#   ./sync_results.sh cgamg                  # only matching name fragment(s)
#   ./sync_results.sh --logs-only 51443560  # just the .cali reports + .log/.out/.err
#   ./sync_results.sh --no-logs             # run dirs only, skip reports/logs
#   ./sync_results.sh --all                 # + fault/volume vtu/vtkhdf (if any) (+ logs)
#   ./sync_results.sh --dry-run             # preview, transfer nothing
#   ./sync_results.sh --dest /Volumes/SSD/seas/bp5_spatial_expanse
#   ./sync_results.sh --host czhao1@login.expanse.sdsc.edu
#   ./sync_results.sh --remote-root /path/to/seas-mfem/miniapps/seas/bp5/spatial_seas_expanse
#   ./sync_results.sh --remote-log-dir /path/to/seas-mfem/miniapps/seas/bp5_spatial/expanse
#
# A JOB ID is the most reliable filter — it appears in every artefact name
# (bp5_seas_*_<jobid>, bp5_spatial_seas_*_<jobid>.log, bp5seas_cg_cali_<jobid>.out).
#
# If --host is not given, $EXPANSE_HOST is used, falling back to the explicit
# czhao1@login.expanse.sdsc.edu (set up an ssh alias `expanse` in ~/.ssh/config
# and pass --host expanse to avoid typing the full login).
#
# Whitelist (default): spatial_seas_*.txt + *.dat + *.pvd
#                      + bp5_spatial_seas_*.log + *.cali* reports + SLURM *.out/*.err
# Whitelist (--all):   + *.vtu + *.vtkhdf  (adds the fault/volume fields, if any)
# Never transferred:   checkpoints (cp_*/*.txt state dumps), .msh meshes, *.sbatch, source.
# =============================================================================

set -u

# Silence the "perl: warning: Setting locale failed" spam (Expanse's login shell
# has no C.UTF-8 locale); ASCII-only artefact names, so plain C is safe.
export LC_ALL=C LANG=C

# Expanse checkout root (the dir that contains miniapps/seas/...).  The run job
# writes station traces under <checkout>/miniapps/seas/bp5/spatial_seas_expanse
# and the .log/.cali reports under <checkout>/miniapps/seas/bp5_spatial/expanse.
REMOTE_CHECKOUT="/expanse/projects/qstore/usc143/qwxdev/apps/expanse/rocky8.8/mfem_seas_versions/seas-mfem"
REMOTE_ROOT="${REMOTE_CHECKOUT}/miniapps/seas/bp5/spatial_seas_expanse"
REMOTE_GLOB="bp5_seas_*"
REMOTE_LOG_DIR_DEFAULT="${REMOTE_CHECKOUT}/miniapps/seas/bp5_spatial/expanse"
LOG_GLOB="bp5_spatial_seas_*.log"
DEFAULT_DEST="$HOME/Downloads/seas-mfem/bp5_spatial_expanse"

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
# (1) Run-dir artefacts (station traces; skipped entirely with --logs-only).
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
        echo "(no station output yet -- job may not have reached the first write; try --logs-only." >&2
        echo " wrong checkout? --remote-root <path>/miniapps/seas/bp5/spatial_seas_expanse ; wrong host? --host ...)" >&2
        [[ -n "$SKIP_LOGS" ]] && exit 1
    else
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

# Build the rsync include set for the run dirs.  Station traces always; the large
# fault/volume fields only with --all.
FILE_INCLUDES=(
    --include='spatial_seas_*.txt'   # SCEC BP5 station traces + global V_max probe
    --include='*.dat'                # any companion .dat traces
    --include='*.pvd'                # ParaView collection files (tiny; if PV was on)
)
WHAT="BP5 station traces (spatial_seas_*.txt) + .dat + PVD (use --all for vtu/vtkhdf)"
if [[ -n "$GRAB_ALL" ]]; then
    FILE_INCLUDES+=(
        --include='*.vtu'
        --include='*.vtkhdf'
        --include='*.pvtu'
    )
    WHAT="all run artefacts (station traces + fault/volume vtu/vtkhdf)"
fi

echo "Remote:  $REMOTE:$REMOTE_ROOT"
echo "Local:   $LOCAL_DEST"
echo "Files:   $WHAT"
echo "Logs:    $([[ -n "$SKIP_LOGS" ]] && echo 'skipped (--no-logs)' || echo "$REMOTE_LOG_DIR -> $LOCAL_DEST/logs/ (.log + .cali reports + SLURM .out/.err)")"
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
# (2) Caliper reports + diagnostic logs (driver .log + .cali* + SLURM .out/.err).
# =============================================================================
if [[ -z "$SKIP_LOGS" ]]; then
    echo "==> reports + logs (${LOG_GLOB} + *.cali* + SLURM .out/.err)"
    mkdir -p "$LOCAL_DEST/logs"
    if [[ ${#JOB_FILTERS[@]} -gt 0 ]]; then
        LOG_INCLUDES=()
        for f in "${JOB_FILTERS[@]}"; do
            LOG_INCLUDES+=( --include="*${f}*.log" \
                            --include="*${f}*.cali" \
                            --include="*${f}*.cali-region-report.txt" \
                            --include="*${f}*.cali-sample-report.txt" \
                            --include="*${f}*.out" \
                            --include="*${f}*.err" )
        done
    else
        LOG_INCLUDES=( --include="${LOG_GLOB}" \
                       --include='*.cali' \
                       --include='*.cali-region-report.txt' \
                       --include='*.cali-sample-report.txt' \
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
