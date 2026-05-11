#!/usr/bin/env bash
# Phase 5 of miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md.
#
# Client-side rsync wrapper for downloading packed fault output from
# Frontera (or any other remote host) to a local workstation.  Uses
# `--partial --partial-dir=.rsync-partial --info=progress2` so a
# dropped network connection resumes from where it left off rather
# than restarting from zero.
#
# Usage:
#   ./download_fault.sh user@frontera:/scratch2/.../packs/job_xxx_fault \
#                       ~/Downloads/job_xxx_fault
#
# The remote path may be either:
#   - A directory containing chunk*.tar.zst files (recursively transferred)
#   - A single packed file (bulk-copied directly)
# The local path is the destination directory (created if missing).

set -euo pipefail

REMOTE="${1:?usage: $0 REMOTE_PATH LOCAL_DIR}"
LOCAL="${2:?usage: $0 REMOTE_PATH LOCAL_DIR}"

mkdir -p "${LOCAL}"

# `--partial-dir` is relative to the destination so partial files
# accumulate in `${LOCAL}/.rsync-partial/`.
RSYNC_OPTS=(
    --archive
    --partial
    --partial-dir=.rsync-partial
    --info=progress2
    --human-readable
    --compress         # zstd files are already compressed but rsync's
                       # own compression handles the small XML headers
)

# Resilient transfer: retry on transient SSH disconnect.
attempt=0
max_attempts=5
until rsync "${RSYNC_OPTS[@]}" "${REMOTE}" "${LOCAL}/"; do
    attempt=$((attempt + 1))
    if (( attempt >= max_attempts )); then
        echo "error: rsync failed after ${max_attempts} attempts" >&2
        exit 1
    fi
    echo "rsync exited non-zero; retrying (${attempt}/${max_attempts}) in 30s..." >&2
    sleep 30
done

echo "done.  contents of ${LOCAL}:"
ls -lh "${LOCAL}" | head -20
