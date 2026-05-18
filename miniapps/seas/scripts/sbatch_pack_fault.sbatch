#!/bin/bash
#SBATCH -J pack_fault           # Job name
#SBATCH -o pack_fault_%j.out    # stdout (%j = job ID)
#SBATCH -e pack_fault_%j.err    # stderr
#SBATCH -p small                # Queue (single-node packing)
#SBATCH -N 1                    # 1 node
#SBATCH -n 1                    # 1 task (Python is single-process)
#SBATCH -t 12:00:00             # wall time — Lustre metadata is the bottleneck
#SBATCH -A EAR20006             # TACC allocation
#
# Phase 5 of miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md.
#
# Pack a legacy per-rank-per-cycle FaultSurface/ directory into chunked
# tar.zst archives suitable for fast `scp` to a workstation.  Designed
# to be run on a Frontera compute node (NOT a login node — the 30-min
# CPU limit will kill the script before it finishes a large run).
#
# Usage:
#   sbatch -A YOUR_ALLOCATION pack_fault.sbatch \
#          /scratch2/.../job_xxx/FaultSurface \
#          /scratch2/.../packs/job_xxx_fault
#
# The first arg is the FaultSurface directory; the second is the
# output basename — chunks land at ${BASENAME}_chunk0001.tar.zst etc.

set -euo pipefail

FAULT_DIR="${1:?usage: sbatch ... pack_fault.sbatch FAULT_DIR OUT_BASENAME}"
OUT_BASENAME="${2:?usage: sbatch ... pack_fault.sbatch FAULT_DIR OUT_BASENAME}"
CHUNK_SIZE="${3:-100}"

export LC_ALL=C
export LANG=C

module load python3 || true   # TACC default; harmless if absent

# Verify the script is reachable.
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PACK_PY="${SCRIPT_DIR}/pack_fault_output.py"
if [ ! -f "${PACK_PY}" ]; then
    echo "error: pack_fault_output.py not found at ${PACK_PY}" >&2
    exit 1
fi

# Make the output directory if it does not exist.
mkdir -p "$(dirname "${OUT_BASENAME}")"

echo "Packing fault dir: ${FAULT_DIR}"
echo "Output basename : ${OUT_BASENAME}"
echo "Chunk size      : ${CHUNK_SIZE} cycles per tar"

python3 "${PACK_PY}" \
    "${FAULT_DIR}" \
    --out-tar "${OUT_BASENAME}" \
    --chunk-size "${CHUNK_SIZE}"

echo "done.  inspect chunks under $(dirname "${OUT_BASENAME}"):"
ls -lah "$(dirname "${OUT_BASENAME}")"/$(basename "${OUT_BASENAME}")_chunk*.tar.* \
    2>/dev/null | head -5 || true
