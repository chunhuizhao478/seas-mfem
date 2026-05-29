#!/bin/bash
#SBATCH -J safs_bigbox
#SBATCH -o safs_bigbox_%j.out
#SBATCH -e safs_bigbox_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 12:00:00
#SBATCH -A EAR20006

# =========================================================================
# Phase 12.3 metric 3 — big-box reference (the gold standard), NO PML.
#
# Runs the SAFS LSW problem on a DEEPER / WIDER mesh whose walls are far enough
# that no reflection contaminates the rupture over [0, ~14 s] (t_reflect grows
# to ~14-17 s).  No PML — over its reflection-free window this is the truth the
# PML run is judged against:
#   pml_reflection_metrics.py --bigbox <pml_on.log> <bigbox.log> --window 14
# Accept the PML if max |V_max^PML - V_max^bigbox| / V_max^bigbox < 5%.
#
# REQUIRES the padded mesh (plan §12.3 req 3 / §12.4):
#   run_z0cut_meshing.py --pad-x 100000 --pad-y 100000 --pad-bottom 60000
# then scp the resulting .msh in and pass it via SAFS_MESH (REQUIRED here — the
# default triq mesh is NOT big enough to be a reflection-free reference).
#
# Overrides: SAFS_MESH (REQUIRED — the big-box .msh), SAFS_CONFIG,
# SAFS_TFINAL (default 16s — a bit past the ~14 s reflection-free window),
# SAFS_CFL, SAFS_RESTART.
# =========================================================================

set -o pipefail
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR}}"
_d="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5; do
    [[ -f "${_d}/miniapps/seas/jobs/safs/_pml_common.sh" ]] && \
        { SEAS_MFEM_ROOT="$(cd "${_d}" && pwd)"; break; }
    _d="${_d}/.."
done
export SEAS_MFEM_ROOT

# Guard: the big-box reference is meaningless on the small default mesh.
if [[ -z "${SAFS_MESH:-}" ]]; then
    echo "ERROR: SAFS_MESH is REQUIRED for the big-box reference — it must point"
    echo "       at the padded mesh (run_z0cut_meshing.py --pad-x/--pad-y 100000"
    echo "       --pad-bottom 60000).  The default triq mesh is too small to be"
    echo "       reflection-free; refusing to run a misleading 'reference'."
    exit 1
fi

source "${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs/_pml_common.sh"

PML_MODE=bigbox

run_pml_driver bigbox_ref \
    --tfinal "${SAFS_TFINAL:-16s}"
exit $?
