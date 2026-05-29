#!/bin/bash
#SBATCH -J safs_pml_on
#SBATCH -o safs_pml_on_%j.out
#SBATCH -e safs_pml_on_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# =========================================================================
# Phase 12.3 metric 1/2 — SAFS LSW reflection study, PML ON.
#
# Runs the SAFS Dc2 LSW problem with the convolutional PML ACTIVE on the five
# non-free-surface walls (x±, y±, z-min; the z=0 free surface stays undamped).
# Pair with spatial_dyn_pml_off_safs.sbatch (identical config, no PML) and feed
# both V_max(t) logs to:
#   pml_reflection_metrics.py --ab <pml_on.log> <pml_off.log>
# Decision rule (plan §12.3 req 1): if the t~28s V_max departure disappears /
# strongly delays WITH PML it was reflection-driven; if it persists it is an
# intrinsic friction/LSW instability.
#
# Thickness defaults to pml_cells(4) * [mesh].lc_far (3000 m for the triq mesh)
# = 12 km (plan §12.0), keeping a 13 km bottom clearance.  Override with
# SAFS_PML_THICKNESS / SAFS_PML_CELLS / SAFS_PML_TARGET_R.  Requires
# [mesh].lc_far_m to be set in the config (or pass SAFS_PML_THICKNESS).
#
# Overrides: SAFS_CONFIG, SAFS_MESH, SAFS_TFINAL (default 40s), SAFS_CFL,
# SAFS_RESTART, SAFS_PML_CELLS, SAFS_PML_TARGET_R, SAFS_PML_THICKNESS.
# =========================================================================

set -o pipefail
# Resolve seas-mfem root, then source the shared PML runner by absolute path
# (sbatch spools the script, so $0/BASH_SOURCE dirname is unreliable).
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR}}"
_d="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5; do
    [[ -f "${_d}/miniapps/seas/jobs/safs/_pml_common.sh" ]] && \
        { SEAS_MFEM_ROOT="$(cd "${_d}" && pwd)"; break; }
    _d="${_d}/.."
done
export SEAS_MFEM_ROOT
source "${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs/_pml_common.sh"

PML_MODE=on

PML_THICK_ARG=""
[[ -n "${SAFS_PML_THICKNESS:-}" ]] && PML_THICK_ARG="--pml-thickness ${SAFS_PML_THICKNESS}"

run_pml_driver pml_on \
    --pml \
    --pml-cells    "${SAFS_PML_CELLS:-4}" \
    --pml-target-R "${SAFS_PML_TARGET_R:-1e-3}" \
    --pml-damp-bottom 1 \
    --pml-damp-top    0 \
    ${PML_THICK_ARG} \
    --tfinal "${SAFS_TFINAL:-40s}"
exit $?
