#!/bin/bash
#SBATCH -J safs_pml_off
#SBATCH -o safs_pml_off_%j.out
#SBATCH -e safs_pml_off_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# =========================================================================
# Phase 12.3 metric 1 — SAFS LSW reflection study, PML OFF (the A/B baseline).
#
# IDENTICAL config / mesh / tfinal to spatial_dyn_pml_on_safs.sbatch but with
# NO PML (no --pml).  This is the reflection-contaminated control; its V_max(t)
# is the "OFF" series in:
#   pml_reflection_metrics.py --ab <pml_on.log> <pml_off.log>
# Because use_pml stays false the driver constructs no PMLLayer, so this run is
# byte-identical to the pre-Phase-12 binary on the same inputs.
#
# Overrides: SAFS_CONFIG, SAFS_MESH, SAFS_TFINAL (default 40s), SAFS_CFL,
# SAFS_RESTART.
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
source "${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs/_pml_common.sh"

PML_MODE=off

run_pml_driver pml_off \
    --tfinal "${SAFS_TFINAL:-40s}"
exit $?
