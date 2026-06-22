#!/bin/bash
# Transfer the gitignored SAFS MFEM sidecars (mesh + stress/velocity h5) to FRONTERA.
# These ~1.6 GB of files are gitignored/large (NOT in the repo), so the safs_frontera
# run jobs abort at their pre-flight unless they are present in the Frontera checkout.
# rsync -R preserves each relative path so the files land in the right subdir under
# the remote repo clone.  Run on your LOCAL machine (Mac):
#   bash miniapps/seas/jobs/safs/safs_frontera/scp_safs_sidecars_frontera.sh
set -euo pipefail

# Silence the macOS-ssh-forwards-LANG / Frontera-Lmod-has-no-C.UTF-8 locale spam.
export LC_ALL=C LANG=C

# Resolve the local repo root from THIS script's location (it lives at
# miniapps/seas/jobs/safs/safs_frontera/ — five levels under the repo root).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../../.." && pwd)"
if [[ ! -d "${REPO_ROOT}/miniapps/seas" ]]; then
    echo "ERROR: cannot locate repo root from ${SCRIPT_DIR} (expected .../miniapps/seas above it)." >&2
    exit 1
fi
cd "${REPO_ROOT}"

# Remote repo clone root on Frontera (the safs_frontera run jobs' SEAS_MFEM_ROOT).
# Edit if your checkout path changes; this must match the run jobs' checkout AND
# the sync_results.sh REMOTE_ROOT.
DEST='zhaochun@frontera.tacc.utexas.edu:/scratch2/10024/zhaochun/seas-project/seas-mfem-spatial-dyn-driver/'

FILES=(
  miniapps/seas/safs/project_7.0_preferred/meshing/results/safv4_deep_500m_opt_safstags.msh
  miniapps/seas/safs/project_7.0_preferred/stress/results/csm_d12_lsw_mfem/safs_stress_csm_d12_mfem.h5
  miniapps/seas/safs/project_7.0_preferred/stress/results/csm_d12_rssrw_mfem/safs_stress_csm_d12_mfem.h5
  miniapps/seas/safs/project_7.0_preferred/velocity/results/multiscale_statewise_cvm/velocity_safs.h5
  miniapps/seas/safs/project_7.0_alternative/experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh
  miniapps/seas/safs/project_7.0_alternative/velocity/results/multiscale_statewise_cvm/velocity_safs.h5
  miniapps/seas/safs/project_7.0_alternative/stress/results/csm_d12_mfem/safs_stress_csm_d12_mfem.h5
)

# Pre-flight: every file must exist locally (they are gitignored, so a fresh
# clone will NOT have them — regenerate via the toolbox scripts first if missing).
missing=0
for f in "${FILES[@]}"; do
    [[ -f "$f" ]] || { echo "MISSING (regenerate first): $f" >&2; missing=1; }
done
[[ "$missing" -eq 0 ]] || { echo "ERROR: one or more sidecars are missing locally." >&2; exit 1; }

echo "Local root:  ${REPO_ROOT}"
echo "Transferring ${#FILES[@]} files ($(du -ch "${FILES[@]}" | tail -1 | cut -f1)) to:"
echo "  ${DEST}"
echo
rsync -avR --progress "${FILES[@]}" "${DEST}"
echo
echo "DONE.  The safs_frontera run jobs' sidecar pre-flight will now pass."
