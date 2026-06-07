#!/bin/bash
#SBATCH -J tpv31_p2_pu_rvsm_split
#SBATCH -o tpv31_p2_pu_rvsm_split_%j.out
#SBATCH -e tpv31_p2_pu_rvsm_split_%j.err
#SBATCH -p normal
#SBATCH -N 16
#SBATCH -n 800
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# SCEC TPV31 RVSMOOTH-SPLIT-MESH control — p2 / ADER-O3 / pure upwind / matrix
#   Riemann, 50 m, NO over-integration, NO resample, 0-15 s.
# Purpose: isolate the MESH effect on the TPV31 LSW blow-up at p2 — the SECOND
# near-fault-asymmetry cure (alternative to the y-mirror hybrid).
# tpv31_50m_rvsmooth_split.msh has the near-fault tetrahedra's volume ratio r_v
# (Eq.18) driven < 1.5 everywhere by node smoothing + a conforming 1->4 tet split
# (smooth_fault_volume_ratio.py -> retriangulate_fault_rv.py).  Compare its
# on-fault sigma_n against the UNSTRUCTURED tpv31_50m.msh (and the hybrid mesh):
# if regularizing the near-fault tet shapes alone tames the drift/blow-up, the
# speckle is mesh-asymmetry-driven (the upwind-flux / Zhang-2023 channel).
# Pairs with the p1 rvsmooth_split control.
#   * p2            : [mesh].order = 2            (config-only; asserted below).
#   * ADER-O3       : --ader-order 3 --time-integrator ader.
#   * pure upwind   : --mixed-flux none --fault-overint 0  (no central flux, no
#                     over-integration; --fault-resample NOT passed).
#   * matrix Riemann: [numerics].interior_flux = "matrix"  (MANDATORY for TPV31's
#                     depth_profile_1d bi-material; config-only, asserted below).
#   Config: tpv31/configs/tpv31_p2.toml.
#
# CROSS-RANK SEAM (np>1, branch fix-cross-rank-material-exchange) ----------
#   This job runs the bi-material Godunov UPWIND seam (mixed_flux=none) on the
#   depth_profile_1d material.  At np>1 ParMETIS can cut the depth profile across
#   a partition seam; the cross-rank material exchange now reads the TRUE peer
#   element's material there (replacing the old local-side stub).  At a seam that
#   crosses a depth discontinuity (2400/5000/10000 m) the upwind seam flux is
#   RE-BASELINED to the correct 2-sided value (== serial / a seam-vertical
#   partition); byte-exact only for partitions that do NOT cut across a depth
#   jump.  ⇒ re-baseline the Frontera gold for this job against a fresh np>1 run
#   on this branch (only cross-depth seams change; the bulk interior is identical).
#
# MESH (no .geo — built by the rv pipeline; stage it or rebuild on a login node):
#   conda activate pythonenv
#   cd miniapps/seas/safs/project_7.0_alternative/experimental_mesh_refinement
#   python smooth_fault_volume_ratio.py --in <repo>/miniapps/seas/tpv31/mesh/tpv31_50m.msh \
#     --out <...>/tpv31/mesh/tpv31_50m_rvsmooth.msh \
#     --fault-phys 101 --rock-phys 1 --boundary-phys 102,103,104 --slide-boundary --rv-ceiling 1.5 --best-effort
#   python retriangulate_fault_rv.py --in <...>/tpv31/mesh/tpv31_50m_rvsmooth.msh \
#     --out <...>/tpv31/mesh/tpv31_50m_rvsmooth_split.msh \
#     --fault-phys 101 --rock-phys 1 --boundary-phys 102,103,104 --rv-ceiling 1.5 --target-rv 1.0
#
# Standalone: build seas_spatial_dyn_driver first (from THIS branch).
# NOTE: p2 dt is smaller than p1; even 48 h won't reach 15 s — checkpoints every
# 2000 steps; resubmit with TPV31_RESTART=<out>/cp to continue.
set -euo pipefail

FAULT_OVERINT=0                                      # no over-integration
MESH="tpv31/mesh/tpv31_50m_rvsmooth_split.msh"
TFINAL="${TPV31_TFINAL:-15}"                          # spec p. 9 = 15 s
CONFIG_TOML="tpv31/configs/tpv31_p2.toml"             # p2, interior_flux="matrix"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
# Pin hypre/2.31.0 at the FRONT (petsc/3.15 + mumps/5.3 prepend a pre-2.31
# libHYPRE that shadows HYPRE_Initialize otherwise — see the noflux TPV31 jobs).
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH:-}"

_d="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "${BASH_SOURCE[0]:-$0}")")}}"
while [[ "${_d}" != "/" && ! -d "${_d}/miniapps/seas" ]]; do _d="$(dirname "${_d}")"; done
cd "${_d}/miniapps/seas"

if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not built.  Build it first:"
    echo "       sbatch miniapps/seas/jobs/spatial_dyn_build.sbatch clean-rebuild"
    exit 1
fi
[[ -f "${CONFIG_TOML}" ]] || { echo "ERROR: ${CONFIG_TOML} not found."; exit 1; }
grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"matrix"' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not interior_flux=\"matrix\" (TPV31 needs the matrix Riemann path)."; exit 4; }
grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*2([[:space:]]|$)' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not [mesh].order=2 (this job is p2)."; exit 4; }

# rvsmooth_split has NO .geo (built by the rv-smooth+split pipeline, ~350 MB).
# It must be staged (scp) or rebuilt on a login node — see the header recipe.
[[ -f "${MESH}" ]] || { echo "ERROR: mesh ${MESH} not found; no .geo to regenerate from.  Stage it (scp) or rebuild via the rv-smooth+split pipeline (see this script's header)."; exit 2; }

OUT="tpv31/out_p2_aderO3_pu_50m_rvsmooth_split_${SLURM_JOB_ID:-local}"
echo "[rvsm-split] TPV31 p2 pure-upwind on RVSMOOTH-SPLIT mesh  mesh=${MESH}  tfinal=${TFINAL}  out=${OUT}"
echo "[rvsm-split] Git: $(git rev-parse --abbrev-ref HEAD) $(git rev-parse --short HEAD)"

RESTART_FLAGS=""
if [[ -n "${TPV31_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV31_RESTART}"
    echo "[rvsm-split] RESTART from ${TPV31_RESTART} (this segment writes ${OUT})."
fi

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --ader-order 3 \
    --time-integrator ader \
    --mixed-flux none \
    --fault-overint ${FAULT_OVERINT} \
    --tfinal "${TFINAL}" \
    --paraview \
    --paraview-fault-vtu \
    --paraview-fault-dt 0.1 \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    --deriv-cache --shared-ck-recursion \
    ${RESTART_FLAGS} \
    --verify-dispatch

echo "[rvsm-split] done (TPV31 p2 pure-upwind, rvsmooth_split mesh).  station traces: ${OUT}/"
echo "[rvsm-split] compare vs unstructured tpv31_50m + hybrid + scec_{eqdyna,seisol} overlays."
