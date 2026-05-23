#!/bin/bash
#SBATCH -J safs_zerodip
#SBATCH -o safs_zerodip_%j.out
#SBATCH -e safs_zerodip_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# DIP-PRESTRESS-ZERO diagnostic for the SAFS blow-up.
# =========================================================================
# Hypothesis under test: the ~32deg dip of the meshed SAFS fault projects
# the regional H&Z stress into a LARGE dip-slip prestress (tau_dip up to
# -34 MPa, comparable to tau_strike ~26 MPa; rake spans -178..+177 deg),
# and that dip-slip drive (via slip<->normal-stress coupling) is what makes
# the rupture blow up.  This run ZEROES the dip pre-stress component
# (tau1_0 -> 0, keeping strike tau2_0 and normal sigma_n0) at the
# projection source (SEAS_ZERO_DIP_PRESTRESS), so BOTH the [derived] gate
# and the dynamics see the zeroed dip.  UNPHYSICAL (discards real stress)
# -- diagnostic only.  Does NOT touch the friction sigma_n term (the
# sigma_n-freeze investigation was reverted; this test is independent).
#
# Read the result (control blows up at t~0.95 s):
#   bounded AND still rupturing (DIAG n_rupturing keeps growing)
#       -> the dip-slip prestress is the driver of the instability.
#   blows up anyway
#       -> the initial dip prestress is NOT the cause (dynamic dip /
#          something else); pivot.
#   bounded but n_rupturing flat (no rupture)
#       -> false stabilization (zeroing dip also killed the event);
#          inconclusive -- the [derived] gate / nucleation needs a relook.
#
# Knobs (override via sbatch --export=ALL,KEY=VAL,...):
#   SAFS_ZERO_DIP=1   (default) zero the dip prestress  | 0 = CONTROL (reproduce blow-up)
#   SAFS_CONFIG=...   default = main projected-stress config (mu_s=0.85, D_c=1.0)
#   SAFS_CFL=0.5      (default) faster than the config's 0.25; control == job 7745103
#   SAFS_TFINAL=2s    (default) past the t~0.95 onset with margin
#   SAFS_MESH=...     default = [mesh].path from config (500 m z0embed)

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

# Resolve seas-mfem root by walking UP until miniapps/seas/Makefile is found.
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR}}"
_walk_up_dir="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5; do
    if [[ -f "${_walk_up_dir}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"
        break
    fi
    _walk_up_dir="${_walk_up_dir}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}"
    exit 1
fi
cd "${SEAS_MFEM_ROOT}"

# ZFP HDF5 filter plugin (libh5zzfp.so), required at runtime for the
# --paraview-fault-zfp-tol below.  build_frontera.sh installs it under
# ${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin (same convention as the
# smoke / resolution sbatches).  Fail fast with a build hint if missing.
H5Z_ZFP_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install"
if [[ -e "${H5Z_ZFP_BUILT_DIR}/plugin/libh5zzfp.so" ]]; then
    export HDF5_PLUGIN_PATH="${H5Z_ZFP_BUILT_DIR}/plugin:${HDF5_PLUGIN_PATH:-}"
else
    echo "ERROR: ZFP plugin not found at"
    echo "         ${H5Z_ZFP_BUILT_DIR}/plugin/libh5zzfp.so"
    echo "       Build it on Frontera with:  bash build_frontera.sh"
    echo "       (or drop --paraview-fault-zfp-tol from the driver invocation)."
    exit 1
fi

grep -q '^MFEM_USE_PETSC *= YES'   config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'   config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'    config/config.mk || { echo "ERROR: no HDF5"; exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: no H5Z_ZFP"; exit 1; }

cd miniapps/seas
# The SEAS_ZERO_DIP_PRESTRESS toggle lives in the header-only
# fault/fault_geometry_safs_templated.inl, which is NOT tracked as a make
# dependency of spatial_dyn_driver.o.  After a `git pull` that only touches
# the .inl, `make` would relink a STALE driver .o without the toggle.  Force
# a recompile of the driver TU so the toggle is actually built in.
rm -f drivers/spatial_dyn_driver.o
make seas_spatial_dyn_driver -j8 || { echo "BUILD FAILED"; exit 1; }

SAFS_ROOT=safs/project_7.0_alternative
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_slip_weakening_safs_projected_stress.toml}"

if [[ -n "${SAFS_MESH:-}" ]]; then
    MESH="${SAFS_MESH}"
else
    MESH=$(sed -n '/^\[mesh\]/,/^\[/p' "${CONFIG_TOML}" \
            | grep -E '^[[:space:]]*path[[:space:]]*=' | head -1 \
            | sed -E 's/^[^=]*=[[:space:]]*"?([^"#]*)"?.*/\1/; s/[[:space:]]+$//')
fi
if [[ -z "${MESH}" || ! -f "${MESH}" ]]; then
    echo "ERROR: SAFS mesh not found (resolved '${MESH}' from ${CONFIG_TOML})."
    echo "       The .msh is gitignored; scp it to the cluster first."
    exit 1
fi

# --- the diagnostic toggle this job exists for ---
ZERO_DIP="${SAFS_ZERO_DIP:-1}"
TAG="control"
if [[ "${ZERO_DIP}" == "1" ]]; then
    export SEAS_ZERO_DIP_PRESTRESS=1   # presence enables it in the driver
    TAG="zerodip"
fi
# rupture-area + onset diagnostic (print-only; confirms the run still ruptures)
export SEAS_DIAG_BLOWUP=1
# Downgrade the R-101 shared-fault DOFData consistency tripwire from a hard
# MFEM_ABORT to a per-check rank-0 diagnostic, so this run proceeds past the
# known curvilinear-fault cross-rank divergence and we can read whether the
# fault-normal canonicalization change (fault_basis.hpp) shrinks worst_rel and
# keeps V_max bounded.  Look for "[R-101 nonfatal] ... worst_rel=... (DIVERGED|ok)".
export SEAS_R101_NONFATAL=1

TFINAL="${SAFS_TFINAL:-2s}"
CFL="${SAFS_CFL:-0.5}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_${TAG}_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_${TAG}_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}"
echo "Config:      ${CONFIG_TOML}"
echo "Zero dip:    ${ZERO_DIP}  (SEAS_ZERO_DIP_PRESTRESS=${SEAS_ZERO_DIP_PRESTRESS:-<unset>})"
echo "R-101 gate:  NON-FATAL (SEAS_R101_NONFATAL=${SEAS_R101_NONFATAL:-<unset>}) -- run proceeds past cross-rank divergence"
echo "ParaView:    ${OUT}/fault.vtkhdf  (hdf5, zfp 1e-12, <=200 snapshots)"
echo "cfl:         ${CFL}    tfinal: ${TFINAL}    (mixed flux = config default = adjacent)"

# ParaView fault.vtkhdf output (matches the smoke / resolution sbatches):
# single .vtkhdf, ZFP-compressed at 1e-12 (slip-rate dynamic-range floor),
# capped at 200 snapshots over the run.  --print-derived still runs the
# well-posedness gate (aborts early on a mis-config, e.g. if zeroing dip
# breaks the nucleation stress drop); V_max + [DIAG] + [R-101 nonfatal] lines
# go to the log.
ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --no-sidecar-material \
    --cfl "${CFL}" \
    --tfinal "${TFINAL}" \
    --paraview \
    --paraview-fault-hdf5 \
    --paraview-fault-zfp-tol 1.0e-12 \
    --paraview-max-snapshots 200 \
    --print-derived \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?
echo "driver_rc=${RC}"
echo "=== V_max + DIAG tail ==="
grep -E "^step |\[DIAG\]|\[DIAG-ONSET\]|\[R-101 nonfatal\]|PASS|FAIL|done\." "${LOG}" | tail -30
echo "=== fault.vtkhdf ==="
ls -la "${OUT}/fault.vtkhdf" 2>/dev/null \
    || echo "  (no fault.vtkhdf written -- run aborted before the first snapshot?)"
exit ${RC}
