#!/bin/bash
#SBATCH -J safs_resDc2_cap
#SBATCH -o safs_resDc2_cap_%j.out
#SBATCH -e safs_resDc2_cap_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# NORMAL-STRESS-CAP (sigma_n strength floor) ARREST TEST — LSW Dc2.
# =========================================================================
# Plan: spatial_dynamic_rupture_sliver_blowup_2026-05-25.md §Phase 4 acceptance.
#
# Question: job 7748818 (LSW Dc2, mixed_flux=adjacent, NO floor) ran a clean
# ~10-12 m/s rupture for ~25 s, then V_max departed the plateau at t~28 s and ran
# away super-exponentially to 5e64 by t~64.6 s.  §2c confirmed the mechanism:
# the sigma_n oscillation crosses TENSILE -> the fault opens -> the LSW strength
# collapses to 0 (max(sigma_n,0)) -> free slip -> positive feedback.
#
# THIS RUN: identical mesh / config / stress / nucleation / ranks / FLUX as
# 7748818 (mixed_flux=adjacent — the config default, NOT passed on the CLI).
# The ONLY change is the **sigma_n strength floor = 10 MPa**, now baked into
#   config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml
#   ([friction].sigma_n_strength_floor_pa = 10.0e6).
# With the floor, a tensile sigma_n no longer drops the strength to 0; the LSW
# strength saturates at mu_eff*10 MPa, so V_abs = (|tau| - mu_eff*10MPa)/eta_s
# stays bounded instead of free-sliding.
#
# ACCEPTANCE (plan §Phase 4):
#   - PASS  : V_max stays bounded (e.g. < 1e3 m/s) through t >= 28-35 s — the
#             free-slip runaway is broken by the floor ALONE (adjacent flux).
#   - FAIL  : V_max still departs at t~28 s — the floor alone is insufficient;
#             the central-flux/sliver sigma_n-OSCILLATION source (not the
#             free-slip) drives the blow-up.  Cross-check against the
#             mixed_flux=none upwind test (..._resDc2_upwind_...sbatch): run the
#             two INDEPENDENTLY so we learn whether floor alone, flux alone, or
#             both are needed.
#
# From t=0 (NOT a restart) on purpose: a restart would carry pre-floor history.
# SEAS_DIAG_BLOWUP=1 so the [DIAG]/[DIAG-ONSET] V_max trace is directly
# comparable to the 7748818 log.
#
# REFLECTION CAVEAT: tfinal 100 s > min_box_dim/cp_max (~6.94 s); the 1st-order
# ABC leaks after ~7 s.  Reflections are present in BOTH this run and 7748818,
# so they cancel as a variable — this isolates the FLOOR effect.  --pml is NOT
# wired (flag only silences the warning).
#
# Overrides: SAFS_CONFIG, SAFS_MESH, SAFS_MIXED_FLUX (default = config adjacent),
# SAFS_TFINAL (default 100s), SAFS_CFL, SAFS_RESTART, SEAS_DIAG_*.

export LC_ALL=C
export LANG=C

# DIAGNOSTICS ON: compare V_max(t) to the 7748818 (no-floor) blow-up log.
export SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-0}"
export SEAS_DIAG_SLIP_VTHR="${SEAS_DIAG_SLIP_VTHR:-10.0}"
# Leave SEAS_NOOPENING OFF — the floor is the mechanism under test, not the
# diagnostic no-opening cap (they are independent; both off here except floor).
export SEAS_NOOPENING="${SEAS_NOOPENING:-0}"
export SEAS_DIAG_FRAME="${SEAS_DIAG_FRAME:-0}"
export SEAS_DIAG_FRAME_XYZ="${SEAS_DIAG_FRAME_XYZ:-365100,3838920,-2488}"
export SEAS_DIAG_FRAME_R="${SEAS_DIAG_FRAME_R:-300}"
export SEAS_DIAG_MACRO="${SEAS_DIAG_MACRO:-0}"
export SEAS_DIAG_BLOWUP="${SEAS_DIAG_BLOWUP:-1}"

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
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR="
    echo "  ${SLURM_SUBMIT_DIR}.  Set SEAS_MFEM_ROOT explicitly via"
    echo "  sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi
cd "${SEAS_MFEM_ROOT}"

# ZFP HDF5 filter plugin (libh5zzfp.so), required at runtime for
# --paraview-*-zfp-tol.
H5Z_ZFP_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install"
if [[ -e "${H5Z_ZFP_BUILT_DIR}/plugin/libh5zzfp.so" ]]; then
    export HDF5_PLUGIN_PATH="${H5Z_ZFP_BUILT_DIR}/plugin:${HDF5_PLUGIN_PATH:-}"
else
    echo "ERROR: ZFP plugin not found at"
    echo "         ${H5Z_ZFP_BUILT_DIR}/plugin/libh5zzfp.so"
    echo "       Build it on Frontera with:  bash build_frontera.sh"
    echo "       (or omit --paraview-*-zfp-tol from the driver invocation)."
    exit 1
fi

grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'   config/config.mk || { echo "ERROR: no HDF5"; exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: no H5Z_ZFP"; exit 1; }

cd miniapps/seas
# This run job does NOT build — the binary must be pre-built (by
# spatial_dyn_build_safs.sbatch, or a login/idev-node `make`).  Building here
# would let concurrently-running jobs race on the same relink (ETXTBSY /
# truncated binary).  Fail fast if it is missing.
if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not found/executable in $(pwd)."
    echo "       Build it FIRST (no in-job build, to avoid the relink race):"
    echo "         sbatch miniapps/seas/jobs/safs/spatial_dyn_build_safs.sbatch"
    echo "       or on a login/idev node:  cd miniapps/seas && make seas_spatial_dyn_driver"
    exit 1
fi

SAFS_ROOT=safs/project_7.0_alternative
# Same config as the 7748818 blow-up (D_c=2.0 / radius=4000 / cfl=0.5) — NOW with
# [friction].sigma_n_strength_floor_pa = 10.0e6.
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml}"

# PRE-FLIGHT: this test is meaningless without the floor — fail loudly if the
# config does not enable it (catches an accidentally-reverted/wrong config).
if ! grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null; then
    echo "ERROR: ${CONFIG_TOML} does not set sigma_n_strength_floor_pa — this is"
    echo "       the NORMAL-STRESS-CAP arrest test and REQUIRES the floor enabled."
    echo "       (It must be a TOP-LEVEL [friction] key, e.g. 10.0e6.)"
    exit 1
fi

# Mesh: read [mesh].path from the config (500 m z0embed); override with SAFS_MESH.
if [[ -n "${SAFS_MESH:-}" ]]; then
    MESH="${SAFS_MESH}"
else
    MESH=$(sed -n '/^\[mesh\]/,/^\[/p' "${CONFIG_TOML}" \
            | grep -E '^[[:space:]]*path[[:space:]]*=' | head -1 \
            | sed -E 's/^[^=]*=[[:space:]]*"?([^"#]*)"?.*/\1/; s/[[:space:]]+$//')
fi
if [[ -z "${MESH}" ]]; then
    echo "ERROR: could not resolve [mesh].path from ${CONFIG_TOML}"
    exit 1
fi
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: SAFS mesh not found:"
    echo "         $(pwd)/${MESH}"
    echo "       The .msh is gitignored and is NOT pulled with the repo."
    echo "       scp it to the path above before submitting this job."
    exit 1
fi

# Flux: default = config's mixed_flux (adjacent) by NOT passing --mixed-flux —
# this is the FLOOR-ONLY experiment (original 7748818 flux + floor).  Override
# with SAFS_MIXED_FLUX={none|adjacent|all_continuous}.
MIXFLUX_ARG=""
if [[ -n "${SAFS_MIXED_FLUX:-}" ]]; then
    case "${SAFS_MIXED_FLUX}" in
        none|adjacent|all_continuous) : ;;
        *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid (none|adjacent|all_continuous)"; exit 1 ;;
    esac
    MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"
fi

# tfinal default 100 s: the full original 7748818 run length (well past the
# t~28 s no-floor onset and the t~64 s runaway).  At ~48 h / 35 s throughput the
# 48 h wall will NOT reach 100 s in one allocation — the job checkpoints
# (--checkpoint-every) and is resumed via SAFS_RESTART.  (Reflection caveat:
# physical only to ~7 s; the V_max(t) trace past t~7 s is still the
# apples-to-apples comparison to 7748818, which had the same reflections.)
TFINAL="${SAFS_TFINAL:-100s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Restart hook (unused for this from-t=0 test; here for parity / 48 h overflow).
RESTART_ARG=""
[[ -n "${SAFS_RESTART:-}" ]] && RESTART_ARG="--restart ${SAFS_RESTART}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_resDc2_cap_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_resDc2_cap_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Queue/size:  normal, 10N x 500r, 48 h"
echo "EXPERIMENT:  sigma_n strength floor (10 MPa) arrest test — FLOOR ONLY"
echo "             (mixed_flux = config default 'adjacent', as in 7748818)"
echo "Compare to:  job 7748818 (no floor, adjacent) blow-up at t~28-64 s"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}"
echo "Config:      ${CONFIG_TOML}   (D_c=2.0 / radius=4000 / floor=10 MPa)"
echo "Diag:        SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP} (V_max/[DIAG-ONSET] trace ON)"
echo "Mixed flux:  ${SAFS_MIXED_FLUX:-<config default = adjacent>}"
echo "tfinal:      ${TFINAL}"
echo "CFL arg:     ${CFL_ARG:-<config default = 0.5>}"

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --no-sidecar-material \
    ${MIXFLUX_ARG} \
    --tfinal "${TFINAL}" \
    ${CFL_ARG} \
    ${RESTART_ARG} \
    --paraview \
    --paraview-fault-hdf5 \
    --paraview-fault-zfp-tol 1.0e-12 \
    --paraview-max-snapshots 200 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

# Confirm the floor was actually active (banner echo) — guards against a build
# that predates the floor or a config that lost it.
if grep -qE 'sigma_n strength floor:[[:space:]]*10[[:space:]]*MPa' "${LOG}"; then
    echo "floor: ACTIVE (banner reports 10 MPa)"
else
    echo "WARNING: did not find 'sigma_n strength floor: 10 MPa' in the banner —"
    echo "         the floor may be DISABLED.  Inspect ${LOG} before trusting the result."
fi

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG} (probably --print-derived ABORTed; "
    echo "  check the [derived] lines)."
    exit ${RC}
fi

echo "driver_rc=${RC} — done.  Inspect V_max(t) in ${LOG}:"
echo "  grep -E '^\\[DIAG\\] step' ${LOG} | awk 'NR%200==1{print \$4,\$5}'"
echo "  PASS (V_max bounded, e.g. <1e3, through t>=28-35 s) -> floor ALONE arrests the runaway."
echo "  FAIL (departs plateau at t~28 s)                    -> floor insufficient; the"
echo "                                                          sigma_n-oscillation source remains"
echo "                                                          (compare the mixed_flux=none upwind test)."
exit ${RC}
