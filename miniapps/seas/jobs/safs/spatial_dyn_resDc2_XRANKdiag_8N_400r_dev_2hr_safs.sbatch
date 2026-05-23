#!/bin/bash
#SBATCH -J safs_resDc2diag
#SBATCH -o safs_resDc2diag_%j.out
#SBATCH -e safs_resDc2diag_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# OBLIQUE (DIP-LOADED) CROSS-RANK BASELINE + Phase-0 seed trace (R-DIP2).
# =========================================================================
# IDENTICAL simulation to spatial_dyn_resolution_Dc2_8N_400r_dev_2hr_safs.sbatch
# (same config / mesh / driver flags — full oblique H&Z prestress, so dip
# traction is PRESENT, unlike spatial_dyn_zerodip_*.sbatch) but with TWO
# runtime env changes:
#
#   (1) SEAS_R101_NONFATAL=1 — DROP the hard R-101 gate.  The hard gate used to
#       abort this run at t~0.455 s on the cross-rank V1 (dip) divergence.  We
#       now downgrade it to a per-check rank-0 log so the run PROCEEDS PAST the
#       divergence.  Two things we want to see, that the abort hid:
#         (a) does the DIP-loaded case stay BOUNDED (V_max peaks then decays,
#             like the strike-only zerodip run, job 7747119) or BLOW UP
#             (V_max -> 1e3, sigma_n -> tension; the dip<->normal runaway)?
#             This is the direct test of PLAN B.2 "the blow-up needs dip."
#         (b) the SPECKLE — the scattered slip-rate elements outside the
#             rupture (the visible cross-rank symptom).  This run REPRODUCES
#             the speckle (ParaView fault.vtkhdf below) as the PRE-FIX
#             BASELINE; the post-reconcile rerun must CLEAR it.  That is the
#             final check on the fix (PLAN Phase 4).
#
#   (2) SEAS_DIAG_XRANK=1 — the env-gated cross-rank trace (Phase 0: localise
#       the ~1e-14 seed).  Default-off in code => only this job pays the cost.
#
# The script extracts the [XRANK] trace and splits it by the (two) ranks that
# own the diverging QP so you can diff them directly:
#
#   diff <per-rank-A>.log <per-rank-B>.log   (or paste/awk step-by-step)
#
# How to read the diff (Phase 0 seed localisation, %.17e full precision):
#   RAWDOF side=+ identical across ranks, but I*_can differ  => source (1)
#                       interpolation (shape1.Loc1 vs shape2.Loc2; inherent)
#   RAWDOF side=+ differ already at the first traced sub-step => source (2)
#                       ghost-exchange not bit-exact (fixable at the source)
#   RAWDOF identical early, then drift apart                  => source (3)
#                       accumulated bulk divergence (inherent)
#   can_t1 differs              => H-C  (canonical frame not bit-identical)
#   tau1_0 / tau1_nuc differ    => H-A  (prestress/nucleation gap)
# Source (1)/(3) => proceed to the reconcile (Phase 2); (2) => fix the
# exchange instead.  See PLAN_shared_fault_reconcile_fix_2026-05-23.md §Phase 0.
#
# Requires the branch-`safs` commit carrying the %.17e + RAWDOF trace extension
# (Phase 0) and the non-aborting VerifySharedFaultDOFDataConsistency mode.

export LC_ALL=C
export LANG=C

# ---- R-DIP2 cross-rank diagnostic flags (ALL baked in here) -------------
# Enable the env-gated [XRANK] trace.  Default-off in the code => only this job
# (and any that export SEAS_DIAG_XRANK) pays any cost; production jobs are
# byte-exact.
export SEAS_DIAG_XRANK=1
# Target the diverging shared fault QP reported by the R-101 abort (job
# 7746601).  Change if a future run reports a different centroid.
export SEAS_DIAG_XRANK_QP="607518,3706359,-4543"
# Capture radius [m] around the target (catches the face's ~3 QPs; the 500 m
# mesh's nearest QPs are well inside this).
export SEAS_DIAG_XRANK_R=200
# Earliest sim-time [s] to print (0 => full history, incl. the bit-identical
# pre-divergence phase so the FIRST divergence step is visible — needed to
# discriminate seed source (2) [differs at the first sub-step] from (3) [drifts
# later]).  With the non-fatal gate this run goes to tfinal (no early abort), so
# the trace volume is larger; raise to e.g. 0.3 if the log gets too big (you
# lose the source-(2)-vs-(3) discrimination but keep the onset).
export SEAS_DIAG_XRANK_TMIN=0.0

# Drop the hard R-101 gate (see header (1)): downgrade the shared-fault DOFData
# tripwire from MFEM_ABORT to a per-check rank-0 log so the oblique run proceeds
# PAST the cross-rank divergence.  Look for "[R-101 nonfatal] ... worst_rel=...
# (DIVERGED|ok)".  This is what lets us see bounded-vs-blowup and the speckle.
export SEAS_R101_NONFATAL=1
# Rupture-area + onset diagnostic (print-only): [DIAG]/[DIAG-ONSET] lines +
# V_max history, so bounded (peaks then decays) vs blow-up is readable from the
# log and we confirm the run still ruptures.
export SEAS_DIAG_BLOWUP=1
# -------------------------------------------------------------------------

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

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'   config/config.mk || { echo "ERROR: no HDF5"; exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: no H5Z_ZFP"; exit 1; }

cd miniapps/seas
make seas_spatial_dyn_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

SAFS_ROOT=safs/project_7.0_alternative
# Default = the resolution-test variant (D_c=2.0 / radius=4000 / cfl=0.5).
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml}"

# Mesh: read [mesh].path from the config (single source of truth = 500 m
# z0embed); override with SAFS_MESH.
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

# Flux: default = config's mixed_flux (adjacent) by NOT passing --mixed-flux.
MIXFLUX_ARG=""
if [[ -n "${SAFS_MIXED_FLUX:-}" ]]; then
    case "${SAFS_MIXED_FLUX}" in
        none|adjacent|all_continuous) : ;;
        *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid (none|adjacent|all_continuous)"; exit 1 ;;
    esac
    MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"
fi

# tfinal default 2 s (covers the divergence + the t~1.0 s blow-up window);
# cfl default = config (0.5).
TFINAL="${SAFS_TFINAL:-2s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_resDc2diag_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_resDc2diag_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}"
echo "Config:      ${CONFIG_TOML}   (RESOLUTION test: D_c=2.0 / radius=4000; OBLIQUE: dip prestress PRESENT)"
echo "R-101 gate:  NON-FATAL (SEAS_R101_NONFATAL=${SEAS_R101_NONFATAL:-<unset>}) -- run proceeds past cross-rank divergence"
echo "Blowup diag: SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP:-<unset>}  (read V_max bounded-vs-blowup + [DIAG-ONSET])"
echo "XRANK diag:  SEAS_DIAG_XRANK=${SEAS_DIAG_XRANK} QP=${SEAS_DIAG_XRANK_QP} R=${SEAS_DIAG_XRANK_R} TMIN=${SEAS_DIAG_XRANK_TMIN}"
echo "ParaView:    ${OUT}/fault.vtkhdf  (speckle baseline -- post-fix rerun must clear it)"
echo "tfinal:      ${TFINAL}"
echo "CFL arg:     ${CFL_ARG:-<config default = 0.5>}"

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --no-sidecar-material \
    ${MIXFLUX_ARG} \
    --tfinal "${TFINAL}" \
    ${CFL_ARG} \
    --paraview \
    --paraview-fault-hdf5 \
    --paraview-fault-zfp-tol 1.0e-12 \
    --paraview-max-snapshots 200 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

# --- Extract the [XRANK] cross-rank trace for offline diffing -------------
# Runs regardless of RC (the R-101 abort is EXPECTED; the trace is the point).
XRANK_LOG="${LOG_DIR}/xrank_resDc2diag_${SLURM_JOB_ID}.log"
grep -h '\[XRANK\]' "${LOG}" > "${XRANK_LOG}" 2>/dev/null || true
NX=$(wc -l < "${XRANK_LOG}" 2>/dev/null || echo 0)
echo "XRANK trace: ${XRANK_LOG}  (${NX} lines)"
if [[ "${NX}" -gt 0 ]]; then
    # Split by the (exactly two) ranks owning the diverging shared QP so the
    # two streams can be diffed directly.
    for r in $(grep -ho 'rank=[0-9][0-9]*' "${XRANK_LOG}" | sort -u); do
        rn=${r#rank=}
        grep -h "${r} " "${XRANK_LOG}" \
            > "${LOG_DIR}/xrank_resDc2diag_${SLURM_JOB_ID}_r${rn}.log"
        echo "  per-rank: ${LOG_DIR}/xrank_resDc2diag_${SLURM_JOB_ID}_r${rn}.log "\
             "($(wc -l < "${LOG_DIR}/xrank_resDc2diag_${SLURM_JOB_ID}_r${rn}.log") lines)"
    done
    echo "  Diff the two per-rank files to localise the seed (Phase 0):"
    echo "    grep 'RAWDOF.*side=+' <rA> vs <rB> -> identical?  (=> source 1/3;"
    echo "      differ at the first sub-step => source 2 ghost-exchange bug)"
    echo "    I*_can (%.17e) differ but RAWDOF identical => source 1 interp."
    echo "    can_t1 => H-C frame;  tau1_0/tau1_nuc => H-A prestress/nuc."
else
    echo "  WARNING: no [XRANK] lines captured.  Either SEAS_DIAG_XRANK was not"
    echo "  exported into the driver env, the run never reached the target QP,"
    echo "  or the QP centroid (${SEAS_DIAG_XRANK_QP}) is wrong for this run"
    echo "  — check the [R-101 nonfatal] centroid in ${LOG}."
fi

# --- bounded-vs-blowup + speckle baseline readout -------------------------
echo "=== V_max + DIAG + R-101 tail (bounded peaks-then-decays vs blow-up) ==="
grep -E "^step |\[DIAG\]|\[DIAG-ONSET\]|\[R-101 nonfatal\]|PASS|FAIL|done\." \
    "${LOG}" | tail -40
echo "=== fault.vtkhdf (speckle baseline for the post-fix comparison) ==="
ls -la "${OUT}/fault.vtkhdf" 2>/dev/null \
    || echo "  (no fault.vtkhdf -- run stopped before the first snapshot?)"

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — with the R-101 gate NON-FATAL this is NOT an R-101"
    echo "  abort.  It is either the --print-derived well-posedness gate (mis-"
    echo "  config) OR a genuine BLOW-UP crash (V_max overflow -> Brent SIGABRT)"
    echo "  — which is itself a RESULT (dip DID blow up).  Read the V_max tail"
    echo "  above: peaks-then-decays => bounded like zerodip; runaway => blow-up."
    exit 0   # diagnostic job: a crash at the divergence is a recorded result.
fi

echo "driver_rc=${RC} (ran to tfinal=${TFINAL} with the gate non-fatal)."
echo "  Inspect ${OUT}/fault.vtkhdf in ParaView: the speckle outside the"
echo "  rupture is the cross-rank symptom (issue B), the pre-fix baseline."
exit 0
