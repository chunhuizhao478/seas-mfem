#!/bin/bash
#SBATCH -J safs_rsDP_triq500_pu_dev
#SBATCH -o safs_rsDP_triq500_pu_dev_%j.out
#SBATCH -e safs_rsDP_triq500_pu_dev_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# SAFS RATE-AND-STATE (aging, depth-varying a(z)/b(z)) + sigma_n CAP, PURE UPWIND
# (mixed_flux=none), on the CHEAP 500 m _triq MESH — SPONTANEOUS-RUPTURE TEST
# (development queue, 8N/400r/2 h).
# =========================================================================
# PURPOSE: test that nucleation works AND that a spontaneous rupture propagates,
# on the 500 m _triq mesh (1.16 M tets) instead of the expensive 250 m
# triqsubdiv mesh (16 M tets, OOM-prone, presplit-only).  The friction is
# DESIGNED so that — at 500 m cells — BOTH:
#   (1) L_b is FULLY RESOLVED: L_b(dip)=2.51 km = 5.0 cells @ 500 m, and
#   (2) the nucleation patch is COMFORTABLY SUPERCRITICAL: 8 km diameter =
#       1.32x the critical nucleation length L_inf (6.06 km), fitting the
#       ~10 km seismogenic depth.
# See the length-scale design in the config header.
#
# This is a SEPARATE case — the Dc=0.05 / 7 km / 250 m dev job
#   spatial_dyn_ratestate_depthprofile_normalcap_pureupwind_triq_8N_400r_dev_2hr_safs.sbatch
# and its config are UNTOUCHED.  The ONLY differences (all in the NEW config /
# NEW CSVs below) are:
#   Dc          0.05 -> 0.10                 (L_b scale: 2.5 -> 5.0 cells @ 500 m)
#   a(z) @patch 0.0170 -> 0.0127 (b-a 0.0091 -> 0.0134; holds L_inf ~6.06 km so the
#                                 8 km patch is supercritical AND fits <10 km depth)
#   nuc radius  3500 -> 4000                 (8 km-DIAMETER patch)
#   mesh        250 m triqsubdiv -> 500 m _triq
#
# NO MIXED FLUX: central/mixed flux is the confirmed SAFS LSW blow-up driver
# (pure-upwind 7751974 STABLE, mixed_flux=adjacent 7751975 ran away).  RS
# regularizes friction but we hold the flux fixed for an apples-to-apples
# comparison.  Analysis:
#   safs/project_7.0_alternative/debug_document/spatial_dynamic_rupture_mixed_flux_analysis_2026-05-26.md
#
# Config (floor enabled):
#   config/spatial_friction_rate_state_safs_projected_stress_depthprofile_Dc010_nuc8km_500m.toml
#   ([friction].sigma_n_strength_floor_pa = 10.0e6; Dc=0.10; a(z) <- param_a_500m_strongvw.csv,
#    (a-b)(z) <- param_a_minus_b_500m_strongvw.csv (VW->VS transition ~11 km);
#    gradual_overstress nucleation, 8 km-DIAMETER patch / 4 km e-fold radius).
#
# CAP CAVEAT (REVIEW R-002): on the RS path the floor touches only the post-solve
# slip-rate DECOMPOSITION strength; the Brent solve setting |V| still uses the
# unfloored |sigma_n|.  Tensile |sigma_n| strengthens RS friction (|V|~0), so the
# cap is largely INERT on RS — RS lacks the LSW tensile free-slip mode.  This run
# is for completeness / mesh comparison, not an active tensile-slip safeguard.
#
# RISK GATE (R-004 / R-020): shared-fault RS psi is 1st-order in time at
# ader_order>=2 and the np>1 RS path is only PARTIALLY validated.  The driver
# prints "Treat parallel RS results as PRELIMINARY" at np>1.
#
# RUN WINDOW (R-025): the SAFS box is reflection-free only for ~6.94 s
# (min_box_dim/cp) and PML is NOT wired.  tfinal defaults to 100 s, so RS
# seismic fields after ~7 s are reflection-CONTAMINATED.  The 500 m mesh's
# larger cells carry a LARGER CFL dt than the 250 m mesh, so the 2 h dev wall
# reaches FURTHER into the window (better for observing propagation); the job
# checkpoints and resumes via SAFS_RESTART.  For physical RS fields use a 6 s
# window (or wire a real PML layer).
#
# Flux: PURE UPWIND (--mixed-flux none) by DEFAULT (header note);
# override with SAFS_MIXED_FLUX={none|adjacent|all_continuous}.

export LC_ALL=C
export LANG=C

# Lean diag (match the slip-weakening jobs): the V_max trace (SEAS_DIAG_BLOWUP)
# is the cap/stability monitor; the verbose per-step traces default OFF.
export SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-0}"
export SEAS_DIAG_BLOWUP="${SEAS_DIAG_BLOWUP:-1}"
export SEAS_DIAG_SLIP_VTHR="${SEAS_DIAG_SLIP_VTHR:-10.0}"
export SEAS_NOOPENING="${SEAS_NOOPENING:-0}"
export SEAS_DIAG_FRAME="${SEAS_DIAG_FRAME:-0}"
export SEAS_DIAG_FRAME_XYZ="${SEAS_DIAG_FRAME_XYZ:-607518,3706359,-4543}"
export SEAS_DIAG_FRAME_R="${SEAS_DIAG_FRAME_R:-300}"
export SEAS_DIAG_MACRO="${SEAS_DIAG_MACRO:-0}"

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
# spatial_dyn_build_safs.sbatch, or a login/idev-node `make`).  Fail fast if missing.
if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not found/executable in $(pwd)."
    echo "       Build it FIRST (no in-job build, to avoid the relink race):"
    echo "         sbatch miniapps/seas/jobs/safs/spatial_dyn_build_safs.sbatch"
    echo "       or on a login/idev node:  cd miniapps/seas && make seas_spatial_dyn_driver"
    exit 1
fi

SAFS_ROOT=safs/project_7.0_alternative
# Default = the Dc=0.10 / 8 km-DIAMETER-nucleation / 500 m-_triq-mesh config
# WITH [friction].sigma_n_strength_floor_pa = 10.0e6.  DEDICATED config; the
# 250 m Dc=0.05 dev config (and the normal-queue parent) are untouched.
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_rate_state_safs_projected_stress_depthprofile_Dc010_nuc8km_500m.toml}"

# PRE-FLIGHT: this is the CAP run — fail loudly if the config does not enable the
# floor.  (For a cap-OFF mesh-isolation run, set SAFS_CONFIG to a no-floor config.)
if ! grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null; then
    echo "ERROR: ${CONFIG_TOML} does not set sigma_n_strength_floor_pa — this CAP"
    echo "       run REQUIRES the floor.  (For a cap-OFF mesh-isolation run, set"
    echo "       SAFS_CONFIG to a config without the floor.)"
    exit 1
fi

# Mesh: DEFAULT = the 500 m shape-preserving _triq mesh; override with SAFS_MESH.
TRIQ_MESH="${SAFS_ROOT}/experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh"

# Mesh INPUT MODE (R-OOM).  The 500 m _triq mesh has ~1.16 M elements; a serial
# read replicates the FULL serial mesh on EVERY rank before partitioning, so
# per-node memory ~ (ranks_per_node x elements).  At 8N/400r (50 ranks/node) the
# 1.16 M-tet replica is ~0.3 GB/rank x 50 = ~15 GB/node — well within the Frontera
# node RAM, so this job DEFAULTS to SERIAL (no presplit needed — the whole point
# of dropping to 500 m from the 16 M-tet 250 m mesh that DID need presplit):
#   SAFS_MESH_MODE=serial (DEFAULT)  --mesh <.msh>; read + in-memory METIS on all
#                                    ranks.  OOM-safe here at 1.16 M tets.
#   SAFS_MESH_MODE=presplit          --mesh <prefix>; per-rank files, no replication.
#                                    Optional; only if you partition the 500 m mesh
#                                    at THIS np first.
SAFS_MESH_MODE="${SAFS_MESH_MODE:-serial}"
case "${SAFS_MESH_MODE}" in
    serial|presplit) : ;;
    *) echo "ERROR: SAFS_MESH_MODE='${SAFS_MESH_MODE}' invalid (serial|presplit)"; exit 1 ;;
esac

# presplit prefix default keyed by task count (only used if SAFS_MESH_MODE=presplit).
PRESPLIT_PREFIX="${SAFS_MESH_PREFIX:-${SAFS_ROOT}/experimental_mesh_refinement/presplit/safs_500m_triq_np${SLURM_NTASKS}}"

if [[ "${SAFS_MESH_MODE}" == "presplit" ]]; then
    MESH="${SAFS_MESH:-${PRESPLIT_PREFIX}}"
    if [[ ! -f "${MESH}.np" ]]; then
        echo "ERROR: pre-split manifest not found: $(pwd)/${MESH}.np"
        echo "       SAFS_MESH_MODE=presplit needs per-rank files built at -n ${SLURM_NTASKS}."
        echo "       The DEFAULT serial mode needs NO presplit — drop SAFS_MESH_MODE to use it."
        exit 1
    fi
else
    MESH="${SAFS_MESH:-${TRIQ_MESH}}"
    if [[ ! -f "${MESH}" ]]; then
        echo "ERROR: triq mesh not found:"
        echo "         $(pwd)/${MESH}"
        echo "       The .msh is gitignored and is NOT pulled with the repo."
        echo "       scp the re-triangulated mesh to the path above before submitting"
        echo "       (it lives under ${SAFS_ROOT}/experimental_mesh_refinement/)."
        exit 1
    fi
fi
MESHMODE_ARG="--mesh-mode ${SAFS_MESH_MODE}"

# Depth-profile CSVs (read at config-load time by the parser).  Fail early here.
for _csv in \
    "${SAFS_ROOT}/friction/rate-and-state/param_a_500m_strongvw.csv" \
    "${SAFS_ROOT}/friction/rate-and-state/param_a_minus_b_500m_strongvw.csv"; do
    if [[ ! -f "${_csv}" ]]; then
        echo "ERROR: depth-profile CSV not found: $(pwd)/${_csv}"
        exit 1
    fi
done

# Flux: PURE UPWIND (--mixed-flux none) by DEFAULT — central/mixed flux is the
# confirmed SAFS LSW blow-up driver; pure upwind carries CFL factor 1.0.
SAFS_MIXED_FLUX="${SAFS_MIXED_FLUX:-none}"
case "${SAFS_MIXED_FLUX}" in
    none|adjacent|all_continuous) : ;;
    *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid (none|adjacent|all_continuous)"; exit 1 ;;
esac
MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"

# tfinal default 100 s.  PAST the ~6.94 s reflection-free window (R-025) — see the
# RUN WINDOW note above.  Checkpoints + SAFS_RESTART for the rest.
TFINAL="${SAFS_TFINAL:-100s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Restart hook — resume past the 2 h dev wall toward tfinal.  Set SAFS_RESTART to
# the PREVIOUS job's checkpoint PREFIX (driver writes "<output-dir>/cp").
RESTART_ARG=""
[[ -n "${SAFS_RESTART:-}" ]] && RESTART_ARG="--restart ${SAFS_RESTART}"

OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_rsDP_triq500_pu_dev_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_rsDP_triq500_pu_dev_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Queue/size:  development, 8N x 400r, 2 h"
echo "EXPERIMENT:  RS depth-profile + sigma_n floor (10 MPa), PURE UPWIND (mixed_flux=none),"
echo "             on the 500 m _triq mesh, Dc=0.10, 8 km-dia SUPERCRITICAL nucleation"
echo "             (L_b=5 cells @500 m; patch 1.32x L_inf) — spontaneous-rupture test"
echo "Output dir:  ${OUT}"
echo "Log file:    ${LOG}"
echo "Mesh:        ${MESH}  (500 m _triq, ~1.16 M tets)"
echo "Mesh mode:   ${SAFS_MESH_MODE}  (serial=read+partition on every rank, OOM-safe at 1.16 M tets; presplit optional)"
echo "Config:      ${CONFIG_TOML}   (RS aging; depth-varying a(z)/b(z); Dc=0.10; 8 km-dia nuc; floor=10 MPa)"
echo "Diag:        SEAS_DIAG_SLIP=${SEAS_DIAG_SLIP} VTHR=${SEAS_DIAG_SLIP_VTHR} SEAS_DIAG_BLOWUP=${SEAS_DIAG_BLOWUP} SEAS_NOOPENING=${SEAS_NOOPENING}"
echo "Diag(frame): SEAS_DIAG_FRAME=${SEAS_DIAG_FRAME} XYZ=${SEAS_DIAG_FRAME_XYZ} R=${SEAS_DIAG_FRAME_R} SEAS_DIAG_MACRO=${SEAS_DIAG_MACRO}"
echo "Mixed flux:  ${SAFS_MIXED_FLUX} (pure upwind; CFL factor 1.0)"
echo "tfinal:      ${TFINAL}  (PAST the R-025 ~6.94 s reflection-free window; 500 m cells reach FURTHER than 250 m in the 2 h wall; SAFS_RESTART for the rest)"
echo "CFL arg:     ${CFL_ARG:-<config default = 0.5>}"
echo "Restart:     ${SAFS_RESTART:-<none — fresh start from t=0>}"

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    ${MESHMODE_ARG} \
    --no-sidecar-material \
    ${MIXFLUX_ARG} \
    --tfinal "${TFINAL}" \
    ${CFL_ARG} \
    ${RESTART_ARG} \
    --paraview \
    --paraview-fault-vtu \
    --paraview-fault-dt 1.0 \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

# Confirm the floor was actually active (banner echo).
if grep -qE 'sigma_n strength floor:[[:space:]]*10[[:space:]]*MPa' "${LOG}"; then
    echo "floor: ACTIVE (banner reports 10 MPa)"
else
    echo "WARNING: did not find 'sigma_n strength floor: 10 MPa' in the banner —"
    echo "         the floor may be DISABLED.  Inspect ${LOG}."
fi

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG} (check the [derived] RS gate lines)."
    exit ${RC}
fi

# Post-run verifier — OPTIONAL (needs h5py, absent on Frontera python3).
module load python3
VERIFIER="${SEAS_MFEM_ROOT}/miniapps/seas/safs/project_7.0_alternative/spatial/code/scripts/verify_spatial_dyn_smoke_safs.py"
CONFIG_TOML_ABS="${SEAS_MFEM_ROOT}/miniapps/seas/${CONFIG_TOML}"
if python3 -c "import h5py" 2>/dev/null; then
    python3 "${VERIFIER}" \
        --fault-vtkhdf "${OUT}/fault.vtkhdf" \
        --toml-config "${CONFIG_TOML_ABS}" \
        --log "${LOG}" \
        > "${LOG_DIR}/verify_summary_rsDP_triq500_pu_dev_${SLURM_JOB_ID}.txt"
    CHECK_RC=$?
    echo "driver_rc=${RC}, verify_rc=${CHECK_RC}"
    exit $(( RC | CHECK_RC ))
else
    echo "driver_rc=${RC}, verify=SKIPPED (no h5py on this host)."
    echo "  Run the verifier on a host that has h5py:"
    echo "    python3 ${VERIFIER} \\"
    echo "        --fault-vtkhdf ${OUT}/fault.vtkhdf \\"
    echo "        --toml-config ${CONFIG_TOML_ABS} \\"
    echo "        --log ${LOG}"
    exit ${RC}
fi
