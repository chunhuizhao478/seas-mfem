#!/bin/bash
#SBATCH -J safs_rs_prf_rst
#SBATCH -o safs_rs_prf_rst_%j.out
#SBATCH -e safs_rs_prf_rst_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=50
#SBATCH --time=02:00:00
#SBATCH --mem=100000M
#SBATCH -A usc143

# =========================================================================
# SAFS RATESTATE v3.0.0 PREFERRED (vw002) — MFEM port with CVM material + MIXED FLUX (adjacent) + RK45.
# *** EXPANSE (SDSC) job — uniform 2N x 50r = 100 ranks / 2 h / 100 GB. ***
# Config: ${SAFS_ROOT}/config/spatial_friction_safs_seisol_v3_0_0_RSSRW_PREFERRED_vw002.toml
# Plan:   safs/project_7.0_alternative/document/PLAN_mixed_flux_seissol_port_2026-06-19.md
# =========================================================================
# RS-SRW vw002 (f_0=0.6, f_w=0.10, Dc=0.10, V_w=0.02->3.0),  multiscale_statewise CVM material (the PRODUCTION material),
# spatially-varying CSM effective stress, mixed_flux=adjacent + time_integrator=rk45
# (central flux needs an RK integrator — UNSTABLE under ADER; rk45 is FSAL so it
# avoids the np>1 endpoint-Mult issue rk4 hit).  CFL 0.5 = STARTING calibration for
# the RK central-flux de-rate (adjacent 0.6x) — sweep via SAFS_CFL.
#
# ---- EXPANSE contract (matches build_expanse.sh) ----------------------------
#   * Launcher : srun --mpi=pmi2 -n $SLURM_NTASKS   (there is NO ibrun on Expanse)
#   * Toolchain: cpu/0.17.3b + gcc/10.2.0 + openmpi/4.1.x  (NOT intel/impi).
#     hypre / metis / PETSc / MUMPS / HDF5 / ZFP are built FROM SOURCE under
#     <root>/extern/*/install by build_expanse.sh — NOT loaded as modules.
#   * Runtime libs: libhdf5.so (sidecar reader) + libzfp.so live in extern/*;
#     added to LD_LIBRARY_PATH below; HDF5_PLUGIN_PATH -> the H5Z-ZFP plugin.
#   * Sizing : compute (exclusive 128-core / 256 GB nodes), 2N x 50r = 100 ranks,
#     --mem=100000M (100 GB/node = 2 GB/rank), 2 h.  This is the SAME resource
#     config the SeisSol ALT/PREFERRED cases used on Expanse.  Output goes UNDER
#     the repo (qstore): <root>/miniapps/seas/safs (override via SAFS_OUTDIR).
#   * Account : EDIT the -A line above, or submit with
#       sbatch --account=<your_expanse_alloc> <this file>
#
# RESTART/CONTINUATION (compute, 2N x 50r = 100 ranks, 2 h): resumes a 100-rank
# checkpoint (RESTART_FROM) toward tfinal=150 s.  Chain as many of these as needed.
#
# Scale up by OVERRIDING the sbatch headers at submit time, e.g.
#   sbatch --nodes=8 --ntasks-per-node=64 --time=48:00:00 <this file>
# (the restart sibling resumes a checkpoint at the SAME rank count it was written).
#
# PRODUCTION CHECKLIST (enforced below):
#   * NO --no-sidecar-material  (production loads the CVM sidecar)
#   * [velocity].use_sidecar=true + the velocity .h5 present under dataset_root
#   * the CSM stress .h5 present (gitignored/large; scp to Expanse scratch)
#
# Overrides: SAFS_CONFIG, SAFS_MESH, SAFS_MIXED_FLUX (default adjacent),
# SAFS_TIME_INTEGRATOR (default rk45), SAFS_TFINAL (default 150s), SAFS_CFL,
# SAFS_DIAG (0/1, default 0 = quiet log), RESTART_FROM (MANDATORY), SAFS_OUTDIR (qstore base).

export LC_ALL=C
export LANG=C
# Per-step [DIAG]/[DIAG-ONSET] monitor (V_max / sigma_n drift).  OFF by default —
# once V_max>10 it prints EVERY step (+ per-DOF onset dumps) and floods the log.
# Set SAFS_DIAG=1 to enable.  (getenv() treats SEAS_DIAG_BLOWUP=0 as "set", so we
# must UNSET it to disable, not set it to 0.)
if [[ "${SAFS_DIAG:-0}" == "1" ]]; then export SEAS_DIAG_BLOWUP=1; else unset SEAS_DIAG_BLOWUP; fi

# --- Expanse modules (MUST match build_expanse.sh: cpu/0.17.3b gcc openmpi) ---
module purge
module load slurm 2>/dev/null || true
module load cpu/0.17.3b
module load gcc/10.2.0
# OpenMPI: first that loads wins (same order as build_expanse.sh OPENMPI_MODULES).
for _mpi in openmpi/4.1.5 openmpi/4.1.3 openmpi/4.1.1 openmpi; do
    if module load "${_mpi}" 2>/dev/null && module is-loaded "${_mpi}" 2>/dev/null; then
        echo "OpenMPI module: ${_mpi}"; break
    fi
done

# Resolve seas-mfem root (walk UP to miniapps/seas/Makefile) BEFORE setting the
# runtime library paths — they live under <root>/extern/*/install.
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR}}"
_walk_up_dir="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5; do
    if [[ -f "${_walk_up_dir}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"; break
    fi
    _walk_up_dir="${_walk_up_dir}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}."
    echo "  Set SEAS_MFEM_ROOT via sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi

# Runtime shared libs from build_expanse.sh's from-source deps.  libhdf5.so is
# REQUIRED (velocity + stress sidecars are HDF5); libzfp.so + the H5Z-ZFP plugin
# are needed only if a ZFP-filtered dataset is read/written (harmless to export
# always).  PETSc/hypre/metis are static; petsc/install/lib covers a shared openblas.
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"

cd "${SEAS_MFEM_ROOT}"

grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
# No in-job build (avoid the concurrent relink race); the binary must be pre-built.
if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not found.  Build it FIRST on Expanse"
    echo "  (SDSC forbids login-node compiles — use an interactive shared session):"
    echo "   bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    exit 1
fi

SAFS_ROOT=safs/project_7.0_alternative
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_safs_seisol_v3_0_0_RSSRW_PREFERRED_vw002.toml}"

# Mesh: DEFAULT = the safv4_deep 500 m mesh (safstags tags 101/102/103/104, non-flat topo); override with SAFS_MESH.
DEFAULT_MESH="safs/project_7.0_preferred/meshing/results/safv4_deep_500m_opt_safstags.msh"
MESH="${SAFS_MESH:-${DEFAULT_MESH}}"
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: mesh not found: $(pwd)/${MESH}  (gitignored/large — scp it first)"; exit 1
fi

# PRODUCTION SIDECAR PRE-FLIGHT — this run uses CVM material + the CSM stress
# field; BOTH .h5 sidecars are gitignored/large and MUST be present (scp'd).
VEL_H5="safs/project_7.0_preferred/velocity/results/multiscale_statewise_cvm/velocity_safs.h5"
STR_H5="safs/project_7.0_preferred/stress/results/csm_d12_rssrw_mfem/safs_stress_csm_d12_mfem.h5"
for _h5 in "${VEL_H5}" "${STR_H5}"; do
    if [[ ! -f "${_h5}" ]]; then
        echo "ERROR: required sidecar missing: $(pwd)/${_h5}"
        echo "       It is gitignored/large; scp it to Expanse before submitting."
        echo "       (Do NOT work around this with --no-sidecar-material — that forces the"
        echo "        constant fallback, which aborts on the matrix path.)"
        exit 1
    fi
done

# Floor pre-flight (v3 PREFERRED ratestate keeps the 10 MPa MFEM-only sigma_n floor).
grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" >/dev/null     || { echo "ERROR: ${CONFIG_TOML} has no sigma_n_strength_floor_pa"; exit 1; }

# Flux: MIXED (adjacent) by DEFAULT — the production choice (central-on-fault-
# adjacent, the sigma_n-leak cure).  Requires RK (set below).
SAFS_MIXED_FLUX="${SAFS_MIXED_FLUX:-adjacent}"
case "${SAFS_MIXED_FLUX}" in
    none|adjacent|all_continuous) : ;;
    *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid"; exit 1 ;;
esac
MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"

# Time integrator: RK45 by DEFAULT — mixed/central flux is UNSTABLE under ADER
# (MatrixMixedFluxUnderAder hard guard); RK is mandatory with mixed flux.  rk45
# (Dormand-Prince) is FSAL, so it skips the endpoint re-eval Mult that rk4 needed
# (the np>1 deadlock fix); use rk45 for mixed flux.
SAFS_TIME_INTEGRATOR="${SAFS_TIME_INTEGRATOR:-rk45}"
case "${SAFS_TIME_INTEGRATOR}" in
    ader|rk4|rk45) : ;;
    *) echo "ERROR: SAFS_TIME_INTEGRATOR='${SAFS_TIME_INTEGRATOR}' invalid (ader|rk4|rk45)"; exit 1 ;;
esac
if [[ "${SAFS_MIXED_FLUX}" != "none" && "${SAFS_TIME_INTEGRATOR}" == "ader" ]]; then
    echo "ERROR: mixed_flux=${SAFS_MIXED_FLUX} + ADER is UNSTABLE (driver aborts). Use rk4/rk45."; exit 1
fi
TI_ARG="--time-integrator ${SAFS_TIME_INTEGRATOR}"

TFINAL="${SAFS_TFINAL:-150s}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}" ]] && CFL_ARG="--cfl ${SAFS_CFL}"

# Restart (MANDATORY for this job): continue RESTART_FROM's checkpoint set.  Set
# RESTART_FROM to the run dir holding cp_checkpoint_r*.txt (one per rank):
#   sbatch --export=ALL,RESTART_FROM=/expanse/lustre/.../run_dir <this file>
# A checkpoint can ONLY resume at the SAME rank count it was written at (the driver
# hard-checks num_ranks; this job is 100 ranks).  The restart OUTPUT dir is distinct
# (job-ID-keyed), so the driver does not abort (restart != output parent).
RESTART_FROM="${RESTART_FROM:-}"
if [[ -z "${RESTART_FROM}" ]]; then
    echo "ERROR: RESTART_FROM is empty.  Set it to the run dir holding cp_checkpoint_r*.txt:"
    echo "         sbatch --export=ALL,RESTART_FROM=/expanse/lustre/.../run_dir <this file>"
    exit 1
fi
RESTART_FROM="${RESTART_FROM%/}"
if [[ ! -f "${RESTART_FROM}/cp_checkpoint_r0.txt" ]]; then
    echo "ERROR: no checkpoint at ${RESTART_FROM}/cp_checkpoint_r0.txt"; exit 1
fi
_n_cp=$(ls "${RESTART_FROM}"/cp_checkpoint_r*.txt 2>/dev/null | wc -l)
if [[ "${_n_cp}" -ne "${SLURM_NTASKS}" ]]; then
    echo "ERROR: found ${_n_cp} checkpoint files in ${RESTART_FROM} (expected ${SLURM_NTASKS} = this job's rank count)."; exit 1
fi
RESTART_ARG="--restart ${RESTART_FROM}/cp"

# Output UNDER the repo on qstore (same convention as the Frontera jobs); the
# per-run dir is job-ID-keyed.  Override the base with SAFS_OUTDIR.
OUT_BASE="${SAFS_OUTDIR:-${SEAS_MFEM_ROOT}/miniapps/seas/safs}"
OUT="${OUT_BASE}/safs_dyn_rs_prf_rst_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs/safs_expanse"
LOG="${LOG_DIR}/spatial_dyn_rs_prf_rst_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Queue/size:  compute, 2N x 50r = ${SLURM_NTASKS} ranks, 100 GB/node, 02:00:00"
echo "EXPERIMENT:  ratestate v3 PREFERRED (vw002), CVM material + CSM stress, MIXED FLUX + RK"
echo "Config:      ${CONFIG_TOML}"
echo "Mesh:        ${MESH}"
echo "Sidecars:    vel=${VEL_H5}  stress=${STR_H5}"
echo "Mixed flux:  ${SAFS_MIXED_FLUX}    Integrator: ${SAFS_TIME_INTEGRATOR}"
echo "tfinal:      ${TFINAL}    CFL: ${CFL_ARG:-<config default 0.5>}"
echo "Restart:     ${RESTART_ARG:-<none — fresh start from t=0>}"
echo "Output dir:  ${OUT}"

# NOTE: production invocation deliberately OMITS --no-sidecar-material.
# stdbuf -oL -eL line-buffers the driver stdout/stderr so ${LOG} updates LIVE
# (output-init messages + the per-100-step [DIAG-SIGN] progress) instead of only
# flushing in ~4 KB blocks — so a stalled run is visible immediately.
srun --mpi=pmi2 -n "${SLURM_NTASKS}" stdbuf -oL -eL ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    ${MIXFLUX_ARG} \
    ${TI_ARG} \
    --tfinal "${TFINAL}" \
    ${CFL_ARG} \
    ${RESTART_ARG} \
    --paraview \
    --paraview-fault-vtu \
    --paraview-free-surface vtu \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

# Banner check (the parser silently ignores unknown flags, so verify the EFFECT):
if grep -qE 'FreeSurf:' "${LOG}"; then
    echo "FreeSurf output: $(grep -m1 -E 'FreeSurf:' "${LOG}")"
else
    echo "WARNING: no 'FreeSurf:' banner line in ${LOG} — free-surface output may be OFF."
fi

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG} (check [derived] lines + any guard ABORT)."; exit ${RC}
fi
echo "driver_rc=${RC} — done.  Inspect ${LOG}: V_max(t) bounded + on-fault sigma_n stable in time."
exit ${RC}
