#!/bin/bash
#SBATCH -J safs_t2_prf
#SBATCH -o safs_t2_prf_%j.out
#SBATCH -e safs_t2_prf_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=100
#SBATCH --time=06:00:00
#SBATCH --mem=200000M
#SBATCH -A lbl107

# =========================================================================
# SAFS v3_4_1 RSSRW PREFERRED THERMAL CASE2 — MFEM port of the SeisSol deck, with
# CVM material + MIXED FLUX (adjacent) + RK.
# *** EXPANSE (SDSC) job — uniform 2N x 50r = 100 ranks / 2 h / 100 GB. ***
# Config: ${SAFS_ROOT}/config/spatial_friction_safs_seisol_v3_4_1_RSSRW_PREFERRED_THERMAL_CASE2.toml
# Plan:   safs/project_7.0_preferred/document/PLAN_thermal_case2_mixedflux_port_2026-07-08.md
# =========================================================================
# Ports ~/Downloads/seisol_quakeworx/safs_seisol_v3_4_1_RSSRW_PREFERRED_THERMAL_CASE2:
#   * FL=103 (RS slip law + strong rate weakening), f_0=0.6, f_w=0.0, Dc=0.10, b=0.019
#   * rs_a / rs_srW SPATIAL, from the SCEC CTM temperature zoning, read through the
#     NEW [friction.rate_state.sidecar] block (the same two baked fields SeisSol
#     reads via ASAGI)
#   * Andersonian C1 k=1.8 effective stress; statewide muscal CVM material
#   * SCEC compact-bell nucleation: 75 MPa, R=2000 m, T_nuc=1 s, at the on-fault
#     10 km hypocentre (609062.8722, 3709528.1324, -10000)
#   * deep19km mesh — its .puml.h5 sibling is sha256-identical to the deck's
#     mesh_preferred.puml.h5, so MFEM and SeisSol run the SAME mesh.
#
# mixed_flux=adjacent + time_integrator=rk4 (central flux needs an RK integrator —
# UNSTABLE under ADER).  DEFAULT rk4 per plan decision D3: rk4 costs 5 wave.Mult per
# step (4 stages + the non-FSAL endpoint re-eval) vs DP45's 6, at the same dt.  The
# old np>1 endpoint-Mult DEADLOCK that once forced rk45 IS FIXED — the endpoint block
# now gates on the TABLEAU, not on `n > 0` (dynamic/rk_time_stepper.hpp:316-321), so
# fault-free ranks no longer skip it.  Set SAFS_TIME_INTEGRATOR=rk45 for the
# conservative fallback.  CFL 0.5 = SeisSol's CFL; the RK central-flux de-rate
# (adjacent 0.6x) is applied internally.  NOTE --cfl-dg-safety is INERT on the RK path.
#
# THIS IS NOT A BIT-FOR-BIT SeisSol REPRODUCTION: SeisSol runs ADER-DG order 4 (p=3)
# with clustered LTS and a pure Godunov flux; this runs p=1 global time stepping with
# a central flux on fault-adjacent faces.  The INPUTS transfer exactly; the
# DISCRETIZATION is the experimental variable.
#
# ARMS (plan decision D1) — select with SAFS_CONFIG:
#   parity   (DEFAULT): ..._THERMAL_CASE2.toml           sigma_n strength floor OMITTED
#   fallback          : ..._THERMAL_CASE2_snfloor10.toml sigma_n strength floor 10 MPa
# Run the PARITY arm first.  Switch only if the first fault snapshots show the
# shallow trace-cell runaway (isolated cells at sigma_n 0.1-5 MPa sliding ~4 m/s).
#
# ---- EXPANSE contract (matches build_expanse.sh) ----------------------------
#   * Launcher : srun --mpi=pmi2 -n $SLURM_NTASKS   (there is NO ibrun on Expanse)
#   * Toolchain: cpu/0.17.3b + gcc/10.2.0 + openmpi/4.1.x  (NOT intel/impi).
#     hypre / metis / PETSc / MUMPS / HDF5 / ZFP are built FROM SOURCE under
#     <root>/extern/*/install by build_expanse.sh — NOT loaded as modules.
#   * Runtime libs: libhdf5.so (sidecar reader) + libzfp.so live in extern/*;
#     added to LD_LIBRARY_PATH below; HDF5_PLUGIN_PATH -> the H5Z-ZFP plugin.
#   * Sizing : compute (exclusive 128-core / 256 GB nodes), 2N x 50r = 100 ranks,
#     --mem=200000M (200 GB/node = 4 GB/rank), 2 h.
#     RAISED from the v3.0.0 job's 100000M.  io/data_field_3d.cpp REPLICATES each
#     sidecar on EVERY rank (full H5Dread into a per-rank vector), and v3_4_1's
#     three sidecars total ~1330 MB/rank -- 66% of a 2 GB/rank budget, before the
#     mesh, wave operator and RK stage vectors.  compute nodes are exclusive and
#     charged whole-node, so the extra memory is free.  The pre-flight below
#     re-checks this against whatever sizing you actually submit with.
#     Output goes UNDER the repo (qstore): <root>/miniapps/seas/safs
#     (override via SAFS_OUTDIR).
#   * Account : EDIT the -A line above, or submit with
#       sbatch --account=<your_expanse_alloc> <this file>
#
# FIRST/SMOKE RUN (compute, 2N x 50r = 100 ranks, 2 h): RK45 + this mesh will NOT
# reach tfinal=150 s in 2 h — the run advances through nucleation/breakout and
# CHECKPOINTS.  ACCEPTANCE = clean start, V_max bounded, on-fault sigma_n stable in
# time (the mixed-flux leak cure), FreeSurf banner active.  Resume toward 150 s with
# the restart sibling (or chain this with SAFS_RESTART).
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
# SAFS_DIAG (0/1, default 0 = quiet log), SAFS_RESTART, SAFS_OUTDIR (qstore base).

export LC_ALL=C
export LANG=C
# Per-step [DIAG]/[DIAG-ONSET] monitor (V_max / sigma_n drift).  OFF by default —
# once V_max>10 it prints EVERY step (+ per-DOF onset dumps) and floods the log.
# Set SAFS_DIAG=1 to enable.  (getenv() treats SEAS_DIAG_BLOWUP=0 as "set", so we
# must UNSET it to disable, not set it to 0.)
if [[ "${SAFS_DIAG:-0}" == "1" ]]; then export SEAS_DIAG_BLOWUP=1; else unset SEAS_DIAG_BLOWUP; fi
# R-101 shared-fault cross-rank consistency tripwire — on the CURVILINEAR SAFS
# fault the per-QP strike frame can flip between the two ranks owning a shared
# face, so this check can abort during nucleation (only ~0.02% of fault DOFs sit
# on shared faces).  SKIPPED by default for SAFS; SAFS_R101=1 re-enables it.
# (Needs the SEAS_R101_SKIP gate in the driver — rebuild to pick it up.)
if [[ "${SAFS_R101:-0}" == "1" ]]; then unset SEAS_R101_SKIP; else export SEAS_R101_SKIP=1; fi

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

SAFS_ROOT=safs/project_7.0_preferred
CONFIG_TOML="${SAFS_CONFIG:-${SAFS_ROOT}/config/spatial_friction_safs_seisol_v3_4_1_RSSRW_PREFERRED_THERMAL_CASE2.toml}"
if [[ ! -f "${CONFIG_TOML}" ]]; then
    echo "ERROR: config not found: $(pwd)/${CONFIG_TOML}"; exit 1
fi

# Mesh: DEFAULT = the deep19km safstags mesh (tags 101/102/103/104, flat top).
# Its .puml.h5 sibling is sha256-identical to the SeisSol deck's mesh_preferred.puml.h5.
DEFAULT_MESH="safs/project_7.0_preferred/meshing_deep19km/results/safv4_deep19km_sub500m_flattop_nwtrim_fixed_opt_safstags.msh"
MESH="${SAFS_MESH:-${DEFAULT_MESH}}"
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: mesh not found: $(pwd)/${MESH}  (gitignored/large — scp it first)"; exit 1
fi

# PRODUCTION SIDECAR PRE-FLIGHT — THREE .h5 sidecars, all gitignored/large.
# The friction sidecar is NEW in v3_4_1 (the CTM temperature zoning).
VEL_H5="safs/project_7.0_preferred/velocity/results/multiscale_statewise_cvm_v3_4_1/velocity_safs.h5"
STR_H5="safs/project_7.0_preferred/stress/results/andersonian_k1p5_freeze500m_mfem/safs_stress_andersonian_k1.5_freeze500m_mfem.h5"
FRI_H5="safs/project_7.0_preferred/friction/results/thermal_case2_mfem/friction_safs.h5"
for _h5 in "${VEL_H5}" "${STR_H5}" "${FRI_H5}"; do
    if [[ ! -f "${_h5}" ]]; then
        echo "ERROR: required sidecar missing: $(pwd)/${_h5}"
        echo "       It is gitignored/large; scp it to Expanse before submitting."
        echo "       Rebuild recipes:"
        echo "         vel : safs/project_7.0_preferred/velocity/code/build_pref_velocity_sidecar.py --in safs_material_cvm.nc --out ${VEL_H5}"
        echo "         str : csm_stress_nc_to_mfem_hdf5.py --pad-m 3000 --in <combined_workflow>/outputs/safs_stress_andersonian_k1.5_freeze500m.nc --out ${STR_H5}"
        echo "         fri : safs/project_7.0_preferred/friction/code/build_pref_friction_sidecar.py --in safs_friction_thermal_case2.nc --out ${FRI_H5}"
        echo "       (Do NOT work around this with --no-sidecar-material — that forces the"
        echo "        constant fallback, which aborts on the matrix path.)"
        exit 1
    fi
done

# The config MUST point at the friction sidecar we just verified, otherwise the
# pre-flight above is checking a file the run never opens.
grep -Fq "${FRI_H5}" "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} does not reference ${FRI_H5}"; exit 1; }

# ---- PER-RANK MEMORY PRE-FLIGHT (v3_4_1) -----------------------------------
# io/data_field_3d.cpp holds the WHOLE field in RAM on EVERY rank
# (data_.assign(Nx*Ny*Nz, 0.0) + a full H5Dread), so all three sidecars are
# REPLICATED per rank, not distributed.  The datasets are uncompressed float64,
# so the on-disk size is an accurate proxy for the resident size.
#
# v3_4_1 raises this from ~780 MB/rank (v3.0.0: 259 MB CVM + 520 MB stress) to
# ~1330 MB/rank, because the statewide CVM cube (452x361x194) replaced the
# smaller v3.0.0 cube (297x196x194) and the friction sidecar is new.  At the
# stock 2N x 50r / --mem=100000M sizing that is 2000 MB/rank -- the sidecars
# alone take ~66% of it, before the mesh, the wave operator and the RK stage
# vectors.  Fail loudly here rather than OOM 20 minutes in.
#
# R-008: this uses ON-DISK size as a proxy for RESIDENT size.  That is exact only
# while the sidecar writers emit UNCOMPRESSED float64 -- they do (sidecar.py and
# csm_stress_nc_to_mfem_hdf5.py pass no `compression=` to h5py).  If a writer ever
# gains gzip/szip, disk size collapses while DataField3D still materialises the
# full float64 array, and this guard would silently UNDER-count.  Compute from the
# dataset shapes (Nx*Ny*Nz*8*n_fields) instead if that day comes.
_sidecar_mb=$(du -m "${VEL_H5}" "${STR_H5}" "${FRI_H5}" | awk '{s+=$1} END {print s}')
_mem_mb="${SLURM_MEM_PER_NODE:-0}"
_tpn="${SLURM_NTASKS_PER_NODE:-1}"
if [[ "${_mem_mb}" -gt 0 && "${_tpn}" -gt 0 ]]; then
    _per_rank_mb=$(( _mem_mb / _tpn ))
    _pct=$(( 100 * _sidecar_mb / _per_rank_mb ))
    echo "Memory:      sidecars ${_sidecar_mb} MB/rank of ${_per_rank_mb} MB/rank budget (${_pct}%)"
    if [[ "${_pct}" -ge 90 ]]; then
        echo "ERROR: the replicated sidecars alone need ${_pct}% of the per-rank memory budget."
        echo "       This job WILL OOM.  Pick one:"
        echo "         (a) fewer ranks per node, more memory each, e.g."
        echo "               sbatch --ntasks-per-node=25 --mem=240000M <this file>"
        echo "         (b) crop the CVM sidecar to the mesh bbox before scp'ing it"
        echo "             (x[253447,711231] y[3532686,3956757] z[-39329,26]) — that alone"
        echo "             removes roughly half of velocity_safs.h5."
        echo "       Override this guard with SAFS_SKIP_MEM_CHECK=1 if you know better."
        [[ "${SAFS_SKIP_MEM_CHECK:-0}" == "1" ]] || exit 1
    elif [[ "${_pct}" -ge 60 ]]; then
        echo "WARNING: replicated sidecars take ${_pct}% of the per-rank memory budget;"
        echo "         the mesh + wave operator + RK stage vectors must fit in the rest."
        echo "         If this job OOMs, use --ntasks-per-node=25 --mem=240000M."
    fi
fi

# Floor pre-flight (plan decision D1).  The PARITY arm OMITS the key (the parser
# REJECTS an explicit negative value, so 'disabled' == 'absent'); the FALLBACK arm
# sets it to 10.0e6.  Both are legal here — echo which arm is running so the log
# records it, and reject a value that is neither.
_floor_line="$(grep -E '^[[:space:]]*sigma_n_strength_floor_pa[[:space:]]*=' "${CONFIG_TOML}" || true)"
if [[ -z "${_floor_line}" ]]; then
    SAFS_ARM="parity (sigma_n strength floor DISABLED — SeisSol parity)"
else
    _floor_val="$(sed -E 's/.*=[[:space:]]*//; s/[[:space:]]*(#.*)?$//' <<< "${_floor_line}")"
    case "${_floor_val}" in
        -*) echo "ERROR: ${CONFIG_TOML} sets sigma_n_strength_floor_pa=${_floor_val}."
            echo "       The parser REJECTS an explicit negative value; OMIT the key to disable."; exit 1 ;;
        *)  SAFS_ARM="fallback (sigma_n strength floor = ${_floor_val} Pa — NOT a SeisSol comparison)" ;;
    esac
fi

# Flux: MIXED (adjacent) by DEFAULT — the production choice (central-on-fault-
# adjacent, the sigma_n-leak cure).  Requires RK (set below).
SAFS_MIXED_FLUX="${SAFS_MIXED_FLUX:-adjacent}"
case "${SAFS_MIXED_FLUX}" in
    none|adjacent|all_continuous) : ;;
    *) echo "ERROR: SAFS_MIXED_FLUX='${SAFS_MIXED_FLUX}' invalid"; exit 1 ;;
esac
MIXFLUX_ARG="--mixed-flux ${SAFS_MIXED_FLUX}"

# Time integrator: RK4 by DEFAULT (plan decision D3) — mixed/central flux is
# UNSTABLE under ADER (MatrixMixedFluxUnderAder hard guard); RK is mandatory with
# mixed flux.  rk4 does 5 wave.Mult per step (4 stages + the non-FSAL endpoint
# re-eval); DP45 does 6, at the same dt.  The np>1 endpoint-Mult deadlock that once
# forced rk45 was fixed by gating that block on the TABLEAU rather than on `n > 0`
# (dynamic/rk_time_stepper.hpp:316-321).  SAFS_TIME_INTEGRATOR=rk45 is the
# conservative fallback.
SAFS_TIME_INTEGRATOR="${SAFS_TIME_INTEGRATOR:-rk4}"
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

# Perf lever: --deriv-cache precomputes the per-element M^-1 K_d mat-vec.  Parity
# is <= 1e-12 (make test-bimaterial-deriv-cache-parity) and it works on the RK path.
# (--face-cache and --shared-ck-recursion are ADER-only; --cfl-dg-safety is inert
# on RK.)  Disable with SAFS_DERIV_CACHE=0.
DERIV_CACHE_ARG=""
[[ "${SAFS_DERIV_CACHE:-1}" == "1" ]] && DERIV_CACHE_ARG="--deriv-cache"

# Restart hook (optional): resume past the wall toward tfinal.  Set SAFS_RESTART to
# the PREVIOUS run's checkpoint PREFIX (driver writes "<output-dir>/cp"); the new
# --output-dir is distinct (job-ID-keyed), so the driver does not abort.
RESTART_ARG=""
[[ -n "${SAFS_RESTART:-}" ]] && RESTART_ARG="--restart ${SAFS_RESTART}"

# Output UNDER the repo on qstore (same convention as the Frontera jobs); the
# per-run dir is job-ID-keyed.  Override the base with SAFS_OUTDIR.
OUT_BASE="${SAFS_OUTDIR:-${SEAS_MFEM_ROOT}/miniapps/seas/safs}"
OUT="${OUT_BASE}/safs_dyn_t2_prf_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs/safs_expanse"
LOG="${LOG_DIR}/spatial_dyn_t2_prf_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Queue/size:  compute, 2N x 50r = ${SLURM_NTASKS} ranks, 200 GB/node, 02:00:00"
echo "EXPERIMENT:  v3_4_1 RSSRW PREFERRED THERMAL CASE2, CVM + Andersonian k1.8 + CTM friction sidecar, MIXED FLUX + RK"
echo "Arm:         ${SAFS_ARM}"
echo "Config:      ${CONFIG_TOML}"
echo "Mesh:        ${MESH}"
echo "Sidecars:    vel=${VEL_H5}"
echo "             str=${STR_H5}"
echo "             fri=${FRI_H5}"
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
    ${DERIV_CACHE_ARG} \
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
