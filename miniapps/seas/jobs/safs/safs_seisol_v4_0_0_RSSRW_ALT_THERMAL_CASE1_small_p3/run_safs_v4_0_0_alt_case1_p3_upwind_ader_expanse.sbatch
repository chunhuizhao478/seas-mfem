#!/bin/bash
#SBATCH --job-name=safs_v400_c1_p3_uw
#SBATCH --account=lbl107
#SBATCH --partition=compute
#SBATCH --constraint="lustre"        # match the SeisSol deck: run lives on /expanse/lustre
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128        # pure-MPI: 128 ranks x 1 core = the WHOLE 128-core node
#SBATCH --cpus-per-task=1            # MFEM has no OpenMP (SeisSol used 8 ranks x 16 threads instead)
#SBATCH --mem=140000M               # MATCHES the SeisSol deck.  Enabled by --sidecar-shared-mem
                                     # (PLAN_sidecar_mpi_shared_memory_2026-07-17): the 3 gridded
                                     # sidecars are ONE MPI-3 shared copy per NODE (~1.4 GB incl.
                                     # rank-0's transient read buffer), not ~1.05 GB x 128 ranks.
                                     # If you must run with SAFS_SHARED_MEM=0 (per-rank replication),
                                     # submit with --mem=249000M instead — the pre-flight enforces this.
#SBATCH -t 15:00:00
#SBATCH -o safs_v400_c1_p3_uw.%j.out
#SBATCH -e safs_v400_c1_p3_uw.%j.out
#SBATCH --export=ALL

# =============================================================================
# SAFS v4_0_0 RSSRW ALT THERMAL CASE1 (coarse "small" twin) -- MFEM port,
#   PURE UPWIND (Godunov) flux + ADER time integration.  *** p3 (order 3). ***
#
# OBJECTIVE: measure MFEM wall-clock vs SeisSol on the SAME coarse mesh
#   (safalt_dl_safgh_2M_walls_mmg.msh == mesh_alt.puml.h5, 1,319,294 tets).
#   SeisSol runs this deck at o4 (=p3), pure Godunov, clustered LTS, on Expanse
#   (2 nodes x 8 MPI x 16 OMP = 256 cores).  This job runs MFEM p3 on the same
#   2 nodes with 256 pure-MPI ranks (128/node = full core parity).  Submit the
#   _p1, _p2, _p3 siblings to bracket SeisSol's p3.
#
# EXPANSE (SDSC) -- srun --mpi=pmi2 (no ibrun); toolchain cpu/0.17.3b + gcc/10.2.0
#   + openmpi/4.1.x, matching build_expanse.sh.  hypre/metis/PETSc/MUMPS/HDF5/ZFP
#   are built FROM SOURCE under <root>/extern/*; their shared libs go on
#   LD_LIBRARY_PATH below.
#
# The binary must be PRE-BUILT (SDSC forbids login-node compiles):
#   bash <root>/build_expanse.sh   (in an interactive shared session)
#
# All inputs are SELF-CONTAINED in this folder as symlinks (mesh + 3 sidecars);
#   the .toml references them by repo-relative path (jobs/safs/<this>/...).  The
#   sidecars are gitignored/large -- rebuild locally (see README) and scp them to
#   their canonical targets under safs/project_7.0_alternative/ on Expanse.
#
# Overrides (env): SAFS_TFINAL (default: config's 150s), SAFS_CFL, SAFS_RESTART,
#   SAFS_OUTDIR (qstore output base), SAFS_DIAG (0/1), SAFS_R101 (0/1),
#   SAFS_DERIV_OPT (default "--deriv-cache --shared-ck-recursion"; set empty to
#   time the un-optimised ADER path),
#   SAFS_SHARED_MEM (default 1 = --sidecar-shared-mem, ONE sidecar copy per node;
#   set 0 for the legacy per-rank replication — then submit with --mem=249000M),
#   SAFS_CALI (default 1 = Caliper per-component time breakdown; 0 = off for a
#   strictly instrumentation-free wall-clock), SAFS_CALI_CONFIG (full CALI_CONFIG
#   override).
#
# CALIPER (per-component time breakdown; cf. the ..._caliper_... reference job):
#   The solver is instrumented with MFEM_PERF_SCOPE regions — on the ADER path:
#     seas::spatial_dyn::step, ComputeADERSubStepStates[AndIntegral] (CK
#     predictor -> ApplySpatialDerivative), seas::spatial_dyn::friction_substep,
#     AdvanceADER (corrector), ComputeVolumeRHS / ComputeFaceFluxRHS /
#     ComputeSharedFaceFluxRHS, ComputeMaxDt.
#   Caliper is driven ENTIRELY by the CALI_CONFIG env var (no driver flag).
#   REQUIRES a Caliper-enabled binary:  USE_CALIPER=YES bash build_expanse.sh
#   (USE_CALIPER defaults to NO -> regions compiled out).  Unlike the dedicated
#   profiling job, a missing libcaliper does NOT abort here — the speed run
#   proceeds unprofiled with a WARNING (this job's primary product is the
#   wall-clock; the breakdown is the diagnosis layer).  In the region report,
#   Max >> Avg across ranks = LOAD IMBALANCE (expected: fault on few ranks).
# Scale/parity: to leave headroom instead of full core parity, submit with e.g.
#   sbatch --ntasks-per-node=100 --mem=240000M <this file>
# =============================================================================

set -o pipefail
export LC_ALL=C
export LANG=C

# --- this run's identity (asserted against the .toml below) ------------------
ORDER=3
ADER_ORDER=4                                  # = ORDER + 1 (SeisSol convention)
CONFIG_BASENAME="spatial_friction_safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_p3.toml"

# --- diagnostics gates (OFF by default; see the ALT upwind_ader sibling) ------
if [[ "${SAFS_DIAG:-0}" == "1" ]]; then export SEAS_DIAG_BLOWUP=1; else unset SEAS_DIAG_BLOWUP; fi
if [[ "${SAFS_R101:-0}" == "1" ]]; then unset SEAS_R101_SKIP; else export SEAS_R101_SKIP=1; fi

# --- Expanse modules (MUST match build_expanse.sh) ---------------------------
module purge
module load slurm 2>/dev/null || true
module load cpu/0.17.3b
module load gcc/10.2.0
for _mpi in openmpi/4.1.5 openmpi/4.1.3 openmpi/4.1.1 openmpi; do
    if module load "${_mpi}" 2>/dev/null && module is-loaded "${_mpi}" 2>/dev/null; then
        echo "OpenMPI module: ${_mpi}"; break
    fi
done

# --- resolve the seas-mfem root (walk UP to miniapps/seas/Makefile) ----------
# This sbatch lives at <root>/miniapps/seas/jobs/safs/<folder>/, so the submit
# dir is 4 levels below <root>.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SUBMIT_DIR}}"
_walk="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5 6; do
    if [[ -f "${_walk}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk}" && pwd)"; break
    fi
    _walk="${_walk}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from ${SUBMIT_DIR}."
    echo "  Set SEAS_MFEM_ROOT via sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi

# --- runtime shared libs (from build_expanse.sh's from-source deps) ----------
# libcaliper.so first (harmless if the build has no Caliper), then hdf5/zfp/petsc.
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"
# HDF5 POSIX file locking FAILS on Lustre (errno 524) and can kill the sidecar
# H5Dread (and any HDF5 output).  Matches the SeisSol deck's runtime env.
export HDF5_USE_FILE_LOCKING=FALSE

cd "${SEAS_MFEM_ROOT}"
grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc in config.mk"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS in config.mk"; exit 1; }

cd miniapps/seas
if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not found.  Build it FIRST on Expanse:"
    echo "   bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    exit 1
fi

# --- config + input pre-flight (all paths resolved from miniapps/seas) --------
JOB_REL="jobs/safs/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_small_p3"
CONFIG_TOML="${JOB_REL}/${CONFIG_BASENAME}"
[[ -f "${CONFIG_TOML}" ]] || { echo "ERROR: config not found: $(pwd)/${CONFIG_TOML}"; exit 1; }

# The .toml MUST declare the order this sbatch was written for (there is no --order
# CLI flag; order is a TOML-only key).  Assert it so a mismatched pairing fails loud.
_cfg_order="$(grep -E '^order[[:space:]]*=' "${CONFIG_TOML}" | head -1 | sed -E 's/[^0-9]*([0-9]+).*/\1/')"
if [[ "${_cfg_order}" != "${ORDER}" ]]; then
    echo "ERROR: ${CONFIG_TOML} sets [mesh].order=${_cfg_order}, but this sbatch is p${ORDER}"
    echo "       (submit run_..._p${_cfg_order}_... instead, or fix the config)."; exit 1
fi
_cfg_ader="$(grep -E '^ader_order[[:space:]]*=' "${CONFIG_TOML}" | head -1 | sed -E 's/[^0-9]*([0-9]+).*/\1/')"
if [[ "${_cfg_ader}" != "${ADER_ORDER}" ]]; then
    echo "ERROR: ${CONFIG_TOML} sets ader_order=${_cfg_ader}, expected ${ADER_ORDER} (=order+1)."; exit 1
fi

# Inputs are folder-local symlinks; -e follows them, so a dangling link (sidecar
# not yet scp'd to its canonical target) fails here rather than 20 min in.
MESH_LN="${JOB_REL}/mesh_alt.msh"
VEL_LN="${JOB_REL}/velocity_safs.h5"
STR_LN="${JOB_REL}/stress_andersonian_k1.65.h5"
FRI_LN="${JOB_REL}/friction_thermal_case1.h5"
for _f in "${MESH_LN}" "${VEL_LN}" "${STR_LN}" "${FRI_LN}"; do
    if [[ ! -e "${_f}" ]]; then
        echo "ERROR: input missing / dangling symlink: $(pwd)/${_f}"
        echo "       -> $(readlink "${_f}" 2>/dev/null)"
        echo "       The sidecars are gitignored/large; rebuild them (see the folder README)"
        echo "       and scp them to their canonical targets under safs/project_7.0_alternative/."
        exit 1
    fi
done

# --- sidecar residency mode + binary support pre-flight ----------------------
# SAFS_SHARED_MEM=1 (default): pass --sidecar-shared-mem, so io/data_field_3d
# holds ONE MPI-3 shared copy of each gridded sidecar per NODE
# (PLAN_sidecar_mpi_shared_memory_2026-07-17.md).  =0: legacy per-rank
# replication (every rank H5Dreads its own full copy).
SAFS_SHARED_MEM="${SAFS_SHARED_MEM:-1}"
SHARED_MEM_ARG=""
if [[ "${SAFS_SHARED_MEM}" == "1" ]]; then
    SHARED_MEM_ARG="--sidecar-shared-mem"
    # CRITICAL: the driver's hand-rolled parser SILENTLY IGNORES unknown
    # flags.  A binary built before the shared-memory feature would ignore
    # --sidecar-shared-mem, replicate per rank, and OOM at this job's
    # --mem=140000M sizing 20 minutes in.  Fail HERE instead: the flag
    # string is a literal in the binary, so a plain binary grep (no
    # `strings` needed on the compute node) proves support.
    if ! LC_ALL=C grep -aq -- "--sidecar-shared-mem" ./seas_spatial_dyn_driver; then
        echo "ERROR: ./seas_spatial_dyn_driver does not know --sidecar-shared-mem"
        echo "       (built before PLAN_sidecar_mpi_shared_memory_2026-07-17)."
        echo "       Rebuild on Expanse first:  bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
        echo "       Or run the legacy mode:    SAFS_SHARED_MEM=0 sbatch --mem=249000M <this file>"
        echo "       Override with SAFS_SKIP_FLAG_CHECK=1 if you know better."
        [[ "${SAFS_SKIP_FLAG_CHECK:-0}" == "1" ]] || exit 1
    fi
fi

# --- memory pre-flight (mode-aware) ------------------------------------------
# On-disk float64 size == resident size (uncompressed datasets).
#   shared mode : the three sidecars are resident ONCE PER NODE (+ the loading
#                 rank's transient read buffer, <= the largest single field);
#                 the guard checks the LEFTOVER per-rank compute budget.
#   legacy mode : sidecars are resident PER RANK (the pre-2026-07-17 math).
_sidecar_mb=$(du -mL "${VEL_LN}" "${STR_LN}" "${FRI_LN}" | awk '{s+=$1} END {print s}')
_mem_mb="${SLURM_MEM_PER_NODE:-0}"
_tpn="${SLURM_NTASKS_PER_NODE:-1}"
if [[ "${_mem_mb}" -gt 0 && "${_tpn}" -gt 0 ]]; then
    if [[ "${SAFS_SHARED_MEM}" == "1" ]]; then
        # Largest single field ~= velocity/3 (3 equal fields in the biggest
        # sidecar) is the transient double buffer on the loading rank.
        _vel_mb=$(du -mL "${VEL_LN}" | awk '{print $1}')
        _node_resident_mb=$(( _sidecar_mb + _vel_mb / 3 ))
        _left_per_rank_mb=$(( (_mem_mb - _node_resident_mb) / _tpn ))
        echo "Memory:      SHARED sidecars ~${_node_resident_mb} MB/NODE of ${_mem_mb} MB/node;"
        echo "             leftover ~${_left_per_rank_mb} MB/rank for mesh + wave operator + ADER"
        if [[ "${_left_per_rank_mb}" -lt 300 ]]; then
            echo "ERROR: <300 MB/rank left after the shared sidecars -> likely OOM."
            echo "       Use fewer ranks per node (e.g. --ntasks-per-node=64) or more --mem."
            echo "       Override this guard with SAFS_SKIP_MEM_CHECK=1 if you know better."
            [[ "${SAFS_SKIP_MEM_CHECK:-0}" == "1" ]] || exit 1
        fi
    else
        _per_rank_mb=$(( _mem_mb / _tpn ))
        _pct=$(( 100 * _sidecar_mb / _per_rank_mb ))
        echo "Memory:      LEGACY replicated sidecars ${_sidecar_mb} MB/rank of ${_per_rank_mb} MB/rank budget (${_pct}%)"
        if [[ "${_pct}" -ge 90 ]]; then
            echo "ERROR: replicated sidecars need ${_pct}% of the per-rank budget -> WILL OOM."
            echo "       SAFS_SHARED_MEM=0 needs the old sizing:"
            echo "         SAFS_SHARED_MEM=0 sbatch --mem=249000M <this file>"
            echo "       (or fewer ranks per node: --ntasks-per-node=64)."
            echo "       Override this guard with SAFS_SKIP_MEM_CHECK=1 if you know better."
            [[ "${SAFS_SKIP_MEM_CHECK:-0}" == "1" ]] || exit 1
        elif [[ "${_pct}" -ge 60 ]]; then
            echo "WARNING: replicated sidecars take ${_pct}% of the per-rank budget; the mesh +"
            echo "         wave operator + ADER work vectors must fit in the rest."
        fi
    fi
fi

# --- numerics + optional ADER perf levers ------------------------------------
# PURE UPWIND + ADER (also the config defaults; passed explicitly + asserted).
# --deriv-cache: per-element cached M^-1 K_d mat-vec (scalar AND matrix path;
#   round-off only).  --shared-ck-recursion: run the ADER CK predictor ONCE per
#   macro-step (ADER-only; byte-identical).  Per the SAFS caliper study these are
#   the biggest wins (CK-fusion ~27% + deriv-cache).  Set SAFS_DERIV_OPT="" to
#   time the un-optimised ADER path.
DERIV_OPT="${SAFS_DERIV_OPT:---deriv-cache --shared-ck-recursion}"

TFINAL_ARG=""
[[ -n "${SAFS_TFINAL:-}" ]] && TFINAL_ARG="--tfinal ${SAFS_TFINAL}"
CFL_ARG=""
[[ -n "${SAFS_CFL:-}"    ]] && CFL_ARG="--cfl ${SAFS_CFL}"
RESTART_ARG=""
[[ -n "${SAFS_RESTART:-}" ]] && RESTART_ARG="--restart ${SAFS_RESTART}"

# --- output dir (under the repo on qstore; per-run job-ID-keyed) --------------
OUT_BASE="${SAFS_OUTDIR:-${SEAS_MFEM_ROOT}/miniapps/seas/safs}"
OUT="${OUT_BASE}/safs_dyn_v4_0_0_alt_c1_p3_uw_${SLURM_JOB_ID}"
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_small_p3"
LOG="${LOG_DIR}/run_p3_upwind_ader_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}"

# --- Caliper per-component time breakdown (SAFS_CALI=1 default) ---------------
# Two channels (same recipe as the ..._caliper_... reference job):
#   runtime-report : per-region table — exclusive + inclusive time, call counts,
#                    per-region MPI time, Min/Avg/Max across ranks (Max >> Avg
#                    = load imbalance).  THE "which part eats the time" answer.
#   spot           : per-rank .cali for Python Hatchet/Thicket / LLNL Spot UI
#                    (which RANKS are slow; deep dive).
# SOFT DEGRADE: if the binary is not linked against libcaliper (built without
# USE_CALIPER=YES), the regions are compiled out — we WARN and run unprofiled
# rather than kill the speed run.  SAFS_CALI=0 disables profiling entirely
# (strictly instrumentation-free wall-clock).
SAFS_CALI="${SAFS_CALI:-1}"
CALI_REGION_TXT="${LOG%.log}.cali-region-report.txt"
CALI_SPOT="${LOG%.log}.cali"
CALI_ACTIVE=0
if [[ "${SAFS_CALI}" == "1" ]]; then
    if ldd ./seas_spatial_dyn_driver 2>/dev/null | grep -qi 'libcaliper'; then
        export CALI_CONFIG="${SAFS_CALI_CONFIG:-runtime-report(calc.inclusive,region.count,profile.mpi,max_column_width=60,output=${CALI_REGION_TXT}),spot(output=${CALI_SPOT})}"
        CALI_ACTIVE=1
    else
        echo "WARNING: SAFS_CALI=1 but ./seas_spatial_dyn_driver is NOT linked against"
        echo "         libcaliper (built with USE_CALIPER=NO) -> the MFEM_PERF_SCOPE"
        echo "         regions are compiled out; running UNPROFILED."
        echo "         For the per-component breakdown rebuild with:"
        echo "           USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
        unset CALI_CONFIG
    fi
else
    unset CALI_CONFIG
fi

echo "Git commit:  $(git rev-parse HEAD 2>/dev/null)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Size:        compute (lustre), ${SLURM_NNODES:-2}N x ${SLURM_NTASKS_PER_NODE:-128}r = ${SLURM_NTASKS} ranks (pure MPI), ${SLURM_MEM_PER_NODE:-140000}M/node, ${SBATCH_TIMELIMIT:-15:00:00}"
echo "Sidecars:    residency mode = $([[ "${SAFS_SHARED_MEM}" == "1" ]] && echo 'SHARED (one copy per node, --sidecar-shared-mem)' || echo 'LEGACY (replicated per rank)')"
if [[ "${CALI_ACTIVE}" == "1" ]]; then
    echo "Caliper:     ON  $(ldd ./seas_spatial_dyn_driver 2>/dev/null | grep -m1 -i libcaliper | awk '{print $1" -> "$3}')"
    echo "CALI_CONFIG: ${CALI_CONFIG}"
else
    echo "Caliper:     OFF ($([[ "${SAFS_CALI}" == "1" ]] && echo 'binary lacks libcaliper' || echo 'SAFS_CALI=0'))"
fi
echo "EXPERIMENT:  v4_0_0 RSSRW ALT THERMAL CASE1 (coarse), p${ORDER}, PURE UPWIND (Godunov) + ADER (order ${ADER_ORDER})"
echo "Config:      ${CONFIG_TOML}   (order ${_cfg_order}, ader_order ${_cfg_ader})"
echo "Inputs:      mesh=$(readlink ${MESH_LN})"
echo "             vel =$(readlink ${VEL_LN})"
echo "             str =$(readlink ${STR_LN})"
echo "             fri =$(readlink ${FRI_LN})"
echo "Numerics:    --mixed-flux none --time-integrator ader --ader-order ${ADER_ORDER}   perf: ${DERIV_OPT:-<none>}"
echo "tfinal:      ${SAFS_TFINAL:-<config 150s>}   CFL: ${SAFS_CFL:-<config 0.5>}"
echo "Output dir:  ${OUT}"
echo "SPEED METRIC: after the run, read the per-step timing from ${LOG}"
echo "  (steps/s and simulated-s per wall-hour); compare to the SeisSol o4 log for"
echo "  the same mesh.  This is a THROUGHPUT comparison, not a bit-for-bit match."

# SLURM 23.02 does not propagate --cpus-per-task to srun; pure-MPI here (no OMP),
# so 1 cpu/rank is correct.  --unbuffered + stdbuf line-buffers the driver log.
srun --mpi=pmi2 -n "${SLURM_NTASKS}" --unbuffered stdbuf -oL -eL ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mixed-flux none \
    --time-integrator ader \
    --ader-order "${ADER_ORDER}" \
    ${SHARED_MEM_ARG} \
    ${DERIV_OPT} \
    ${TFINAL_ARG} \
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

if grep -qE 'FreeSurf:' "${LOG}"; then
    echo "FreeSurf output: $(grep -m1 -E 'FreeSurf:' "${LOG}")"
else
    echo "WARNING: no 'FreeSurf:' banner in ${LOG} -- free-surface output may be OFF."
fi
# Shared-memory banner check (the parser ignores unknown flags, so verify the
# EFFECT): the driver prints this line from rank 0 when the mode is active.
if [[ "${SAFS_SHARED_MEM}" == "1" ]]; then
    if grep -q 'MPI-3 shared-memory' "${LOG}"; then
        echo "Sidecar shared-mem: $(grep -m1 'MPI-3 shared-memory' "${LOG}")"
    else
        echo "WARNING: no '[sidecar] ... MPI-3 shared-memory' banner in ${LOG} --"
        echo "         the binary may predate the feature; sidecars were likely"
        echo "         REPLICATED per rank (OOM risk at this --mem sizing)."
    fi
fi

# --- Caliper outputs ----------------------------------------------------------
if [[ "${CALI_ACTIVE}" == "1" ]]; then
    echo ""
    echo "=== Caliper per-component breakdown ==="
    if [[ -f "${CALI_REGION_TXT}" ]]; then
        echo "region report -> ${CALI_REGION_TXT}   (Max>>Avg across ranks = imbalance)"
        echo "--- instrumented regions (head) ---"
        head -40 "${CALI_REGION_TXT}"
    else
        echo "WARNING: ${CALI_REGION_TXT} not written — Caliper channel failed; CALI messages:"
        grep -iE 'cali|caliper' "${LOG}" | head
    fi
    [[ -f "${CALI_SPOT}" ]] && echo "spot/.cali    -> ${CALI_SPOT}   (Hatchet/Thicket or Spot UI: per-rank detail)"
fi
if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} -- see ${LOG} (check [derived] lines + any guard ABORT)."; exit ${RC}
fi
echo "driver_rc=${RC} -- done.  Inspect ${LOG}: V_max(t) bounded, on-fault sigma_n stable,"
echo "  and the per-step wall-time for the MFEM-vs-SeisSol speed comparison."
exit ${RC}
