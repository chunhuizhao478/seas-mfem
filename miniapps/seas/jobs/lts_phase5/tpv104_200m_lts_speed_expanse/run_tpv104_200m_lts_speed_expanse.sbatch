#!/bin/bash
#SBATCH --job-name=tpv104_200m_lts_speed
#SBATCH --account=lbl107
#SBATCH --partition=compute
#SBATCH --constraint="lustre"        # run lives on /expanse/lustre (match the SeisSol deck)
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128        # pure-MPI: 128 ranks x 1 core = the whole node
#SBATCH --cpus-per-task=1            # MFEM has no OpenMP
#SBATCH --mem=140000M               # TPV104 has NO sidecars; matches the proven SAFS sizing
#SBATCH -t 24:00:00
#SBATCH -o tpv104_200m_lts_speed.%j.out
#SBATCH -e tpv104_200m_lts_speed.%j.out
#SBATCH --export=ALL

# =============================================================================
# LTS PHASE 5 — REALIZED-SPEEDUP + BENCHMARK deck.  SCEC TPV104, 200 m production
#   mesh (2.46 M tets), p3 / ADER-O4, pure upwind (Godunov) scalar Riemann.
#   *** EXPANSE (SDSC) — srun --mpi=pmi2, toolchain cpu/0.17.3b + gcc/10.2.0
#       + openmpi/4.1.x (matches build_expanse.sh). ***
#
# WHAT THIS MEASURES ----------------------------------------------------------
#   Two runs that differ in EXACTLY ONE config line ([numerics].lts):
#     OFF   : lts="off"    global time stepping (GTS) baseline
#     RATE2 : lts="rate2"  clustered rate-2 LTS (bulk + fault-half interleave)
#   Both use IDENTICAL physics, order, flux, CFL safety, and perf levers, so the
#   ONLY difference is the time-stepping scheme.  From the run logs we read:
#     (1) MFEM-LTS vs MFEM-GTS  sim-seconds / wall-hour  -> the REALIZED speedup
#         (the number the plan projects at ~38.7x ideal, >=15x acceptance gate);
#     (2) the --lts-report PREDICTED harmonic speedup + the [lts-partition]
#         achieved per-cluster / fault-work imbalance (the realization-fraction
#         attribution — bulk ~1.16x, fault ~2.34x measured locally at np=256);
#     (3) MFEM-LTS vs the SeisSol o4 log on the SAME mesh (throughput);
#     (4) MFEM-LTS on-fault stations vs the SCEC TPV104 reference (physics).
#
# WHY 200 m / p3 / ADER-O4 ----------------------------------------------------
#   * 200 m is the SCEC TPV104 production resolution; the 1000 m coarse mesh is
#     too small for 256 ranks (fine clusters have fewer cells than ranks -> the
#     partition imbalance is a pure over-decomposition artifact, not physics).
#   * order 3 / ader_order 4 is order-matched to SeisSol o4 -> the honest,
#     accuracy-matched speed comparison (the plan's ~1.2-1.5x-faster headline;
#     a p1 run would be naturally faster still but is a different accuracy).
#   * perf levers ON for BOTH runs: --deriv-cache (cached M^-1 K_d mat-vec) and
#     --shared-ck-recursion (ADER CK predictor once per macro-step).  This is a
#     THROUGHPUT comparison, not a bit-exact gate, so the levers belong on.
#
# SELF-CONTAINED: both decks + the 200 m mesh symlink live in this folder.  The
#   200 m TPV104 mesh is a committed repo file (tpv104/mesh/tpv104_200m.msh,
#   ~126 MB); TPV104 has NO sidecars (analytic material) so nothing is scp'd.
#
# The binary must be PRE-BUILT (SDSC forbids login-node compiles):
#   bash <root>/build_expanse.sh            (bit-exact / no Caliper)
#   USE_CALIPER=YES bash <root>/build_expanse.sh   (for the per-component report)
#   The binary MUST carry LTS Phase 0-4b (clustering + partition + bulk seams +
#   fault-half interleave + the [lts-partition] imbalance diagnostic).
#
# Overrides (env):
#   P5_TFINAL_GTS   (default 2.0s)  GTS baseline window — only needs a steady
#                                   per-step rate; the full 12 s GTS run is ~6x
#                                   slower than LTS and unnecessary for the rate.
#   P5_TFINAL_LTS   (default 12.0s) LTS window — full TPV104 rupture (benchmark
#                                   stations + rupture speed) AND steady rate.
#   P5_RUN_GTS / P5_RUN_LTS (default 1/1)   toggle either run.
#   P5_DERIV_OPT    (default "--deriv-cache --shared-ck-recursion"; "" = off).
#   P5_CALI         (default 1)     Caliper per-component report (soft-degrades
#                                   if the binary lacks libcaliper).
#   P5_CALI_CONFIG  full CALI_CONFIG override.
#   P5_OUT_BASE     (default <folder>/runs)   output root (job-ID-keyed).
#   P5_FAULT_INTERLEAVE (default 1) SEAS_LTS_FAULT_INTERLEAVE for the RATE2 run.
# =============================================================================

set -o pipefail
export LC_ALL=C
export LANG=C

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
# This sbatch lives at <root>/miniapps/seas/jobs/lts_phase5/<folder>/.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SUBMIT_DIR}}"
_walk="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5 6 7; do
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

# --- runtime shared libs (caliper first — harmless if absent; then hdf5/zfp/petsc)
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"
# HDF5 POSIX file locking FAILS on Lustre (errno 524); matches the SAFS decks.
export HDF5_USE_FILE_LOCKING=FALSE

cd "${SEAS_MFEM_ROOT}"
grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc in config.mk"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS in config.mk"; exit 1; }

cd miniapps/seas
DRIVER=./seas_spatial_dyn_driver
if [[ ! -x "${DRIVER}" ]]; then
    echo "ERROR: ${DRIVER} not found.  Build it FIRST on Expanse (interactive shared session):"
    echo "   bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    echo "   The binary must carry LTS Phase 0-4b (clustering + partition + fault-half interleave)."
    exit 1
fi

# --- inputs (all resolved from miniapps/seas) --------------------------------
JOB_REL="jobs/lts_phase5/tpv104_200m_lts_speed_expanse"
CFG_OFF="${JOB_REL}/tpv104_200m_lts_off.toml"
CFG_RATE2="${JOB_REL}/tpv104_200m_lts_rate2.toml"
MESH_LN="${JOB_REL}/mesh_tpv104_200m.msh"     # the decks reference this via [mesh].path

for _f in "${CFG_OFF}" "${CFG_RATE2}"; do
    [[ -f "${_f}" ]] || { echo "ERROR: required config missing: $(pwd)/${_f}"; exit 1; }
done
# The in-folder mesh symlink and the mesh itself are *.msh (gitignored), so a
# fresh checkout may lack the symlink even when the mesh is present.  Re-create
# it from the known relative target (resolves to miniapps/seas/tpv104/mesh/).
if [[ ! -e "${MESH_LN}" ]]; then
    ln -sf ../../../tpv104/mesh/tpv104_200m.msh "${MESH_LN}" 2>/dev/null || true
fi
if [[ ! -e "${MESH_LN}" ]]; then
    echo "ERROR: mesh missing / dangling symlink: $(pwd)/${MESH_LN} -> $(readlink "${MESH_LN}" 2>/dev/null)"
    echo "       The 200 m TPV104 mesh (tpv104/mesh/tpv104_200m.msh, ~126 MB) is *.msh"
    echo "       and therefore GITIGNORED — it does NOT arrive via git clone/pull.  Stage"
    echo "       it onto Expanse (scp/rsync into tpv104/mesh/) before submitting."
    exit 1
fi

# The two decks MUST differ in EXACTLY one line — the [numerics].lts value.
DIFF_OUT="$(diff "${CFG_OFF}" "${CFG_RATE2}" || true)"
NDIFF="$(printf '%s\n' "${DIFF_OUT}" | grep -cE '^[<>]' || true)"
if [[ "${NDIFF}" != "2" ]]; then
    echo "ERROR: the two decks differ in more than the lts line (${NDIFF} changed"
    echo "       lines; expected 2 = one '<' + one '>').  diff:"
    printf '%s\n' "${DIFF_OUT}"
    exit 2
fi
if ! printf '%s\n' "${DIFF_OUT}" | grep -qE '^[<>].*lts[[:space:]]*='; then
    echo "ERROR: the one differing line is not the lts line — refusing to run."
    printf '%s\n' "${DIFF_OUT}"
    exit 2
fi
echo "[phase5] deck diff OK — the two decks differ only in [numerics].lts."

# --- knobs -------------------------------------------------------------------
TFINAL_GTS="${P5_TFINAL_GTS:-2.0s}"
TFINAL_LTS="${P5_TFINAL_LTS:-12.0s}"
RUN_GTS="${P5_RUN_GTS:-1}"
RUN_LTS="${P5_RUN_LTS:-1}"
# Kernel-plan Phase 0: --face-cache added to the default (it was OFF in the
# original Phase-5 runs — the ~15% Vandermonde-LU + ~6% CalcShape interior-face
# tax was paid in full).  Composes with the deriv levers (only exclusivity is
# UsePrecomputedFaceFluxes).  Override P5_DERIV_OPT to change.
DERIV_OPT="${P5_DERIV_OPT:---deriv-cache --shared-ck-recursion --face-cache}"
FAULT_INTERLEAVE="${P5_FAULT_INTERLEAVE:-1}"
ADER_ORDER=4

# --- output dirs (self-contained under the folder; job-ID-keyed) -------------
JOB_ABS="${SEAS_MFEM_ROOT}/miniapps/seas/${JOB_REL}"
OUT_BASE="${P5_OUT_BASE:-${JOB_ABS}/runs}"
OUT_ROOT="${OUT_BASE}/${SLURM_JOB_ID:-local}"
OUT_OFF="${OUT_ROOT}/gts_off"
OUT_RATE2="${OUT_ROOT}/lts_rate2"
LOG_OFF="${OUT_ROOT}/run_gts_off.log"
LOG_RATE2="${OUT_ROOT}/run_lts_rate2.log"
LOG_REPORT="${OUT_ROOT}/lts_report.log"
mkdir -p "${OUT_OFF}" "${OUT_RATE2}"

echo "Git commit:  $(git rev-parse HEAD 2>/dev/null)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Size:        compute (lustre), ${SLURM_NNODES:-2}N x ${SLURM_NTASKS_PER_NODE:-128}r = ${SLURM_NTASKS} ranks (pure MPI), ${SLURM_MEM_PER_NODE:-240000}M/node, ${SBATCH_TIMELIMIT:-24:00:00}"
echo "EXPERIMENT:  TPV104 200 m, p3 / ADER-O4, PURE UPWIND — LTS(rate2) vs GTS(off)"
echo "Decks:       OFF=${CFG_OFF}"
echo "             RATE2=${CFG_RATE2}   (differ only in [numerics].lts)"
echo "Mesh:        $(readlink "${MESH_LN}")"
echo "Numerics:    --mixed-flux none --time-integrator ader --ader-order ${ADER_ORDER} --cfl-dg-safety 1.0   perf: ${DERIV_OPT:-<none>}"
echo "dt NOTE:     Courant cfl=0.5 in both decks — MATCHED to SeisSol's CFL=0.5 (identical h"
echo "             formula + order factor 7 => identical dt, same step count per sim-s)."
echo "             Runs before commit $(git rev-parse --short HEAD 2>/dev/null) (52332649 GTS, 52340525 LTS) used cfl=0.25:"
echo "             their LTS-vs-GTS RATIO is valid; their sim-s/h convert x2 for cross-code use."
echo "Windows:     GTS tfinal=${TFINAL_GTS} (steady rate), LTS tfinal=${TFINAL_LTS} (full rupture)"
echo "Fault ILV:   RATE2 SEAS_LTS_FAULT_INTERLEAVE=${FAULT_INTERLEAVE}"
echo "Output root: ${OUT_ROOT}"
echo "SPEED METRIC: read simulated-s / wall-hour from each run log; the realized"
echo "  speedup = (LTS sim-s/hr) / (GTS sim-s/hr).  Compare LTS to the SeisSol o4"
echo "  log for the same mesh.  THROUGHPUT comparison, not a bit-for-bit match."

# -----------------------------------------------------------------------------
# 0. --lts-report PRE-STEP (np=1): PREDICTED harmonic speedup + the hypothetical
#    ${SLURM_NTASKS}-rank partition imbalance (the realization-fraction predictor)
#    WITHOUT paying the full run.  Report-and-exit (--lts-report implies dry-run).
# -----------------------------------------------------------------------------
echo ""
echo "=== [phase5] --lts-report (predicted speedup + ${SLURM_NTASKS}-rank partition imbalance) ==="
srun --mpi=pmi2 -n 1 --unbuffered stdbuf -oL -eL "${DRIVER}" \
    --config "${CFG_RATE2}" \
    --mixed-flux none --time-integrator ader --ader-order "${ADER_ORDER}" \
    --lts-report --lts-report-nparts "${SLURM_NTASKS}" \
    > "${LOG_REPORT}" 2>&1 || true
grep -E "lts-report|lts-partition|imbalance|per-cluster|speedup|histogram|Nc =" "${LOG_REPORT}" || \
    { echo "WARNING: --lts-report produced no [lts-report] lines; see ${LOG_REPORT}"; tail -5 "${LOG_REPORT}"; }

# --- Caliper per-component breakdown (P5_CALI=1 default; soft-degrade) --------
P5_CALI="${P5_CALI:-1}"
CALI_ACTIVE=0
if [[ "${P5_CALI}" == "1" ]]; then
    if ldd "${DRIVER}" 2>/dev/null | grep -qi 'libcaliper'; then
        CALI_ACTIVE=1
        echo "Caliper:     ON  $(ldd "${DRIVER}" 2>/dev/null | grep -m1 -i libcaliper | awk '{print $1" -> "$3}')"
    else
        echo "Caliper:     OFF (binary lacks libcaliper; rebuild USE_CALIPER=YES for the per-component report)"
    fi
else
    echo "Caliper:     OFF (P5_CALI=0)"
fi

# run_one: one timed srun over all ranks.  $6 = extra --export var (fault ILV;
#          comma-free so it is safe on srun --export), $7 = per-run Caliper
#          output stem (empty => no Caliper).  CALI_CONFIG contains commas/parens,
#          so it is EXPORTED to the shell env (propagated by --export=ALL) rather
#          than passed on --export (which is comma-separated and would corrupt it).
run_one () {
    local cfg="$1" out="$2" log="$3" tfinal="$4" label="$5" export_extra="$6" cali_stem="$7"
    local export_arg="ALL"
    [[ -n "${export_extra}" ]] && export_arg="ALL,${export_extra}"
    if [[ "${CALI_ACTIVE}" == "1" && -n "${cali_stem}" ]]; then
        export CALI_CONFIG="${P5_CALI_CONFIG:-runtime-report(calc.inclusive,region.count,profile.mpi,max_column_width=60,output=${cali_stem}.cali-region-report.txt),spot(output=${cali_stem}.cali)}"
    else
        unset CALI_CONFIG
    fi
    echo ""
    echo "[phase5] === run ${label}: cfg=${cfg} tfinal=${tfinal} out=${out} ==="
    srun --mpi=pmi2 -n "${SLURM_NTASKS}" --export="${export_arg}" --unbuffered stdbuf -oL -eL "${DRIVER}" \
        --config "${cfg}" \
        --mixed-flux none \
        --time-integrator ader \
        --ader-order "${ADER_ORDER}" \
        --cfl-dg-safety 1.0 \
        ${DERIV_OPT} \
        --tfinal "${tfinal}" \
        --paraview \
        --paraview-fault-vtu \
        --paraview-max-snapshots 300 \
        --print-derived \
        --checkpoint-every 4000 \
        --output-dir "${out}" \
        > "${log}" 2>&1
    return $?
}

RC_OFF=0; RC_RATE2=0

if [[ "${RUN_GTS}" == "1" ]]; then
    run_one "${CFG_OFF}" "${OUT_OFF}" "${LOG_OFF}" "${TFINAL_GTS}" \
            "GTS (lts=off, baseline)" "" "${LOG_OFF%.log}"
    RC_OFF=$?
    echo "[phase5] GTS rc=${RC_OFF}"
fi

if [[ "${RUN_LTS}" == "1" ]]; then
    run_one "${CFG_RATE2}" "${OUT_RATE2}" "${LOG_RATE2}" "${TFINAL_LTS}" \
            "LTS (lts=rate2, clustered)" "SEAS_LTS_FAULT_INTERLEAVE=${FAULT_INTERLEAVE}" \
            "${LOG_RATE2%.log}"
    RC_RATE2=$?
    echo "[phase5] LTS rc=${RC_RATE2}"
fi

# -----------------------------------------------------------------------------
# Post-run verification + speed extraction.
# -----------------------------------------------------------------------------
post_checks () {
    local log="$1" label="$2"
    [[ -f "${log}" ]] || { echo "  (${label}: no log)"; return; }
    echo "--- ${label}: ${log} ---"
    # cfl_dg_safety EFFECT: the [time] line must report safety = 1 (SeisSol dt).
    if grep -q 'cfl_dg_safety = 1 ' "${log}"; then
        echo "  dt: $(grep -m1 -E '^\[time\] dt ' "${log}")  (cfl_dg_safety=1)"
    elif grep -q 'cfl_dg_safety' "${log}"; then
        echo "  ERROR: run did NOT use cfl_dg_safety=1.0 -> timing is inflated:"
        grep -m1 'cfl_dg_safety' "${log}"
    fi
    # LTS activation banners (RATE2 only): partition inject + imbalance.
    grep -m1 -E '\[lts\] LTS-aware partition injected' "${log}" | sed 's/^/  /' || true
    grep -E '\[lts-partition\]' "${log}" | sed 's/^/  /' || true
    # Best-effort speed summary — the driver prints per-step / throughput timing.
    grep -iE 'sim.?s.*(hour|hr)|steps?/s|wall.?time|throughput|per-step' "${log}" | tail -4 | sed 's/^/  /' || true
}

echo ""
echo "=== [phase5] speed + activation summary ==="
[[ "${RUN_GTS}" == "1" ]] && post_checks "${LOG_OFF}"   "GTS (off)"
[[ "${RUN_LTS}" == "1" ]] && post_checks "${LOG_RATE2}" "LTS (rate2)"

# --- Caliper outputs ---------------------------------------------------------
if [[ "${CALI_ACTIVE}" == "1" ]]; then
    for stem in "${LOG_OFF%.log}" "${LOG_RATE2%.log}"; do
        rpt="${stem}.cali-region-report.txt"
        if [[ -f "${rpt}" ]]; then
            echo ""
            echo "=== Caliper region report: ${rpt} (Max>>Avg across ranks = imbalance) ==="
            head -40 "${rpt}"
        fi
    done
fi

echo ""
if [[ ${RC_OFF} -ne 0 || ${RC_RATE2} -ne 0 ]]; then
    echo "ERROR: a run failed (GTS rc=${RC_OFF}, LTS rc=${RC_RATE2}).  Inspect the logs"
    echo "       under ${OUT_ROOT} for [derived] lines / guard ABORTs."
    exit 3
fi
echo "[phase5] done.  Next: compute (LTS sim-s/hr)/(GTS sim-s/hr) = realized speedup,"
echo "  compare LTS throughput to the SeisSol o4 log, and diff the LTS on-fault"
echo "  stations against the SCEC TPV104 reference (see README.md)."
exit 0
