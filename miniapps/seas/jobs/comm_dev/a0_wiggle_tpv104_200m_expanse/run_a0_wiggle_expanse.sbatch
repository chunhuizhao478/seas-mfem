#!/bin/bash
#SBATCH --job-name=a0_wiggle_tpv104_200m
#SBATCH --account=lbl107
#SBATCH --partition=compute
#SBATCH --constraint="lustre"        # match the Phase-5 / SeisSol decks
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128        # pure-MPI: 256 ranks total — MUST match baseline job 52344266
#SBATCH --cpus-per-task=1
#SBATCH --mem=140000M                # TPV104 has NO sidecars
#SBATCH -t 4:00:00                   # two short (0.6 s) LTS legs + --lts-report + init
#SBATCH -o a0_wiggle_tpv104_200m.%j.out
#SBATCH -e a0_wiggle_tpv104_200m.%j.out
#SBATCH --export=ALL

# =============================================================================
# TRACK-A  A0 — WIGGLE A/B.  SCEC TPV104, 200 m production mesh (2.46 M tets),
#   p3 / ADER-O4, pure upwind (Godunov).  *** EXPANSE (SDSC). ***
#
# THE EXPERIMENT (the highest-value one-line test in the plan) -----------------
#   Two clustered-LTS runs that differ in EXACTLY ONE config line
#   ([numerics].lts_wiggle):
#     RATE2     : lts_wiggle unset  -> default "scan"  -> lambda=0.63  (baseline)
#     WIGGLEOFF : lts_wiggle="off"                     -> lambda=1
#   lts_nc_cap defaults to 6 in BOTH decks and the lambda=1 nc_cap merge is
#   verified to fire, so BOTH cluster to Nc=6.  Rounds/sync = 125 is IDENTICAL;
#   the sync LENGTH T_s is proportional to lambda, so the ONLY thing that changes
#   is exchange rounds per SIM-SECOND, which falls exactly 1/0.63 = 1.587x.
#
# WHAT A0 DECIDES -------------------------------------------------------------
#   (a) MEASURES the wiggle lever (rounds/sim-s down 1.587x — near-certain by
#       construction; confirmed here as a sanity check on Nc and the sync count).
#   (b) TESTS THE ASSUMPTION THE WHOLE COMMUNICATION TRACK RESTS ON: whether the
#       exposed MPI_Waitall time scales with the exchange-round rate at all.  NO
#       run has ever varied the round rate; the comm plan's payoff was fitted to
#       a hypothetical.  If exposed wait (s/sim-s) falls ~1.587x from RATE2 to
#       WIGGLEOFF, the round-merge levers (A1-A3) are real; if it barely moves,
#       the wait is dominated by per-round latency / load skew and A1-A3 are
#       worth far less than believed — and we learn it in one job, not weeks.
#   lambda=1 is also SeisSol's own configuration on this benchmark.
#
# SCOPE / WHAT A0 DOES NOT DO -------------------------------------------------
#   A0 measures the BULK seam-exchange wire scaling, which is a STEADY-STATE
#   property, so a short PRE-rupture window (tfinal=0.6 s; nucleation is at
#   t_0=1.0 s) is sufficient and cheap.  A0 does NOT measure fault/friction skew
#   (that needs the rupture window — see the plan's known-unknowns 4 & 5).  The
#   earlier 24% skew proxy was taken at tfinal=0.5 s, pre-nucleation, and is NOT
#   a bound on fault skew; A0 does not attempt to fix that.
#
# WHY THIS IS A CLEAN A/B -----------------------------------------------------
#   * SAME job, SAME 256 ranks, SAME flags, SAME (short) window, SAME binary.
#   * ParaView OUTPUT IS OFF for both legs: the B0 measurement was contaminated
#     by ParaView + perf; A0 removes that.  Caliper runtime-report is ON for both
#     (it is A0's instrument) and, being identical in both legs, cancels in the
#     RATIO wait(lambda=0.63)/wait(lambda=1) that A0 actually reads.
#   * A0 is a RATIO test — anything identical in both legs cancels.
#
# The binary must be PRE-BUILT WITH CALIPER (SDSC forbids login-node compiles)
# — without libcaliper there is NO MPI_Waitall number and A0 answers nothing:
#   USE_CALIPER=YES bash <root>/build_expanse.sh
#   The binary MUST carry LTS Phase 0-4b AND commit 884349a (the LTS step
#   Caliper region — without it MPI_Waitall reports as %-of-program, not %-of-step).
#
# Overrides (env):
#   A0_TFINAL   (default 0.6s)  short pre-rupture window; enough for ~50-80 syncs.
#   A0_NODES-equivalent is fixed at 256 ranks (do NOT change — comm share is
#               rank-count dependent and must match the 52344266 baseline).
#   A0_DERIV_OPT (default "--deriv-cache --shared-ck-recursion --face-cache").
#               face-cache is LTS-inert (proven: it only reaches the GTS
#               corrector) so it does not perturb this LTS A/B; kept for parity
#               with the Phase-5 decks.
#   A0_RUN_BASELINE / A0_RUN_WIGGLEOFF (default 1/1)  toggle either leg.
#   A0_OUT_BASE (default <folder>/runs)   output root (job-ID-keyed).
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
# This sbatch lives at <root>/miniapps/seas/jobs/comm_dev/<folder>/.
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

# --- runtime shared libs (caliper first — REQUIRED here; then hdf5/zfp/petsc) -
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"
export HDF5_USE_FILE_LOCKING=FALSE

cd "${SEAS_MFEM_ROOT}"
grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc in config.mk"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS in config.mk"; exit 1; }

cd miniapps/seas
DRIVER=./seas_spatial_dyn_driver
if [[ ! -x "${DRIVER}" ]]; then
    echo "ERROR: ${DRIVER} not found.  Build it FIRST on Expanse (interactive shared session):"
    echo "   USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    exit 1
fi

# --- inputs ------------------------------------------------------------------
JOB_REL="jobs/comm_dev/a0_wiggle_tpv104_200m_expanse"
P5_REL="jobs/lts_phase5/tpv104_200m_lts_speed_expanse"
CFG_BASE="${P5_REL}/tpv104_200m_lts_rate2.toml"       # lambda=0.63 (default "scan")
CFG_WOFF="${JOB_REL}/tpv104_200m_lts_wiggleoff.toml"  # lambda=1   (lts_wiggle="off")
MESH_LN="${P5_REL}/mesh_tpv104_200m.msh"

for _f in "${CFG_BASE}" "${CFG_WOFF}"; do
    [[ -f "${_f}" ]] || { echo "ERROR: required config missing: $(pwd)/${_f}"; exit 1; }
done
# Re-create the mesh symlink if a fresh checkout dropped it (mesh is gitignored).
if [[ ! -e "${MESH_LN}" ]]; then
    ln -sf ../../../tpv104/mesh/tpv104_200m.msh "${MESH_LN}" 2>/dev/null || true
fi
if [[ ! -e "${MESH_LN}" ]]; then
    echo "ERROR: mesh missing / dangling symlink: $(pwd)/${MESH_LN}"
    echo "       The 200 m TPV104 mesh (tpv104/mesh/tpv104_200m.msh, ~126 MB) is *.msh"
    echo "       and GITIGNORED — stage it onto Expanse before submitting."
    exit 1
fi

# The two decks MUST differ in EXACTLY one NON-COMMENT line — the added
# [numerics].lts_wiggle value.  The two decks also carry different HEADER COMMENT
# blocks (each describes its own purpose), so comment-only '<'/'>' lines are
# expected and ignored; only the non-comment body must be a single-line add.
DIFF_OUT="$(diff "${CFG_BASE}" "${CFG_WOFF}" || true)"
ADD_NC="$(printf '%s\n' "${DIFF_OUT}" | grep -E '^>' | grep -vE '^>[[:space:]]*#')"
DEL_NC="$(printf '%s\n' "${DIFF_OUT}" | grep -E '^<' | grep -vE '^<[[:space:]]*#')"
NADD_NC="$(printf '%s\n' "${ADD_NC}" | grep -cE '.')"
NDEL_NC="$(printf '%s\n' "${DEL_NC}" | grep -cE '.')"
# The single added non-comment line must itself be the lts_wiggle assignment.
ADD_IS_WIGGLE="$(printf '%s\n' "${ADD_NC}" | grep -cE '^>[[:space:]]*lts_wiggle[[:space:]]*=')"
if [[ "${NADD_NC}" != "1" || "${NDEL_NC}" != "0" || "${ADD_IS_WIGGLE}" != "1" ]]; then
    echo "ERROR: the two decks differ by more than the single lts_wiggle line."
    echo "       non-comment added=${NADD_NC} (want 1), non-comment removed=${NDEL_NC} (want 0),"
    echo "       added line is lts_wiggle=${ADD_IS_WIGGLE} (want 1).  non-comment diff:"
    printf '  ADD: %s\n' "${ADD_NC}"
    printf '  DEL: %s\n' "${DEL_NC}"
    exit 2
fi
echo "[a0] deck diff OK — WIGGLEOFF adds exactly one non-comment line: lts_wiggle=\"off\"."

# --- knobs -------------------------------------------------------------------
TFINAL="${A0_TFINAL:-0.6s}"
RUN_BASE="${A0_RUN_BASELINE:-1}"
RUN_WOFF="${A0_RUN_WIGGLEOFF:-1}"
DERIV_OPT="${A0_DERIV_OPT:---deriv-cache --shared-ck-recursion --face-cache}"
ADER_ORDER=4
FAULT_INTERLEAVE=1

# --- output dirs -------------------------------------------------------------
JOB_ABS="${SEAS_MFEM_ROOT}/miniapps/seas/${JOB_REL}"
OUT_BASE="${A0_OUT_BASE:-${JOB_ABS}/runs}"
OUT_ROOT="${OUT_BASE}/${SLURM_JOB_ID:-local}"
OUT_M063="${OUT_ROOT}/lambda063_baseline"
OUT_M100="${OUT_ROOT}/lambda100_wiggleoff"
LOG_M063="${OUT_ROOT}/run_lambda063_baseline.log"
LOG_M100="${OUT_ROOT}/run_lambda100_wiggleoff.log"
mkdir -p "${OUT_M063}" "${OUT_M100}"

echo "Git commit:  $(git rev-parse HEAD 2>/dev/null)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Size:        compute (lustre), ${SLURM_NNODES:-2}N x ${SLURM_NTASKS_PER_NODE:-128}r = ${SLURM_NTASKS} ranks (pure MPI)"
echo "EXPERIMENT:  TPV104 200 m A0 wiggle A/B — RATE2(lambda=0.63) vs WIGGLEOFF(lambda=1), both LTS"
echo "Decks:       BASE=${CFG_BASE}"
echo "             WOFF=${CFG_WOFF}   (differ only in lts_wiggle)"
echo "Window:      tfinal=${TFINAL} (pre-rupture; ~50-80 syncs) — measures BULK wire scaling only"
echo "ParaView:    OFF for both legs (clean MPI_Waitall; no IO contamination)"
echo "Output root: ${OUT_ROOT}"

# --- Caliper is REQUIRED (A0's instrument) -----------------------------------
if ! ldd "${DRIVER}" 2>/dev/null | grep -qi 'libcaliper'; then
    echo "ERROR: the binary lacks libcaliper — A0 CANNOT measure exposed MPI_Waitall,"
    echo "       which is its entire purpose.  Rebuild: USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    exit 4
fi
echo "Caliper:     ON  $(ldd "${DRIVER}" 2>/dev/null | grep -m1 -i libcaliper | awk '{print $1" -> "$3}')"
# profile.mpi gives per-MPI-call inclusive time (MPI_Waitall = the exposed wait).
CALI_BASE="runtime-report(calc.inclusive,region.count,profile.mpi,max_column_width=60"

# -----------------------------------------------------------------------------
# 0. --lts-report PRE-STEP (np=1): prints BOTH clusterings (RAW lambda=1/Nc=7 and
#    PRODUCTION scan lambda=0.63/Nc=6).  NOTE: its RAW mode uses nc_cap=0 => Nc=7,
#    which is NOT the runtime lambda=1 config (nc_cap=6 => Nc=6).  So the true
#    runtime round counts come from the [lts] ENABLED banners in the timed legs
#    below, not from this report.  Kept as an independent cross-check only.
# -----------------------------------------------------------------------------
echo ""
echo "=== [a0] --lts-report (cross-check clustering; timed legs are authoritative) ==="
LOG_REPORT="${OUT_ROOT}/lts_report.log"
srun --mpi=pmi2 -n 1 --unbuffered stdbuf -oL -eL "${DRIVER}" \
    --config "${CFG_BASE}" \
    --mixed-flux none --time-integrator ader --ader-order "${ADER_ORDER}" \
    --lts-report --lts-report-nparts "${SLURM_NTASKS}" \
    > "${LOG_REPORT}" 2>&1 || true
grep -E "lts-report|Nc =|histogram|speedup" "${LOG_REPORT}" || \
    { echo "WARNING: --lts-report produced no lines; see ${LOG_REPORT}"; tail -5 "${LOG_REPORT}"; }

# run_leg: one timed srun over all ranks with Caliper on and NO ParaView.
run_leg () {
    local cfg="$1" out="$2" log="$3" label="$4" cali_stem="$5"
    export CALI_CONFIG="${CALI_BASE},output=${cali_stem}.cali-region-report.txt),spot(output=${cali_stem}.cali)"
    echo ""
    echo "[a0] === run ${label}: cfg=${cfg} tfinal=${TFINAL} out=${out} ==="
    srun --mpi=pmi2 -n "${SLURM_NTASKS}" --export="ALL,SEAS_LTS_FAULT_INTERLEAVE=${FAULT_INTERLEAVE}" \
        --unbuffered stdbuf -oL -eL "${DRIVER}" \
        --config "${cfg}" \
        --mixed-flux none \
        --time-integrator ader \
        --ader-order "${ADER_ORDER}" \
        --cfl-dg-safety 1.0 \
        ${DERIV_OPT} \
        --tfinal "${TFINAL}" \
        --print-derived \
        --output-dir "${out}" \
        > "${log}" 2>&1
    return $?
}

RC_M063=0; RC_M100=0
if [[ "${RUN_BASE}" == "1" ]]; then
    run_leg "${CFG_BASE}" "${OUT_M063}" "${LOG_M063}" "BASELINE lambda=0.63 (scan)" "${LOG_M063%.log}"
    RC_M063=$?; echo "[a0] baseline rc=${RC_M063}"
fi
if [[ "${RUN_WOFF}" == "1" ]]; then
    run_leg "${CFG_WOFF}" "${OUT_M100}" "${LOG_M100}" "WIGGLEOFF lambda=1" "${LOG_M100%.log}"
    RC_M100=$?; echo "[a0] wiggleoff rc=${RC_M100}"
fi

# -----------------------------------------------------------------------------
# A0 RESULT EXTRACTION — the A/B the whole job exists to produce.
# -----------------------------------------------------------------------------
# Pull, from each leg: the [lts] ENABLED banner (Nc + lambda), and the
# MPI_Waitall inclusive Avg from the Caliper profile.mpi report.  Then print the
# two ratios that decide A0:
#   R_rounds = rounds/sim-s(0.63) / rounds/sim-s(1.0)   (expect ~1.587x)
#   R_wait   = exposed_wait(0.63) / exposed_wait(1.0)   (THE test: ~1.587x => wait
#              scales with rounds; ~1.0 => it does not, and A1-A3 are devalued)
extract_lts_banner () {  # $1 = log
    grep -m1 -E '^\[lts\] ENABLED' "$1" | sed 's/^/  /' || echo "  (no [lts] ENABLED banner)"
}
# Caliper runtime-report is FIXED-WIDTH with spaces in region names, and its
# column set/order depends on the CALI_CONFIG (here: Min, Max, Avg time/rank,
# then Time %, plus a region Count).  Rather than guess a column index (fragile
# and unvalidated off-machine), surface the HEADER + the raw MPI_Waitall row(s)
# verbatim so the numbers are self-identifying, and compute the ratio from the
# "Avg time/rank" column located BY THE HEADER, falling back to hand-analysis.
report_header () {  # $1 = report ; echoes the column-header line
    [[ -f "$1" ]] || return 0
    grep -m1 -iE 'time/rank|Avg|Time[[:space:]]*%|Path' "$1"
}
waitall_rows () {   # $1 = report ; echoes verbatim MPI_Waitall row(s)
    [[ -f "$1" ]] || return 0
    grep -iE 'MPI_Waitall' "$1"
}
# Avg extraction: Caliper runtime-report's time columns are, in fixed order,
# Min / Max / Avg time/rank / Time %  (verified against
# document/code_optimization_dev/profiling_analysis_2026-06-24.md, where the
# MPI_Waitall row reads 0.42 / 42.10 / 14.76 / 1.77 = Min / Max / Avg / Time%).
# The region name "MPI_Waitall" contains no digits, so ALL numeric tokens on that
# row are the data columns and Avg is the 3rd numeric token.  A trailing
# region.count column (if present) comes after Time% and does not shift Avg.
# A report contains SEVERAL MPI_Waitall rows at different nesting depths (job
# 52416603 had a negligible 11-call row AND the real 20336-call row worth 52.6 %
# of wall).  Taking the FIRST match picked the wrong one and produced a wrong
# verdict.  Select the row with the LARGEST Avg instead.  $2 selects which column
# to echo: 3=Avg (default), 1=Min, 2=Max — Min/Max give the rank skew.
waitall_col () {    # $1 = report, $2 = column index (1=Min 2=Max 3=Avg); echoes value or empty
    local rpt="$1" col="${2:-3}"; [[ -f "${rpt}" ]] || return 0
    awk -v col="${col}" 'tolower($0) ~ /mpi_waitall/ {
        n=0; s=$0;
        while (match(s, /[0-9]+\.?[0-9]*/)) {
            n++; tok[n]=substr(s, RSTART, RLENGTH)+0; s=substr(s, RSTART+RLENGTH);
        }
        if (n >= 3 && tok[3] > bestavg) { bestavg=tok[3]; for(i=1;i<=3;i++) best[i]=tok[i] }
    } END { if (bestavg > 0) printf "%.6f\n", best[col] }' "${rpt}"
}
waitall_avg () { waitall_col "$1" 3; }

echo ""
echo "=== [a0] RESULT — wiggle A/B (the round-scaling test) ==="
RPT_M063="${LOG_M063%.log}.cali-region-report.txt"
RPT_M100="${LOG_M100%.log}.cali-region-report.txt"
echo "--- BASELINE lambda=0.63 ---"; extract_lts_banner "${LOG_M063}"
echo "--- WIGGLEOFF lambda=1.0 ---"; extract_lts_banner "${LOG_M100}"
echo ""
echo "--- raw MPI_Waitall rows (READ THESE; the ratio below is best-effort) ---"
echo "  [header] $(report_header "${RPT_M063}")"
echo "  lambda=0.63:"; waitall_rows "${RPT_M063}" | sed 's/^/    /'
echo "  lambda=1.00:"; waitall_rows "${RPT_M100}" | sed 's/^/    /'
# Physics sanity + the run's own step totals (the honest end-to-end wall).
echo ""
echo "--- physics / completion sanity ---"
for pair in "0.63:${LOG_M063}" "1.00:${LOG_M100}"; do
    lam="${pair%%:*}"; lg="${pair#*:}"
    printf "  lambda=%s: %s\n" "${lam}" "$(grep -m1 -E '^\[lts\] done:' "${lg}" 2>/dev/null || echo '(no [lts] done line — run may not have completed)')"
    grep -m1 -iE 'nan|inf|ABORT|diverg' "${lg}" 2>/dev/null | sed 's/^/    WARNING: /' || true
done
echo "--- end-to-end wall (Caliper seas::spatial_dyn::step, Avg time/rank) ---"
for pair in "0.63:${RPT_M063}" "1.00:${RPT_M100}"; do
    lam="${pair%%:*}"; rp="${pair#*:}"
    printf "  lambda=%s: %s\n" "${lam}" "$(grep -m1 -E 'seas::spatial_dyn::step' "${rp}" 2>/dev/null | sed 's/  */ /g' || echo '(no step region)')"
done

W063="$(waitall_avg "${RPT_M063}")"
W100="$(waitall_avg "${RPT_M100}")"
# Rank skew through the rupture window (retires the pre-nucleation 24% proxy).
for pair in "0.63:${RPT_M063}" "1.00:${RPT_M100}"; do
    lam="${pair%%:*}"; rp="${pair#*:}"
    mn="$(waitall_col "${rp}" 1)"; mx="$(waitall_col "${rp}" 2)"; av="$(waitall_col "${rp}" 3)"
    if [[ -n "${mn}" && -n "${av}" ]]; then
        awk -v l="${lam}" -v mn="${mn}" -v mx="${mx}" -v av="${av}" 'BEGIN{
            printf "  lambda=%s MPI_Waitall  Min %.1f / Avg %.1f / Max %.1f s   skew (Avg-Min)/Avg = %.1f%%   Max/Avg = %.2f\n",
                   l, mn, av, mx, 100*(av-mn)/av, mx/av }'
    fi
done
echo ""
echo "Best-effort MPI_Waitall Avg (largest-Avg row; VERIFY against the rows above):  lambda=0.63 -> ${W063:-?} s ;  lambda=1.0 -> ${W100:-?} s"
if [[ -n "${W063}" && -n "${W100}" ]] && awk -v b="${W100}" 'BEGIN{exit !(b>0)}'; then
    awk -v a="${W063}" -v b="${W100}" 'BEGIN{
        d=a/b;
        printf "R_wait = wait(0.63)/wait(1.0) = %.3f  (expect ~1.587 IF exposed wait scales with the round rate)\n", d;
        if (d >= 1.45)      print "VERDICT (best-effort): exposed wait SCALES with round rate -> A1-A3 merge/overlap are real. Track A GO.";
        else if (d >= 1.15) print "VERDICT (best-effort): PARTIAL scaling -> A1-A3 worth less than 1.587x; re-size the payoff before committing.";
        else                print "VERDICT (best-effort): exposed wait does NOT track rounds -> latency/skew-bound; A1-A3 devalued, revisit Track A.";
    }'
    echo "  ^ CONFIRM this ratio by hand from the raw rows before acting on it."
else
    echo "VERDICT: UNDETERMINED from auto-parse — read the raw MPI_Waitall rows above and compute R_wait by hand."
fi
echo "NOTE: full Caliper reports (all MPI calls, Min/Max/Avg spread) are at:"
echo "      ${RPT_M063}"
echo "      ${RPT_M100}"

echo ""
if [[ ${RC_M063} -ne 0 || ${RC_M100} -ne 0 ]]; then
    echo "ERROR: a leg failed (baseline rc=${RC_M063}, wiggleoff rc=${RC_M100}).  Inspect logs under ${OUT_ROOT}."
    exit 3
fi
echo "[a0] done.  Read R_wait above; it is the go/no-go for the Track-A merge levers."
