#!/bin/bash
#SBATCH --job-name=a1_merge_tpv104_200m
#SBATCH --account=lbl107
#SBATCH --partition=compute
#SBATCH --constraint="lustre"        # match the Phase-5 / SeisSol decks
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128        # pure-MPI: 256 ranks total — MUST match baseline job 52344266
#SBATCH --cpus-per-task=1
#SBATCH --mem=140000M                # TPV104 has NO sidecars
#SBATCH -t 6:00:00                   # two short (0.6 s) LTS legs + --lts-report + init
#SBATCH -o a1_merge_tpv104_200m.%j.out
#SBATCH -e a1_merge_tpv104_200m.%j.out
#SBATCH --export=ALL

# =============================================================================
# TRACK-A  A1 -- PER-TICK EXCHANGE MERGE A/B.  SCEC TPV104, 200 m production mesh
#   (2.46 M tets), p3 / ADER-O4, pure upwind.  *** EXPANSE (SDSC). ***
#
# THE EXPERIMENT ---------------------------------------------------------------
#   Both legs are IDENTICAL in every way except ONE CLI FLAG:
#     OFF : per-correct seam exchanges  -> 125 rounds/sync at Nc=6 (today)
#     ON  : --lts-merge-exchanges       ->  64 rounds/sync (the A1 merge)
#   Deck is the lambda=1 wiggle-off config (A0's measured 1.370x winner), so A1 is
#   measured ON TOP of A0 rather than tangled with it.
#
# FOUR LEGS, IN TWO PHASES.  The two questions need OPPOSITE output settings, so
# they are run separately rather than compromised into one:
#   PHASE 1 (bitwise, tfinal=A1_TF_BIT, ParaView ON):  produces fault output BYTES
#       to diff.  A1 claims BITWISE IDENTITY, so this is pass/fail, no tolerance.
#       Short window -- identity either holds from the first sync or it does not.
#   PHASE 2 (timing, tfinal=A1_TF_TIME, ParaView OFF): clean MPI_Waitall, no IO
#       contamination (the same cleanliness A0 established; B0 was polluted by
#       ParaView + perf and that is why it is not a trusted denominator).
#
#   *** READ PHASE 1 FIRST.  A speed number from a path that is not bit-identical
#   *** is worth nothing.  The job prints the bitwise verdict before the timing.
#
# LOCAL EVIDENCE ALREADY IN HAND (so this run is a SCALE check, not first contact):
#   both LTS steppers gated ON-vs-OFF by memcmp at 2/3/4 clusters --
#   tests/unit/test_lts_a1_merge_np2.cpp (bulk) 138/207/276 checks all passed and
#   test_lts_a1_merge_fault_np2.cpp (fault-interleave, the stepper THIS deck runs)
#   138/207/276 all passed; collectives 38->30, 102->62, 230->126.
#   What the local fixture CANNOT cover, and this run must: np=256, the LTS-aware
#   partition, the real 6-cluster layout, and rank seams at production scale.
#
# EXPECTED: rounds/rank fall 125->64 = 1.95x; exposed wait falls ~1.8x on A0's
#   measured super-proportional law; end-to-end ~1.24x (3123 -> ~2512 s/sim-s).
#   ALL THREE ARE PROJECTIONS.  The measurement decides.
#
# Overrides (env):
#   A1_TF_BIT   (default 0.05s) bitwise window -- a few syncs is enough.
#   A1_TF_TIME  (default 0.6s)  timing window  -- matches A0's clean baseline.
#   A1_RUN_BITWISE / A1_RUN_TIMING (default 1/1)  toggle either phase.
#   A1_DERIV_OPT (default "--deriv-cache --shared-ck-recursion --face-cache").
#   A1_OUT_BASE (default <folder>/runs)
# =============================================================================
========================================================================

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
JOB_REL="jobs/comm_dev/a1_merge_tpv104_200m_expanse"
P5_REL="jobs/lts_phase5/tpv104_200m_lts_speed_expanse"
CFG_BASE="${P5_REL}/tpv104_200m_lts_rate2.toml"       # lambda=0.63 (default "scan")
A0_REL="jobs/comm_dev/a0_wiggle_tpv104_200m_expanse"
CFG="${A0_REL}/tpv104_200m_lts_wiggleoff.toml"        # lambda=1 (A0 winner)
MESH_LN="${P5_REL}/mesh_tpv104_200m.msh"

for _f in "${CFG}"; do
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


# --- knobs -------------------------------------------------------------------
TF_BIT="${A1_TF_BIT:-0.05s}"
TF_TIME="${A1_TF_TIME:-0.6s}"
RUN_BIT="${A1_RUN_BITWISE:-1}"
RUN_TIME="${A1_RUN_TIMING:-1}"
DERIV_OPT="${A1_DERIV_OPT:---deriv-cache --shared-ck-recursion --face-cache}"
ADER_ORDER=4
FAULT_INTERLEAVE=1

JOB_ABS="${SEAS_MFEM_ROOT}/miniapps/seas/${JOB_REL}"
OUT_ROOT="${A1_OUT_BASE:-${JOB_ABS}/runs}/${SLURM_JOB_ID:-local}"
mkdir -p "${OUT_ROOT}"

echo "Git commit:  $(git rev-parse HEAD 2>/dev/null)"
echo "Size:        ${SLURM_NNODES:-2}N x ${SLURM_NTASKS_PER_NODE:-128}r = ${SLURM_NTASKS} ranks"
echo "EXPERIMENT:  A1 per-tick exchange merge -- OFF vs ON, one CLI flag apart"
echo "Deck:        ${CFG}  (lambda=1, A0's measured winner)"
echo "Phase 1:     BITWISE  tfinal=${TF_BIT}  ParaView ON  (pass/fail, no tolerance)"
echo "Phase 2:     TIMING   tfinal=${TF_TIME} ParaView OFF (clean MPI_Waitall)"
echo "Output root: ${OUT_ROOT}"

if ! ldd "${DRIVER}" 2>/dev/null | grep -qi 'libcaliper'; then
    echo "ERROR: binary lacks libcaliper -- no MPI_Waitall, phase 2 is meaningless."
    echo "       Rebuild: sbatch jobs/comm_dev/a0_wiggle_tpv104_200m_expanse/build_driver_expanse.sbatch"
    exit 4
fi
CALI_BASE="runtime-report(calc.inclusive,region.count,profile.mpi,max_column_width=60"

# run_leg <merge:0|1> <phase:bit|time> <tfinal> <outdir> <log> <cali_stem>
run_leg () {
    local merge="$1" phase="$2" tf="$3" out="$4" log="$5" stem="$6"
    local extra=(); local pv=()
    [[ "${merge}" == "1" ]] && extra+=(--lts-merge-exchanges)
    if [[ "${phase}" == "bit" ]]; then
        # Deterministic, low-cadence fault output: the BYTES we diff.  Identical
        # settings in both legs, so anything it adds cancels in the comparison.
        pv=(--paraview --paraview-fault-vtu --paraview-max-snapshots 4)
        unset CALI_CONFIG
    else
        pv=()   # no IO: clean wait measurement
        export CALI_CONFIG="${CALI_BASE},output=${stem}.cali-region-report.txt)"
    fi
    mkdir -p "${out}"
    echo ""
    echo "[a1] === ${phase} leg, merge=${merge}, tfinal=${tf} ==="
    srun --mpi=pmi2 -n "${SLURM_NTASKS}" \
        --export="ALL,SEAS_LTS_FAULT_INTERLEAVE=${FAULT_INTERLEAVE}" \
        --unbuffered stdbuf -oL -eL "${DRIVER}" \
        --config "${CFG}" --mixed-flux none --time-integrator ader \
        --ader-order "${ADER_ORDER}" --cfl-dg-safety 1.0 ${DERIV_OPT} \
        "${extra[@]}" "${pv[@]}" \
        --tfinal "${tf}" --print-derived --output-dir "${out}" \
        > "${log}" 2>&1
    return $?
}

RC=0
# =============================== PHASE 1: BITWISE =============================
BIT_VERDICT="SKIPPED"
if [[ "${RUN_BIT}" == "1" ]]; then
    B_OFF="${OUT_ROOT}/bit_off"; B_ON="${OUT_ROOT}/bit_on"
    run_leg 0 bit "${TF_BIT}" "${B_OFF}" "${OUT_ROOT}/bit_off.log" "" || RC=$?
    run_leg 1 bit "${TF_BIT}" "${B_ON}"  "${OUT_ROOT}/bit_on.log"  "" || RC=$?
    echo ""
    echo "=== [a1] PHASE 1 RESULT -- BITWISE IDENTITY (read this FIRST) ==="
    for L in bit_off bit_on; do
        printf "  %-8s %s\n" "${L}" \
          "$(grep -m1 -E '^\[lts\] done:' "${OUT_ROOT}/${L}.log" 2>/dev/null || echo '(no [lts] done line -- leg did not finish)')"
    done
    # Byte-compare every fault output file produced by both legs.
    NDIFF=0; NSAME=0; NMISS=0
    while IFS= read -r f; do
        rel="${f#${B_OFF}/}"
        if [[ ! -f "${B_ON}/${rel}" ]]; then NMISS=$((NMISS+1)); continue; fi
        if cmp -s "${f}" "${B_ON}/${rel}"; then NSAME=$((NSAME+1)); else
            NDIFF=$((NDIFF+1)); echo "  DIFFERS: ${rel}"; fi
    done < <(find "${B_OFF}" -type f \( -name '*.vtu' -o -name '*.pvd' -o -name '*.vtkhdf' \) 2>/dev/null | sort)
    echo "  fault-output files: ${NSAME} identical, ${NDIFF} differing, ${NMISS} missing-on-ON"
    # V_max to full printed precision is an independent, IO-free cross-check.
    V_OFF="$(grep -m1 -oE 'V_max_global = [0-9.eE+-]+' "${OUT_ROOT}/bit_off.log" 2>/dev/null | awk '{print $3}')"
    V_ON="$(grep -m1 -oE 'V_max_global = [0-9.eE+-]+' "${OUT_ROOT}/bit_on.log"  2>/dev/null | awk '{print $3}')"
    echo "  V_max_global: OFF=${V_OFF:-?}  ON=${V_ON:-?}"
    if [[ "${NDIFF}" == "0" && "${NMISS}" == "0" && "${NSAME}" -gt 0 && -n "${V_OFF}" && "${V_OFF}" == "${V_ON}" ]]; then
        BIT_VERDICT="PASS"
        echo "  VERDICT: BITWISE IDENTICAL at np=${SLURM_NTASKS} -- A1 is correct at production scale."
    else
        BIT_VERDICT="FAIL"
        echo "  VERDICT: *** NOT BITWISE IDENTICAL -- DO NOT USE THE TIMING BELOW. ***"
        echo "           A1 claims bit-exactness; a difference is a BUG, not a tolerance issue."
    fi
fi

# =============================== PHASE 2: TIMING ==============================
if [[ "${RUN_TIME}" == "1" ]]; then
    T_OFF="${OUT_ROOT}/time_off"; T_ON="${OUT_ROOT}/time_on"
    run_leg 0 time "${TF_TIME}" "${T_OFF}" "${OUT_ROOT}/time_off.log" "${OUT_ROOT}/time_off" || RC=$?
    run_leg 1 time "${TF_TIME}" "${T_ON}"  "${OUT_ROOT}/time_on.log"  "${OUT_ROOT}/time_on"  || RC=$?

    # Largest-Avg MPI_Waitall row (a report holds several; taking the first picks a
    # negligible one -- that bug produced a wrong verdict in job 52416603).
    wait_col () {
        local rpt="$1" col="${2:-3}"; [[ -f "${rpt}" ]] || return 0
        awk -v col="${col}" 'tolower($0) ~ /mpi_waitall/ {
            n=0; s=$0;
            while (match(s,/[0-9]+\.?[0-9]*/)) { n++; tok[n]=substr(s,RSTART,RLENGTH)+0; s=substr(s,RSTART+RLENGTH); }
            if (n>=3 && tok[3]>best) { best=tok[3]; for(i=1;i<=3;i++) b[i]=tok[i] }
        } END { if (best>0) printf "%.6f\n", b[col] }' "${rpt}"
    }
    step_row () { grep -m1 -E 'seas::spatial_dyn::step' "$1" 2>/dev/null | sed 's/  */ /g'; }
    R_OFF="${OUT_ROOT}/time_off.cali-region-report.txt"
    R_ON="${OUT_ROOT}/time_on.cali-region-report.txt"

    echo ""
    echo "=== [a1] PHASE 2 RESULT -- TIMING ==="
    echo "  bitwise verdict carried forward: ${BIT_VERDICT}"
    echo "  --- raw MPI_Waitall rows (verify the ratio by hand) ---"
    for p in "OFF:${R_OFF}" "ON:${R_ON}"; do
        echo "  ${p%%:*}:"; grep -iE 'MPI_Waitall' "${p#*:}" 2>/dev/null | sed 's/^/    /'
    done
    echo "  --- end-to-end wall (seas::spatial_dyn::step) ---"
    echo "  OFF: $(step_row "${R_OFF}")"
    echo "  ON : $(step_row "${R_ON}")"
    WOFF="$(wait_col "${R_OFF}" 3)"; WON="$(wait_col "${R_ON}" 3)"
    SOFF="$(step_row "${R_OFF}" | awk '{print $4}')"; SON="$(step_row "${R_ON}" | awk '{print $4}')"
    echo "  exchange rounds/rank: OFF $(grep -oE '[0-9]+' <<<"$(grep -iE 'MPI_Waitall' "${R_OFF}" 2>/dev/null | tail -1 | awk '{print $(NF-3)}')" | head -1)  ON $(grep -oE '[0-9]+' <<<"$(grep -iE 'MPI_Waitall' "${R_ON}" 2>/dev/null | tail -1 | awk '{print $(NF-3)}')" | head -1)"
    if [[ -n "${WOFF}" && -n "${WON}" ]]; then
        awk -v a="${WOFF}" -v b="${WON}" -v so="${SOFF}" -v sn="${SON}" -v v="${BIT_VERDICT}" 'BEGIN{
            if (b>0) printf "  R_wait = wait(OFF)/wait(ON) = %.3f   (predicted ~1.8 on A0'"'"'s law)\n", a/b;
            if (sn>0 && so>0) printf "  end-to-end = %.3fx   (%.0f -> %.0f s wall; predicted ~1.24x)\n", so/sn, so, sn;
            if (v != "PASS") print "  *** TIMING IS NOT ADMISSIBLE: phase 1 did not PASS. ***";
        }'
    else
        echo "  UNDETERMINED -- read the raw rows above."
    fi
fi

echo ""
echo "[a1] done.  Phase 1 (bitwise) gates Phase 2 (timing).  Reports under ${OUT_ROOT}"
exit ${RC}
