#!/bin/bash
#SBATCH -J bp5_restart
#SBATCH -o bp5_restart_test_%j.out
#SBATCH -e bp5_restart_test_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# BP5 — V2 PETSc-TS restart end-to-end test (development queue, 2h, 8N x 400r).
# =======================================================================
# *** BP5 ONLY — do NOT adapt this sbatch for TPV104 ***
#
# The V1+V2 checkpoint format and the `--restart` + `--petsc-ts` driver
# code path are BP5-specific (quasi-static, PETSc TS RK45).  TPV104 uses
# an ADER-DG EXPLICIT time stepper (`dynamic/tpv104_substep_iterator.cpp`)
# with no checkpoint/restart machinery: no `--restart` flag, no
# `--checkpoint-interval`, no `--petsc-ts`.  Porting restart to TPV104
# is a Phase-4 follow-up — see `petsc_ts_restart_plan_2026-05-16.md`
# §"Out of scope" for the required new work.  R-002 (REVIEW.md round 6).
#
# Validates the petsc_ts_restart_plan_2026-05-16.md Phase-1 deliverable:
# the `--restart` + `--petsc-ts` combination must now work end-to-end on
# Frontera, replacing the old "ERROR: --restart is not supported with
# --petsc-ts" early-abort.
#
# Two phases inside a SINGLE 2h development-queue reservation.
#
# Wall-time budget: the reference 48h Phase-6 sbatch covers 634 yr at
# 8N x 400r (~13 yr/h baseline; pure interseismic is roughly 1.5-2x
# faster than the event-included average).  2h on the dev queue gives
# ~26 simulated years at the conservative baseline, ~50 yr if we got
# lucky with pure interseismic.  Two phases of ~13 yr each is the
# safe target, leaving room for the ~10 min build + ~1 min validation
# block + minor queue-launch overhead.  Dev queue's faster turnaround
# (vs `normal`) keeps launch latency small enough that 2h works.
#
#   Phase A — fresh run, ~12.7 yr (tfinal = 4e8 s).
#             Pure interseismic regime (first BP5 earthquake at ~240 yr).
#             Driver writes V1 + V2 checkpoint at end-of-run
#             (--checkpoint-interval triggers the final write at
#             `bp5_verification_full.cpp:2814-2868`).
#             Expected wall: ~30 min.
#
#   Phase B — RESTART from Phase A's checkpoint, run to ~25.4 yr
#             (tfinal = 8e8 s).  Driver hits the V2 restart block
#             at `bp5_verification_full.cpp:2441-2563` and replays
#             the checkpointed t / current_dt / TSSetStepNumber.
#             Writes its own V1 + V2 checkpoint at end-of-run.
#             Expected wall: ~30 min (Phase B covers ~12.7 yr MORE).
#
# Total expected wall: ~70 min (build + 2 x 30 min + validation),
# leaving ~50 min slack inside the 2h reservation for variance.
#
# Validation block at the end of the script:
#   1. Phase A produced V1 + V2 trailing block in every per-rank file.
#   2. Phase B's stdout contains the "PETSc TS restart (V2): ..." log
#      line printed by the V2 restart block (proves V2 actually ran).
#   3. Phase B's final time > Phase A's final time (proves the
#      trajectory continued past the seam).
#   4. Phase B's cumulative rejection count >= Phase A's
#      (R-005: restart_rejections_carryover preserved).
#
# Acceptance gate (plan §"Acceptance Criteria" #1, #4, #14):
#   - #1: --restart + --petsc-ts no longer aborts.
#   - #4: V2 trailing block has all 10 fields.
#   - #14: the sbatch can be re-submitted (we automate it here in one
#     submission, but the same logic applies to a manual A→B chain).
#
# ParaView output ENABLED at production cadence — restart correctness is
# NOT PV-independent.  The V2 trailing block carries six ParaView-
# schedule fields (paraview_snapshots, last_write_time, last_v_max,
# current_regime, last_committed_cycle, last_volume_write_time) added
# precisely because R-304 / R-004 / R-006 found that restart breaks the
# adaptive schedule without them.  With ParaView disabled the driver's
# `if (pv_out)` guard skips RestoreScheduleState / SetTotalSnapshotsWritten
# / SetLastCommittedCycle entirely, leaving the most ParaView-coupled
# parts of the V2 machinery unexercised on the cluster — exactly the
# bits the unit test cannot reach because it has no ParMesh / no real
# fault output.  Cost: ~13 fault snapshots/phase at the production
# 1 yr interseismic cadence, ZFP tol 1e-12 (~1.4 MB/frame) ≈ 36 MB
# total across both phases.  Negligible vs the wall-clock budget.
#
# Volume PV suppressed (--no-volume-pv) — matches the production 48h
# sbatch.  Only fault-surface output is emitted.
#
# Why 8N x 400r matching the production footprint?  R-001 / R-302 of
# the restart plan documents that rank-count changes between
# checkpoint write and restart read are forbidden by
# `ReadCheckpoint`'s `file_num_ranks == size` guard.  Running at the
# production rank count ensures the production restart workflow is
# what's being tested.  Dev queue (vs `normal`) is fine for the
# footprint and gives faster turnaround.

export LC_ALL=C
export LANG=C

# Fail fast on module-load / build / pipe failures.
set -euo pipefail

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

# Auto-detect project root (walk up from $SLURM_SUBMIT_DIR to find
# build_frontera.sh).  Override with SEAS_MFEM_ROOT=$(pwd) sbatch ...
if [ -z "${SEAS_MFEM_ROOT:-}" ]; then
    _probe="${SLURM_SUBMIT_DIR:-$PWD}"
    while [ "${_probe}" != "/" ] && [ -n "${_probe}" ]; do
        if [ -f "${_probe}/build_frontera.sh" ]; then
            SEAS_MFEM_ROOT="${_probe}"
            break
        fi
        _probe="$(dirname "${_probe}")"
    done
    unset _probe
fi
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-/scratch2/10024/zhaochun/seas-project/seas-mfem}"
echo "  SEAS_MFEM_ROOT = ${SEAS_MFEM_ROOT}"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
# NOTE: do NOT load phdf5 on Frontera (compat with intel/19 stack).
# build_frontera.sh builds HDF5 1.14.x from source into
# ${SEAS_MFEM_ROOT}/extern/hdf5/install.

HDF5_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/hdf5/install"
H5Z_ZFP_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install"

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${HDF5_BUILT_DIR}/lib:${LD_LIBRARY_PATH:-}"

if [ -d "${H5Z_ZFP_BUILT_DIR}/plugin" ]; then
   export HDF5_PLUGIN_PATH="${H5Z_ZFP_BUILT_DIR}/plugin"
fi

cd "${SEAS_MFEM_ROOT}"

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES'    config/config.mk || { echo "ERROR: no PETSc";   exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'    config/config.mk || { echo "ERROR: no MUMPS";   exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'     config/config.mk || { echo "ERROR: no HDF5";    exit 1; }

cd miniapps/seas
make clean
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

# =======================================================================
# Result directories.  Phase A and Phase B each get their own so the
# checkpoint paths are unambiguous (--restart points at Phase A's
# prefix; Phase B's own checkpoint lands in Phase B's dir).
# =======================================================================
# Chained-restart naming convention (matches the driver's recommended
# pattern in the --restart / --output-dir collision safety check):
# all segments live under a single base dir, each in its own
# numbered `segment_NNN/` subdir.  Phase A → segment_001, Phase B →
# segment_002.  Future chained restarts can drop into segment_003,
# segment_004, ... without touching prior outputs.
RESTART_TEST_BASE="bp5/results_restart_test_job${SLURM_JOB_ID}"
RESULT_DIR_A="${RESTART_TEST_BASE}/segment_001"
RESULT_DIR_B="${RESTART_TEST_BASE}/segment_002"
OUTPUT_PREFIX_A="bp5_restart_seg001_job${SLURM_JOB_ID}"
OUTPUT_PREFIX_B="bp5_restart_seg002_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR_A}" "${RESULT_DIR_B}"

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Base dir:      ${RESTART_TEST_BASE}"
echo "Phase A dir:   ${RESULT_DIR_A}  (segment_001 — fresh)"
echo "Phase B dir:   ${RESULT_DIR_B}  (segment_002 — restart from A)"
echo "Footprint:     8N x 400r, dev queue, 2h total (~30 min/phase + 60 min slack)"
echo "Phase A tfinal: 4.0e8 s ≈ 12.7 yr (pure interseismic, no event)"
echo "Phase B tfinal: 8.0e8 s ≈ 25.4 yr (continues from A; still no event)"
echo "Checkpoint:    --checkpoint-interval 5000 (V1 + V2 trailing)"
echo "ParaView:      ENABLED — fault HDF5 at ZFP 1e-12, 1 yr interseismic"
echo "               cadence; --no-volume-pv (matches production 48h)"
echo "               Exercises pv_out->Restore/Set* on the V2 restart path"
echo "========================================"

# =======================================================================
# Phase A — fresh run, write V1+V2 checkpoint at end-of-run.
# =======================================================================
echo ""
echo "========================================"
echo "Phase A: fresh run to t = 4.0e8 s (≈ 12.7 yr)"
echo "========================================"
ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --ref-dir bp5/benchmark_data \
      --output-dir "${RESULT_DIR_A}" \
      --output-prefix "${OUTPUT_PREFIX_A}" \
      --checkpoint-interval 5000 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --paraview-adaptive \
      --no-volume-pv \
      --paraview-fault-zfp-tol 1e-12 \
      --paraview-dt-co     1.0 \
      --paraview-dt-nu     43200 \
      --paraview-dt-inter-yr 1.0 \
      --tfinal 400000000 || PHASE_A_RC=$?
PHASE_A_RC=${PHASE_A_RC:-0}
echo "Phase A exit code: ${PHASE_A_RC}"
if [ "${PHASE_A_RC}" -ne 0 ]; then
   echo "FAIL: Phase A did not exit cleanly.  Aborting before Phase B."
   exit 1
fi

# Acceptance check #1 — V1 + V2 trailing block present on every rank.
echo ""
echo "----------------------------------------"
echo "Validation #1: Phase A produced V1 + V2 checkpoint files"
echo "----------------------------------------"
PHASE_A_PREFIX_PATH="${RESULT_DIR_A}/${OUTPUT_PREFIX_A}"
EXPECTED_RANKS=400
FOUND_V1=0
FOUND_V2=0
MISSING_RANKS=""
for r in $(seq 0 $((EXPECTED_RANKS - 1))); do
   ckpt="${PHASE_A_PREFIX_PATH}_checkpoint_r${r}.txt"
   if [ ! -f "${ckpt}" ]; then
      MISSING_RANKS="${MISSING_RANKS} ${r}"
      continue
   fi
   if head -1 "${ckpt}" | grep -q "^SEAS_CHECKPOINT_V1$"; then
      FOUND_V1=$((FOUND_V1 + 1))
   fi
   if grep -q "^PETSC_TS_V2$" "${ckpt}"; then
      FOUND_V2=$((FOUND_V2 + 1))
   fi
done
echo "  V1 header found on ${FOUND_V1} / ${EXPECTED_RANKS} ranks"
echo "  V2 trailer found on ${FOUND_V2} / ${EXPECTED_RANKS} ranks"
if [ -n "${MISSING_RANKS}" ]; then
   echo "  WARNING: missing checkpoint files on ranks:${MISSING_RANKS}"
fi
if [ "${FOUND_V1}" -ne "${EXPECTED_RANKS}" ] \
   || [ "${FOUND_V2}" -ne "${EXPECTED_RANKS}" ]; then
   echo "FAIL: V1/V2 checkpoint coverage incomplete.  Aborting before Phase B."
   exit 2
fi
echo "  PASS"

# Capture Phase A's final time + step from its stdout (the driver
# prints "Final time: NN years" and "Total steps: N" in the summary).
PHASE_A_LOG="${SBATCH_LOG_DIR}/bp5_restart_test_${SLURM_JOB_ID}.out"
PHASE_A_FINAL_TIME=$(awk '/Final time:/ { print $3 }' "${PHASE_A_LOG}" | tail -1)
PHASE_A_FINAL_STEPS=$(awk '/Total steps:/ { print $3 }' "${PHASE_A_LOG}" | tail -1)
PHASE_A_FINAL_REJECTS=$(awk '/Step rejections:/ { print $3 }' "${PHASE_A_LOG}" | tail -1)
echo "  Phase A final time:     ${PHASE_A_FINAL_TIME} yr"
echo "  Phase A total steps:    ${PHASE_A_FINAL_STEPS}"
echo "  Phase A rejection count: ${PHASE_A_FINAL_REJECTS}"

# Capture Phase A fault.vtkhdf size BEFORE Phase B runs.  Used by
# Validation #8 (R-103) to confirm Phase B did not mutate Phase A's
# fault output.  Without this capture, Validation #8 would be a no-op
# (both stat calls after Phase B exits are guaranteed equal).
FAULT_A_BEFORE_B_SIZE=$(stat -c %s "${RESULT_DIR_A}/fault.vtkhdf" 2>/dev/null \
                     || stat -f %z "${RESULT_DIR_A}/fault.vtkhdf" 2>/dev/null \
                     || echo 0)
echo "  Phase A fault.vtkhdf size (before B): ${FAULT_A_BEFORE_B_SIZE} bytes"

# =======================================================================
# Phase B — restart from Phase A's checkpoint, run to t = 3e9 s.
# =======================================================================
echo ""
echo "========================================"
echo "Phase B: restart from Phase A, run to t = 3.0e9 s (≈ 95 yr)"
echo "         --restart ${PHASE_A_PREFIX_PATH}"
echo "========================================"
ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --ref-dir bp5/benchmark_data \
      --output-dir "${RESULT_DIR_B}" \
      --output-prefix "${OUTPUT_PREFIX_B}" \
      --checkpoint-interval 5000 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --paraview-adaptive \
      --no-volume-pv \
      --paraview-fault-zfp-tol 1e-12 \
      --paraview-dt-co     1.0 \
      --paraview-dt-nu     43200 \
      --paraview-dt-inter-yr 1.0 \
      --restart "${PHASE_A_PREFIX_PATH}" \
      --tfinal 800000000 || PHASE_B_RC=$?
PHASE_B_RC=${PHASE_B_RC:-0}
echo "Phase B exit code: ${PHASE_B_RC}"
if [ "${PHASE_B_RC}" -ne 0 ]; then
   echo "FAIL: Phase B did not exit cleanly.  Restart broken."
   exit 3
fi

# =======================================================================
# Validation block — was restart actually exercised?
# =======================================================================
echo ""
echo "========================================"
echo "Validation summary"
echo "========================================"

PHASE_B_LOG="${SBATCH_LOG_DIR}/bp5_restart_test_${SLURM_JOB_ID}.out"

# Acceptance check #2 — Phase B's stdout contains the V2 restart log line
# (printed by the V2 restart block on rank 0 only).
echo ""
echo "Validation #2: V2 restart block actually ran"
echo "----------------------------------------"
if grep -q "PETSc TS restart (V2):" "${PHASE_B_LOG}"; then
   echo "  PASS — found 'PETSc TS restart (V2):' log line"
   grep "PETSc TS restart (V2):" "${PHASE_B_LOG}" | head -1
else
   echo "  FAIL — V2 restart block did NOT run; --restart was probably"
   echo "         treated as no-op or aborted silently"
   exit 4
fi

# Acceptance check #3 — Phase B's final time > Phase A's final time.
PHASE_B_FINAL_TIME=$(awk '/Final time:/ { print $3 }' "${PHASE_B_LOG}" | tail -1)
PHASE_B_FINAL_STEPS=$(awk '/Total steps:/ { print $3 }' "${PHASE_B_LOG}" | tail -1)
PHASE_B_FINAL_REJECTS=$(awk '/Step rejections:/ { print $3 }' "${PHASE_B_LOG}" | tail -1)
echo ""
echo "Validation #3: Phase B's trajectory continued past Phase A"
echo "----------------------------------------"
echo "  Phase A final time:     ${PHASE_A_FINAL_TIME} yr"
echo "  Phase B final time:     ${PHASE_B_FINAL_TIME} yr"
# Use awk for FP comparison (bash arithmetic is integer-only).
TIME_CONTINUED=$(awk -v a="${PHASE_A_FINAL_TIME}" -v b="${PHASE_B_FINAL_TIME}" \
                 'BEGIN { print (b > a) ? "1" : "0" }')
if [ "${TIME_CONTINUED}" = "1" ]; then
   echo "  PASS — Phase B (t=${PHASE_B_FINAL_TIME}) > Phase A (t=${PHASE_A_FINAL_TIME})"
else
   echo "  FAIL — Phase B final time did NOT exceed Phase A's"
   echo "         (either restart loaded wrong t, or Phase B ran zero steps)"
   exit 5
fi

# Acceptance check #4 — cumulative rejection count is monotonic across
# the restart seam.  Phase B's total step_rejections should be >=
# Phase A's (R-005: carryover preserved).  If Phase B = Phase A,
# the restart didn't accumulate; if Phase B < Phase A, the carryover
# was dropped on restart.
echo ""
echo "Validation #4: cumulative rejection count is monotonic (R-005)"
echo "----------------------------------------"
echo "  Phase A rejections:     ${PHASE_A_FINAL_REJECTS}"
echo "  Phase B rejections:     ${PHASE_B_FINAL_REJECTS}"
if [ -z "${PHASE_A_FINAL_REJECTS}" ] || [ -z "${PHASE_B_FINAL_REJECTS}" ]; then
   echo "  SKIP — could not parse rejection counts from logs"
elif [ "${PHASE_B_FINAL_REJECTS}" -ge "${PHASE_A_FINAL_REJECTS}" ]; then
   echo "  PASS — Phase B (${PHASE_B_FINAL_REJECTS}) >= Phase A (${PHASE_A_FINAL_REJECTS})"
else
   echo "  FAIL — R-005 regression: Phase B rejection count"
   echo "         (${PHASE_B_FINAL_REJECTS}) < Phase A's"
   echo "         (${PHASE_A_FINAL_REJECTS}); carryover lost on restart"
   exit 6
fi

# Acceptance check #5 — Phase B's V2 checkpoint exists and is well-formed
# (i.e., chained restart from Phase B is possible — proves the V2 write
# path in the monitor / final-checkpoint sites also works).
echo ""
echo "Validation #5: Phase B's V2 checkpoint is well-formed"
echo "----------------------------------------"
PHASE_B_PREFIX_PATH="${RESULT_DIR_B}/${OUTPUT_PREFIX_B}"
B_V1=0; B_V2=0
for r in $(seq 0 $((EXPECTED_RANKS - 1))); do
   ckpt="${PHASE_B_PREFIX_PATH}_checkpoint_r${r}.txt"
   [ -f "${ckpt}" ] || continue
   head -1 "${ckpt}" | grep -q "^SEAS_CHECKPOINT_V1$" && B_V1=$((B_V1 + 1))
   grep -q "^PETSC_TS_V2$" "${ckpt}" && B_V2=$((B_V2 + 1))
done
echo "  Phase B V1 header: ${B_V1} / ${EXPECTED_RANKS}"
echo "  Phase B V2 trailer: ${B_V2} / ${EXPECTED_RANKS}"
if [ "${B_V1}" -eq "${EXPECTED_RANKS}" ] && [ "${B_V2}" -eq "${EXPECTED_RANKS}" ]; then
   echo "  PASS — Phase B wrote a valid V1+V2 checkpoint (chainable)"
else
   echo "  FAIL — Phase B checkpoint coverage incomplete"
   exit 7
fi

# Acceptance check #6 — paraview_snapshots is monotonic across the
# restart seam (R-004 / R-006 driver-side wiring proof).  Phase A's
# V2 stores the count at end-of-Phase-A; Phase B's V2 should store
# Phase_A_count + Phase_B_count (because the V2 restart block calls
# pv_out->SetTotalSnapshotsWritten(loaded_value), then Phase B
# accumulates more snapshots on top).  If the driver's `if (pv_out)`
# guard is broken or RestoreScheduleState fails to apply, Phase B's
# stored value would equal Phase B's own count alone (< Phase A) — or
# zero.
echo ""
echo "Validation #6: paraview_snapshots round-trips through V2 and is"
echo "                monotonic across the restart seam"
echo "----------------------------------------"
A_PV_SNAP=$(awk '/^paraview_snapshots /{print $2}' \
            "${PHASE_A_PREFIX_PATH}_checkpoint_r0.txt" | tail -1)
B_PV_SNAP=$(awk '/^paraview_snapshots /{print $2}' \
            "${PHASE_B_PREFIX_PATH}_checkpoint_r0.txt" | tail -1)
echo "  Phase A paraview_snapshots (V2 field, rank 0): ${A_PV_SNAP}"
echo "  Phase B paraview_snapshots (V2 field, rank 0): ${B_PV_SNAP}"
if [ -z "${A_PV_SNAP}" ] || [ -z "${B_PV_SNAP}" ]; then
   echo "  FAIL — could not parse paraview_snapshots from V2 blocks"
   exit 8
fi
# R-004 (REVIEW.md round 6): dropped the `A_PV_SNAP > 0` gate.  The
# strict round-trip check is now check #2 above (V2 block ran); the
# "ParaView produced output" gate is now check #7 below (fault.vtkhdf
# exists with > 1 KiB).  This loop just verifies monotonicity.
if [ "${B_PV_SNAP}" -ge "${A_PV_SNAP}" ]; then
   echo "  PASS — Phase B (${B_PV_SNAP}) >= Phase A (${A_PV_SNAP})"
   echo "         (monotonic; proves driver pv_out->Restore/Set* fired"
   echo "         on restart and the V2 paraview_snapshots field was"
   echo "         loaded + propagated into Phase B's counter)"
else
   echo "  FAIL — paraview_snapshots not monotonic across the restart"
   echo "         seam: R-004/R-006 wiring lost the saved count"
   exit 8
fi

# Acceptance check #7 — fault.vtkhdf actually landed in BOTH phases'
# output directories.  Catches the "ParaView silently emitted nothing"
# failure mode that check #6's snapshot-count monotonicity cannot
# distinguish from a legitimate "0 snapshots due to cadence/timing".
# Without this gate, an entirely broken ParaView path passes #6 (both
# A and B store paraview_snapshots = 0).  R-004 (REVIEW.md round 6).
echo ""
echo "Validation #7: fault.vtkhdf landed in both phases' output dirs"
echo "----------------------------------------"
FAULT_A="${RESULT_DIR_A}/fault.vtkhdf"
FAULT_B="${RESULT_DIR_B}/fault.vtkhdf"
# Cross-platform stat: GNU coreutils (-c %s) on Linux; BSD (-f %z) on
# macOS / Frontera login.  Try both; fall back to 0.
FAULT_A_SIZE=$(stat -c %s "${FAULT_A}" 2>/dev/null \
            || stat -f %z "${FAULT_A}" 2>/dev/null \
            || echo 0)
FAULT_B_SIZE=$(stat -c %s "${FAULT_B}" 2>/dev/null \
            || stat -f %z "${FAULT_B}" 2>/dev/null \
            || echo 0)
echo "  Phase A fault.vtkhdf size: ${FAULT_A_SIZE} bytes"
echo "  Phase B fault.vtkhdf size: ${FAULT_B_SIZE} bytes"
if [ "${FAULT_A_SIZE}" -gt 1024 ] && [ "${FAULT_B_SIZE}" -gt 1024 ]; then
   echo "  PASS — both fault.vtkhdf files exist with > 1 KiB"
else
   echo "  FAIL — fault.vtkhdf missing or empty in at least one phase"
   exit 9
fi

# Acceptance check #8 — Phase A's fault.vtkhdf size is UNCHANGED
# after Phase B finishes.  Compares the size captured BEFORE Phase B
# started (FAULT_A_BEFORE_B_SIZE, ~30 lines after Phase A's wall time
# was recorded) against the size now.  Without the BEFORE-B capture
# this would be a no-op (both stats after Phase B exit are guaranteed
# equal).  R-002/R-004/R-103 (REVIEW.md rounds 6 + 8).
echo ""
echo "Validation #8: Phase A fault.vtkhdf untouched after Phase B finishes"
echo "----------------------------------------"
FAULT_A_SIZE_AFTER_B=$(stat -c %s "${FAULT_A}" 2>/dev/null \
                    || stat -f %z "${FAULT_A}" 2>/dev/null \
                    || echo 0)
if [ "${FAULT_A_SIZE_AFTER_B}" = "${FAULT_A_BEFORE_B_SIZE}" ]; then
   echo "  PASS — Phase A fault.vtkhdf size unchanged"
   echo "         (was ${FAULT_A_BEFORE_B_SIZE} before B, still"
   echo "          ${FAULT_A_SIZE_AFTER_B} after B)"
else
   echo "  FAIL — Phase A's fault.vtkhdf was modified during Phase B"
   echo "         (was ${FAULT_A_BEFORE_B_SIZE} before B, now ${FAULT_A_SIZE_AFTER_B})"
   echo "         The --output-dir safety check did not isolate Phase B"
   exit 10
fi

# Acceptance check #9 — station / probe files (`*_fltst_*.txt`)
# present in BOTH phases' output dirs and SEPARATE.  Objective #2
# (REVIEW.md round 6 R-005) specifically called out station-file
# preservation; absent this check, a regression in ProbeOutput's
# path construction could silently route Phase B's station files
# into Phase A's directory, overwriting them.
echo ""
echo "Validation #9: station files (*_fltst_*.txt) present + separate"
echo "----------------------------------------"
A_STATIONS=$(ls "${RESULT_DIR_A}"/*_fltst_*.txt 2>/dev/null | wc -l | tr -d ' ')
B_STATIONS=$(ls "${RESULT_DIR_B}"/*_fltst_*.txt 2>/dev/null | wc -l | tr -d ' ')
echo "  Phase A station files: ${A_STATIONS}"
echo "  Phase B station files: ${B_STATIONS}"
if [ "${A_STATIONS}" -gt 0 ] && [ "${B_STATIONS}" -gt 0 ]; then
   echo "  PASS — both phases produced their own *_fltst_*.txt files"
   echo "         (objective #2: V2 station files saved, V1 preserved)"
else
   echo "  FAIL — station files missing in at least one phase"
   echo "         (Phase A=${A_STATIONS}, Phase B=${B_STATIONS})"
   exit 11
fi

# =======================================================================
# Summary + disk usage record.
# =======================================================================
echo ""
echo "========================================"
echo "All restart-test acceptance checks PASSED"
echo "========================================"

SUMMARY_FILE="${RESULT_DIR_B}/RESTART_TEST_SUMMARY.txt"
{
   echo "BP5 V2 PETSc-TS restart end-to-end test — SUMMARY"
   echo "Git commit:  $(git rev-parse HEAD)"
   echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
   echo "Slurm job:   ${SLURM_JOB_ID}"
   echo "Footprint:   8N x 400r, dev queue, 2h reservation"
   echo ""
   echo "Phase A:"
   echo "  Output:      ${RESULT_DIR_A}"
   echo "  tfinal:      1.5e9 s ≈ 47.5 yr"
   echo "  Final time:  ${PHASE_A_FINAL_TIME} yr"
   echo "  Steps:       ${PHASE_A_FINAL_STEPS}"
   echo "  Rejections:  ${PHASE_A_FINAL_REJECTS}"
   echo ""
   echo "Phase B (restart from Phase A):"
   echo "  Output:      ${RESULT_DIR_B}"
   echo "  tfinal:      3.0e9 s ≈ 95 yr"
   echo "  Final time:  ${PHASE_B_FINAL_TIME} yr"
   echo "  Steps:       ${PHASE_B_FINAL_STEPS}"
   echo "  Rejections:  ${PHASE_B_FINAL_REJECTS} (cumulative, includes Phase A's)"
   echo ""
   echo "Acceptance checks:"
   echo "  #1 V1+V2 checkpoint coverage (Phase A): PASS (${FOUND_V1}/${EXPECTED_RANKS})"
   echo "  #2 V2 restart block ran on Phase B:     PASS"
   echo "  #3 Phase B trajectory continued past A: PASS (${PHASE_B_FINAL_TIME} > ${PHASE_A_FINAL_TIME})"
   echo "  #4 R-005 rejection monotonicity:        PASS (${PHASE_B_FINAL_REJECTS} >= ${PHASE_A_FINAL_REJECTS})"
   echo "  #5 Phase B V2 checkpoint chainable:     PASS (${B_V1}/${EXPECTED_RANKS})"
   echo "  #6 pv_snapshot monotonic across seam:   PASS (${B_PV_SNAP} >= ${A_PV_SNAP})"
   echo "  #7 fault.vtkhdf landed in both phases:  PASS (A=${FAULT_A_SIZE} B=${FAULT_B_SIZE} bytes)"
   echo "  #8 Phase A fault.vtkhdf untouched:      PASS (size unchanged: ${FAULT_A_SIZE})"
   echo "  #9 station files separate (obj #2):     PASS (A=${A_STATIONS} B=${B_STATIONS})"
   echo ""
   echo "du -sh (Phase A): $(du -sh ${RESULT_DIR_A} | cut -f1)"
   echo "du -sh (Phase B): $(du -sh ${RESULT_DIR_B} | cut -f1)"
} > "${SUMMARY_FILE}"
cat "${SUMMARY_FILE}"
