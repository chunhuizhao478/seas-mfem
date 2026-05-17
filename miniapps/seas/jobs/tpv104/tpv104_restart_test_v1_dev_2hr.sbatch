#!/bin/bash
#SBATCH -J tpv104_restart_v1
#SBATCH -o tpv104_restart_v1_%j.out
#SBATCH -e tpv104_restart_v1_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 — V1 restart end-to-end test (dev queue, 2h, 8N x 400r).
# =======================================================================
# Validates the Phase-4 TPV104 restart port: --restart +
# --checkpoint-interval now work end-to-end, mirroring the BP5 V1
# restart pattern.  Format details: `io/tpv104_checkpoint.hpp`
# (TPV104_CHECKPOINT_V1 magic tag distinct from BP5's
# SEAS_CHECKPOINT_V1; per-rank text files using the same
# {prefix}_checkpoint_r{rank}.txt naming as BP5).
#
# Two phases inside a SINGLE 2h dev-queue reservation:
#
#   Phase A — fresh run 0 → 2.5 s (~1 hr wall at 8N x 400r, dynamic
#             rupture).  Driver writes V1 checkpoint at end-of-run
#             (always when --checkpoint-interval > 0).
#
#   Phase B — RESTART from Phase A's checkpoint, runs 2.5 → 5.0 s
#             (~1 hr wall — driver loads Q + dof_data dynamic fields
#             + scalars, then resumes the time loop from
#             `step = restart_step` instead of 0).  Writes its own
#             V1 checkpoint at end-of-run (chainable).
#
# Rate: 2.5 s simulated time per 1 hr wall.  Total wall ~2 hr fits
# the dev-queue reservation exactly.
#
# V3 deferral (out of scope for this V1 port; tracked in
# `debug_document/paraview_output_debug_document/petsc_ts_restart_plan_2026-05-16.md`
# §"Out of scope"):
#   - SECONDARY ParaView collection's schedule state.  V2 (BP5)
#     carries ONE collection's regime + last_write_time; TPV104 has
#     TWO (pv_out + pv_bulk_out).  On restart the driver emits
#     `WARNING: --restart resumes the PRIMARY ParaView collection's
#     schedule state ... V3 will carry both collections.`  Each phase
#     still writes a fresh fault.vtkhdf + volume.vtkhdf +
#     ParaView_bulk/volume.vtkhdf in its own --output-dir; only the
#     pv_bulk_out CADENCE memory is lost across the seam.
#   - substep-iterator internal state (only matters on
#     --fault-iterator substep; one-shot AdvanceADER is stateless
#     between steps).
#   - DOFData STATIC fields (impedances, a, Dc, background prestress,
#     LSW params) — re-initialised by InitializeFaultDOFs_TPV104 on
#     restart.  Checkpoint preserves only DYNAMIC fields (psi,
#     slip_rate, V1/V2, slip1/slip2, tau*_nuc, sigma_n_nuc).

export LC_ALL=C
export LANG=C

# Fail fast on module-load / build / pipe failures.
set -euo pipefail

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

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
# NOTE: do NOT load phdf5 on Frontera.  build_frontera.sh builds HDF5
# 1.14.x from source into ${SEAS_MFEM_ROOT}/extern/hdf5/install.

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
grep -q '^MFEM_USE_H5Z_ZFP *= YES'  config/config.mk || { echo "ERROR: no H5Z-ZFP"; exit 1; }

cd miniapps/seas

# Mesh regen (TPV104 200m mesh is too large to commit; .geo is tracked).
MESH_SRC="tpv104/mesh/tpv104_200m.geo"
MESH_OUT="tpv104/mesh/tpv104_200m.msh"
GMSH_BIN="${SEAS_MFEM_ROOT}/extern/gmsh/bin/gmsh"
if [ ! -s "${MESH_OUT}" ] || [ "${FORCE_MESHGEN:-0}" = "1" ]; then
    if [ ! -x "${GMSH_BIN}" ]; then
        echo "ERROR: ${GMSH_BIN} not found.  Run build_frontera.sh first."
        exit 1
    fi
    if [ ! -s "${MESH_SRC}" ]; then
        echo "ERROR: ${MESH_SRC} missing from the checkout."
        exit 1
    fi
    echo "=== Regenerating ${MESH_OUT} from ${MESH_SRC} ==="
    "${GMSH_BIN}" -3 "${MESH_SRC}" -o "${MESH_OUT}" \
        || { echo "MESH GENERATION FAILED"; exit 1; }
fi

rm -f seas_tpv104_driver drivers/tpv104_driver.o
make seas_tpv104_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

# =======================================================================
# Chained-restart naming convention (matches the BP5 sbatch +
# bp5_verification_full.cpp safety-check hint): single base dir +
# segment_NNN subdirs.  segment_002 is RESTARTED from segment_001's
# end-of-run V1 checkpoint via `--restart ${RESULT_DIR_A}/${OUTPUT_PREFIX_A}`.
# =======================================================================
RESTART_TEST_BASE="tpv104/results_restart_test_job${SLURM_JOB_ID}"
RESULT_DIR_A="${RESTART_TEST_BASE}/segment_001"
RESULT_DIR_B="${RESTART_TEST_BASE}/segment_002"
OUTPUT_PREFIX_A="tpv104_restart_seg001_job${SLURM_JOB_ID}"
OUTPUT_PREFIX_B="tpv104_restart_seg002_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR_A}" "${RESULT_DIR_B}"

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Base dir:      ${RESTART_TEST_BASE}"
echo "Phase A dir:   ${RESULT_DIR_A}  (segment_001 — fresh, 0 → 2.5 s)"
echo "Phase B dir:   ${RESULT_DIR_B}  (segment_002 — RESTART, 2.5 → 5.0 s)"
echo "Footprint:     8N x 400r, dev queue, 2h total"
echo "Phase A tfinal: 2.5 s  (~1 hr wall, fresh)"
echo "Phase B tfinal: 5.0 s  (~1 hr wall — driver loads Phase A's V1"
echo "                        checkpoint and resumes from t=2.5 s)"
echo "Rate:          2.5 s / 1 hr  →  total wall ~2 hr matches reservation"
echo "ParaView:      ENABLED.  PRIMARY collection schedule state"
echo "                        (last_write_time, regime) carries across"
echo "                        the seam (V2 BP5 schema).  SECONDARY"
echo "                        ParaView_bulk collection re-inits in its"
echo "                        default regime (V3 deferral).  Each phase"
echo "                        writes a fresh fault.vtkhdf + volume.vtkhdf"
echo "                        + ParaView_bulk/volume.vtkhdf in its own"
echo "                        --output-dir."
echo "Restart:       V1 checkpoint (TPV104_CHECKPOINT_V1) — preserves Q,"
echo "                dof_data dynamic fields, t, dt, step.  See"
echo "                io/tpv104_checkpoint.hpp for the on-disk schema."
echo "========================================"

# =======================================================================
# Phase A — fresh TPV104 run.  Driver writes the V1 checkpoint at
# end-of-run (always when --checkpoint-interval > 0) plus produces
# fault.vtkhdf + volume.vtkhdf + ParaView_bulk/volume.vtkhdf +
# *_station_*.dat probe files inside ${RESULT_DIR_A}.
# =======================================================================
echo ""
echo "========================================"
echo "Phase A: fresh TPV104 run, 0 → 2.5 s  (~1 hr wall)"
echo "========================================"
ibrun ./seas_tpv104_driver \
      --mesh tpv104/mesh/tpv104_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 2.5 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law aging \
      --fault-iterator one-shot \
      --mixed-flux adjacent \
      --output-dir "${RESULT_DIR_A}" \
      --output-prefix "${OUTPUT_PREFIX_A}" \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.1 \
      --paraview-bulk-dt 0.1 \
      --paraview-volume-zfp-tol 1e-3 \
      --paraview-bulk-zfp-tol 1e-3 \
      --paraview-fault-zfp-tol 1e-12 \
      --paraview-max-snapshots 50 \
      --checkpoint-interval 100 \
      || PHASE_A_RC=$?
PHASE_A_RC=${PHASE_A_RC:-0}
echo "Phase A exit code: ${PHASE_A_RC}"
if [ "${PHASE_A_RC}" -ne 0 ]; then
   echo "FAIL: Phase A did not exit cleanly.  Aborting before Phase B."
   exit 1
fi

# Capture Phase A station-file fingerprints BEFORE Phase B runs.
# Used by Validation #7 (R-004) to confirm Phase B did not mutate
# Phase A's outputs.  Stored in a temp file so the validation block
# can compare after Phase B finishes.
PHASE_A_STATION_HASHES_BEFORE_B="${RESTART_TEST_BASE}/phaseA_station_hashes_before_B.md5"
( cd "${RESULT_DIR_A}" && \
  ls *_station_*.dat 2>/dev/null | sort | xargs -r md5sum \
) > "${PHASE_A_STATION_HASHES_BEFORE_B}" 2>/dev/null || true
echo "Captured $(wc -l < "${PHASE_A_STATION_HASHES_BEFORE_B}") Phase A"
echo "  station-file md5sums before Phase B starts (R-004 unchanged check)."

# Capture Phase A fault.vtkhdf size BEFORE Phase B runs.  Used by
# Validation #3 (R-103) to confirm Phase B did not mutate Phase A's
# fault output.  Without this capture, Validation #3 would be a no-op
# (both stat calls after Phase B exits are guaranteed equal).
FAULT_A_BEFORE_B_SIZE=$(stat -c %s "${RESULT_DIR_A}/fault.vtkhdf" 2>/dev/null \
                     || stat -f %z "${RESULT_DIR_A}/fault.vtkhdf" 2>/dev/null \
                     || echo 0)
echo "Captured Phase A fault.vtkhdf size: ${FAULT_A_BEFORE_B_SIZE} bytes"
echo "  (before Phase B starts; R-103 unchanged check)."

# =======================================================================
# Phase B — RESTART from Phase A's V1 checkpoint.  Driver loads
# Q + dof_data dynamic fields + (t, dt, step) from
# {RESULT_DIR_A}/{OUTPUT_PREFIX_A}_checkpoint_r{rank}.txt, then
# resumes the ADER time loop from `step = restart_step`.  tfinal=5.0
# is ABSOLUTE (BP5 convention) so Phase B advances 2.5 s past
# Phase A's tfinal, completing in ~1 hr wall (same rate as A).
# Output to a SEPARATE dir (segment_002) enforced by the driver's
# safety check.
# =======================================================================
echo ""
echo "========================================"
echo "Phase B: TPV104 RESTART, 2.5 → 5.0 s  (~1 hr wall)"
echo "         --restart ${RESULT_DIR_A}/${OUTPUT_PREFIX_A}"
echo "========================================"
ibrun ./seas_tpv104_driver \
      --mesh tpv104/mesh/tpv104_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 5.0 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law aging \
      --fault-iterator one-shot \
      --mixed-flux adjacent \
      --output-dir "${RESULT_DIR_B}" \
      --output-prefix "${OUTPUT_PREFIX_B}" \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.1 \
      --paraview-bulk-dt 0.1 \
      --paraview-volume-zfp-tol 1e-3 \
      --paraview-bulk-zfp-tol 1e-3 \
      --paraview-fault-zfp-tol 1e-12 \
      --paraview-max-snapshots 50 \
      --checkpoint-interval 100 \
      --restart "${RESULT_DIR_A}/${OUTPUT_PREFIX_A}" \
      || PHASE_B_RC=$?
PHASE_B_RC=${PHASE_B_RC:-0}
echo "Phase B exit code: ${PHASE_B_RC}"
if [ "${PHASE_B_RC}" -ne 0 ]; then
   echo "FAIL: Phase B did not exit cleanly."
   exit 2
fi

# =======================================================================
# Validation block — verifies TPV104 V1 restart end-to-end.
#
# Expected outcome (Phase-4 V1 port landed): Phase B's --restart loads
# Phase A's checkpoint, prints "TPV104 restart loaded:" on rank 0, and
# resumes the time loop from t ≈ 2.5 s.  Each phase's outputs land in
# its own --output-dir (no clobber).
#
# Acceptance checks (exit codes 3-9):
#   #1 Phase A fault.vtkhdf landed (>1 KiB)
#   #2 Phase B fault.vtkhdf landed in its own dir (>1 KiB)
#   #3 Phase A fault.vtkhdf untouched after Phase B (size compare)
#   #4 Phase B's "TPV104 restart loaded:" log line present
#   #5 Phase B station files start at t ≈ 2.5 s
#   #6 Both phases produced their own station files
#   #7 Phase A station files untouched after Phase B (md5sum compare)
# =======================================================================
echo ""
echo "========================================"
echo "Validation summary (TPV104 V1 restart end-to-end)"
echo "========================================"

PHASE_LOG="${SBATCH_LOG_DIR}/tpv104_restart_v1_${SLURM_JOB_ID}.out"

# Acceptance check #1 — Phase A produced fault.vtkhdf.  Sanity check
# that the driver is functioning at all and Phase A's output layout
# matches what Phase B's --restart load expects.
echo ""
echo "Validation #1: Phase A produced fault.vtkhdf"
echo "----------------------------------------"
FAULT_A="${RESULT_DIR_A}/fault.vtkhdf"
FAULT_A_SIZE=$(stat -c %s "${FAULT_A}" 2>/dev/null \
            || stat -f %z "${FAULT_A}" 2>/dev/null \
            || echo 0)
echo "  Phase A fault.vtkhdf size: ${FAULT_A_SIZE} bytes"
if [ "${FAULT_A_SIZE}" -gt 1024 ]; then
   echo "  PASS — Phase A wrote fault.vtkhdf (> 1 KiB)"
else
   echo "  FAIL — Phase A fault.vtkhdf missing or empty"
   exit 3
fi

# Acceptance check #2 — Phase B produced its own fault.vtkhdf in a
# DIFFERENT directory.  The --output-dir separation is enforced by
# the canonical-path collision check in tpv104_driver.cpp:495-579;
# this validation confirms Phase B's outputs landed in segment_002.
echo ""
echo "Validation #2: Phase B produced fault.vtkhdf in a separate dir"
echo "----------------------------------------"
FAULT_B="${RESULT_DIR_B}/fault.vtkhdf"
FAULT_B_SIZE=$(stat -c %s "${FAULT_B}" 2>/dev/null \
            || stat -f %z "${FAULT_B}" 2>/dev/null \
            || echo 0)
echo "  Phase B fault.vtkhdf size: ${FAULT_B_SIZE} bytes"
if [ "${FAULT_B_SIZE}" -gt 1024 ]; then
   echo "  PASS — Phase B wrote fault.vtkhdf (> 1 KiB) in segment_002"
else
   echo "  FAIL — Phase B fault.vtkhdf missing or empty"
   exit 4
fi

# Acceptance check #3 — Phase A's fault.vtkhdf size unchanged across
# Phase B's run.  Confirms the --output-dir separation + driver
# safety check kept Phase B's writes out of Phase A's directory.
# Size captured BEFORE Phase B started (line 217-219) and re-captured
# after — mirrors Validation #7's md5sum-based station-file check.
echo ""
echo "Validation #3: Phase A fault.vtkhdf untouched after Phase B finishes"
echo "----------------------------------------"
FAULT_A_SIZE_AFTER_B=$(stat -c %s "${FAULT_A}" 2>/dev/null \
                    || stat -f %z "${FAULT_A}" 2>/dev/null \
                    || echo 0)
if [ "${FAULT_A_SIZE_AFTER_B}" = "${FAULT_A_BEFORE_B_SIZE}" ]; then
   echo "  PASS — Phase A fault.vtkhdf size unchanged"
   echo "         (was ${FAULT_A_BEFORE_B_SIZE} before B, still"
   echo "          ${FAULT_A_SIZE_AFTER_B} after B)"
else
   echo "  FAIL — Phase A fault.vtkhdf was modified during Phase B"
   echo "         (was ${FAULT_A_BEFORE_B_SIZE} before B, now ${FAULT_A_SIZE_AFTER_B})"
   exit 5
fi

# Acceptance check #4 — proves Phase B's V1 restart actually fired.
# Phase B's stdout MUST contain the "TPV104 restart loaded:" log line
# printed by the driver's restart block on rank 0 (proves the V1 file
# was loaded and the time loop actually resumed from t = restart_t,
# not started from t = 0).
echo ""
echo "Validation #4: Phase B's V1 restart block actually ran"
echo "                (looks for 'TPV104 restart loaded:' on rank 0)"
echo "----------------------------------------"
if grep -qE "TPV104 restart loaded:" "${PHASE_LOG}" 2>/dev/null; then
   echo "  PASS — Phase B log contains 'TPV104 restart loaded:' line"
   grep "TPV104 restart loaded:" "${PHASE_LOG}" | head -1
else
   echo "  FAIL — 'TPV104 restart loaded:' line missing from Phase B"
   echo "         stdout.  Either the V1 file did not exist, or the"
   echo "         driver's restart block silently skipped it.  Check:"
   echo "          - Phase A wrote the per-rank checkpoint files"
   echo "          - drivers/tpv104_driver.cpp still calls"
   echo "            ReadTpv104Checkpoint when --restart is set"
   echo "          - The safety check did not abort Phase B (would"
   echo "            be 'Continuing would clobber' instead)"
   exit 6
fi

# Acceptance check #5 — Phase B's *_station_*.dat station files start
# at t ≈ 2.5 s (Phase A's tfinal), NOT t=0.  Proves the V1 restart
# correctly loaded t and the time loop is integrating from the saved
# state.  Threshold 2.4 (rather than 2.5) gives a small tolerance
# margin for the first post-restart probe-write that might fall a
# few output_dt's after the restart point.  TPV104StationWriter
# names files <prefix>_station_<name>.dat (dynamic/tpv104_setup.hpp:400);
# the BP5 glob *_fltst_*.txt does NOT match.
echo ""
echo "Validation #5: Phase B station files start near t = 2.5 s"
echo "                (proves V1 restart loaded the time correctly)"
echo "----------------------------------------"
B_STATIONS_FIRST=$(ls "${RESULT_DIR_B}"/*_station_*.dat 2>/dev/null | head -1)
if [ -z "${B_STATIONS_FIRST}" ]; then
   echo "  FAIL — no *_station_*.dat files in ${RESULT_DIR_B}.  Phase B"
   echo "         restart did not produce station outputs.  Check:"
   echo "          - TPV104StationWriter wired in (dynamic/tpv104_setup.hpp)"
   echo "          - --output-dir / --output-prefix passed through to it"
   exit 7
fi
FIRST_T=$(awk '/^[[:space:]]*[0-9-]/ { print $1; exit }' \
          "${B_STATIONS_FIRST}" 2>/dev/null || echo "")
echo "  First t in ${B_STATIONS_FIRST}: ${FIRST_T:-(unparseable)}"
if [ -z "${FIRST_T}" ]; then
   echo "  FAIL — could not parse first t from Phase B station file"
   exit 7
fi
RESTART_TIME_OK=$(awk -v t="${FIRST_T}" \
                 'BEGIN { print (t >= 2.4) ? "1" : "0" }')
if [ "${RESTART_TIME_OK}" = "1" ]; then
   echo "  PASS — Phase B station file starts at t=${FIRST_T}"
   echo "         (>= 2.4, proves restart loaded t ≈ Phase A tfinal)"
else
   echo "  FAIL — Phase B station file starts at t=${FIRST_T},"
   echo "         expected >= 2.4 (Phase A tfinal was 2.5).  The"
   echo "         restart block didn't load t correctly, OR the"
   echo "         driver re-initialised station_writer from t=0"
   echo "         on restart.  Check ReadTpv104Checkpoint's t"
   echo "         output and the station_writer init path."
   exit 7
fi

# Acceptance check #6 — station files in both phases' --output-dirs.
# Mirrors BP5 check #9.  Both counts MUST be > 0; the prior INFO
# fall-through (when both were 0) silently passed the check on
# regressions where the station writer disappeared entirely.
echo ""
echo "Validation #6: station files present in both phases' dirs"
echo "----------------------------------------"
A_STATIONS=$(ls "${RESULT_DIR_A}"/*_station_*.dat 2>/dev/null | wc -l | tr -d ' ')
B_STATIONS=$(ls "${RESULT_DIR_B}"/*_station_*.dat 2>/dev/null | wc -l | tr -d ' ')
echo "  Phase A station files: ${A_STATIONS}"
echo "  Phase B station files: ${B_STATIONS}"
if [ "${A_STATIONS}" -gt 0 ] && [ "${B_STATIONS}" -gt 0 ]; then
   echo "  PASS — both phases produced their own *_station_*.dat files"
   echo "         (A=${A_STATIONS}, B=${B_STATIONS})"
else
   echo "  FAIL — at least one phase produced no *_station_*.dat files"
   echo "         (A=${A_STATIONS}, B=${B_STATIONS})."
   echo "         Either the station writer regressed or the glob"
   echo "         convention drifted (was BP5's *_fltst_*.txt prior"
   echo "         to R-004 fix; TPV104 writes *_station_*.dat)."
   exit 8
fi

# Acceptance check #7 — Phase A station files unchanged after Phase B.
# Re-hash the same files we hashed before Phase B; the md5sums must
# match.  Mirror of validation #3 for fault.vtkhdf (size compare).
echo ""
echo "Validation #7: Phase A station files md5sum unchanged after Phase B"
echo "----------------------------------------"
PHASE_A_STATION_HASHES_AFTER_B="${RESTART_TEST_BASE}/phaseA_station_hashes_after_B.md5"
( cd "${RESULT_DIR_A}" && \
  ls *_station_*.dat 2>/dev/null | sort | xargs -r md5sum \
) > "${PHASE_A_STATION_HASHES_AFTER_B}" 2>/dev/null || true
if [ ! -s "${PHASE_A_STATION_HASHES_BEFORE_B}" ]; then
   echo "  FAIL — pre-B md5sum file is empty; Phase A produced no"
   echo "         station files (check Validation #6)."
   exit 9
fi
if diff -q "${PHASE_A_STATION_HASHES_BEFORE_B}" \
           "${PHASE_A_STATION_HASHES_AFTER_B}" > /dev/null 2>&1; then
   echo "  PASS — Phase A station files byte-identical before/after Phase B"
   echo "         ($(wc -l < "${PHASE_A_STATION_HASHES_BEFORE_B}") files)"
else
   echo "  FAIL — Phase A station files were mutated during Phase B."
   echo "         Diff:"
   diff "${PHASE_A_STATION_HASHES_BEFORE_B}" \
        "${PHASE_A_STATION_HASHES_AFTER_B}" | head -20 | sed 's/^/         /'
   echo "         This means --output-dir isolation between phases is"
   echo "         broken — Phase B clobbered Phase A's writes despite"
   echo "         the safety check.  Investigate the canonical-path"
   echo "         comparison in tpv104_driver.cpp:495-579."
   exit 9
fi

# =======================================================================
# Summary + disk usage record.
# =======================================================================
echo ""
echo "========================================"
echo "All TPV104 V1 restart acceptance checks PASSED"
echo "========================================"

SUMMARY_FILE="${RESULT_DIR_B}/RESTART_TEST_SUMMARY.txt"
{
   echo "TPV104 V1 restart end-to-end test — SUMMARY"
   echo "Git commit:  $(git rev-parse HEAD)"
   echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
   echo "Slurm job:   ${SLURM_JOB_ID}"
   echo "Footprint:   8N x 400r, dev queue, 2h reservation"
   echo ""
   echo "Phase A: fresh 0 → 2.5 s"
   echo "  Output: ${RESULT_DIR_A}"
   echo "  fault.vtkhdf: ${FAULT_A_SIZE} bytes"
   echo "  station files: ${A_STATIONS}"
   echo ""
   echo "Phase B: RESTART 2.5 → 5.0 s (--restart ${RESULT_DIR_A}/${OUTPUT_PREFIX_A})"
   echo "  Output: ${RESULT_DIR_B}"
   echo "  fault.vtkhdf: ${FAULT_B_SIZE} bytes"
   echo "  station files: ${B_STATIONS}"
   echo ""
   echo "Acceptance checks:"
   echo "  #1 Phase A fault.vtkhdf landed:           PASS (${FAULT_A_SIZE} bytes)"
   echo "  #2 Phase B fault.vtkhdf in own dir:       PASS (${FAULT_B_SIZE} bytes)"
   echo "  #3 Phase A fault.vtkhdf untouched:        PASS"
   echo "  #4 V1 restart block fired on Phase B:     PASS"
   echo "  #5 Phase B station files start near 2.5:  PASS"
   echo "  #6 station files present + separate:      PASS (A=${A_STATIONS} B=${B_STATIONS})"
   echo "  #7 Phase A station md5 unchanged by B:    PASS"
   echo ""
   echo "V3 deferred: per-collection ParaView schedule state (carry the"
   echo "             secondary pv_bulk_out collection's regime +"
   echo "             last_write_time across the seam).  Today V2"
   echo "             carries only the primary pv_out collection."
   echo ""
   echo "du -sh (Phase A): $(du -sh ${RESULT_DIR_A} | cut -f1)"
   echo "du -sh (Phase B): $(du -sh ${RESULT_DIR_B} | cut -f1)"
} > "${SUMMARY_FILE}"
cat "${SUMMARY_FILE}"
