#!/bin/bash
# check_dyn_gate.sh — the 200-step Nitsche gate of the GPU dynamic-rupture
# plan (D0 onwards), the counterpart of check_gate.sh.
#
# A gate run is
#   seas_nitsche_dyn_driver --config DECK --mesh MESH --order P \
#       --max-steps 200 --output-dir RUN_DIR --prefix RUN_PREFIX  > run.log 2>&1
# with the deck's stations written every step, and this script decides
# whether it passed.  The run passes only if
#   1. scripts/nitsche_gate_check.py passes at TOL (the same non-empty station
#      set as the golden, every line parsed, finite values, increasing times,
#      matching row count and times, every station field's error finite and
#      <= TOL);
#   2. the run exited with code 0;
#   3. its log has the driver's "=== done:" line.
#
# Usage:
#   check_dyn_gate.sh RUN_DIR RUN_PREFIX GOLDEN_DIR GOLDEN_PREFIX TOL RUN_LOG RUN_EXIT_CODE [N_STATIONS]
# Run from miniapps/seas (or set SEAS_DIR).  N_STATIONS, when given, is the
# number of stations the golden must hold.  Exit 0 = PASS, 1 = FAIL,
# 2 = usage error.  Writes RUN_DIR/RUN_PREFIX_gate.json.

set -u
if [ "$#" -ne 7 ] && [ "$#" -ne 8 ]; then
   sed -n '2,22p' "$0" >&2
   exit 2
fi
RUN_DIR="$1"; RUN_PREFIX="$2"; GOLDEN_DIR="$3"; GOLDEN_PREFIX="$4"
TOL="$5"; RUN_LOG="$6"; RUN_RC="$7"; N_STATIONS="${8:-}"
SEAS_DIR="${SEAS_DIR:-$(pwd)}"
CHECKER="${SEAS_DIR}/scripts/nitsche_gate_check.py"

fail=0
echo "=== check_dyn_gate: ${RUN_DIR}/${RUN_PREFIX} vs ${GOLDEN_DIR}/${GOLDEN_PREFIX}, TOL ${TOL} ==="

expect=()
[ -n "${N_STATIONS}" ] && expect=(--expect-stations "${N_STATIONS}")
if ! python3 "${CHECKER}" --run-dir "${RUN_DIR}" --run-prefix "${RUN_PREFIX}" \
        --golden-dir "${GOLDEN_DIR}" --golden-prefix "${GOLDEN_PREFIX}" \
        --tol "${TOL}" "${expect[@]}" --json "${RUN_DIR}/${RUN_PREFIX}_gate.json"; then
   echo "check_dyn_gate: FAIL (file validation or errors)"
   fail=1
fi

if [ "${RUN_RC}" != "0" ]; then
   echo "check_dyn_gate: FAIL (run exit code ${RUN_RC})"
   fail=1
fi

if [ ! -f "${RUN_LOG}" ]; then
   echo "check_dyn_gate: FAIL (no run log ${RUN_LOG})"
   fail=1
else
   grep -q "^=== done: " "${RUN_LOG}" \
      || { echo "check_dyn_gate: FAIL (log lacks the '=== done:' line)"; fail=1; }
fi

if [ "${fail}" -eq 0 ]; then
   echo "check_dyn_gate: PASS"
fi
exit "${fail}"
