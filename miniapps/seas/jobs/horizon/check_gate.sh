#!/bin/bash
# check_gate.sh — the 50-step BP5 gate of the GPU plan (G0 onwards).
#
# A gate run is
#   seas_bp5_full <v63 flags> --order P --max-steps 50 --write-every-step \
#       --ref-dir <golden dir> --ref-prefix <golden prefix> \
#       --regression-tolerance TOL  > run.log 2>&1
# and this script decides whether it passed.  The built-in comparator alone
# is not enough (it skips a missing reference, drops lines it cannot parse,
# and passes a NaN error), so the files are validated first by
# bp5/bp5_gate_check.py, which also recomputes every error with the
# comparator's formula.  The run passes only if
#   1. bp5_gate_check.py passes at TOL (full station set and global file in
#      run and golden, every line parsed, finite values, increasing times,
#      matching start time and row count, every error finite and <= TOL,
#      the end time's relative difference included, matching event count);
#   2. the run exited with code 0;
#   3. its log says "Stations compared: 10 / 10", "REGRESSION CHECK: PASS"
#      and "Event count: N [ok]".
#
# Usage:
#   check_gate.sh RUN_DIR RUN_PREFIX GOLDEN_DIR GOLDEN_PREFIX TOL RUN_LOG RUN_EXIT_CODE
# Run from miniapps/seas (or set SEAS_DIR).  Exit 0 = PASS, 1 = FAIL,
# 2 = usage error.  Writes RUN_DIR/RUN_PREFIX_gate.json.

set -u
if [ "$#" -ne 7 ]; then
   sed -n '2,23p' "$0" >&2
   exit 2
fi
RUN_DIR="$1"; RUN_PREFIX="$2"; GOLDEN_DIR="$3"; GOLDEN_PREFIX="$4"
TOL="$5"; RUN_LOG="$6"; RUN_RC="$7"
SEAS_DIR="${SEAS_DIR:-$(pwd)}"
CHECKER="${SEAS_DIR}/bp5/bp5_gate_check.py"

fail=0
echo "=== check_gate: ${RUN_DIR}/${RUN_PREFIX} vs ${GOLDEN_DIR}/${GOLDEN_PREFIX}, TOL ${TOL} ==="

if ! python3 "${CHECKER}" --run-dir "${RUN_DIR}" --run-prefix "${RUN_PREFIX}" \
        --golden-dir "${GOLDEN_DIR}" --golden-prefix "${GOLDEN_PREFIX}" \
        --tol "${TOL}" --json "${RUN_DIR}/${RUN_PREFIX}_gate.json"; then
   echo "check_gate: FAIL (file validation or errors)"
   fail=1
fi

if [ "${RUN_RC}" != "0" ]; then
   echo "check_gate: FAIL (run exit code ${RUN_RC})"
   fail=1
fi

if [ ! -f "${RUN_LOG}" ]; then
   echo "check_gate: FAIL (no run log ${RUN_LOG})"
   fail=1
else
   grep -q "Stations compared: 10 / 10" "${RUN_LOG}" \
      || { echo "check_gate: FAIL (log lacks 'Stations compared: 10 / 10')"; fail=1; }
   grep -q "REGRESSION CHECK: PASS" "${RUN_LOG}" \
      || { echo "check_gate: FAIL (log lacks 'REGRESSION CHECK: PASS')"; fail=1; }
   grep -Eq "Event count: [0-9]+ \[ok\]" "${RUN_LOG}" \
      || { echo "check_gate: FAIL (log lacks a matching event count)"; fail=1; }
fi

if [ "${fail}" -eq 0 ]; then
   echo "check_gate: PASS"
fi
exit "${fail}"
