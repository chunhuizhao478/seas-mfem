#!/bin/bash
# bp5_restart_compare.sh — REFERENCE comparison between the BP5 V2
# restart pair (Phase A + Phase B) and the single-shot reference run.
#
# Usage:
#   scripts/bp5_restart_compare.sh <pair_base> <ref_base>
#
# Where:
#   <pair_base> = output dir of the pair sbatch, e.g.
#                 bp5/results_restart_test_job<pair_jobid>
#                 (must contain segment_001/ and segment_002/)
#   <ref_base>  = output dir of the reference sbatch, e.g.
#                 bp5/results_reference_job<ref_jobid>
#                 (must contain single_shot/)
#
# Reads the final station row from each run at the central
# fault station `fltst_strk+00dp+10` (X=0 along-strike, depth=10km
# in the seismogenic zone — active state evolution).
#
# BP5 station file column layout (bp5_benchmark_output.hpp:427-431):
#   $1 t(s)  $2 slip_strike(m)  $3 slip_dip(m)
#   $4 log10(V_strike)(m/s)  $5 log10(V_dip)(m/s)
#   $6 tau_strike(MPa)  $7 tau_dip(MPa)  $8 log10(state)(s)
#
# Compares slip_strike, log10(V_strike), log10(state) between the
# pair's Phase B final row and the reference's final row.  Passes
# if abs_tol OR rel_tol are met per variable.  This is the gold-
# standard restart-correctness check for BP5 V2 — a state-corruption
# bug in the V2 checkpoint (e.g., field-order swap in the trailing
# block) would produce divergence at t ≈ 3 yr even though the pair
# sbatch's own #1-#9 file-isolation validations all PASS.
#
# Tolerance design (initial values; tighten if real Frontera runs
# show tighter agreement):
#   slip_strike:     abs 1e-3 m, rel 1%
#   log10(V_strike): abs 0.05 (in log10 space — ~12% relative on V)
#   log10(state):    abs 0.01 (in log10 space — ~2.3% relative on state)
#
# Exit codes:
#   0  PASS — pair's t≈3yr state matches reference within tolerance
#   1  USAGE error or missing input dirs
#   2  station file missing in pair or reference
#   3  REFERENCE divergence beyond tolerance on at least one column
#
# Designed to run on the login node after BOTH sbatches finish; no
# MPI, no module loads, no Slurm needed.  Pure bash + awk.

set -euo pipefail

if [ "$#" -ne 2 ]; then
   echo "Usage: $0 <pair_base> <ref_base>" >&2
   echo "" >&2
   echo "  <pair_base>  Pair sbatch output dir (contains segment_001 +" >&2
   echo "               segment_002).  Typically:" >&2
   echo "                 bp5/results_restart_test_job<pair_jobid>" >&2
   echo "  <ref_base>   Reference sbatch output dir (contains" >&2
   echo "               single_shot).  Typically:" >&2
   echo "                 bp5/results_reference_job<ref_jobid>" >&2
   exit 1
fi

PAIR_BASE=$1
REF_BASE=$2

if [ ! -d "${PAIR_BASE}" ]; then
   echo "ERROR: pair_base directory does not exist: ${PAIR_BASE}" >&2
   exit 1
fi
if [ ! -d "${REF_BASE}" ]; then
   echo "ERROR: ref_base directory does not exist: ${REF_BASE}" >&2
   exit 1
fi

# Station file: fltst_strk+00dp+10 (central, 10 km depth, in
# seismogenic zone).  BP5 uses .txt extension.
SEAM_STATION="fltst_strk+00dp+10"

PAIR_B_DIR="${PAIR_BASE}/segment_002"
REF_DIR="${REF_BASE}/single_shot"

PAIR_B_SEAM=$(ls "${PAIR_B_DIR}"/*_"${SEAM_STATION}".txt 2>/dev/null \
              | head -1)
REF_SEAM=$(ls "${REF_DIR}"/*_"${SEAM_STATION}".txt 2>/dev/null \
           | head -1)

if [ -z "${PAIR_B_SEAM}" ] || [ ! -s "${PAIR_B_SEAM}" ]; then
   echo "ERROR: pair Phase B seam station file missing or empty:" >&2
   echo "       ${PAIR_B_DIR}/*_${SEAM_STATION}.txt" >&2
   exit 2
fi
if [ -z "${REF_SEAM}" ] || [ ! -s "${REF_SEAM}" ]; then
   echo "ERROR: reference seam station file missing or empty:" >&2
   echo "       ${REF_DIR}/*_${SEAM_STATION}.txt" >&2
   exit 2
fi

echo "================================================================"
echo "BP5 V2 restart REFERENCE comparison"
echo "  Pair Phase B file: ${PAIR_B_SEAM}"
echo "  Reference file:    ${REF_SEAM}"
echo "  Station:           ${SEAM_STATION}"
echo "================================================================"

# Helpers — mirror the TPV104 compare script + pair sbatch's helpers.
extract_row() {
   local file=$1 which=$2
   if [ "$which" = "first" ]; then
      awk '/^[[:space:]]*[0-9-]/ { print; exit }' "$file"
   else
      awk '/^[[:space:]]*[0-9-]/ { last=$0 } END { print last }' "$file"
   fi
}

col() {
   echo "$1" | awk -v c="$2" '{ print $c }'
}

compare_val() {
   local name=$1 a=$2 b=$3 abs_tol=$4 rel_tol=$5
   awk -v name="$name" -v a="$a" -v b="$b" \
       -v abs_tol="$abs_tol" -v rel_tol="$rel_tol" '
   BEGIN {
      diff = (a > b) ? a - b : b - a;
      mag_a = (a >= 0) ? a : -a;
      mag_b = (b >= 0) ? b : -b;
      mag = (mag_a > mag_b) ? mag_a : mag_b;
      rel = (mag > 1e-30) ? (diff / mag) : 0;
      if (diff <= abs_tol || rel <= rel_tol) {
         printf "  PASS  %-22s  pair=%+.6e  ref=%+.6e  diff=%.3e  rel=%.3e\n",
                name, a, b, diff, rel;
         exit 0;
      } else {
         printf "  FAIL  %-22s  pair=%+.6e  ref=%+.6e  diff=%.3e  rel=%.3e (abs_tol=%.1e rel_tol=%.1e)\n",
                name, a, b, diff, rel, abs_tol, rel_tol;
         exit 1;
      }
   }'
}

PAIR_B_LAST=$(extract_row "${PAIR_B_SEAM}" last)
REF_LAST=$(extract_row "${REF_SEAM}" last)

echo ""
echo "Last data rows:"
echo "  pair B: ${PAIR_B_LAST}"
echo "  ref:    ${REF_LAST}"
echo ""

# Sanity: both endpoints should be at the same t (≈ 3 yr = 9.467e7 s).
PAIR_T=$(col "${PAIR_B_LAST}" 1)
REF_T=$(col "${REF_LAST}" 1)
echo "Endpoint times (s):"
echo "  pair t = ${PAIR_T}"
echo "  ref  t = ${REF_T}"
# BP5 output cadence is adaptive (paraview-dt-inter-yr 1.0 → ~yr in
# interseismic), so endpoints may differ by ~yr ≈ 3.16e7 s.  Warn
# if they differ by more than 2 yr.
T_MATCH_OK=$(awk -v a="${PAIR_T}" -v b="${REF_T}" \
             'BEGIN { diff = (a > b) ? a - b : b - a;
                      print (diff <= 6.3e7) ? "1" : "0" }')
if [ "${T_MATCH_OK}" != "1" ]; then
   echo "  WARN — endpoint times differ by more than 2 yr."
   echo "         Comparison may be confounded by trajectory phase"
   echo "         offset.  Check both sbatches' --tfinal."
fi

echo ""
echo "Column-by-column comparison (tolerances per debug doc):"
ref_fail=0

# slip_strike — integrated quantity, smooth.  Tight tolerance.
compare_val "slip_strike[t≈3yr]" "$(col "${PAIR_B_LAST}" 2)" \
            "$(col "${REF_LAST}" 2)" 1e-3 0.01 || ref_fail=1

# log10(V_strike) — slip-rate in log space; absolute tolerance in
# log10 units (~0.05 = 12% relative on V_strike itself).
compare_val "log10(V_strike)" "$(col "${PAIR_B_LAST}" 4)" \
            "$(col "${REF_LAST}" 4)" 0.05 0.05 || ref_fail=1

# log10(state) — state variable, directly checkpointed.  Tightest
# tolerance; 0.01 log10 units ≈ 2.3% on state itself.
compare_val "log10(state)" "$(col "${PAIR_B_LAST}" 8)" \
            "$(col "${REF_LAST}" 8)" 0.01 0.01 || ref_fail=1

echo ""
echo "================================================================"
if [ "${ref_fail}" -eq 0 ]; then
   echo "REFERENCE COMPARISON: PASS"
   echo "  Pair's Phase B final state matches the single-shot reference"
   echo "  within tolerance at station ${SEAM_STATION}."
   echo "  V2 restart preserves state-trajectory correctness."
   echo "================================================================"
   exit 0
else
   echo "REFERENCE COMPARISON: FAIL"
   echo "  At least one column at station ${SEAM_STATION} diverges"
   echo "  between the pair (Phase A + Phase B with V2 restart) and the"
   echo "  single-shot reference.  Restart load may have corrupted"
   echo "  state OR introduced silent drift."
   echo ""
   echo "  Next steps:"
   echo "  - Try other central stations (fltst_strk+00dp+00,"
   echo "    fltst_strk+00dp+22) to confirm divergence is global."
   echo "  - Inspect V2 trailing-block round-trip in unit tests:"
   echo "    make seas_test_bp5_petsc_ts_restart && ./seas_test_bp5_petsc_ts_restart"
   echo "  - See debug_document/paraview_output_debug_document/"
   echo "      tpv104_restart_state_verification_2026-05-17.md"
   echo "    §\"Tolerance design\" + §\"Future tightening\" for the"
   echo "    diagnostic menu (binary checkpoint, deterministic reductions,"
   echo "    h5diff cross-check on fault.vtkhdf snapshots)."
   echo "================================================================"
   exit 3
fi
