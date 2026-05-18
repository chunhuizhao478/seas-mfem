#!/bin/bash
# tpv104_restart_compare.sh — REFERENCE comparison between the TPV104
# restart pair (Phase A + Phase B) and the single-shot reference run.
#
# Usage:
#   scripts/tpv104_restart_compare.sh <pair_base> <ref_base>
#
# Where:
#   <pair_base> = output dir of the pair sbatch, e.g.
#                 tpv104/results_pair_test_job123456
#                 (must contain segment_001/ and segment_002/)
#   <ref_base>  = output dir of the reference sbatch, e.g.
#                 tpv104/results_reference_job123457
#                 (must contain single_shot/)
#
# Reads the final station row from each run at the central nucleation
# station x2_0_x3_7.5 (TPV104StationWriter column layout:
#   $1 t  $2 h-slip  $3 h-slip-rate  $4 h-shear-stress
#   $5 v-slip  $6 v-slip-rate  $7 v-shear-stress  $8 n-stress  $9 psi).
#
# Compares h-slip, h-slip-rate, and psi between the pair's Phase B
# t=5.0 row and the reference's t=5.0 row.  Passes if both abs_tol
# OR rel_tol are met per variable (whichever is looser passes the
# pair).  This is the gold-standard restart-correctness check; a
# state-corruption bug in the V1 checkpoint (e.g., field-order swap)
# would produce divergence at t=5.0 even though the pair sbatch's
# own #1-#9 file-isolation validations all PASS.
#
# Tolerance design lives in
# debug_document/paraview_output_debug_document/
#   tpv104_restart_state_verification_2026-05-17.md
# (§"Tolerance design" — initial values; tighten/loosen there).
#
# Exit codes:
#   0  PASS — pair's t=5.0 state matches reference within tolerance
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
   echo "                 tpv104/results_pair_test_job<pair_jobid>" >&2
   echo "  <ref_base>   Reference sbatch output dir (contains" >&2
   echo "               single_shot).  Typically:" >&2
   echo "                 tpv104/results_reference_job<ref_jobid>" >&2
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

# Locate the seam station file in each run.  TPV104 station files
# are named <prefix>_station_<station>.dat; both sbatches encode
# the slurm job id in the prefix, so we glob.
SEAM_STATION="x2_0_x3_7.5"

PAIR_B_DIR="${PAIR_BASE}/segment_002"
REF_DIR="${REF_BASE}/single_shot"

PAIR_B_SEAM=$(ls "${PAIR_B_DIR}"/*_station_"${SEAM_STATION}".dat 2>/dev/null \
              | head -1)
REF_SEAM=$(ls "${REF_DIR}"/*_station_"${SEAM_STATION}".dat 2>/dev/null \
           | head -1)

if [ -z "${PAIR_B_SEAM}" ] || [ ! -s "${PAIR_B_SEAM}" ]; then
   echo "ERROR: pair Phase B seam station file missing or empty:" >&2
   echo "       ${PAIR_B_DIR}/*_station_${SEAM_STATION}.dat" >&2
   exit 2
fi
if [ -z "${REF_SEAM}" ] || [ ! -s "${REF_SEAM}" ]; then
   echo "ERROR: reference seam station file missing or empty:" >&2
   echo "       ${REF_DIR}/*_station_${SEAM_STATION}.dat" >&2
   exit 2
fi

echo "================================================================"
echo "TPV104 restart REFERENCE comparison"
echo "  Pair Phase B file: ${PAIR_B_SEAM}"
echo "  Reference file:    ${REF_SEAM}"
echo "  Station:           ${SEAM_STATION}"
echo "================================================================"

# Helpers — mirror the pair sbatch's compare_val / extract_row / col
# implementations so the validation logic stays in sync.
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
         printf "  PASS  %-20s  pair=%+.6e  ref=%+.6e  diff=%.3e  rel=%.3e\n",
                name, a, b, diff, rel;
         exit 0;
      } else {
         printf "  FAIL  %-20s  pair=%+.6e  ref=%+.6e  diff=%.3e  rel=%.3e (abs_tol=%.1e rel_tol=%.1e)\n",
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

# Sanity: both endpoints should be at the same t (≈ 5.0 s).
PAIR_T=$(col "${PAIR_B_LAST}" 1)
REF_T=$(col "${REF_LAST}" 1)
echo "Endpoint times:"
echo "  pair t = ${PAIR_T}"
echo "  ref  t = ${REF_T}"
T_MATCH_OK=$(awk -v a="${PAIR_T}" -v b="${REF_T}" \
             'BEGIN { diff = (a > b) ? a - b : b - a;
                      print (diff <= 0.06) ? "1" : "0" }')
if [ "${T_MATCH_OK}" != "1" ]; then
   echo "  WARN — endpoint times differ by more than one output_dt (0.05)."
   echo "         Comparison may be confounded.  If diff is large,"
   echo "         re-check both sbatches' --tfinal and --output-dt."
fi

echo ""
echo "Column-by-column comparison (tolerances per debug doc):"
ref_fail=0

# h-slip — integrated quantity, smooth.  Tight tolerance.
compare_val "h-slip[t=5.0]" "$(col "${PAIR_B_LAST}" 2)" \
            "$(col "${REF_LAST}" 2)" 1e-3 0.01 || ref_fail=1

# h-slip-rate — evolves fast during rupture; looser rel_tol.
compare_val "h-slip-rate"   "$(col "${PAIR_B_LAST}" 3)" \
            "$(col "${REF_LAST}" 3)" 1e-3 0.05 || ref_fail=1

# psi — directly checkpointed; tightest tolerance (round-trip canary).
compare_val "psi"           "$(col "${PAIR_B_LAST}" 9)" \
            "$(col "${REF_LAST}" 9)" 1e-6 0.001 || ref_fail=1

echo ""
echo "================================================================"
if [ "${ref_fail}" -eq 0 ]; then
   echo "REFERENCE COMPARISON: PASS"
   echo "  Pair's Phase B final state matches the single-shot reference"
   echo "  within tolerance at station ${SEAM_STATION}."
   echo "  Restart preserves state-trajectory correctness."
   echo "================================================================"
   exit 0
else
   echo "REFERENCE COMPARISON: FAIL"
   echo "  At least one column at station ${SEAM_STATION} diverges"
   echo "  between the pair (Phase A + Phase B with restart) and the"
   echo "  single-shot reference.  Restart load may have corrupted"
   echo "  state OR introduced silent drift."
   echo ""
   echo "  Next steps:"
   echo "  - Re-check the comparison at additional stations (try"
   echo "    x2_0_x3_3, x2_0_x3_12 for ranges of nucleation activity)."
   echo "  - Inspect Q vs dof_data round-trip in unit tests:"
   echo "    make seas_test_tpv104_checkpoint && ./seas_test_tpv104_checkpoint"
   echo "  - See debug_document/paraview_output_debug_document/"
   echo "      tpv104_restart_state_verification_2026-05-17.md"
   echo "    §\"Tolerance design\" + §\"Future tightening\" for the"
   echo "    diagnostic menu (binary checkpoint, deterministic reductions,"
   echo "    h5diff cross-check on fault.vtkhdf snapshots)."
   echo "================================================================"
   exit 3
fi
