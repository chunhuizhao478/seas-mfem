#!/bin/bash
#SBATCH -J tpv104_dthalf
#SBATCH -o tpv104_dthalf_%j.out
#SBATCH -e tpv104_dthalf_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 dt-halving diagnostic — R7-001(b) cadence audit
# ======================================================
#
# Purpose
# -------
# Determine whether the macro-step ψ-cadence (R7-001 option (b),
# documented in drivers/tpv104_driver.cpp:1023-1037 + the I-17 entry of
# debug_document/tpv104_debug_document/research_seissol_tpv104_2026-04-24.md)
# is the dominant source of the 1.7-2.7× peak-V over-shoot observed
# against the SeisSol reference at the propagating-front stations.
#
# Hypothesis (per tpv104_psi_cadence_review_2026-04-24.md §C-1):
# the macro-vs-per-sub-step ψ deviation is bounded by O(dt_macro²) (Simpson-
# rule / ADER-consistency argument applied to dψ/dt = -(V/L)(ψ - ψ_ss(V))).
# So a 5.7× reduction in dt_macro should reduce the peak-V over-shoot
# by ≈ 32× IF the cadence gap is dominant.
#
# Run parameters
# --------------
#   --dt 5e-5         # 5.7× smaller than the t3 production dt = 2.85e-4 s
#   --tfinal 1.5      # spans nucleation (t < 1 s) + early breakout
#   --ader-order 2    # matches the t3 production run for apples-to-apples
#   --verify-dispatch # so [dispatch] lines land in the log; Phase-3 probe-
#                     # diff reads these to confirm the production code
#                     # path is unchanged (Brent + one-shot ADER + macro-
#                     # step ψ + macro-step nucleation).
#
# Cost
# ----
#   tfinal / dt = 1.5 / 5e-5 = 30,000 macro steps.  At ~30 k macro
#   steps/h on 400 ranks (extrapolating from the persistent-nuc 400r
#   t=12 production observed cadence in jobs/tpv104/README.md), this
#   fits in the 2 h development queue.  Submit with -p development;
#   the production -p normal queue is reserved for the full 12 s run.
#
# Acceptance criteria — applied post-run by the analyst
# ------------------------------------------------------
# Compare the peak |V_strike| at each of the 9 fault stations between
# this run and the production t3 run, and against the SeisSol reference
# at the same wall-clock (use tpv104/scripts/probe_diff.py + the
# visualize_results.py script for plots).  Apply the table from
# tpv104_psi_cadence_review_2026-04-24.md §C-1:
#
#   over-shoot reduction   |  conclusion
#   ≥ 25×                  |  macro-dt cadence is dominant — wire iterator
#   4× – 25×               |  cadence is one of several errors; also audit
#                          |  nucleation cadence (I-04) and friction
#                          |  dispatch (I-15)
#   < 4×                   |  cadence is NOT dominant — investigate
#                          |  elsewhere; iterator wiring will not fix it
#
# Pre-requisites — same as tpv104_t3_mesh_coupled.sbatch:
#   1. tpv104_build_banner_sanity.sbatch  (PASS)
#   2. tpv104_mesh_build.sbatch           (produced tpv104_200m.msh)
#   3. tpv104_init_stations.sbatch        (PASS — P2_D ψ_ini gate)
#
# Reads after the run
# -------------------
# The post-run summary block writes RESULT.txt with:
#   - the 9 station files' row counts and NaN-freeness (gate, like t3)
#   - the actual dispatched dispatch summary lines (so a future audit
#     reading this run's logs sees that the [dispatch] line confirmed
#     Brent + one-shot ADER + macro-step ψ — i.e. the SAME code path as
#     the baseline production t3 run, only with smaller dt).
#
# This script does NOT touch any [C2] file.  The only change vs the t3
# template is `--dt`, `--tfinal`, output prefix, and the `--verify-dispatch`
# flag.  The build, dispatched code path, and post-run gates are
# byte-identical to the production run.

export LC_ALL=C
export LANG=C

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas

# Integration-readiness gate — same as t3 template: require a passing
# init record so the dt-halving run starts from a verified ψ_ini state.
PREV_INIT_RESULT=$(ls -dt tpv104/results_init_stations_job* 2>/dev/null | head -1)
if [ -z "${PREV_INIT_RESULT}" ] \
   || ! grep -q "^STATUS: PASS" "${PREV_INIT_RESULT}/RESULT.txt" 2>/dev/null; then
   echo "ERROR: no passing tpv104_init_stations.sbatch result found."
   echo "       Run that job first and verify STATUS: PASS before proceeding."
   exit 1
fi
echo "Found init PASS at: ${PREV_INIT_RESULT}"

RESULT_DIR="tpv104/results_dthalf_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

# Force a clean rebuild of every TPV104 unit so we know this run's
# binaries match the current source tree.  Same scope as the t3 template.
rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o
make seas_tpv104_driver -j8 || { echo "BUILD FAILED"; exit 1; }

if [ ! -s tpv104/mesh/tpv104_200m.msh ]; then
   echo "ERROR: tpv104/mesh/tpv104_200m.msh missing."
   exit 1
fi

echo "========================================"
echo "Git commit:   $(git rev-parse HEAD)"
echo "Output dir:   ${RESULT_DIR}"
echo "Target:       tfinal = 1.5 s, dt = 5e-5 s (30,000 macro steps)"
echo "Mesh:         tpv104_200m.msh at 400 ranks"
echo "Diagnostic:   R7-001(b) cadence audit (dt-halving)"
echo "========================================"

# `--verify-dispatch` emits one [dispatch] line per rank disclosing the
# actual production-path solver/iterator/friction-law dispatch.  The
# post-run summary parses these for cadence-audit attribution.
ibrun ./seas_tpv104_driver \
      --mesh tpv104/mesh/tpv104_200m.msh \
      --tfinal 1.5 \
      --dt 5e-5 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law slip-srw \
      --fault-iterator substep \
      --verify-dispatch \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv104_dthalf

# Post-run sanity: same gate as the t3 production run + a [dispatch]
# parity check that confirms this run's code path matches the t3
# baseline (same Brent + one-shot ADER + macro-step ψ + macro-step
# nucleation; only dt differs).
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
STATIONS_NONEMPTY=0
STATIONS_NANFREE=0
STATIONS_TOTAL=0
for f in "${RESULT_DIR}"/tpv104_dthalf_station_x2_*.dat; do
   STATIONS_TOTAL=$((STATIONS_TOTAL + 1))
   if [ $(wc -l < "${f}") -ge 50 ]; then
      STATIONS_NONEMPTY=$((STATIONS_NONEMPTY + 1))
   fi
   if ! grep -qiE "nan|inf" "${f}"; then
      STATIONS_NANFREE=$((STATIONS_NANFREE + 1))
   fi
done

# Extract the rank-0 [dispatch] line from the stdout log (the most
# recent .out file in SBATCH_LOG_DIR matching this job).  If the
# production code path silently changed (CLI typo, default flip, etc.)
# the gate fails loudly.
DISPATCH_LINE=$(grep -m1 '\[dispatch\] rank=0' \
                "${SBATCH_LOG_DIR}/tpv104_dthalf_${SLURM_JOB_ID}.out" \
                2>/dev/null || echo "MISSING")
DISPATCH_OK=0
if echo "${DISPATCH_LINE}" | grep -q 'solver=brent' \
   && echo "${DISPATCH_LINE}" | grep -q 'iterator=oneshot' \
   && echo "${DISPATCH_LINE}" | grep -q 'fric-law=slip-srw'; then
   DISPATCH_OK=1
fi

# Peak |V_strike| at the hypocenter (x2_0_x3_7.5) — the headline
# acceptance number.  Column 3 of the SCEC trace layout is h-slip-rate.
HYPO_FILE="${RESULT_DIR}/tpv104_dthalf_station_x2_0_x3_7.5.dat"
PEAK_V_HYPO="N/A"
if [ -s "${HYPO_FILE}" ]; then
   PEAK_V_HYPO=$(awk 'NF>=9 && $1+0==$1 {v=$3<0?-$3:$3; if(v>m)m=v} END{print m}' \
                 "${HYPO_FILE}")
fi

{
   echo "TPV104 dt-halving diagnostic run"
   echo "dt                          : 5e-5 s"
   echo "tfinal                      : 1.5 s"
   echo "macro steps                 : 30000"
   echo ""
   echo "fault stations total        : ${STATIONS_TOTAL}"
   echo "fault stations >= 50 rows   : ${STATIONS_NONEMPTY}"
   echo "fault stations NaN/Inf-free : ${STATIONS_NANFREE}"
   echo ""
   echo "[dispatch] rank=0 line      : ${DISPATCH_LINE}"
   if [ "${DISPATCH_OK}" -eq 1 ]; then
      echo "dispatch parity vs t3       : PASS (Brent + oneshot + slip-srw)"
   else
      echo "dispatch parity vs t3       : FAIL — code path drifted from baseline"
   fi
   echo ""
   echo "peak |V_strike| at hypo     : ${PEAK_V_HYPO} m/s"
   echo "  (compare to t3 baseline + SeisSol reference; see"
   echo "   tpv104_psi_cadence_review_2026-04-24.md §C-1 acceptance table)"
   echo ""
   if [ "${STATIONS_TOTAL}" -eq 9 ] \
      && [ "${STATIONS_NONEMPTY}" -eq 9 ] \
      && [ "${STATIONS_NANFREE}" -eq 9 ] \
      && [ "${DISPATCH_OK}" -eq 1 ]; then
      echo "STATUS: PASS"
   else
      echo "STATUS: FAIL"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
