#!/bin/bash
#SBATCH -J tpv104_fwfix
#SBATCH -o tpv104_fwfix_%j.out
#SBATCH -e tpv104_fwfix_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 f_w fix validation — t = 2 s, 200 m, 400 ranks
# =====================================================
#
# Purpose
# -------
# Validate the f_w = 0.1 → 0.2 fix landed in `config/tpv104_params.hpp`
# (commit corresponding to tpv104_v_overshoot_divergence_2026-04-24.md
# follow-up).  The SCEC TPV104 spec table in
# `benchmark_document/SCEC_validation_slip_law.pdf` page 2 lists
# f_w = 0.2; MFEM had 0.1 (half), which produced systematic 1.7-2.7×
# V_strike over-shoot vs the SeisSol reference traces in
# `benchmark_data/seisol/`.
#
# Hypothesis: with the corrected f_w, the V_strike peak at the
# hypocenter (t ≈ 1.0-1.1 s) should drop from MFEM's t12 baseline of
# 9.45 m/s toward SeisSol's reference of 7.39 m/s (1.28× over-shoot).
# Off-hypocenter propagating-front stations show a larger 1.2-2.7×
# over-shoot in the t12 baseline; if cleanly attributable to f_w, the
# corrected run should bring all 9 stations within ~10% of SeisSol's
# values (per measurements just before commit):
#   x2_0_x3_7.5     MFEM=9.45  SeisSol=7.39   1.28×
#   x2_-12_x3_3    MFEM=18.70  SeisSol=9.93   1.88×
#   x2_+12_x3_3    MFEM=18.14  SeisSol=9.52   1.91×
#   x2_-12_x3_12   MFEM=16.11  SeisSol=8.37   1.92×
#   x2_+12_x3_12   MFEM=14.89  SeisSol=8.28   1.80×
#   x2_-9_x3_7.5   MFEM=13.83  SeisSol=6.29   2.20×
#   x2_+9_x3_7.5   MFEM=14.65  SeisSol=6.33   2.31×
#   x2_0_x3_3       MFEM=13.97  SeisSol=5.14   2.72×
#   x2_0_x3_12      MFEM=13.67  SeisSol=11.40  1.20×
#
# Mesh / dispatch parity with the t12 dev baseline
# ------------------------------------------------
#   - tpv104_200m.msh (same mesh as t12 baseline, so the mesh-asymmetry
#     contribution to σ_n leak is HELD CONSTANT — only f_w changes).
#   - 8 nodes × 400 ranks (same).
#   - --ader-order 2, --order 1, --cfl 0.5 (same).
#   - --friction-solver newton-stable, --fault-iterator substep,
#     --fric-law slip-srw  (banner-only; production path uses Brent +
#     one-shot ADER + macro-step ψ per R7-001(b)).
#
# The ONLY production-path change between this run and the t12 baseline
# is the f_w constant inside the rebuilt seas_tpv104_driver.
#
# Run economics
# -------------
#   tfinal = 2.0 s; expected wall ~30 min at 400 ranks (extrapolating
#   from the t12 dev cadence ~6 s/macro-step → ~3000 macro-steps in
#   t = 2 s ≈ 5 h, but station+driver-only output (no paraview) is
#   substantially cheaper than the t12 dev with paraview — should
#   land closer to 1 h wall).  2 h dev wall is comfortable margin.
#
# Acceptance criterion (post-run, applied by hand)
# ------------------------------------------------
# Compute peak |V_strike| at the hypocenter from
# `${RESULT_DIR}/tpv104_fwfix_station_x2_0_x3_7.5.dat` column 3 and
# compare against SeisSol's 7.39 m/s reference:
#
#   peak_V_hypo > 9.0 m/s   →  unchanged from t12 baseline (9.45);
#                              f_w fix did not collapse over-shoot,
#                              investigate τ_2_trial instrumentation
#                              + dt-halving + mesh asymmetry next.
#   8.0 < peak_V ≤ 9.0      →  partial fix; f_w bug accounts for some
#                              of the over-shoot but other
#                              contributors remain.  Also check
#                              off-hypocenter stations (where t12
#                              baseline shows 1.8-2.7× over-shoot vs
#                              hypocenter's 1.28× — if those drop
#                              substantially, f_w was the dominant
#                              contributor at those stations).
#   peak_V ≤ 8.0 m/s        →  full fix at hypocenter (within ~10% of
#                              SeisSol's 7.39).  Confirm at
#                              off-hypocenter stations; expect them
#                              also to drop into the SeisSol-reference
#                              ±10% band.
#
# Pre-requisites (verified at startup)
# ------------------------------------
#   1. tpv104/mesh/tpv104_200m.msh exists.
#   2. tpv104_init_stations.sbatch PASS — confirms ψ_ini gate at
#      station QPs matches the SeisSol reference 5.6359184e-01.
#      (P2_D ψ_ini gate is independent of f_w and remains valid.)
#
# R7-001 disclosure
# -----------------
# Production-path runtime is Brent + one-shot ADER + macro-step ψ
# regardless of CLI flags.  --friction-solver / --fault-iterator
# values are parsed but banner-only.  See jobs/tpv104/README.md and
# debug_document/tpv104_debug_document/research_seissol_tpv104_2026-04-24.md.

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

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas

# Verify the rebuilt code carries the f_w = 0.2 fix.
F_W_LINE=$(grep -E '^\s*static\s+constexpr\s+real_t\s+f_w' config/tpv104_params.hpp | head -1)
if ! echo "${F_W_LINE}" | grep -q '0\.2'; then
   echo "ERROR: config/tpv104_params.hpp f_w is NOT 0.2 — fix did not land."
   echo "       got: ${F_W_LINE}"
   exit 1
fi
echo "f_w confirmed: ${F_W_LINE}"

PREV_INIT_RESULT=$(ls -dt tpv104/results_init_stations_job* 2>/dev/null | head -1)
if [ -z "${PREV_INIT_RESULT}" ] \
   || ! grep -q "^STATUS: PASS" "${PREV_INIT_RESULT}/RESULT.txt" 2>/dev/null; then
   echo "ERROR: no passing tpv104_init_stations.sbatch result found."
   exit 1
fi
echo "Found init PASS at: ${PREV_INIT_RESULT}"

RESULT_DIR="tpv104/results_fwfix_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o
make seas_tpv104_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

if [ ! -s tpv104/mesh/tpv104_200m.msh ]; then
   echo "ERROR: tpv104/mesh/tpv104_200m.msh missing."
   exit 1
fi

echo "========================================"
echo "Git commit:   $(git rev-parse HEAD)"
echo "Git branch:   $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:   ${RESULT_DIR}"
echo "Target:       tfinal = 2.0 s on tpv104_200m.msh at 400 ranks"
echo "Goal:         validate f_w = 0.2 fix; expect peak V_hypo ≤ 3.8 m/s"
echo "========================================"

ibrun ./seas_tpv104_driver \
      --mesh tpv104/mesh/tpv104_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 2.0 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law slip-srw \
      --fault-iterator substep \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv104_fwfix \
      --output-dt 0.05 \
      --verify-dispatch

# Post-run quick-check: extract peak |V_strike| at the hypocenter
# station (x2_0_x3_7.5) and report next to SeisSol's reference value.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv104_fwfix_station_x2_0_x3_7.5.dat"

PEAK_V_HYPO="N/A"
PEAK_T_HYPO="N/A"
if [ -s "${HYPO_FILE}" ]; then
   read PEAK_V_HYPO PEAK_T_HYPO < <(awk '!/^#/ && NF>=3 {
      v=$3; if (v<0) v=-v;
      if (v>m) { m=v; tm=$1 }
   } END { printf "%g %g\n", m, tm }' "${HYPO_FILE}")
fi

# SeisSol reference peak at the hypocenter (computed from
# benchmark_data/seisol/tpv104_seisol_x2_0_x3_7.5.txt col 3, max).
SEISSOL_PEAK_V=7.39
SEISSOL_PEAK_T=1.07

# Verdict per acceptance table.  Bands:
#   PASS    : peak ≤ 8.0 m/s  (within ~10% of SeisSol 7.39)
#   PARTIAL : 8.0 < peak ≤ 9.0 m/s  (some improvement vs t12 baseline 9.45)
#   FAIL    : peak > 9.0 m/s  (essentially unchanged from t12 baseline)
VERDICT="UNKNOWN"
if [ "${PEAK_V_HYPO}" != "N/A" ]; then
   if awk "BEGIN { exit !(${PEAK_V_HYPO} <= 8.0) }"; then
      VERDICT="PASS — f_w fix collapses over-shoot to within ~10% of SeisSol"
   elif awk "BEGIN { exit !(${PEAK_V_HYPO} <= 9.0) }"; then
      VERDICT="PARTIAL — f_w fix improves but other contributors remain"
   else
      VERDICT="FAIL — f_w fix did not collapse over-shoot; investigate τ_2_trial"
   fi
fi

# Also report off-hypocenter peaks (over-shoot was 1.8-2.7× there in
# t12 baseline — if those drop close to SeisSol values too, f_w was the
# dominant contributor across the rupture front, not just at the hypo).
declare -A SEISSOL_PEAKS=(
   ["x2_0_x3_7.5"]="7.39"
   ["x2_-12_x3_3"]="9.93"   ["x2_12_x3_3"]="9.52"
   ["x2_-12_x3_12"]="8.37"  ["x2_12_x3_12"]="8.28"
   ["x2_-9_x3_7.5"]="6.29"  ["x2_9_x3_7.5"]="6.33"
   ["x2_0_x3_3"]="5.14"     ["x2_0_x3_12"]="11.40"
)
OFF_HYPO_TABLE=""
for tag in x2_0_x3_3 x2_0_x3_12 x2_-9_x3_7.5 x2_9_x3_7.5 \
           x2_-12_x3_3 x2_12_x3_3 x2_-12_x3_12 x2_12_x3_12; do
   SF="${RESULT_DIR}/tpv104_fwfix_station_${tag}.dat"
   if [ -s "${SF}" ]; then
      MV=$(awk '!/^#/ && NF>=3 && $1+0==$1 {v=$3;if(v<0)v=-v;if(v>m)m=v} END{print m+0}' "${SF}")
      SV="${SEISSOL_PEAKS[$tag]}"
      RAT=$(awk "BEGIN{printf \"%.2f\", $MV/$SV}")
      OFF_HYPO_TABLE+=$(printf "  %-18s  MFEM=%6.3f  SeisSol=%5.2f  ratio=%s\n" "$tag" "$MV" "$SV" "$RAT")
      OFF_HYPO_TABLE+=$'\n'
   fi
done

# Dispatch parity check (rank-0 [dispatch] line from --verify-dispatch).
DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv104_fwfix_${SLURM_JOB_ID}.out"
DISPATCH_LINES=$(grep '\[dispatch\] rank=0' "${DISPATCH_LOG}" 2>/dev/null \
                 || echo "MISSING")
DISPATCH_OK=0
if echo "${DISPATCH_LINES}" | grep -q 'solver=brent' \
   && echo "${DISPATCH_LINES}" | grep -q 'iterator=oneshot' \
   && echo "${DISPATCH_LINES}" | grep -q 'fric-law=slip-srw'; then
   DISPATCH_OK=1
fi

{
   echo "TPV104 f_w fix validation run"
   echo "Mesh:                  tpv104_200m"
   echo "tfinal:                2.0 s"
   echo "f_w (compiled):        $(echo "${F_W_LINE}" | grep -oE '[01]\.[0-9]+' | head -1)"
   echo ""
   echo "peak |V_strike| at hypo  : ${PEAK_V_HYPO} m/s at t = ${PEAK_T_HYPO} s"
   echo "SeisSol reference (hypo) : ${SEISSOL_PEAK_V} m/s at t = ${SEISSOL_PEAK_T} s"
   echo "t12 baseline (hypo)      : 9.45 m/s (over-shoot 1.28× pre-fix)"
   echo "verdict                  : ${VERDICT}"
   echo ""
   echo "Off-hypocenter peak comparison (lower ratio = better fix):"
   echo "${OFF_HYPO_TABLE}"
   echo "[dispatch] rank=0 lines  : ${DISPATCH_LINES}"
   if [ "${DISPATCH_OK}" -eq 1 ]; then
      echo "dispatch parity vs t12   : PASS (Brent + oneshot + slip-srw)"
   else
      echo "dispatch parity vs t12   : FAIL — code path drifted"
   fi
   echo ""
   if echo "${VERDICT}" | grep -q "^PASS" && [ "${DISPATCH_OK}" -eq 1 ]; then
      echo "STATUS: PASS"
   else
      echo "STATUS: FAIL"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
