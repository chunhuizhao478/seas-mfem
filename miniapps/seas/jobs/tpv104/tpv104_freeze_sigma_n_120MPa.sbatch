#!/bin/bash
#SBATCH -J tpv104_frzsn
#SBATCH -o tpv104_frzsn_%j.out
#SBATCH -e tpv104_frzsn_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 σ_n freeze diagnostic — t = 2 s, 200 m, 400 ranks
# ========================================================
#
# Purpose
# -------
# Test whether holding the fault-normal traction at the constant
# background value (120 MPa, SCEC TPV104 spec) — i.e. disabling the
# trial-σ_n update inside the friction Riemann solve — improves the
# rupture-front V_strike profile vs the SeisSol reference.
#
# Background
# ----------
# Post-f_w-fix, the hypocenter peak |V_strike| collapsed to within
# ~0.35% of SeisSol (see tpv104_fw_fix_validation.sbatch verdict).
# Three secondary symptoms remain at the hypocenter:
#   (a) σ_n perturbation ±0.5 MPa during rupture transit
#       (SeisSol shows flat 120 MPa)
#   (b) V_dip leakage ±0.03 m/s
#   (c) slip_dip accumulates ~3.5 mm by t = 2 s
# (a) is upstream of (b)/(c) via the friction strength
#     F_strength = |σ_n_total| · f(ψ, |V|),
# so freezing σ_n_total at the background should remove (a) entirely
# and tells us how much of (b)/(c) is downstream of the σ_n leak vs
# independently driven.
#
# How the freeze is wired
# -----------------------
# Single env-var gate added to dynamic/fault_face_flux.cpp::
# ComputeStageState (the chokepoint reached by both Evaluate and
# EvaluateADER, hence by every TPV104 RK4-stage and ADER-corrector
# fault Riemann solve):
#
#   SEAS_TPV104_FREEZE_SIGMA_N=120e6   →  s.sigma_n_trial = 120 MPa
#   SEAS_TPV104_FREEZE_SIGMA_N=1       →  120 MPa (default)
#   unset / =0                          →  baseline (no override)
#
# Sign convention: TPV104Params::sigma_n is positive-compression in
# the fault-local frame, so 120e6 corresponds to 120 MPa compression.
# The override propagates through CompleteFromTrial → CompleteFromTheta
# → CompleteFromVabs → BuildImposedState → WriteBackState, so the
# friction solver, slip-rate denominator, imposed-state SXX, and
# DOFData.sigma_n_corr station output ALL see the frozen value.
#
# Locally: 25/25 unit tests pass with env unset (byte-identical
# baseline); locked_fault gate-2 worst |sigma_n_corr| pins to
# 1.200e+08 Pa with =120e6 (was 2.083e+08 baseline drift).
#
# Mesh / dispatch parity with the fwfix run
# -----------------------------------------
#   - tpv104_200m.msh (same as fwfix and t12 baseline → mesh asymmetry
#     contribution to σ_n leak is HELD CONSTANT; only the freeze gate
#     changes between this run and fwfix).
#   - 8 nodes × 400 ranks (same).
#   - --ader-order 2, --order 1, --cfl 0.5 (same).
#   - --tfinal 2.0 (same window, direct trace overlay vs fwfix).
#   - --friction-solver newton-stable, --fric-law slip-srw,
#     --fault-iterator substep (banner-only; production path is
#     Brent + one-shot ADER + macro-step ψ per R7-001(b)).
#
# The ONLY production-path change between this run and the fwfix run
# is the SEAS_TPV104_FREEZE_SIGMA_N env var.
#
# Run economics
# -------------
#   tfinal = 2.0 s; expected wall ~30 min at 400 ranks (matches fwfix
#   wall — the freeze adds one branch+strtod per ComputeStageState
#   call, negligible).  2 h dev wall is comfortable margin.
#
# Acceptance criterion (post-run, applied by hand)
# ------------------------------------------------
# Per-station peak |V_strike| comparison vs the SeisSol reference and
# vs the fwfix baseline:
#
#   freeze improves over fwfix
#     → σ_n leak was a contributor to the residual off-hypocenter
#       over-shoot; investigate root cause of σ_n leak (mesh non-mirror
#       vs imposed-state back-flux) next.
#   freeze unchanged vs fwfix
#     → σ_n leak is downstream of the residual (b)/(c) drivers, not a
#       cause; revert freeze and look elsewhere (V_dip / slip_dip
#       channels).
#   freeze worse than fwfix
#     → σ_n was self-correcting; freezing breaks the friction
#       feedback that was helping.  Revert.
#
# Pre-requisites (verified at startup)
# ------------------------------------
#   1. tpv104/mesh/tpv104_200m.msh exists.
#   2. tpv104_init_stations.sbatch PASS — confirms ψ_ini gate (the
#      freeze does not touch initialization, so this remains valid).
#   3. dynamic/fault_face_flux.cpp carries the freeze gate (verified
#      below by grep for the env-var name).
#
# R7-001 disclosure
# -----------------
# Production-path runtime is Brent + one-shot ADER + macro-step ψ
# regardless of CLI flags.  --friction-solver / --fault-iterator
# values are parsed but banner-only.

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

# Verify the freeze gate landed in the source tree.  If this grep
# misses, the source on Frontera predates the diagnostic and the run
# would silently behave as baseline.
if ! grep -q 'SEAS_TPV104_FREEZE_SIGMA_N' dynamic/fault_face_flux.cpp; then
   echo "ERROR: dynamic/fault_face_flux.cpp does not carry the SEAS_TPV104_FREEZE_SIGMA_N"
   echo "       gate.  git pull on Frontera before submitting."
   exit 1
fi
echo "freeze gate confirmed in dynamic/fault_face_flux.cpp"

# Confirm f_w fix is still in place (we want this to compose with f_w
# = 0.2, not silently regress).
F_W_LINE=$(grep -E '^\s*static\s+constexpr\s+real_t\s+f_w' config/tpv104_params.hpp | head -1)
if ! echo "${F_W_LINE}" | grep -q '0\.2'; then
   echo "ERROR: config/tpv104_params.hpp f_w is NOT 0.2 — pull latest source."
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

RESULT_DIR="tpv104/results_freezesn_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

# Force rebuild so the new ComputeStageState carrying the freeze gate
# is linked into the driver — the gate is plain runtime code (no
# preprocessor flag), so a default `make` rebuild is sufficient.
rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o
make seas_tpv104_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

if [ ! -s tpv104/mesh/tpv104_200m.msh ]; then
   echo "ERROR: tpv104/mesh/tpv104_200m.msh missing."
   exit 1
fi

# === FREEZE ENABLE =====================================================
# 120 MPa = TPV104 SCEC background normal stress (positive compression).
export SEAS_TPV104_FREEZE_SIGMA_N=120e6
# =======================================================================

echo "========================================"
echo "Git commit:   $(git rev-parse HEAD)"
echo "Git branch:   $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:   ${RESULT_DIR}"
echo "Target:       tfinal = 2.0 s on tpv104_200m.msh at 400 ranks"
echo "Diagnostic:   SEAS_TPV104_FREEZE_SIGMA_N = ${SEAS_TPV104_FREEZE_SIGMA_N}"
echo "              (sigma_n_trial pinned at 120 MPa in ComputeStageState)"
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
      --output-prefix tpv104_frzsn \
      --output-dt 0.05 \
      --verify-dispatch

# Post-run quick-check: extract peak |V_strike| at every TPV104 station
# and report next to the SeisSol reference and the fwfix baseline.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv104_frzsn_station_x2_0_x3_7.5.dat"

PEAK_V_HYPO="N/A"
PEAK_T_HYPO="N/A"
if [ -s "${HYPO_FILE}" ]; then
   read PEAK_V_HYPO PEAK_T_HYPO < <(awk '!/^#/ && NF>=3 {
      v=$3; if (v<0) v=-v;
      if (v>m) { m=v; tm=$1 }
   } END { printf "%g %g\n", m, tm }' "${HYPO_FILE}")
fi

# σ_n_corr flatness check at hypocenter — direct verification the
# freeze actually pinned σ_n.  Column 8 of the station file is
# sigma_n_corr (per tpv104_driver.cpp StationWriter format).
SIGMA_N_RANGE_HYPO="N/A"
if [ -s "${HYPO_FILE}" ]; then
   SIGMA_N_RANGE_HYPO=$(awk '!/^#/ && NF>=8 && $1+0==$1 {
      s=$8;
      if (NR==1 || s<smin) smin=s;
      if (NR==1 || s>smax) smax=s;
   } END {
      if (smin=="" || smax=="") { print "N/A"; exit }
      printf "[min=%.4e max=%.4e span=%.4e]", smin, smax, (smax-smin)
   }' "${HYPO_FILE}")
fi

# SeisSol reference and fwfix-baseline peak |V_strike| per station,
# from benchmark_data/seisol/ and tpv104/results_fwfix_job7677251.
declare -A SEISSOL_PEAKS=(
   ["x2_0_x3_7.5"]="7.39"
   ["x2_-12_x3_3"]="9.93"   ["x2_12_x3_3"]="9.52"
   ["x2_-12_x3_12"]="8.37"  ["x2_12_x3_12"]="8.28"
   ["x2_-9_x3_7.5"]="6.29"  ["x2_9_x3_7.5"]="6.33"
   ["x2_0_x3_3"]="5.14"     ["x2_0_x3_12"]="11.40"
)

OFF_HYPO_TABLE=""
for tag in x2_0_x3_7.5 x2_0_x3_3 x2_0_x3_12 \
           x2_-9_x3_7.5 x2_9_x3_7.5 \
           x2_-12_x3_3 x2_12_x3_3 \
           x2_-12_x3_12 x2_12_x3_12; do
   SF="${RESULT_DIR}/tpv104_frzsn_station_${tag}.dat"
   if [ -s "${SF}" ]; then
      MV=$(awk '!/^#/ && NF>=3 && $1+0==$1 {v=$3;if(v<0)v=-v;if(v>m)m=v} END{print m+0}' "${SF}")
      SV="${SEISSOL_PEAKS[$tag]}"
      RAT=$(awk "BEGIN{printf \"%.2f\", $MV/$SV}")
      OFF_HYPO_TABLE+=$(printf "  %-18s  MFEM=%6.3f  SeisSol=%5.2f  ratio=%s\n" "$tag" "$MV" "$SV" "$RAT")
      OFF_HYPO_TABLE+=$'\n'
   fi
done

# Dispatch parity check (rank-0 [dispatch] line from --verify-dispatch).
DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv104_frzsn_${SLURM_JOB_ID}.out"
DISPATCH_LINES=$(grep '\[dispatch\] rank=0' "${DISPATCH_LOG}" 2>/dev/null \
                 || echo "MISSING")
DISPATCH_OK=0
if echo "${DISPATCH_LINES}" | grep -q 'solver=brent' \
   && echo "${DISPATCH_LINES}" | grep -q 'iterator=oneshot' \
   && echo "${DISPATCH_LINES}" | grep -q 'fric-law=slip-srw'; then
   DISPATCH_OK=1
fi

# Acceptance: freeze should reproduce ψ-equilibrium and not blow up.
# The σ_n flatness check gives a quantitative confirmation the override
# fired at the hypocenter station.  The 1 kPa span tolerance absorbs
# fp accumulation in the Riemann velocity-jump term, which still rides
# on dynamic Q[VX] even when sigma_n_corr is pinned.
SIGMA_FLAT_OK=0
if [ "${SIGMA_N_RANGE_HYPO}" != "N/A" ]; then
   SPAN=$(echo "${SIGMA_N_RANGE_HYPO}" | grep -oE 'span=[0-9eE.+-]+' | sed 's/span=//')
   if [ -n "${SPAN}" ] && awk "BEGIN { exit !(${SPAN} <= 1e3) }"; then
      SIGMA_FLAT_OK=1
   fi
fi

{
   echo "TPV104 σ_n freeze diagnostic run"
   echo "Mesh:                  tpv104_200m"
   echo "tfinal:                2.0 s"
   echo "freeze:                SEAS_TPV104_FREEZE_SIGMA_N=${SEAS_TPV104_FREEZE_SIGMA_N}"
   echo "f_w (compiled):        $(echo "${F_W_LINE}" | grep -oE '[01]\.[0-9]+' | head -1)"
   echo ""
   echo "peak |V_strike| at hypo  : ${PEAK_V_HYPO} m/s at t = ${PEAK_T_HYPO} s"
   echo "SeisSol reference (hypo) : 7.39 m/s at t = 1.07 s"
   echo "fwfix baseline (hypo)    : 7.42 m/s (~0.35% over SeisSol)"
   echo ""
   echo "sigma_n_corr at hypo     : ${SIGMA_N_RANGE_HYPO}"
   if [ "${SIGMA_FLAT_OK}" -eq 1 ]; then
      echo "  freeze fired           : PASS (span ≤ 1 kPa)"
   else
      echo "  freeze fired           : FAIL — gate did not pin σ_n; check env var propagation through ibrun"
   fi
   echo ""
   echo "Per-station peak |V_strike| comparison (lower ratio = better fit to SeisSol):"
   echo "${OFF_HYPO_TABLE}"
   echo "[dispatch] rank=0 lines  : ${DISPATCH_LINES}"
   if [ "${DISPATCH_OK}" -eq 1 ]; then
      echo "dispatch parity vs fwfix : PASS (Brent + oneshot + slip-srw)"
   else
      echo "dispatch parity vs fwfix : FAIL — code path drifted"
   fi
   echo ""
   echo "Interpretation guide:"
   echo "  per-station ratios ≤ fwfix ratios  → σ_n leak was a contributor;"
   echo "                                       investigate σ_n leak root cause next."
   echo "  per-station ratios ≈ fwfix ratios  → σ_n leak is downstream / not causal;"
   echo "                                       revert freeze, look at V_dip / slip_dip."
   echo "  per-station ratios > fwfix ratios  → freezing breaks helpful friction"
   echo "                                       feedback; revert freeze."
   echo ""
   if [ "${SIGMA_FLAT_OK}" -eq 1 ] && [ "${DISPATCH_OK}" -eq 1 ] \
      && [ "${PEAK_V_HYPO}" != "N/A" ]; then
      echo "STATUS: PASS"
   else
      echo "STATUS: FAIL"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
