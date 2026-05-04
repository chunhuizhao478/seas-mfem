#!/bin/bash
#SBATCH -J tpv104_mfadj_p2_O2
#SBATCH -o tpv104_mfadj_p2_O2_%j.out
#SBATCH -e tpv104_mfadj_p2_O2_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# TPV104 mixed-flux ADER-O × ADJACENT discriminator —
# MIXED FLUX ADJACENT, P=2 + ADER-O2.
# =======================================================================
# Discriminator job for the spurious-rupture instability observed at
#   tpv104_mixed_flux_adjacent_200m_p2_O3_normal job 7680797
# (point-source hot spots during rupture propagation, t ≈ 6.5 s on 200 m
# mesh + 500 ranks; Brent NaN spam at shared-fault QPs; cross-confirmed
# on TPV102 P=2 + ADER-O3 + adjacent).  The companion clean run
#   tpv104_mixed_flux_none_200m_p2_O3_normal job 7680798
# (same P=2, same ADER-O3, mixed-flux NONE) ruled out:
#   - ADER-O3 alone as the destabilizer
#   - MPI rank-boundary topology at np=500
# Cross-friction-law evidence (TPV104 SRW + TPV102 ageing-law RS both
# blow up at the same configuration) ruled out friction-solver-specific
# explanations (Newton-asinh vs Brent vs constitutive law).
#
# Two candidates remain (REVIEW.md DRDG3D comparison
# document/system_dev/drdg3d_mixed_flux_comparison_2026-04-28.md):
#
#   (1) ADER-O3's cubic-in-time Cauchy–Kowalewski predictor injects
#       high-frequency content into the central-flux RHS that the
#       non-dissipative Mixed-Flux 2 dispatch cannot damp.  RK54 (the
#       DRDG3D reference) does not have an analogous polynomial-in-time
#       reconstruction → DRDG3D's TPV5 results are stable at
#       p=2 + RK54 + mixed-flux ADJACENT.
#   (2) Mixed-flux ADJACENT itself is unstable at p ≥ 2 (independent of
#       ADER order), via either:
#         (a) the wider element-based face set vs DRDG3D's node-based
#             face set, or
#         (b) the central-flux dispatch's underdamping on tetrahedral
#             cells at p ≥ 2 with the (v, σ) non-conservative form.
#
# This job pins (1) vs (2):
#   - If TPV104 P=2 + ADER-O2 + adjacent runs CLEAN to t = 12 s →
#     Candidate 1 (ADER-O3 × ADJACENT interaction) is the root cause.
#     The fix lives in the ADER predictor / dispatch interface.
#   - If TPV104 P=2 + ADER-O2 + adjacent ALSO BLOWS UP →
#     Candidate 2 (mixed-flux ADJACENT itself unstable at p ≥ 2) is the
#     root cause.  The fix lives in the central-flux dispatch (face set
#     tightening, alpha-blend regulariser, or modal filter).
#
# Either outcome rules out half the candidate space and tells us where
# to direct the next round of work.
#
# Diff vs the broken sbatch
# (tpv104_mixed_flux_adjacent_200m_p2_O3_normal.sbatch):
#   --ader-order 3 → 2
#   --output-prefix tpv104_mfadj_p2_O3 → tpv104_mfadj_p2_O2
#   RESULT_DIR / DISPATCH_LOG / job-name suffix p2_O3 → p2_O2
#   SEAS_FORCE_MIXED_FLUX_ADER_O_GT2 → unset (R-1204 guard does NOT
#     fire at ADER-O ≤ 2; this is the validated combination).
#   R-1407 verdict block removed (no override warning to expect).
#
# All other parameters identical to the broken sbatch:
#   - Same mesh tpv104_200m.msh
#   - Same P=2, same --mixed-flux adjacent
#   - Same CFL 0.5, same tfinal 12.0
#   - Same partition (default ParMETIS, no fault-locality flag)
#   - Same ParaView cadence
#   - Same module loads + LD_LIBRARY_PATH
#   - Same 10 N × 500 r normal-queue 48 h envelope
#
# Wall-time budget
# ----------------
# Normal queue 48 h envelope.  ADER-O2 has ~30% lower per-step cost
# than ADER-O3 (cubic predictor → quadratic predictor) but a slightly
# tighter dt (R-1502 mixed-flux CFL × Pelties O2 stability bound) gives
# a similar number of macro steps as the O3 baseline.  Net effect:
# wall-time ≈ 0.7 × 11.8 h ≈ 8.3 h ≪ 48 h.

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

F_W_LINE=$(grep -E '^\s*static\s+constexpr\s+real_t\s+f_w' config/tpv104_params.hpp | head -1)
if ! echo "${F_W_LINE}" | grep -q '0\.2'; then
   echo "ERROR: f_w is NOT 0.2"; exit 1
fi
echo "f_w confirmed: ${F_W_LINE}"

if [ ! -s tpv104/mesh/tpv104_200m.msh ]; then
   echo "ERROR: tpv104/mesh/tpv104_200m.msh missing."
   exit 1
fi

RESULT_DIR="tpv104/results_mixed_flux_adjacent_p2_O2_normal_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o dynamic/wave_operator.o
make seas_tpv104_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

unset SEAS_TPV104_FREEZE_SIGMA_N
unset SEAS_FORCE_MIXED_FLUX_SUBSTEP_MPI
# R-1204: --mixed-flux adjacent is VALIDATED at ADER-O2.  No override
# needed and the env var is intentionally NOT set so any environment
# leak would be caught (the driver guard would abort cleanly).
unset SEAS_FORCE_MIXED_FLUX_ADER_O_GT2

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          tpv104_200m.msh"
echo "Polynomial:    --order 2 + --ader-order 2 (Pelties O2; R-1204 validated)"
echo "Iterator:      one-shot (R7-001 dispatched; --fault-iterator one-shot)"
echo "Mixed flux:    adjacent (Mixed-Flux 2, Zhang 2023)"
echo "R-1204 ovrd:   UNSET (not needed at ADER-O ≤ 2)"
echo "Queue:         normal (10N x 500r, 48h)"
echo "tfinal:        12.0 s (full production tfinal target)"
echo "ParaView:      surface 0.5 s, bulk 0.5 s"
echo "Discriminator: pins ADER-O3 × ADJACENT (Cand. 1) vs"
echo "               ADJACENT-itself-unstable-at-p≥2 (Cand. 2)"
echo "========================================"

ibrun ./seas_tpv104_driver \
      --mesh tpv104/mesh/tpv104_200m.msh \
      --mesh-scale 1 \
      --order 2 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 12.0 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law slip-srw \
      --fault-iterator one-shot \
      --mixed-flux adjacent \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv104_mfadj_p2_O2 \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.5 \
      --paraview-bulk-dt 0.5 \
      --verify-dispatch

# Post-run summary — same template as the broken O3 sbatch, minus the
# R-1407 override-trail block (no warning expected at ADER-O2).
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv104_mfadj_p2_O2_station_x2_0_x3_7.5.dat"

PEAK_V_HYPO="N/A"
PEAK_T_HYPO="N/A"
SIGMA_N_PERTURB_HYPO="N/A"
SLIP_DIP_END="N/A"
T_END="N/A"
N_ROWS="N/A"
if [ -s "${HYPO_FILE}" ]; then
   read PEAK_V_HYPO PEAK_T_HYPO < <(awk '!/^#/ && NF>=3 {
      v=$3; if (v<0) v=-v;
      if (v>m) { m=v; tm=$1 }
   } END { printf "%g %g\n", m, tm }' "${HYPO_FILE}")
   SIGMA_N_PERTURB_HYPO=$(awk '!/^#/ && NF>=8 && $1+0==$1 {
      d=$8 - 120e6; if (d<0) d=-d;
      if (d>m) m=d
   } END { printf "%.4e", m+0 }' "${HYPO_FILE}")
   SLIP_DIP_END=$(awk '!/^#/ && NF>=5 && $1+0==$1 {a=$5} END{printf "%.4e", a+0}' "${HYPO_FILE}")
   T_END=$(awk '!/^#/ && NF>=1 && $1+0==$1 {tt=$1} END{printf "%.4f", tt+0}' "${HYPO_FILE}")
   N_ROWS=$(awk '!/^#/ && NF>=1 && $1+0==$1' "${HYPO_FILE}" | wc -l | tr -d ' ')
fi

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv104_mfadj_p2_O2_${SLURM_JOB_ID}.out"
DISPATCH_LINES=$(grep '\[dispatch\] rank=0' "${DISPATCH_LOG}" 2>/dev/null \
                 || echo "MISSING")
MIXED_FLUX_BANNER=$(grep -E '\[mixed-flux\] mode=' "${DISPATCH_LOG}" 2>/dev/null \
                    | head -1 || echo "MISSING")
MIXED_FLUX_RANK_LINES=$(grep -c '\[mixed-flux\] rank=' "${DISPATCH_LOG}" 2>/dev/null \
                        || echo 0)

# Engagement gate.
MIXED_FLUX_VERDICT="UNKNOWN"
if echo "${MIXED_FLUX_BANNER}" | grep -q "mode=adjacent"; then
   GSUM=$(echo "${MIXED_FLUX_BANNER}" \
          | grep -oE '_global=([0-9]+)' | head -1 \
          | grep -oE '[0-9]+')
   if [ -z "${GSUM}" ]; then
      GSUM=$(echo "${MIXED_FLUX_BANNER}" \
             | grep -oE 'global=([0-9]+)' | head -1 \
             | grep -oE '[0-9]+')
   fi
   if [ -n "${GSUM}" ] && [ "${GSUM}" -gt 0 ]; then
      MIXED_FLUX_VERDICT="ENGAGED — |central_set|_global=${GSUM}"
   else
      MIXED_FLUX_VERDICT="MISSING — banner present but global count is 0/missing"
   fi
else
   MIXED_FLUX_VERDICT="NOT-ENGAGED — no [mixed-flux] mode=adjacent line in log"
fi

# tfinal completion gate — same thresholds as the O3 broken job.
TFINAL_VERDICT="UNKNOWN"
if [ "${T_END}" != "N/A" ]; then
   if awk "BEGIN { exit !(${T_END} >= 11.95) }"; then
      TFINAL_VERDICT="PASS — reached tfinal"
   elif awk "BEGIN { exit !(${T_END} >= 8.0) }"; then
      TFINAL_VERDICT="PARTIAL — tfinal ≥ 8 s (post-rupture)"
   else
      TFINAL_VERDICT="WEAK — tfinal < 8 s (likely the same instability)"
   fi
fi

# Discriminator gate: distinguishes Candidate 1 (ADER-O3 × ADJACENT) vs
# Candidate 2 (ADJACENT-itself-unstable-at-p≥2) by comparing the
# hypocenter peak |V_strike| against a sane physical bound.  The broken
# O3 run reached |V| > 1e11 m/s by t ≈ 7.8 s; the clean none-O3 run
# stays below ~10 m/s throughout.  Use 100 m/s as the discriminator
# threshold (well above any physical TPV104 hypocenter slip rate yet
# well below the O3-broken-run scale).
DISCRIMINATOR_VERDICT="UNKNOWN"
if [ "${PEAK_V_HYPO}" != "N/A" ]; then
   if awk "BEGIN { exit !(${PEAK_V_HYPO} <= 100.0) }"; then
      DISCRIMINATOR_VERDICT="CAND-1 — ADER-O3 × ADJACENT is the destabilizer (this O2 run is clean)"
   else
      DISCRIMINATOR_VERDICT="CAND-2 — ADJACENT itself is unstable at p≥2 (peak |V| > 100 m/s, same pathology as O3 job 7680797)"
   fi
fi

{
   echo "TPV104 mixed-flux ADJACENT (P=2 + ADER-O2 + --mixed-flux adjacent) NORMAL [DISCRIMINATOR]"
   echo "Mesh:                    tpv104_200m.msh"
   echo "Polynomial order:        --order 2 + --ader-order 2 (R-1204 validated)"
   echo "Iterator:                one-shot"
   echo "Mixed flux:              adjacent (Mixed-Flux 2, Zhang 2023)"
   echo "R-1204 override:         UNSET (not needed at ADER-O ≤ 2)"
   echo "Queue / resources:       normal, 10N x 500r, 48 h"
   echo "tfinal target:           12.0 s"
   echo "tfinal achieved (hypo):  ${T_END} s   (${TFINAL_VERDICT})"
   echo "hypo station rows:       ${N_ROWS}"
   echo ""
   echo "peak |V_strike| at hypo  : ${PEAK_V_HYPO} m/s at t = ${PEAK_T_HYPO} s"
   echo "peak |σ_n − 120 MPa| @ h : ${SIGMA_N_PERTURB_HYPO} Pa"
   echo "slip_dip at hypo (t_end) : ${SLIP_DIP_END} m"
   echo ""
   echo "Mixed-flux banner        : ${MIXED_FLUX_BANNER}"
   echo "Mixed-flux engagement    : ${MIXED_FLUX_VERDICT}"
   echo "Mixed-flux rank lines    : ${MIXED_FLUX_RANK_LINES}"
   echo "[dispatch] rank=0 lines  : ${DISPATCH_LINES}"
   echo ""
   echo "DISCRIMINATOR VERDICT    : ${DISCRIMINATOR_VERDICT}"
   echo ""
   if [ "${PEAK_V_HYPO}" != "N/A" ] \
      && [ "${SIGMA_N_PERTURB_HYPO}" != "N/A" ] \
      && echo "${MIXED_FLUX_VERDICT}" | grep -q "^ENGAGED"; then
      if echo "${DISCRIMINATOR_VERDICT}" | grep -q "^CAND-1"; then
         echo "STATUS: PASS — clean run; ADER-O3 × ADJACENT identified as root cause"
      elif echo "${DISCRIMINATOR_VERDICT}" | grep -q "^CAND-2"; then
         echo "STATUS: FAIL-AS-EXPECTED — ADJACENT itself unstable; central-flux dispatch is root cause"
      else
         echo "STATUS: INCONCLUSIVE — see verdicts above"
      fi
   else
      echo "STATUS: FAIL — see verdicts above"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
