#!/bin/bash
#SBATCH -J tpv205_mfnone_p2_O3
#SBATCH -o tpv205_mfnone_p2_O3_%j.out
#SBATCH -e tpv205_mfnone_p2_O3_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# TPV205 (SCEC TPV5) baseline — MIXED FLUX NONE, P=2 + ADER-O3.
# =======================================================================
# Clean-control / baseline pair to
#   tpv205_mixed_flux_adjacent_200m_p2_O3_normal.sbatch
# Same rebuilt TPV205 driver after the R-001 round-5 fixes (REVIEW.md
# 2026-04-27 final review):
#   - DOFData::lsw_mu_s/lsw_mu_d/lsw_d_c LSW-native fields (Phase 1)
#   - FaultFaceFlux::EvaluateADER_LSW (no slip accumulation; iterator
#     owns slip evolution exclusively — REVIEW R-001)
#   - WaveOperator::SetFaultFrictionLaw(LSW) wired in driver
#   - Two LSW dispatch branches in wave_operator.inl (interior +
#     R-1600/R-1601 shared-fault fallback)
#
# Goal
# ----
# Verify that TPV205 (linear slip-weakening, SCEC TPV5) at P=2 + ADER-O3
# without mixed-flux ADJACENT is STABLE through tfinal=12 s.  This is
# the third leg of the cross-friction-law diagnosis matrix:
#   - TPV104 SRW + adjacent + O3 (job 7680797): BLOWS UP
#   - TPV104 SRW + none     + O3 (job 7680798): CLEAN
#   - TPV102 ageing-law RS + adjacent + O3:     BLOWS UP (user 2026-04-28)
#   - TPV102 ageing-law RS + none + O3:         pending (companion sbatch)
#   - TPV205 LSW + adjacent + O3:               pending
#   - TPV205 LSW + none + O3:                   THIS JOB
#
# If TPV205 LSW + none completes cleanly → the failure mode is purely
# bulk-side (mixed-flux ADJACENT × ADER-O3 in wave_operator.inl), and
# the LSW closed-form physics is a CLEAN fault solver under the upwind-
# everywhere baseline.  Combined with the LSW-version of the adjacent
# job, the matrix becomes friction-law-independent on the failure axis.
#
# If TPV205 LSW + none ALSO blows up → there is a TPV205-specific
# regression in the LSW dispatch path independent of mixed-flux.  This
# would be a R-001 / R-016 round-5 regression and require triage of
# the LSW-dispatch wiring (FaultFrictionLaw::LSW reach into the wave
# operator, EvaluateADER_LSW correctness, R-1600 shared-fault fallback).
#
# Wall-time budget
# ----------------
# Normal queue 48 h envelope.  At P=2 + ADER-O3 + upwind everywhere on
# the 200 m mesh at 500 ranks, no mixed-flux CFL reduction → ~2550 macro
# steps to t=12 s.  At ~14 s/step → ~9.9 h wall.  48 h envelope holds
# with margin for the bulk ParaView I/O.
#
# Output cadence (matches the broken adjacent companion for clean diff)
# ---------------------------------------------------------------------
# --output-dt 0.05                  # 16 SCEC fault stations every 0.05 s
# --paraview --paraview-dt 0.5      # fault-surface VTU/PVD every 0.5 s
# --paraview-bulk-dt 0.5            # bulk volume VTU/PVD every 0.5 s
#
# R-1204 / mixed-flux ADER-O note
# -------------------------------
# `--mixed-flux none` does not engage central flux on any face, so the
# R-1204 driver guard does not fire and no override env var is needed.
# `SEAS_FORCE_MIXED_FLUX_ADER_O_GT2` is intentionally UNSET so any
# environment leak would be caught.
#
# R-1503 / mixed-flux + substep iterator note
# -------------------------------------------
# Job uses `--fault-iterator one-shot` (O=1 quadrature inside the
# Tpv205SubStepIterator).  R-1503 is silent and no override is needed.
# All TPV205 dispatch — interior + shared-fault — flows through
# AdvanceADERWithSubStep regardless of the --fault-iterator value (LSW
# is not implemented in the legacy rate-and-state wave.AdvanceADER
# fault branch).

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

if [ ! -s tpv205/mesh/tpv2053d_200m.msh ]; then
   echo "ERROR: tpv205/mesh/tpv2053d_200m.msh missing.  Run gmsh on tpv205/mesh/tpv2053d_200m.geo first."
   exit 1
fi

RESULT_DIR="tpv205/results_mixed_flux_none_p2_O3_normal_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv205_driver drivers/tpv205_driver.o \
      dynamic/tpv205_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o dynamic/wave_operator.o
make seas_tpv205_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

# `--mixed-flux none` does not require the R-1204 override.  Leave the
# env var UNSET so any environment leak would not silently affect this job.
unset SEAS_FORCE_MIXED_FLUX_ADER_O_GT2
unset SEAS_FORCE_MIXED_FLUX_SUBSTEP_MPI

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          tpv2053d_200m.msh"
echo "Polynomial:    --order 2 + --ader-order 3"
echo "Iterator:      one-shot (Tpv205SubStepIterator with O=1 quadrature)"
echo "Mixed flux:    none (upwind everywhere; baseline / clean control)"
echo "Friction law:  lsw (linear slip-weakening; closed-form, SCEC TPV5)"
echo "Queue:         normal (10N x 500r, 48h)"
echo "tfinal:        12.0 s (full SCEC TPV5 §III window)"
echo "ParaView:      surface 0.5 s, bulk 0.5 s"
echo "R-1204 ovrd:   UNSET (not needed at --mixed-flux none)"
echo "========================================"

ibrun ./seas_tpv205_driver \
      --mesh tpv205/mesh/tpv2053d_200m.msh \
      --mesh-scale 1 \
      --order 2 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 12.0 \
      --ader-order 3 \
      --friction-solver newton-stable \
      --fric-law lsw \
      --fault-iterator one-shot \
      --mixed-flux none \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv205_mfnone_p2_O3 \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.5 \
      --paraview-bulk-dt 0.5 \
      --verify-dispatch

# Post-run summary — mirrors the TPV104 mixed_flux_none normal-queue
# baseline.  No mixed-flux engagement gate (mode=none → no [mixed-flux]
# banner expected); no R-1407 override-trail block (no [WARNING] line
# expected at --mixed-flux none).  R-016 LSW dispatch verdict retained
# (the LSW closed-form must dispatch correctly even without mixed-flux).
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv205_mfnone_p2_O3_station_x2_0_x3_7.5.dat"

PEAK_V_HYPO="N/A"
PEAK_T_HYPO="N/A"
SIGMA_N_PERTURB_HYPO="N/A"
SLIP_STRIKE_END="N/A"
T_END="N/A"
N_ROWS="N/A"
if [ -s "${HYPO_FILE}" ]; then
   # Column layout (tpv205_setup.hpp::TPV205StationWriter::WriteStep):
   #   1: time   2: h-slip(strike=slip2)   3: h-slip-rate(V2)
   #   4: h-shear-stress(tau2_corr)        5: v-slip(dip=slip1)
   #   6: v-slip-rate(V1)                  7: v-shear-stress(tau1_corr)
   #   8: n-stress(sigma_n_corr)           9: mu_eff(δ)
   read PEAK_V_HYPO PEAK_T_HYPO < <(awk '!/^#/ && NF>=3 {
      v=$3; if (v<0) v=-v;
      if (v>m) { m=v; tm=$1 }
   } END { printf "%g %g\n", m, tm }' "${HYPO_FILE}")
   SIGMA_N_PERTURB_HYPO=$(awk '!/^#/ && NF>=8 && $1+0==$1 {
      d=$8 - 120e6; if (d<0) d=-d;
      if (d>m) m=d
   } END { printf "%.4e", m+0 }' "${HYPO_FILE}")
   SLIP_STRIKE_END=$(awk '!/^#/ && NF>=2 && $1+0==$1 {a=$2} END{printf "%.4e", a+0}' "${HYPO_FILE}")
   T_END=$(awk '!/^#/ && NF>=1 && $1+0==$1 {tt=$1} END{printf "%.4f", tt+0}' "${HYPO_FILE}")
   N_ROWS=$(awk '!/^#/ && NF>=1 && $1+0==$1' "${HYPO_FILE}" | wc -l | tr -d ' ')
fi

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv205_mfnone_p2_O3_${SLURM_JOB_ID}.out"
DISPATCH_LINES=$(grep '\[dispatch\] rank=0' "${DISPATCH_LOG}" 2>/dev/null \
                 || echo "MISSING")
MIXED_FLUX_BANNER=$(grep -E '\[mixed-flux\] mode=' "${DISPATCH_LOG}" 2>/dev/null \
                    | head -1 || echo "MISSING")

# REVIEW R-016 sanity: the dispatch banner must report the LSW closed-
# form solver, not the rate-and-state Brent path.  A REGRESSION verdict
# means the SetFaultFrictionLaw(LSW) wiring did not reach the wave
# operator's fault dispatch and the run silently produced rate-and-state
# physics on the LSW DOFData — which would invalidate this baseline as a
# clean control.
LSW_DISPATCH_VERDICT="UNKNOWN"
if grep -q "friction_solver_actual=lsw-closed-form" "${DISPATCH_LOG}" 2>/dev/null; then
   LSW_DISPATCH_VERDICT="ENGAGED — friction_solver_actual=lsw-closed-form"
elif grep -q "friction_solver_actual=brent" "${DISPATCH_LOG}" 2>/dev/null; then
   LSW_DISPATCH_VERDICT="REGRESSION — banner reports brent (R-016 fix not in effect)"
else
   LSW_DISPATCH_VERDICT="UNKNOWN — no friction_solver_actual line in log"
fi

# tfinal completion gate — normal queue should reach t=12 s.
TFINAL_VERDICT="UNKNOWN"
if [ "${T_END}" != "N/A" ]; then
   if awk "BEGIN { exit !(${T_END} >= 11.95) }"; then
      TFINAL_VERDICT="PASS — reached tfinal"
   elif awk "BEGIN { exit !(${T_END} >= 8.0) }"; then
      TFINAL_VERDICT="PARTIAL — tfinal ≥ 8 s (post-rupture)"
   else
      TFINAL_VERDICT="WEAK — tfinal < 8 s"
   fi
fi

# Sanity gate: peak |V_strike| should stay in physical range.  TPV205
# analytic LSW initial slip rate at the nucleation patch is 0.0779 m/s;
# rupture peaks at a few m/s as the front sweeps the fault.  Anything
# ≥ 100 m/s means this baseline run also went unstable — would
# invalidate the comparison and indicate a deeper issue.
SANE_V_VERDICT="UNKNOWN"
if [ "${PEAK_V_HYPO}" != "N/A" ]; then
   if awk "BEGIN { exit !(${PEAK_V_HYPO} <= 100.0) }"; then
      SANE_V_VERDICT="PASS — peak |V| in physical range"
   else
      SANE_V_VERDICT="FAIL — peak |V| > 100 m/s; baseline ALSO blew up"
   fi
fi

{
   echo "TPV205 mixed-flux baseline (P=2 + ADER-O3 + --mixed-flux none) NORMAL"
   echo "Mesh:                    tpv2053d_200m.msh"
   echo "Polynomial order:        --order 2 + --ader-order 3"
   echo "Iterator:                one-shot"
   echo "Mixed flux:              none (upwind everywhere; clean control)"
   echo "Friction law:            lsw (linear slip-weakening; SCEC TPV5)"
   echo "Queue / resources:       normal, 10N x 500r, 48 h"
   echo "tfinal target:           12.0 s"
   echo "tfinal achieved (hypo):  ${T_END} s   (${TFINAL_VERDICT})"
   echo "hypo station rows:       ${N_ROWS}"
   echo ""
   echo "peak |V_strike| at hypo  : ${PEAK_V_HYPO} m/s at t = ${PEAK_T_HYPO} s   (${SANE_V_VERDICT})"
   echo "peak |σ_n − 120 MPa| @ h : ${SIGMA_N_PERTURB_HYPO} Pa"
   echo "slip_strike at hypo (end): ${SLIP_STRIKE_END} m"
   echo ""
   echo "Mixed-flux banner        : ${MIXED_FLUX_BANNER}"
   echo "[dispatch] rank=0 lines  : ${DISPATCH_LINES}"
   echo "LSW dispatch (R-016)     : ${LSW_DISPATCH_VERDICT}"
   echo ""
   if [ "${PEAK_V_HYPO}" != "N/A" ] \
      && [ "${SIGMA_N_PERTURB_HYPO}" != "N/A" ] \
      && echo "${TFINAL_VERDICT}" | grep -qE "^PASS|^PARTIAL" \
      && echo "${SANE_V_VERDICT}" | grep -q "^PASS" \
      && echo "${LSW_DISPATCH_VERDICT}" | grep -q "^ENGAGED"; then
      echo "STATUS: PASS — clean control + LSW closed-form active; diff vs --mixed-flux adjacent companion"
   else
      echo "STATUS: FAIL — see verdicts above"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
