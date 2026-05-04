#!/bin/bash
#SBATCH -J tpv102_mfnone_p2_O3
#SBATCH -o tpv102_mfnone_p2_O3_%j.out
#SBATCH -e tpv102_mfnone_p2_O3_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# TPV102 mixed-flux baseline — MIXED FLUX NONE, P=2 + ADER-O3.
# =======================================================================
# Clean-control / baseline pair to the broken
#   tpv102_mixed_flux_adjacent_200m_p2_O3_normal.sbatch
# (user 2026-04-28: TPV102 P=2 + ADER-O3 + --mixed-flux adjacent shows
# the same spurious-rupture instability as TPV104 SRW job 7680797 —
# point-source hot spots during rupture propagation, not at the rupture-
# area boundary).  Mirrors the TPV104 mixed-flux NONE clean-control
# pattern (job 7680798) which completed cleanly at the same P=2 + O=3
# with --mixed-flux none.
#
# Goal
# ----
# Verify that TPV102 (regularised RS, ageing law) at P=2 + ADER-O3
# without mixed-flux ADJACENT is STABLE through tfinal=12 s.  Together
# with the broken adjacent companion, this isolates `--mixed-flux
# adjacent` as the responsible factor for TPV102 too — extending the
# cross-friction-law evidence (TPV104 SRW + TPV102 ageing-law RS both
# fail under adjacent, both pass under none) and ruling out friction-
# solver-specific explanations for the failure.
#
# Goal pair
# ---------
# (a) TPV102 hypocenter normal-stress variation through rupture
#     (σ_n perturbation),
# (b) TPV102 hypocenter dip-slip drift at t = 12 s,
# (c) compare against the broken --mixed-flux adjacent companion to
#     quantify the central-flux dispatch's pollution at TPV102 spec
#     parameters (τ_ini=75 MPa, σ_n=120 MPa, V_ini=1e-12 m/s,
#     Δτ₀=25 MPa, R=3 km, T=1 s).
#
# Wall-time budget
# ----------------
# Normal queue 48 h envelope.  At P=2 + ADER-O3 + upwind everywhere on
# the 200 m mesh at 500 ranks, no mixed-flux CFL reduction → ~2550
# macro steps to t=12 s.  At ~14 s/step → ~9.9 h wall.  48 h envelope
# holds with ample margin for the bulk ParaView I/O.
#
# Output cadence (matches the broken adjacent companion for clean diff)
# ---------------------------------------------------------------------
# --output-dt 0.05                  # station files every 0.05 s
# --paraview --paraview-dt 0.5      # fault-surface VTU/PVD every 0.5 s
# --paraview-bulk-dt 0.5            # bulk volume VTU/PVD every 0.5 s
#
# R-1204 / mixed-flux ADER-O note
# -------------------------------
# `--mixed-flux none` does not engage central flux on any face, so the
# R-1204 driver guard does not fire and no override env var is needed.
# `SEAS_FORCE_MIXED_FLUX_ADER_O_GT2` is intentionally UNSET so any
# environment leak would be caught.

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

if [ ! -s tpv102/mesh/tpv102_200m.msh ]; then
   echo "ERROR: tpv102/mesh/tpv102_200m.msh missing."
   exit 1
fi

RESULT_DIR="tpv102/results_mixed_flux_none_p2_O3_normal_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv102_driver drivers/tpv102_driver.o \
      dynamic/tpv102_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o dynamic/wave_operator.o
make seas_tpv102_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

# `--mixed-flux none` does not require the R-1204 override.  Leave the
# env var UNSET so any environment leak (e.g. from a previous mixed-flux
# adjacent run on the same node) would not silently affect this job.
unset SEAS_FORCE_MIXED_FLUX_ADER_O_GT2
unset SEAS_FORCE_MIXED_FLUX_SUBSTEP_MPI

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          tpv102_200m.msh"
echo "Polynomial:    --order 2 + --ader-order 3"
echo "Iterator:      one-shot (legacy wave.AdvanceADER dispatch)"
echo "Mixed flux:    none (upwind everywhere; baseline / clean control)"
echo "Friction law:  aging (SCEC TPV101/102, dθ/dt = 1 − Vθ/L)"
echo "Queue:         normal (10N x 500r, 48h)"
echo "tfinal:        12.0 s (full SCEC TPV101/102 window)"
echo "ParaView:      surface 0.5 s, bulk 0.5 s"
echo "R-1204 ovrd:   UNSET (not needed at --mixed-flux none)"
echo "========================================"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 2 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 12.0 \
      --ader-order 3 \
      --friction-solver newton-stable \
      --fric-law aging \
      --fault-iterator one-shot \
      --mixed-flux none \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv102_mfnone_p2_O3 \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.5 \
      --paraview-bulk-dt 0.5 \
      --verify-dispatch

# Post-run summary — mirrors the TPV104 mixed_flux_none normal-queue
# baseline.  Quick-check pattern: hypocenter peak |V_strike|, peak
# |σ_n − 120 MPa|, slip2(t_end).  No mixed-flux engagement gate (mode=none
# → no [mixed-flux] banner expected); no R-1407 override-trail block
# (no [WARNING] line expected at --mixed-flux none).
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv102_mfnone_p2_O3_station_flt_0_7.5.dat"

PEAK_V_HYPO="N/A"
PEAK_T_HYPO="N/A"
SIGMA_N_PERTURB_HYPO="N/A"
SLIP_STRIKE_END="N/A"
T_END="N/A"
N_ROWS="N/A"
if [ -s "${HYPO_FILE}" ]; then
   # Column layout (tpv102_setup.hpp::TPV102StationWriter):
   #   1: time   2: slip1(dip)  3: slip2(strike)
   #   4: V1     5: V2          6: tau1   7: tau2   8: sigma_n  9: log10_theta
   read PEAK_V_HYPO PEAK_T_HYPO < <(awk '!/^#/ && NF>=5 {
      v=$5; if (v<0) v=-v;
      if (v>m) { m=v; tm=$1 }
   } END { printf "%g %g\n", m, tm }' "${HYPO_FILE}")
   SIGMA_N_PERTURB_HYPO=$(awk '!/^#/ && NF>=8 && $1+0==$1 {
      d=$8 - 120e6; if (d<0) d=-d;
      if (d>m) m=d
   } END { printf "%.4e", m+0 }' "${HYPO_FILE}")
   SLIP_STRIKE_END=$(awk '!/^#/ && NF>=3 && $1+0==$1 {a=$3} END{printf "%.4e", a+0}' "${HYPO_FILE}")
   T_END=$(awk '!/^#/ && NF>=1 && $1+0==$1 {tt=$1} END{printf "%.4f", tt+0}' "${HYPO_FILE}")
   N_ROWS=$(awk '!/^#/ && NF>=1 && $1+0==$1' "${HYPO_FILE}" | wc -l | tr -d ' ')
fi

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv102_mfnone_p2_O3_${SLURM_JOB_ID}.out"
DISPATCH_LINES=$(grep '\[dispatch\] rank=0' "${DISPATCH_LOG}" 2>/dev/null \
                 || echo "MISSING")
MIXED_FLUX_BANNER=$(grep -E '\[mixed-flux\] mode=' "${DISPATCH_LOG}" 2>/dev/null \
                    | head -1 || echo "MISSING")

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

# Sanity gate: peak |V_strike| should stay in physical range (≲ 30 m/s
# at the TPV102 hypocenter under the SCEC spec).  Anything ≥ 100 m/s
# means this baseline run also went unstable, which would invalidate
# the comparison and indicate a deeper issue than mixed-flux dispatch.
SANE_V_VERDICT="UNKNOWN"
if [ "${PEAK_V_HYPO}" != "N/A" ]; then
   if awk "BEGIN { exit !(${PEAK_V_HYPO} <= 100.0) }"; then
      SANE_V_VERDICT="PASS — peak |V| in physical range"
   else
      SANE_V_VERDICT="FAIL — peak |V| > 100 m/s; baseline ALSO blew up"
   fi
fi

{
   echo "TPV102 mixed-flux baseline (P=2 + ADER-O3 + --mixed-flux none) NORMAL"
   echo "Mesh:                    tpv102_200m.msh"
   echo "Polynomial order:        --order 2 + --ader-order 3"
   echo "Iterator:                one-shot"
   echo "Mixed flux:              none (upwind everywhere; clean control)"
   echo "Friction law:            aging (SCEC TPV101/102)"
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
   echo ""
   if [ "${PEAK_V_HYPO}" != "N/A" ] \
      && [ "${SIGMA_N_PERTURB_HYPO}" != "N/A" ] \
      && echo "${TFINAL_VERDICT}" | grep -qE "^PASS|^PARTIAL" \
      && echo "${SANE_V_VERDICT}" | grep -q "^PASS"; then
      echo "STATUS: PASS — clean control; diff vs --mixed-flux adjacent companion"
   else
      echo "STATUS: FAIL — see verdicts above"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
