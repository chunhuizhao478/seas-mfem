#!/bin/bash
#SBATCH -J tpv104_mfnone_p2_O3
#SBATCH -o tpv104_mfnone_p2_O3_%j.out
#SBATCH -e tpv104_mfnone_p2_O3_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# TPV104 mixed-flux comparison — BASELINE (no mixed flux), P=2 + ADER-O3.
# =======================================================================
# Pair 1 of 2 in the P=2 production sweep on tpv104_200m.msh:
#   [1] this job:  --mixed-flux none      (upwind everywhere, baseline)
#   [2] companion: --mixed-flux adjacent  (Mixed-Flux 2 per Zhang 2023)
#
# Both jobs share identical parameters (mesh, P, ADER-O, tfinal, CFL, BC,
# ParaView cadence, queue, allocation) so the only physics-relevant
# difference is the central-flux dispatch on fault-adjacent non-fault
# interior faces.  Pairs with the round-1 dev-queue P=1+O=2 doublet to
# span the Pelties precision ladder under both upwind and Mixed-Flux 2
# dispatch.
#
# Goal
# ----
# At the recommended Pelties O3 (P=2 + ADER-O3) on the production
# tpv104_200m.msh, quantify Mixed-Flux 2's effect on:
#   (a) hypocenter normal-stress variation through rupture (σ_n perturbation),
#   (b) hypocenter dip-slip drift at t = 12 s (the spurious-pollution channel),
# vs. the upwind baseline.  P=2 reduces the precision-floor by ~8× per the
# symmirror reference (ADER_ITERATOR_SYM1000_RESULTS.md), so the dip-drift
# magnitude should land in the ~1e-2 m regime here vs ~1e-1 m at P=1.
#
# Wall-time budget
# ----------------
# Normal queue 48 h envelope.  At 500 ranks on the 200 m mesh, P=2+O=3:
# ~3× per-step cost vs P=1+O=2 (denser shape tables, larger ADER predictor)
# and ~1.7× more steps (CFL ~ h_min / cp_eff scales with element diameter,
# only ~1.7× tighter at P=2 since dt is set by element h not by P directly
# in our --cfl 0.5 path).  Extrapolating from the t12_dev P=1+O=2 baseline
# (~6 s/step on 400 ranks → ~5 s/step on 500 ranks), P=2 should land near
# ~15 s/step.  ~2550 macro steps to t=12s → ~10.6 h wall.  48 h envelope
# is intentionally generous (3.5x headroom) so the run completes through
# t=12 s including any rupture-front slowdown in the post-nucleation phase.
#
# Output cadence
# --------------
# --output-dt 0.05                  # station files every 0.05 s
# --paraview --paraview-dt 0.5      # fault-surface VTU/PVD every 0.5 s
# --paraview-bulk-dt 0.5            # bulk volume VTU/PVD every 0.5 s
# (12 s tfinal × 2 frames/s = 24 surface + 24 bulk frames; bulk volume on
# 200 m mesh at P=2 is ~3-5 GB/frame, total ~75-120 GB on /scratch2.)
#
# R-1204 / mixed-flux ADER-O note
# -------------------------------
# This is the BASELINE (--mixed-flux none).  No R-1204 guard concern —
# the guard fires only when --mixed-flux != none AND --ader-order > 2.
# The companion --mixed-flux adjacent job at O=3 will need to set the
# SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1 override env var (see that sbatch).

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

RESULT_DIR="tpv104/results_mixed_flux_none_p2_O3_normal_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o dynamic/wave_operator.o
make seas_tpv104_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

unset SEAS_TPV104_FREEZE_SIGMA_N
unset SEAS_FORCE_MIXED_FLUX_ADER_O_GT2
unset SEAS_FORCE_MIXED_FLUX_SUBSTEP_MPI

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          tpv104_200m.msh"
echo "Polynomial:    --order 2 + --ader-order 3 (canonical Pelties O3)"
echo "Iterator:      one-shot (R7-001 dispatched; --fault-iterator one-shot)"
echo "Mixed flux:    none (upwind everywhere, baseline)"
echo "Queue:         normal (10N x 500r, 48h)"
echo "tfinal:        12.0 s (full production tfinal target)"
echo "ParaView:      surface 0.5 s, bulk 0.5 s"
echo "========================================"

ibrun ./seas_tpv104_driver \
      --mesh tpv104/mesh/tpv104_200m.msh \
      --mesh-scale 1 \
      --order 2 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 12.0 \
      --ader-order 3 \
      --friction-solver newton-stable \
      --fric-law slip-srw \
      --fault-iterator one-shot \
      --mixed-flux none \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv104_mfnone_p2_O3 \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.5 \
      --paraview-bulk-dt 0.5 \
      --verify-dispatch

# Post-run summary — same template as the dev-queue P=1 pair.  Normal
# queue with 48 h envelope: tfinal=12 s should complete cleanly so the
# t_end ≈ 12 s gate matters here (unlike the dev-queue 2 h truncation).
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv104_mfnone_p2_O3_station_x2_0_x3_7.5.dat"

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

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv104_mfnone_p2_O3_${SLURM_JOB_ID}.out"
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

{
   echo "TPV104 mixed-flux baseline (P=2 + ADER-O3 + --mixed-flux none) NORMAL"
   echo "Mesh:                    tpv104_200m.msh"
   echo "Polynomial order:        --order 2 + --ader-order 3"
   echo "Iterator:                one-shot"
   echo "Mixed flux:              none (upwind everywhere)"
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
   echo "[dispatch] rank=0 lines  : ${DISPATCH_LINES}"
   echo ""
   if [ "${PEAK_V_HYPO}" != "N/A" ] \
      && [ "${SIGMA_N_PERTURB_HYPO}" != "N/A" ] \
      && echo "${TFINAL_VERDICT}" | grep -qE "^PASS|^PARTIAL"; then
      echo "STATUS: PASS (companion --mixed-flux adjacent job to be diffed)"
   else
      echo "STATUS: FAIL — see verdict above"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
