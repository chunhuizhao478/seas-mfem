#!/bin/bash
#SBATCH -J tpv104_mfadj_p1_O2_dev
#SBATCH -o tpv104_mfadj_p1_O2_dev_%j.out
#SBATCH -e tpv104_mfadj_p1_O2_dev_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 mixed-flux comparison — MIXED FLUX ADJACENT, P=1 + ADER-O2.
# =======================================================================
# Pair 2 of 2 in the mixed-flux production sweep on tpv104_200m.msh:
#   [1] companion: --mixed-flux none      (upwind everywhere, baseline)
#   [2] this job:  --mixed-flux adjacent  (Mixed-Flux 2 per Zhang 2023)
#
# Both jobs share identical parameters (mesh, P, ADER-O, tfinal, CFL,
# BC, ParaView cadence, queue, allocation) so the only physics-relevant
# difference is the central-flux dispatch on fault-adjacent non-fault
# interior faces.
#
# Goal
# ----
# Quantify the effect of Zhang 2023 Mixed-Flux 2 on:
#   (a) hypocenter normal-stress variation through rupture (σ_n perturbation),
#   (b) hypocenter dip-slip drift at t = 12 s (the spurious-pollution channel),
# vs. the upwind baseline at production resolution (200 m mesh, P=1, t=12 s).
#
# Wall-time budget
# ----------------
# Dev queue 2 h envelope.  At 400 ranks on the 200 m mesh, P=1+O=2:
# extrapolating from tpv104_t12_dev_paraview.sbatch (~6 s/step → ~1200
# steps in 2 h → tfinal_actual ≈ 5.6 s) the run will be SLURM-killed
# before t=12 s.  R-1502 reduces the effective CFL by 0.9× when
# --mixed-flux adjacent is engaged, so this job will run a bit slower
# than the baseline (≈1080 steps in 2 h → tfinal_actual ≈ 5.0 s).
#
# Output cadence (must match the baseline for a clean diff)
# ---------------------------------------------------------
# --output-dt 0.05                  # station files every 0.05 s
# --paraview --paraview-dt 0.5      # fault-surface VTU/PVD every 0.5 s
# --paraview-bulk-dt 0.5            # bulk volume VTU/PVD every 0.5 s
#
# R-1204 / mixed-flux ADER-O note
# -------------------------------
# The R-1204 driver guard fires only when --mixed-flux != none AND
# --ader-order > 2.  This job uses ADER-O 2 (the in-spec, validated
# combination), so NO override env var is needed.
# The companion R-1409 driver banner will print:
#     [mixed-flux] mode=adjacent  |central_set|_global=N  ...
# which the post-run summary parses.

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

RESULT_DIR="tpv104/results_mixed_flux_adjacent_p1_O2_dev_job${SLURM_JOB_ID}"
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
echo "Polynomial:    --order 1 + --ader-order 2 (canonical Pelties O2)"
echo "Iterator:      one-shot (R7-001 dispatched; --fault-iterator one-shot)"
echo "Mixed flux:    adjacent (Mixed-Flux 2, Zhang 2023)"
echo "tfinal:        12.0 s (will be SLURM-killed at 2 h dev-queue limit)"
echo "ParaView:      surface 0.5 s, bulk 0.5 s"
echo "========================================"

ibrun ./seas_tpv104_driver \
      --mesh tpv104/mesh/tpv104_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 12.0 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law slip-srw \
      --fault-iterator one-shot \
      --mixed-flux adjacent \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv104_mfadj_p1_O2 \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.5 \
      --paraview-bulk-dt 0.5 \
      --verify-dispatch

# Post-run summary — same template as the baseline companion sbatch.
# Dev-queue run is intentionally SLURM-killed before tfinal=12 s, so
# the "≥50 rows / 9 fault stations" gate does not apply.  Quick-check
# pattern: hypocenter peak |V_strike|, peak |σ_n − 120 MPa|, and
# slip_dip(t=last_recorded), plus the [mixed-flux] banner so the
# central-flux dispatch is verified to have engaged.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv104_mfadj_p1_O2_station_x2_0_x3_7.5.dat"

PEAK_V_HYPO="N/A"
PEAK_T_HYPO="N/A"
SIGMA_N_PERTURB_HYPO="N/A"
SLIP_DIP_END="N/A"
T_END="N/A"
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
fi

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv104_mfadj_p1_O2_dev_${SLURM_JOB_ID}.out"
DISPATCH_LINES=$(grep '\[dispatch\] rank=0' "${DISPATCH_LOG}" 2>/dev/null \
                 || echo "MISSING")
MIXED_FLUX_BANNER=$(grep -E '\[mixed-flux\] mode=' "${DISPATCH_LOG}" 2>/dev/null \
                    | head -1 || echo "MISSING")
MIXED_FLUX_RANK_LINES=$(grep -c '\[mixed-flux\] rank=' "${DISPATCH_LOG}" 2>/dev/null \
                        || echo 0)

# Sanity gate: --mixed-flux adjacent should produce a non-empty
# central set on the production mesh.  If the global count is 0, the
# dispatch did not engage and the run is invalid for the comparison.
MIXED_FLUX_VERDICT="UNKNOWN"
if echo "${MIXED_FLUX_BANNER}" | grep -q "mode=adjacent"; then
   GSUM=$(echo "${MIXED_FLUX_BANNER}" \
          | grep -oE 'global=([0-9]+)' | head -1 \
          | grep -oE '[0-9]+')
   if [ -z "${GSUM}" ]; then
      GSUM=$(echo "${MIXED_FLUX_BANNER}" \
             | grep -oE '\|central_set\|_global=([0-9]+)' | head -1 \
             | grep -oE '[0-9]+$')
   fi
   if [ -n "${GSUM}" ] && [ "${GSUM}" -gt 0 ]; then
      MIXED_FLUX_VERDICT="ENGAGED — |central_set|_global=${GSUM}"
   else
      MIXED_FLUX_VERDICT="MISSING — banner present but global count is 0/missing"
   fi
else
   MIXED_FLUX_VERDICT="NOT-ENGAGED — no [mixed-flux] mode=adjacent line in log"
fi

{
   echo "TPV104 mixed-flux ADJACENT (P=1 + ADER-O2 + --mixed-flux adjacent) DEV-QUEUE"
   echo "Mesh:                    tpv104_200m.msh"
   echo "Polynomial order:        --order 1 + --ader-order 2"
   echo "Iterator:                one-shot"
   echo "Mixed flux:              adjacent (Mixed-Flux 2, Zhang 2023)"
   echo "tfinal target:           12.0 s (dev-queue 2 h SLURM-kill expected before tfinal)"
   echo "tfinal achieved (hypo):  ${T_END} s"
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
   if [ "${PEAK_V_HYPO}" != "N/A" ] \
      && [ "${SIGMA_N_PERTURB_HYPO}" != "N/A" ] \
      && echo "${MIXED_FLUX_VERDICT}" | grep -q "^ENGAGED"; then
      echo "STATUS: PASS — central dispatch engaged; diff with baseline companion job"
   else
      echo "STATUS: FAIL — see verdict above"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
