#!/bin/bash
#SBATCH -J tpv102_mfadj_p2_O3
#SBATCH -o tpv102_mfadj_p2_O3_%j.out
#SBATCH -e tpv102_mfadj_p2_O3_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# TPV102 production run — MIXED FLUX ADJACENT, P=2 + ADER-O3 one-shot.
# =======================================================================
# Production-scale follow-up to the TPV102 mixed-flux dev-queue sbatch
# (tpv102_mixed_flux_adjacent_200m_p1_O2_dev.sbatch).  Same rebuilt
# TPV102 driver (TPV104 code flow + AgingLawPsi swap; commits adding
# tpv102_substep_iterator, tpv102_nucleation, tpv102_friction_solver,
# drivers/tpv102_driver.cpp), bumped from P=1+O=2 to P=2+O=3 with the
# normal-queue 48 h envelope so tfinal=12 s completes.
#
# Baseline ADER-O3 + mixed-flux ADJACENT (Mixed-Flux 2 per Zhang 2023)
# at the SCEC TPV101/102 ageing-law spec parameters:
#   τ_ini=75 MPa, σ_n=120 MPa, V_ini=1e-12 m/s, Δτ₀=25 MPa, R=3 km, T=1 s.
#
# Goal
# ----
# Produce the canonical TPV102 traces for the 9 SCEC fault-station +
# 6 free-surface-station outputs at production resolution
# (tpv102_200m.msh, P=2 + ADER-O3 + mixed-flux adjacent), comparable to
# the TPV104 mixed-flux normal-queue pair (commit 4189c9c).
#
# Wall-time budget
# ----------------
# Normal queue 48 h envelope.  At P=2 + ADER-O3 + mixed-flux on the
# 200 m mesh at 500 ranks, dt scales as P=2 → ~1/3 of the P=1 dt and
# ADER-O3 has higher per-step cost.  Budget allows tfinal=12 s
# completion with margin for the bulk ParaView I/O.
#
# Output cadence
# --------------
# --output-dt 0.05                  # station files every 0.05 s
# --paraview --paraview-dt 0.5      # fault-surface VTU/PVD every 0.5 s
# --paraview-bulk-dt 0.5            # bulk volume VTU/PVD every 0.5 s
#
# R-1204 / mixed-flux ADER-O note
# -------------------------------
# Mixed-flux dispatch was validated at ADER-O2 only.  This job uses
# O=3, so the R-1204 driver guard fires.  The override env var
# SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1 below acknowledges this is an
# UNTESTED combination (cite R-1204/R-1407 if a regression appears).
#
# R-1503 / mixed-flux + substep iterator note
# -------------------------------------------
# Driver guard would fire for `--fault-iterator substep` + `--mixed-flux
# != none`.  This job uses `--fault-iterator one-shot`, so R-1503 is
# silent and no override is needed.

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

RESULT_DIR="tpv102/results_mixed_flux_adjacent_p2_O3_normal_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv102_driver drivers/tpv102_driver.o \
      dynamic/tpv102_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o dynamic/wave_operator.o
make seas_tpv102_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

# R-1204: mixed-flux ADJACENT was only validated at ADER-O2.  This run
# uses ADER-O3, so the driver's R-1204 guard would abort without the
# override.  Acknowledge the UNTESTED combination per R-1407.
export SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1
unset SEAS_FORCE_MIXED_FLUX_SUBSTEP_MPI

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          tpv102_200m.msh"
echo "Polynomial:    --order 2 + --ader-order 3"
echo "Iterator:      one-shot (legacy wave.AdvanceADER dispatch)"
echo "Mixed flux:    adjacent (Mixed-Flux 2, Zhang et al. 2023)"
echo "Friction law:  aging (SCEC TPV101/102, dθ/dt = 1 − Vθ/L)"
echo "Queue:         normal (10N x 500r, 48h)"
echo "tfinal:        12.0 s (full SCEC TPV101/102 window)"
echo "ParaView:      surface 0.5 s, bulk 0.5 s"
echo "Override:      SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1 (R-1204 / R-1407)"
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
      --mixed-flux adjacent \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv102_mfadj_p2_O3 \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.5 \
      --paraview-bulk-dt 0.5 \
      --verify-dispatch

# Post-run summary — same template as the TPV102 mixed-flux dev-queue
# sbatch (tpv102_mixed_flux_adjacent_200m_p1_O2_dev.sbatch).  Quick-check
# pattern: hypocenter peak |V_strike|, peak |σ_n − 120 MPa|, slip2(t_end),
# plus the [mixed-flux] banner.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv102_mfadj_p2_O3_station_flt_0_7.5.dat"

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

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv102_mfadj_p2_O3_${SLURM_JOB_ID}.out"
DISPATCH_LINES=$(grep '\[dispatch\] rank=0' "${DISPATCH_LOG}" 2>/dev/null \
                 || echo "MISSING")
MIXED_FLUX_BANNER=$(grep -E '\[mixed-flux\] mode=' "${DISPATCH_LOG}" 2>/dev/null \
                    | head -1 || echo "MISSING")
MIXED_FLUX_RANK_LINES=$(grep -c '\[mixed-flux\] rank=' "${DISPATCH_LOG}" 2>/dev/null \
                        || echo 0)

# Sanity gate: --mixed-flux adjacent should produce a non-empty central
# set on the production mesh.  If the global count is 0 the dispatch did
# not engage and the run is invalid for the comparison.
MIXED_FLUX_VERDICT="UNKNOWN"
if echo "${MIXED_FLUX_BANNER}" | grep -q "mode=adjacent"; then
   GSUM=$(echo "${MIXED_FLUX_BANNER}" \
          | grep -oE '\|central_set\|_global=([0-9]+)' | head -1 \
          | grep -oE '[0-9]+$')
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

{
   echo "TPV102 mixed-flux ADJACENT (P=2 + ADER-O3 + --mixed-flux adjacent) NORMAL-QUEUE"
   echo "Mesh:                    tpv102_200m.msh"
   echo "Polynomial order:        --order 2 + --ader-order 3"
   echo "Iterator:                one-shot"
   echo "Mixed flux:              adjacent (Mixed-Flux 2, Zhang et al. 2023)"
   echo "Friction law:            aging (SCEC TPV101/102)"
   echo "tfinal target:           12.0 s"
   echo "tfinal achieved (hypo):  ${T_END} s"
   echo "hypo station rows:       ${N_ROWS}"
   echo ""
   echo "peak |V_strike| at hypo  : ${PEAK_V_HYPO} m/s at t = ${PEAK_T_HYPO} s"
   echo "peak |σ_n − 120 MPa| @ h : ${SIGMA_N_PERTURB_HYPO} Pa"
   echo "slip_strike at hypo (end): ${SLIP_STRIKE_END} m"
   echo ""
   echo "Mixed-flux banner        : ${MIXED_FLUX_BANNER}"
   echo "Mixed-flux engagement    : ${MIXED_FLUX_VERDICT}"
   echo "Mixed-flux rank lines    : ${MIXED_FLUX_RANK_LINES}"
   echo "[dispatch] rank=0 lines  : ${DISPATCH_LINES}"
   echo ""
   if [ "${PEAK_V_HYPO}" != "N/A" ] \
      && [ "${SIGMA_N_PERTURB_HYPO}" != "N/A" ] \
      && echo "${MIXED_FLUX_VERDICT}" | grep -q "^ENGAGED"; then
      echo "STATUS: PASS — central dispatch engaged at production scale"
   else
      echo "STATUS: FAIL — see verdict above"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
