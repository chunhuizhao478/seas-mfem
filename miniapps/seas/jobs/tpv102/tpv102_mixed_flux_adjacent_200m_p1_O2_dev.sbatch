#!/bin/bash
#SBATCH -J tpv102_mfadj_p1_O2_dev
#SBATCH -o tpv102_mfadj_p1_O2_dev_%j.out
#SBATCH -e tpv102_mfadj_p1_O2_dev_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV102 mixed-flux test — MIXED FLUX ADJACENT, P=1 + ADER-O2 one-shot.
# =======================================================================
# First mixed-flux smoke run on the rebuilt TPV102 driver (TPV104 code
# flow + AgingLawPsi swap; commits adding tpv102_substep_iterator,
# tpv102_nucleation, tpv102_friction_solver, drivers/tpv102_driver.cpp).
#
# Baseline ADER-O2 + mixed-flux ADJACENT (Mixed-Flux 2 per Zhang 2023)
# at the SCEC TPV101/102 ageing-law spec parameters:
#   τ_ini=75 MPa, σ_n=120 MPa, V_ini=1e-12 m/s, Δτ₀=25 MPa, R=3 km, T=1 s.
#
# Goal
# ----
# Verify the new TPV102 driver runs cleanly under MPI at production scale
# on tpv102_200m.msh, with the mixed-flux dispatch engaged, and produces
# the 9 SCEC fault-station + 6 free-surface-station outputs.
#
# Wall-time budget
# ----------------
# Dev queue 2 h envelope.  Same 200 m mesh as the TPV104 mixed-flux pair
# but with the much smaller V_ini (1e-12 vs 1e-16) the early loading is
# slower and rupture nucleates ~t=1 s once the Δτ ramp completes.  At
# 400 ranks on the 200 m mesh, P=1+O=2 + R-1502 mixed-flux CFL factor
# 0.9× → ~5.5 s/step, ~1080 steps in 2 h → tfinal_actual ≈ 5.0 s,
# which covers nucleation through early rupture.
#
# Output cadence
# --------------
# --output-dt 0.05                  # station files every 0.05 s
# --paraview --paraview-dt 0.5      # fault-surface VTU/PVD every 0.5 s
# --paraview-bulk-dt 0.5            # bulk volume VTU/PVD every 0.5 s
#
# R-1204 / mixed-flux ADER-O note
# -------------------------------
# Mixed-flux dispatch is validated at ADER-O2 only.  This job uses O=2,
# so the R-1204 driver guard does NOT fire and no override env var is
# needed.  The companion banner will print:
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

if [ ! -s tpv102/mesh/tpv102_200m.msh ]; then
   echo "ERROR: tpv102/mesh/tpv102_200m.msh missing."
   exit 1
fi

RESULT_DIR="tpv102/results_mixed_flux_adjacent_p1_O2_dev_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv102_driver drivers/tpv102_driver.o \
      dynamic/tpv102_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o dynamic/wave_operator.o
make seas_tpv102_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

unset SEAS_FORCE_MIXED_FLUX_ADER_O_GT2
unset SEAS_FORCE_MIXED_FLUX_SUBSTEP_MPI

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          tpv102_200m.msh"
echo "Polynomial:    --order 1 + --ader-order 2 (Pelties O2)"
echo "Iterator:      one-shot (legacy wave.AdvanceADER dispatch)"
echo "Mixed flux:    adjacent (Mixed-Flux 2, Zhang et al. 2023)"
echo "Friction law:  aging (SCEC TPV101/102, dθ/dt = 1 − Vθ/L)"
echo "Queue:         development (8N x 400r, 2h)"
echo "tfinal:        12.0 s (dev-queue 2 h SLURM-kill expected before tfinal)"
echo "ParaView:      surface 0.5 s, bulk 0.5 s"
echo "========================================"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 12.0 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law aging \
      --fault-iterator one-shot \
      --mixed-flux adjacent \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv102_mfadj_p1_O2 \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.5 \
      --paraview-bulk-dt 0.5 \
      --verify-dispatch

# Post-run summary — same template as the TPV104 mixed-flux dev pair.
# Dev-queue 2 h envelope, so tfinal=12 s is intentionally SLURM-killable
# before completion.  Quick-check pattern: hypocenter peak |V_strike|,
# peak |σ_n − 120 MPa|, slip2(t_end), plus the [mixed-flux] banner.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv102_mfadj_p1_O2_station_flt_0_7.5.dat"

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

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv102_mfadj_p1_O2_dev_${SLURM_JOB_ID}.out"
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
   echo "TPV102 mixed-flux ADJACENT (P=1 + ADER-O2 + --mixed-flux adjacent) DEV-QUEUE"
   echo "Mesh:                    tpv102_200m.msh"
   echo "Polynomial order:        --order 1 + --ader-order 2"
   echo "Iterator:                one-shot"
   echo "Mixed flux:              adjacent (Mixed-Flux 2, Zhang et al. 2023)"
   echo "Friction law:            aging (SCEC TPV101/102)"
   echo "tfinal target:           12.0 s (dev-queue 2 h SLURM-kill expected before tfinal)"
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
      echo "STATUS: PASS — central dispatch engaged; rebuild verified for TPV102"
   else
      echo "STATUS: FAIL — see verdict above"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
