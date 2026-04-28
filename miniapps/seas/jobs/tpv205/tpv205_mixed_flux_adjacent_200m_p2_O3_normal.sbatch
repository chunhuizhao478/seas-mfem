#!/bin/bash
#SBATCH -J tpv205_mfadj_p2_O3
#SBATCH -o tpv205_mfadj_p2_O3_%j.out
#SBATCH -e tpv205_mfadj_p2_O3_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# TPV205 (SCEC TPV5) production run — MIXED FLUX ADJACENT,
# P=2 + ADER-O3 one-shot.
# =======================================================================
# First production-scale TPV205 mixed-flux run on the rebuilt LSW
# dispatch (Phase 1-4 of tpv205_lsw_native_fields_plan_2026-04-27):
#   - DOFData::lsw_mu_s/lsw_mu_d/lsw_d_c LSW-native fields
#   - FaultFaceFlux::EvaluateADER_LSW
#   - WaveOperator::SetFaultFrictionLaw(LSW) wired in driver
#   - Two LSW dispatch branches in wave_operator.inl (interior + R-1600
#     shared-fault fallback) so np > 1 produces correct rupture physics
#     across MPI rank boundaries.
#
# Baseline ADER-O3 + mixed-flux ADJACENT (Mixed-Flux 2 per Zhang 2023)
# at the SCEC TPV5 spec parameters:
#   τ_back=70 MPa, τ_nuc=81.6 MPa, τ_left=78 MPa, τ_right=62 MPa,
#   σ_n=120 MPa, μ_s=0.677 (rupture area) / 10000 (barrier),
#   μ_d=0.525, d_c=0.40 m, V_ini=0 m/s.  No time-varying nucleation —
#   the τ_strike(x, z) field is at its full patch values from t = 0
#   and rupture initiates spontaneously from the 81.6 MPa central
#   patch (μ_s · σ_n = 81.24 MPa < 81.6 MPa).
#
# Goal
# ----
# Produce the canonical TPV205 traces for the 16 SCEC fault stations
# (kStationsTPV205 in config/tpv205_params.hpp), the fault-surface
# ParaView output, and the bulk-volume ParaView output for comparison
# against the DRDG3D reference under tpv205/benchmark_data/DRDG3D_*.
#
# Wall-time budget
# ----------------
# Normal queue 48 h envelope.  At P=2 + ADER-O3 + mixed-flux on the
# 200 m mesh at 500 ranks, dt scales as P=2 → ~1/3 of the P=1 dt and
# ADER-O3 has higher per-step cost.  Assume effective ~12 s/step;
# tfinal = 12.0 s typically requires ~2400 macro steps; budget allows
# completion with margin for the bulk ParaView I/O.
#
# Output cadence
# --------------
# --output-dt 0.01                  # station files every 10 ms (matches
#                                     SCEC trace cadence)
# --paraview --paraview-dt 0.05     # fault-surface VTU/PVD every 50 ms
# --paraview-bulk-dt 0.10           # bulk volume VTU/PVD every 100 ms
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
# The driver's R-1503 guard would fire for `--fault-iterator substep`
# + `--mixed-flux != none` at any np.  This job uses
# `--fault-iterator one-shot` so R-1503 is silent and no override is
# needed.

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

RESULT_DIR="tpv205/results_mixed_flux_adjacent_p2_O3_normal_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv205_driver drivers/tpv205_driver.o \
      dynamic/tpv205_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o dynamic/wave_operator.o
make seas_tpv205_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

# R-1204: mixed-flux ADJACENT was only validated at ADER-O2.  This run
# uses ADER-O3, so the driver's R-1204 guard would abort without the
# override.  Acknowledge the UNTESTED combination per R-1407.
export SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1
unset SEAS_FORCE_MIXED_FLUX_SUBSTEP_MPI

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          tpv2053d_200m.msh"
echo "Polynomial:    --order 2 + --ader-order 3"
echo "Iterator:      one-shot (legacy AdvanceADERWithSubStep with O=1 quadrature)"
echo "Mixed flux:    adjacent (Mixed-Flux 2, Zhang et al. 2023)"
echo "Friction law:  lsw (linear slip-weakening; closed-form, SCEC TPV5)"
echo "Queue:         normal (10N x 500r, 48h)"
echo "tfinal:        12.0 s (full SCEC TPV5 §III window)"
echo "ParaView:      surface 0.05 s, bulk 0.10 s"
echo "Override:      SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1 (R-1204 / R-1407)"
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
      --mixed-flux adjacent \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv205_mfadj_p2_O3 \
      --output-dt 0.01 \
      --paraview \
      --paraview-dt 0.05 \
      --paraview-bulk-dt 0.10 \
      --verify-dispatch

# Post-run summary — same template as the TPV102/TPV104 mixed-flux
# dev-queue jobs.  Quick-check pattern: hypocenter-station peak |V_strike|,
# peak |σ_n − 120 MPa|, slip2(t_end), plus the [mixed-flux] banner and
# the [dispatch] rank=0 LSW banner (REVIEW R-004 / R-016).
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv205_mfadj_p2_O3_station_x2_0_x3_7.5.dat"

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

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv205_mfadj_p2_O3_${SLURM_JOB_ID}.out"
DISPATCH_LINES=$(grep '\[dispatch\] rank=0' "${DISPATCH_LOG}" 2>/dev/null \
                 || echo "MISSING")
MIXED_FLUX_BANNER=$(grep -E '\[mixed-flux\] mode=' "${DISPATCH_LOG}" 2>/dev/null \
                    | head -1 || echo "MISSING")
MIXED_FLUX_RANK_LINES=$(grep -c '\[mixed-flux\] rank=' "${DISPATCH_LOG}" 2>/dev/null \
                        || echo 0)

# Sanity gate: --mixed-flux adjacent should produce a non-empty central
# set on the production mesh.  If the global count is 0 the dispatch did
# not engage and the run is invalid.
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

# REVIEW R-016 sanity: the dispatch banner must report the LSW
# closed-form solver, not the rate-and-state Brent path.
LSW_DISPATCH_VERDICT="UNKNOWN"
if grep -q "friction_solver_actual=lsw-closed-form" "${DISPATCH_LOG}" 2>/dev/null; then
   LSW_DISPATCH_VERDICT="ENGAGED — friction_solver_actual=lsw-closed-form"
elif grep -q "friction_solver_actual=brent" "${DISPATCH_LOG}" 2>/dev/null; then
   LSW_DISPATCH_VERDICT="REGRESSION — banner reports brent (R-016 fix not in effect)"
else
   LSW_DISPATCH_VERDICT="UNKNOWN — no friction_solver_actual line in log"
fi

{
   echo "TPV205 mixed-flux ADJACENT (P=2 + ADER-O3 + --mixed-flux adjacent) NORMAL-QUEUE"
   echo "Mesh:                    tpv2053d_200m.msh"
   echo "Polynomial order:        --order 2 + --ader-order 3"
   echo "Iterator:                one-shot"
   echo "Mixed flux:              adjacent (Mixed-Flux 2, Zhang et al. 2023)"
   echo "Friction law:            lsw (linear slip-weakening; SCEC TPV5)"
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
   echo "LSW dispatch (R-016)     : ${LSW_DISPATCH_VERDICT}"
   echo ""
   if [ "${PEAK_V_HYPO}" != "N/A" ] \
      && [ "${SIGMA_N_PERTURB_HYPO}" != "N/A" ] \
      && echo "${MIXED_FLUX_VERDICT}" | grep -q "^ENGAGED" \
      && echo "${LSW_DISPATCH_VERDICT}" | grep -q "^ENGAGED"; then
      echo "STATUS: PASS — central dispatch engaged + LSW closed-form active"
   else
      echo "STATUS: FAIL — see verdict above"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
