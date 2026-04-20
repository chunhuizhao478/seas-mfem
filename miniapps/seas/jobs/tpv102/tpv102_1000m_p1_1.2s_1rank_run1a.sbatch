#!/bin/bash
#SBATCH -J tpv102_run1a_1r
#SBATCH -o tpv102_run1a_1r_%j.out
#SBATCH -e tpv102_run1a_1r_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV102 v8.0.0 debug plan — Phase 1 Run 1a.
#
# Purpose: test the **interior-fault branch in isolation** at np=1.
# With a single rank there are NO shared fault faces (no partition
# seam crosses the fault); only the interior-face branch of
# `ComputeFaceFluxRHS` runs.  If the v8.0.0 pinned-V_max symptom
# (rupture stalled at hypocenter, all off-hypo stations at V_ini)
# reproduces here, the bug is in the interior branch, Evaluate,
# or accumulation — Case A in the plan's decision table.  If
# propagation is healthy here but fails at N-rank, the bug lives
# in the shared-fault branch — Case B.
#
# Mesh: tpv102_1000m.msh.  A serial run on tpv102_200m.msh would
# take ~100+ hours wall-clock and is infeasible even on a Frontera
# compute node.  Using the 1000m mesh is a deliberate numerical
# trade — the fault geometry, material, boundary layout, and
# SCEC-spec rate-and-state parameters are IDENTICAL across mesh
# sizes; only the cohesive-zone resolution changes.  At h=1000 m,
# Λ_dyn/h ≈ 0.16 — under-resolved but sufficient to detect a
# global symptom (pinned V_max) or its absence.
#
# WALL-CLOCK: 1000 m / p=1 / serial — dt_cfl ≈ 1.47e-3 s (from
# Phase 0 local probe), so 817 RK4 steps for tfinal=1.2 s.  At
# ~0.4 s/step on a single Cascade Lake core, ~5–6 min.  Padding
# generously for I/O + init to 2 hr so the dev queue accepts.
# SU COST: ~2 SU.
#
# Usage:
#   sbatch jobs/tpv102/tpv102_1000m_p1_1.2s_1rank_run1a.sbatch

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
rm -f drivers/tpv102_driver.o dynamic/wave_operator.o dynamic/fault_face_flux.o
make seas_tpv102_driver -j"${JOBS}" SEAS_EXTRA_CPPFLAGS="" \
   || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_1000m_p1_1.2s_run1a_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Phase 1 Run 1a: np=1 interior-fault branch only, tfinal=1.2 s"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_1000m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 1.2 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --no-domain-pv \
      --debug-qnorm

# Phase 1 success criterion: score how many of the 9 SCEC fault
# stations have V2 > 1e-6 m/s at the final time.  Station file
# names are `flt_0_3.dat`, `flt_0_7.5.dat`, `flt_0_12.dat`,
# `flt_9_7.5.dat`, `flt_12_3.dat`, `flt_12_12.dat`,
# `flt_n9_7.5.dat`, `flt_n12_3.dat`, `flt_n12_12.dat`.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv102_run1a_1r_${SLURM_JOB_ID}.out"
if [ ! -f "${JOB_LOG}" ]; then
   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}'). Hint: SLURM_SUBMIT_DIR='${SLURM_SUBMIT_DIR:-<unset>}'." > "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

SLIP_COUNT=0
SLIPPED_STATIONS=""
SILENT_STATIONS=""
for sta in flt_0_3 flt_0_7.5 flt_0_12 flt_9_7.5 flt_12_3 flt_12_12 \
           flt_n9_7.5 flt_n12_3 flt_n12_12; do
   STA_FILE="${RESULT_DIR}/${sta}.dat"
   if [ ! -f "${STA_FILE}" ]; then
      SILENT_STATIONS="${SILENT_STATIONS} ${sta}(missing)"
      continue
   fi
   # Take the last data line and read V2.  Station file columns
   # (tpv102_setup.hpp:346): time(1) slip1(2) slip2(3) V1(4) V2(5)
   # tau1(6) tau2(7) sigma_n(8) log10_theta(9).  Consider "slipped"
   # iff |V2| > 1e-6 (= 1e-12 after squaring for awk comparison).
   V2=$(awk '!/^#/ && NF>=5 {v=$5} END {print (v=="" ? "0" : v)}' "${STA_FILE}")
   if awk "BEGIN {exit !( ($V2)*($V2) > 1e-12 )}"; then
      SLIP_COUNT=$((SLIP_COUNT + 1))
      SLIPPED_STATIONS="${SLIPPED_STATIONS} ${sta}"
   else
      SILENT_STATIONS="${SILENT_STATIONS} ${sta}"
   fi
done

if grep -q "MFEM abort\|MFEM_ABORT" "${JOB_LOG}"; then
   echo "FAIL: MFEM abort during Run 1a" > "${RESULT_FILE}"
else
   echo "Run 1a propagation score: ${SLIP_COUNT} / 9 stations slipped" > "${RESULT_FILE}"
   echo "  Slipped:${SLIPPED_STATIONS}" >> "${RESULT_FILE}"
   echo "  Silent :${SILENT_STATIONS}" >> "${RESULT_FILE}"
   if [ "${SLIP_COUNT}" -eq 0 ]; then
      echo "  Decision: 0/9 — Case A (interior / Evaluate / accumulation)" >> "${RESULT_FILE}"
   elif [ "${SLIP_COUNT}" -lt 4 ]; then
      echo "  Decision: 1-3/9 partial — Case B-partial" >> "${RESULT_FILE}"
   else
      echo "  Decision: 4-9/9 — healthy at this scale" >> "${RESULT_FILE}"
   fi
fi
cat "${RESULT_FILE}"
