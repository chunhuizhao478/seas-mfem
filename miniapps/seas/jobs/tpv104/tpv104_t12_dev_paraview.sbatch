#!/bin/bash
#SBATCH -J tpv104_t12_dev
#SBATCH -o tpv104_t12_dev_%j.out
#SBATCH -e tpv104_t12_dev_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 dev-queue partial-run with ParaView output for visual sanity.
#
# Purpose: dispatch the full t = 12 s production target to the dev queue
# under a 2-hour wall.  The job WILL be killed by SLURM before t = 12 s
# finishes — this is intentional.  The point is to verify on a short
# turnaround that:
#   (a) banner + dispatch consistency holds at 400 ranks,
#   (b) the ADER one-shot loop runs without NaN past nucleation
#       (t > 1 s) and into early rupture propagation (1 s < t < ~2 s),
#   (c) ParaView output lands and renders in the same pattern as the
#       TPV102 dev-queue runs (jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch),
#   (d) station + fault-surface VTU files are produced at the expected
#       cadence so the Phase-3 probe-diff harness can be exercised on
#       partial output before committing to the 5-h normal-queue run.
#
# WALL-CLOCK ESTIMATE:
#   At 400 ranks on the 200 m mesh, ~4.7 ms per macro-step → ~2550
#   steps to reach t = 12 s.  Throughput (extrapolated from TPV102
#   400-rank dev runs at ~0.1 s/step): ~6 s/step on 200 m.  In 2 h
#   that's ~1200 steps → tfinal_actual ≈ 5.6 s.  Plenty for nucleation
#   + breakout + early propagation visual check.
# SU COST: ~800 SU (8 nodes × 2 h × ~50 SU/node-h).
#
# Output cadence:
#   --output-dt 0.1                  # 9 fault stations + 6 surface stations every 0.1 s
#   --paraview --paraview-dt 0.1     # fault-surface VTU/PVD every 0.1 s
#   --paraview-bulk-dt 1.0           # volume velocity + sigma_yy/xy/xz every 1.0 s
#                                    # (bulk volume is GB per frame; 1 s cadence
#                                    # gives ~12 frames over tfinal=12 s, total
#                                    # ~25-60 GB on /scratch2).
#
# Pre-requisites:
#   1. tpv104_build_banner_sanity.sbatch — PASS (banner + dispatch).
#   2. tpv104 .msh files in miniapps/seas/tpv104/mesh/ (built locally
#      and scp'd, OR via tpv104_mesh_build.sbatch).
#   3. tpv104_init_stations.sbatch — PASS (P2_D ψ_ini gate).
#
# R7-001 disclosure (banner-vs-dispatch): runtime is Brent + one-shot
# + macro-step ψ.  --friction-solver / --fault-iterator flags are NOT
# dispatched (banner-only).  See jobs/tpv104/README.md.

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

# Integration-readiness gate — fail fast if the prior init job did not PASS.
PREV_INIT_RESULT=$(ls -dt tpv104/results_init_stations_job* 2>/dev/null | head -1)
if [ -z "${PREV_INIT_RESULT}" ] \
   || ! grep -q "^STATUS: PASS" "${PREV_INIT_RESULT}/RESULT.txt" 2>/dev/null; then
   echo "ERROR: no passing tpv104_init_stations.sbatch result found."
   echo "       Run that job first and verify STATUS: PASS before proceeding."
   exit 1
fi
echo "Found init PASS at: ${PREV_INIT_RESULT}"

RESULT_DIR="tpv104/results_t12_dev_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o
make seas_tpv104_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

if [ ! -s tpv104/mesh/tpv104_200m.msh ]; then
   echo "ERROR: tpv104/mesh/tpv104_200m.msh missing."
   exit 1
fi

echo "========================================"
echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Target: tfinal = 12.0 s on tpv104_200m.msh at 400 ranks (DEV QUEUE — will be killed at 2 h)"
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
      --fault-iterator substep \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv104 \
      --output-dt 0.1 \
      --paraview \
      --paraview-dt 0.1 \
      --paraview-bulk-dt 1.0

# Post-run check: dev-queue run is intentionally killed before tfinal,
# so the production sbatch's "≥50 rows / 9 fault stations" gate does
# not apply.  Quick-check decision pattern mirrors
# jobs/tpv102/tpv102_200m_p1_2.0s_400rank_vmax_halt.sbatch:
#   (a) MFEM abort short-circuit (FAIL fast on solver/R-101 trip),
#   (b) per-station slip-progression score (how many of the 9 fault
#       stations have non-trivial strike-slip rate, the TPV104 mode-II
#       diagnostic),
#   (c) V_max timeline (first / peak / final lines from the job log) —
#       confirms nucleation drove V_max above V_ini=1e-16 and rupture
#       propagated rather than pinning,
#   (d) ParaView frame count + station NaN-free check.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
JOB_LOG="${SBATCH_LOG_DIR}/tpv104_t12_dev_${SLURM_JOB_ID}.out"
JOB_ERR="${SBATCH_LOG_DIR}/tpv104_t12_dev_${SLURM_JOB_ID}.err"

# (a) MFEM abort short-circuit.
if [ -f "${JOB_ERR}" ] && grep -q "MFEM abort\|MFEM_ABORT\|R-101 shared-fault DOFData" "${JOB_ERR}"; then
   ABORT_LINE=$(grep -E "MFEM abort|MFEM_ABORT|R-101 shared-fault DOFData" "${JOB_ERR}" | head -1)
   echo "FAIL: MFEM abort during dev-queue run" > "${RESULT_FILE}"
   echo "      ${ABORT_LINE}" >> "${RESULT_FILE}"
   cat "${RESULT_FILE}"
   exit 0
fi

# (b) Per-station slip-progression score.  TPV104 station column layout
# (per TPV104StationWriter, plan §4.4):
#   col 1 = t       col 2 = h-slip      col 3 = h-slip-rate (strike rate)
#   col 4 = h-stress col 5 = v-slip      col 6 = v-slip-rate (dip rate)
#   col 7 = v-stress col 8 = n-stress    col 9 = psi
# TPV104 is pure strike-slip → score on col 3 (h-rate).  A station has
# "slipped" if any row has |h-rate|^2 > 1e-12 → |h-rate| > ~1e-6 m/s.
SLIP_COUNT=0
SLIPPED_STATIONS=""
SILENT_STATIONS=""
for sta in tpv104_station_x2_0_x3_3 tpv104_station_x2_0_x3_7.5 tpv104_station_x2_0_x3_12 \
           tpv104_station_x2_9_x3_7.5 tpv104_station_x2_12_x3_3 tpv104_station_x2_12_x3_12 \
           tpv104_station_x2_-9_x3_7.5 tpv104_station_x2_-12_x3_3 tpv104_station_x2_-12_x3_12; do
   STA_FILE="${RESULT_DIR}/${sta}.dat"
   if [ ! -f "${STA_FILE}" ]; then
      SILENT_STATIONS="${SILENT_STATIONS} ${sta}(missing)"
      continue
   fi
   # Maximum |h-rate| across all rows (col 3).  Returns 0 if no data.
   V_PEAK=$(awk '!/^#/ && NF>=3 { v=$3; if (v<0) v=-v; if (v>m) m=v } END { print (m=="" ? "0" : m) }' "${STA_FILE}")
   if awk "BEGIN {exit !( ($V_PEAK)*($V_PEAK) > 1e-12 )}"; then
      SLIP_COUNT=$((SLIP_COUNT + 1))
      SLIPPED_STATIONS="${SLIPPED_STATIONS} ${sta}"
   else
      SILENT_STATIONS="${SILENT_STATIONS} ${sta}"
   fi
done

# (c) V_max timeline from job log.
V_MAX_INIT=""
V_MAX_PEAK=""
V_MAX_FINAL=""
if [ -f "${JOB_LOG}" ]; then
   V_MAX_INIT=$(grep "V_max" "${JOB_LOG}" | head -1)
   V_MAX_FINAL=$(grep "V_max" "${JOB_LOG}" | tail -1)
   V_MAX_PEAK=$(grep "V_max" "${JOB_LOG}" \
                | awk '{for(i=1;i<=NF;i++) if ($i ~ /^[0-9.e+-]+$/) v=$i; print v}' \
                | sort -g | tail -1)
fi
LAST_STEP=$([ -f "${JOB_LOG}" ] && grep "^Step " "${JOB_LOG}" | tail -1 || echo "")

# (d) Frame counts + NaN sweep.
STATIONS_TOTAL=0
STATIONS_NANFREE=0
for f in "${RESULT_DIR}"/tpv104_station_x2_*.dat; do
   STATIONS_TOTAL=$((STATIONS_TOTAL + 1))
   if ! grep -qiE "nan|inf" "${f}"; then
      STATIONS_NANFREE=$((STATIONS_NANFREE + 1))
   fi
done
PV_FRAMES=$(find "${RESULT_DIR}/ParaView" -name "Cycle*.vtu" 2>/dev/null | wc -l || echo 0)
PV_BULK_FRAMES=$(find "${RESULT_DIR}/ParaView_bulk" -name "Cycle*.vtu" 2>/dev/null | wc -l || echo 0)
FAULT_VTU=$(find "${RESULT_DIR}" -name "tpv104_fault_t*.vtu" 2>/dev/null | wc -l || echo 0)

{
   echo "TPV104 t=12.0 s DEV-QUEUE partial run (intentionally killed at 2h)"
   echo "  Slip score : ${SLIP_COUNT} / 9 fault stations have |h-rate| > 1e-6 m/s"
   echo "  Slipped    :${SLIPPED_STATIONS}"
   echo "  Silent     :${SILENT_STATIONS}"
   echo "  fault stations total       : ${STATIONS_TOTAL}"
   echo "  fault stations NaN/Inf-free: ${STATIONS_NANFREE}"
   echo "  ParaView/ volume frames    : ${PV_FRAMES}"
   echo "  ParaView_bulk/ frames      : ${PV_BULK_FRAMES}"
   echo "  fault-surface VTU frames   : ${FAULT_VTU}"
   [ -n "${LAST_STEP}"    ] && echo "  Last step line : ${LAST_STEP}"
   [ -n "${V_MAX_INIT}"   ] && echo "  First V_max line: ${V_MAX_INIT}"
   [ -n "${V_MAX_PEAK}"   ] && echo "  Peak  V_max value: ${V_MAX_PEAK}"
   [ -n "${V_MAX_FINAL}"  ] && echo "  Final V_max line: ${V_MAX_FINAL}"
   # Decision:
   #   - 0/9 slipped AND V_max stays at V_ini=1e-16: ADER one-shot did
   #     NOT propagate rupture — escalate (probably nucleation off, or
   #     wave operator silent on shared faces, or P2_D wired wrong).
   #   - >=6/9 slipped AND V_max peak > 0.1 m/s: rupture broke through
   #     hypocenter, propagated past partition seams — proceed to the
   #     5-h normal-queue run.
   #   - Intermediate: partial — investigate per-station before the 5h.
   if [ "${STATIONS_TOTAL}" -ne 9 ] || [ "${STATIONS_NANFREE}" -ne 9 ] \
      || [ "${FAULT_VTU}" -lt 5 ]; then
      echo "STATUS: FAIL — output incomplete or NaN/Inf in stations"
   elif [ "${SLIP_COUNT}" -eq 0 ]; then
      echo "STATUS: FAIL — 0/9 slipped; rupture did not break through hypocenter"
      echo "  Decision: do NOT submit the 5-h normal-queue run; diagnose first"
   elif [ "${SLIP_COUNT}" -ge 6 ]; then
      echo "STATUS: PASS — ${SLIP_COUNT}/9 slipped; rupture propagated past seams"
      echo "  Decision: cleared to submit jobs/tpv104/tpv104_t3_mesh_coupled.sbatch"
   else
      echo "STATUS: PARTIAL — ${SLIP_COUNT}/9 slipped; investigate before 5-h run"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
