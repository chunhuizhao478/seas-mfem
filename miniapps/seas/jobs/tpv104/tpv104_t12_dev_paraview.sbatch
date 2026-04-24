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
#   --paraview-bulk-dt 0.5           # volume velocity + sigma_yy/xy/xz every 0.5 s
#                                    # (bulk volume is large; coarser cadence).
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
      --paraview-bulk-dt 0.5

# Post-run check: dev-queue run is intentionally killed before tfinal,
# so the production sbatch's "≥ 50 rows / 9 fault stations" gate does
# not apply.  Instead verify:
#   (a) at least one fault station has ≥ 5 rows (early rupture visible),
#   (b) no NaN/Inf in any station file,
#   (c) ParaView/ directory exists and contains at least 5 fault-surface
#       VTU frames (covers the t < 0.5 s nucleation interval at min).
RESULT_FILE="${RESULT_DIR}/RESULT.txt"

STATIONS_TOTAL=0
STATIONS_WITH_DATA=0
STATIONS_NANFREE=0
for f in "${RESULT_DIR}"/tpv104_station_x2_*.dat; do
   STATIONS_TOTAL=$((STATIONS_TOTAL + 1))
   ROWS=$(awk '!/^#/' "${f}" | wc -l)
   if [ "${ROWS}" -ge 5 ]; then
      STATIONS_WITH_DATA=$((STATIONS_WITH_DATA + 1))
   fi
   if ! grep -qiE "nan|inf" "${f}"; then
      STATIONS_NANFREE=$((STATIONS_NANFREE + 1))
   fi
done

PV_DIR="${RESULT_DIR}/ParaView"
PV_FRAMES=$(find "${PV_DIR}" -name "Cycle*.vtu" 2>/dev/null | wc -l || echo 0)
PV_BULK_DIR="${RESULT_DIR}/ParaView_bulk"
PV_BULK_FRAMES=$(find "${PV_BULK_DIR}" -name "Cycle*.vtu" 2>/dev/null | wc -l || echo 0)
FAULT_VTU=$(find "${RESULT_DIR}" -name "tpv104_fault_t*.vtu" 2>/dev/null | wc -l || echo 0)

{
   echo "TPV104 t=12.0 s DEV-QUEUE partial run (intentionally killed at 2h)"
   echo "fault stations total       : ${STATIONS_TOTAL}"
   echo "fault stations with ≥5 rows: ${STATIONS_WITH_DATA}"
   echo "fault stations NaN/Inf-free: ${STATIONS_NANFREE}"
   echo "ParaView/ volume frames    : ${PV_FRAMES}"
   echo "ParaView_bulk/ frames      : ${PV_BULK_FRAMES}"
   echo "fault-surface VTU frames   : ${FAULT_VTU}"
   if [ "${STATIONS_TOTAL}" -eq 9 ] \
      && [ "${STATIONS_WITH_DATA}" -ge 1 ] \
      && [ "${STATIONS_NANFREE}" -eq 9 ] \
      && [ "${FAULT_VTU}" -ge 5 ]; then
      echo "STATUS: PASS (partial output as expected for dev-queue cap)"
   else
      echo "STATUS: FAIL"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
