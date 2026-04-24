#!/bin/bash
#SBATCH -J tpv104_init_stn
#SBATCH -o tpv104_init_stn_%j.out
#SBATCH -e tpv104_init_stn_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 50
#SBATCH -t 00:20:00
#SBATCH -A EAR20006

# TPV104 Phase-2 init-time station sanity on Frontera — gate P2_D.
#
# Purpose: verify at tfinal=0 (init-only) that every one of the nine
# SCEC TPV104 fault stations reports ψ_ini = 5.6359184e-01 to 1e-8.
# This is plan §4.6 Phase-2 P2_D gate for the TPV104 driver.  Init-
# only means tfinal=0 and the ADER time loop is skipped entirely, so
# the R7-001 solver-dispatch caveat below is IRRELEVANT for this run:
# the check reads from InitializeFaultDOFs_TPV104's ψ_ini output
# before any solver is invoked.
#
# Mesh: tpv104_200m.msh (production canonical mesh per plan §4.1).
# Time integrator: ADER-O2 (one-shot via wave.AdvanceADER; sub-step
#                  iterator NOT wired — R7-001 open).
# Fault iterator : one-shot.  --fault-iterator=substep below is
#                  parsed for banner/smoke-test parity but has NO
#                  effect on the dispatched loop.  For tfinal=0 this
#                  makes no difference: the time loop doesn't run.
# Friction solver: Brent, hard-coded via EvaluateADERTotal.
#                  --friction-solver=newton-stable below is parsed
#                  but has NO effect on the dispatched solver.
#                  Again, irrelevant at tfinal=0.
# Friction law   : slip-SRW (ψ-space, free-function
#                  UpdateStateAnalyticSlipLawSRW, macro-step cadence).
#
# WALL-CLOCK: ~5 min for init + 1 time step that dumps station t=0 rows.
# SU COST: ~0.5 SU.
#
# Pre-requisites (run first in this order):
#   1. `sbatch tpv104_build_banner_sanity.sbatch` — banner + dry-run
#   2. `sbatch tpv104_mesh_build.sbatch`          — produces .msh files
#
# Mesh-coupled init path — driver loads the .msh, sets up WaveOperator,
# initialises DOFData via InitializeFaultDOFs_TPV104, writes the t=0
# station row, exits.  tfinal=0 skips the ADER loop entirely; ~5 min
# even at 200 m because no time integration runs.

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

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas

RESULT_DIR="tpv104/results_init_stations_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o
make seas_tpv104_driver -j8 || { echo "BUILD FAILED"; exit 1; }

# Verify the .msh exists.
if [ ! -s tpv104/mesh/tpv104_200m.msh ]; then
   echo "ERROR: tpv104/mesh/tpv104_200m.msh missing."
   echo "       Run tpv104_mesh_build.sbatch first."
   exit 1
fi

echo "========================================"
echo "Git commit: $(git rev-parse HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "========================================"

ibrun ./seas_tpv104_driver \
      --mesh tpv104/mesh/tpv104_200m.msh \
      --tfinal 0.0 \
      --dt 1.0e-3 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law slip-srw \
      --fault-iterator substep \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv104

# Gate P2_D: ψ_ini at every FAULT station == 5.6359184e-01 ± 1e-8.
# The driver writes TWO families of station files:
#   - tpv104_station_x2_<y>_x3_<z>.dat  — 9 fault stations (ψ in col 9)
#   - tpv104_station_surf_<...>.dat      — 6 surface stations (no ψ col)
# This gate targets only the 9 fault stations.  Station file layout
# (per TPV104StationWriter, plan §4.4):
#   column 1 : t (== 0 for init)
#   column 9 : ψ (plan §5.1 SCEC trace format)
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
STATIONS_OK=0
STATIONS_TOTAL=0
MAX_ABS_ERR=0
for f in "${RESULT_DIR}"/tpv104_station_x2_*.dat; do
   if [ ! -s "${f}" ]; then continue; fi
   STATIONS_TOTAL=$((STATIONS_TOTAL + 1))
   # First non-comment row = t=0 snapshot.
   PSI=$(awk '!/^#/ && NF>=9 {print $9; exit}' "${f}")
   if [ -z "${PSI}" ]; then continue; fi
   ERR=$(awk -v p="${PSI}" 'BEGIN{ e = p - 0.56359184; if (e<0) e = -e; printf "%.3e\n", e }')
   if awk -v e="${ERR}" 'BEGIN{ exit !(e < 1.0e-8) }'; then
      STATIONS_OK=$((STATIONS_OK + 1))
   fi
   # Track worst error across all stations.
   MAX_ABS_ERR=$(awk -v m="${MAX_ABS_ERR}" -v e="${ERR}" \
                  'BEGIN{ printf "%.3e\n", (e > m ? e : m) }')
done

{
   echo "P2_D gate — ψ_ini station check"
   echo "stations OK     : ${STATIONS_OK}"
   echo "stations total  : ${STATIONS_TOTAL}"
   echo "max |ψ - ψ_ini| : ${MAX_ABS_ERR}"
   if [ "${STATIONS_TOTAL}" -eq 9 ] && [ "${STATIONS_OK}" -eq 9 ]; then
      echo "STATUS: PASS"
   else
      echo "STATUS: FAIL"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
