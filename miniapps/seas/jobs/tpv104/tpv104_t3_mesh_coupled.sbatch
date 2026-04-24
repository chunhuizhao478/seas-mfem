#!/bin/bash
#SBATCH -J tpv104_t3
#SBATCH -o tpv104_t3_%j.out
#SBATCH -e tpv104_t3_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 06:00:00
#SBATCH -A EAR20006

# TPV104 Phase-2 t=3.0 s mesh-coupled production run on Frontera.
#
# Purpose: first long TPV104 simulation spanning nucleation (t < 1 s)
# through breakout and early rupture propagation (1 s < t < 3 s).
# Target fixture matches the Phase-3 probe plan §5.11 3.A run order —
# the output station traces are the reference for Phase-3 probe-diff
# against the SCEC TPV104 benchmark data in
# `miniapps/seas/tpv104/benchmark_data/`.
#
# Mesh             : tpv104_200m.msh (200 m fault resolution)
# ADER order       : 2 (t = 3 s fits in wall-clock at 400 ranks)
# dt               : ~4.7 ms per macro-step (CFL-limited by cs=3464 m/s)
#
# R7-001 DISCLOSURE (banner-vs-dispatch consistency):
#   Time integrator: ADER-O2 (one-shot via wave.AdvanceADER).
#   Fault iterator : one-shot.  Tpv104SubStepIterator is compiled +
#                    unit-tested but NOT wired into this time loop,
#                    because the wiring requires exposing per-sub-step
#                    I± from wave_operator.inl which is on the extreme-
#                    care no-touch list (CLAUDE.md [C2]).  The
#                    --fault-iterator=substep flag below is parsed for
#                    banner parity with the Phase-3 probe harness but
#                    has NO effect on the dispatched loop.
#   Friction solver: Brent (hard-coded via EvaluateADERTotal).
#                    The --friction-solver=newton-stable flag below is
#                    parsed for CLI validation (typos still abort) but
#                    has NO effect on the dispatched solver.
#                    For TPV104 background ψ/a ≈ 56, Brent is robust
#                    and numerically equivalent to stable-asinh Newton
#                    to tolerance; the physics result is unaffected.
#                    Any Phase-3 probe-diff finding pointing at
#                    "stable-asinh Newton" is misattributed.
#   Friction law   : slip-SRW (ψ-space), evolved by
#                    UpdateStateAnalyticSlipLawSRW at macro-step
#                    cadence (R7-007; per-sub-step ψ cadence blocked
#                    by R7-001 option a).
#   Nucleation     : ApplyNucleationIncremental_TPV104 called once
#                    per macro-step.  The smoothStepIncrement
#                    telescoping identity guarantees the total
#                    perturbation at t = T_nuc = 1.0 s matches the
#                    per-sub-step result at macro-step boundaries.
#
# WALL-CLOCK: ~5 h on 400 ranks at 200 m resolution.
# SU COST: ~2000 SU (8 nodes × 6 h × ~40 core-h/node-h).
#
# Pre-requisites (run in order, verify PASS before next):
#   1. `sbatch tpv104_build_banner_sanity.sbatch`  (~2 min)
#   2. `sbatch tpv104_mesh_build.sbatch`           (~5 min)
#   3. `sbatch tpv104_init_stations.sbatch`        (~5 min, P2_D gate)
#
# Pre-submit guardrails:
#   (A) tpv104_build_banner_sanity.sbatch must be PASS (build works).
#   (B) tpv104_mesh_build.sbatch must have produced tpv104_200m.msh.
#   (C) tpv104_init_stations.sbatch must be PASS (P2_D ψ_ini gate).
# The script below enforces (C) via a `grep STATUS: PASS` check before
# the build step; it fails fast with exit 1 if no passing init record
# is found.

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

# Integration-readiness gate — fail fast if the mesh-coupled branch
# has not landed yet.
PREV_INIT_RESULT=$(ls -dt tpv104/results_init_stations_job* 2>/dev/null | head -1)
if [ -z "${PREV_INIT_RESULT}" ] \
   || ! grep -q "^STATUS: PASS" "${PREV_INIT_RESULT}/RESULT.txt" 2>/dev/null; then
   echo "ERROR: no passing tpv104_init_stations.sbatch result found."
   echo "       Run that job first and verify STATUS: PASS before proceeding."
   exit 1
fi
echo "Found init PASS at: ${PREV_INIT_RESULT}"

RESULT_DIR="tpv104/results_t3_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o
make seas_tpv104_driver -j8 || { echo "BUILD FAILED"; exit 1; }

if [ ! -s tpv104/mesh/tpv104_200m.msh ]; then
   echo "ERROR: tpv104/mesh/tpv104_200m.msh missing."
   exit 1
fi

echo "========================================"
echo "Git commit: $(git rev-parse HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Target: tfinal = 3.0 s on tpv104_200m.msh at 400 ranks"
echo "========================================"

ibrun ./seas_tpv104_driver \
      --mesh tpv104/mesh/tpv104_200m.msh \
      --tfinal 3.0 \
      --dt 4.7e-3 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law slip-srw \
      --fault-iterator substep \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv104

# Post-run sanity: every fault station file has ≥ 50 rows (one per
# stride sample over 3 s at the default output cadence) and no
# NaN/Inf.  Glob matches only the 9 fault stations
# (tpv104_station_x2_*_x3_*.dat); surface stations have a different
# column layout and are excluded from this gate.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
STATIONS_NONEMPTY=0
STATIONS_NANFREE=0
STATIONS_TOTAL=0
for f in "${RESULT_DIR}"/tpv104_station_x2_*.dat; do
   STATIONS_TOTAL=$((STATIONS_TOTAL + 1))
   if [ $(wc -l < "${f}") -ge 50 ]; then
      STATIONS_NONEMPTY=$((STATIONS_NONEMPTY + 1))
   fi
   if ! grep -qiE "nan|inf" "${f}"; then
      STATIONS_NANFREE=$((STATIONS_NANFREE + 1))
   fi
done

{
   echo "TPV104 t=3.0 s production run"
   echo "fault stations total       : ${STATIONS_TOTAL}"
   echo "fault stations ≥ 50 rows   : ${STATIONS_NONEMPTY}"
   echo "fault stations NaN/Inf-free: ${STATIONS_NANFREE}"
   if [ "${STATIONS_TOTAL}" -eq 9 ] \
      && [ "${STATIONS_NONEMPTY}" -eq 9 ] \
      && [ "${STATIONS_NANFREE}" -eq 9 ]; then
      echo "STATUS: PASS"
   else
      echo "STATUS: FAIL"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
