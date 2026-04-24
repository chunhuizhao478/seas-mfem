#!/bin/bash
#SBATCH -J tpv104_banner
#SBATCH -o tpv104_banner_%j.out
#SBATCH -e tpv104_banner_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:05:00
#SBATCH -A EAR20006

# TPV104 Phase-2 build + banner sanity on Frontera.
#
# Purpose: cheapest possible verification that:
#   (1) the TPV104 driver compiles on Frontera with current HEAD,
#   (2) the banner prints the R7-001-honest disclosure strings
#       (Brent hard-coded via EvaluateADERTotal; sub-step iterator
#       not wired; --friction-solver / --fault-iterator banner-only),
#   (3) --dry-run --verify-dispatch tri-consistency passes (banner
#       claims match the [dispatch] actual tags),
#   (4) --dry-run wiring path completes without NaN.
#
# NOT a physics run — no mesh is loaded.  Pre-requisite for the
# init + t=3.0 s mesh-coupled runs below.
#
# WALL-CLOCK: ~2 min on 1 development node.
# SU COST: ~0.04 SU.

export LC_ALL=C
export LANG=C

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

# Module list copied VERBATIM from the working TPV102 init sbatch
# (`tpv102_200m_p1_0.01s_50rank_init.sbatch`) — per
# `feedback_sbatch_modules.md`: do not drop modules, even ones that
# seem unrelated.  fftw3 is pulled in transitively by petsc; dropping
# it has been seen to break runtime linker resolution.
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

echo "========================================"
echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Working dir: $(pwd)"
echo "========================================"

# Force full rebuild — TPV104 changes frequently during Phase 2.
rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o
make seas_tpv104_driver -j"${JOBS}" \
   || { echo "BUILD FAILED"; exit 1; }

echo "========================================"
echo "[banner + dry-run check]"
echo "========================================"

# Run with --verify-dispatch so the tri-consistency lines print too.
./seas_tpv104_driver --dry-run --verify-dispatch 2>&1 | tee tpv104_banner_output.txt

# Success criteria (machine-checkable, R7-001-honest + R8-004 per-rank):
#   - banner contains "Friction solver: Brent (hard-coded via EvaluateADERTotal"
#   - banner contains "Fault iterator: one-shot (sub-step iterator NOT wired"
#   - dispatch line  "[dispatch] rank=0 friction_solver_actual=brent"
#   - dispatch line  "[dispatch] rank=0 fault_iterator_actual=oneshot"
#   - dispatch line  "[dispatch] rank=0 friction_law_actual=slip-srw"
#   - "[dry-run] OK." at end
#   - no "NaN" or "FAIL" in output
RESULT_FILE="tpv104_banner_result.txt"

if grep -q "Friction solver: Brent (hard-coded via EvaluateADERTotal" tpv104_banner_output.txt \
   && grep -q "Fault iterator: one-shot (sub-step iterator NOT wired" tpv104_banner_output.txt \
   && grep -q "\[dispatch\] rank=0 friction_solver_actual=brent" tpv104_banner_output.txt \
   && grep -q "\[dispatch\] rank=0 fault_iterator_actual=oneshot" tpv104_banner_output.txt \
   && grep -q "\[dispatch\] rank=0 friction_law_actual=slip-srw" tpv104_banner_output.txt \
   && grep -q "\[dry-run\] OK\." tpv104_banner_output.txt \
   && ! grep -qiE "NaN|FAIL|MFEM_ABORT" tpv104_banner_output.txt; then
   echo "PASS: TPV104 driver banner + dry-run sanity complete" > "${RESULT_FILE}"
   grep -E "Friction solver:|Fault iterator:|Time integrator:|Friction law:|\[dispatch\]|\[dry-run\]" \
        tpv104_banner_output.txt >> "${RESULT_FILE}"
else
   echo "FAIL: banner check or dry-run output missing expected lines" > "${RESULT_FILE}"
   echo "--- banner output ---" >> "${RESULT_FILE}"
   cat tpv104_banner_output.txt >> "${RESULT_FILE}"
fi

echo "========================================"
cat "${RESULT_FILE}"
echo "========================================"
