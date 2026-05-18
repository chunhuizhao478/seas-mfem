#!/bin/bash
#SBATCH -J bp5_ref_3yr
#SBATCH -o bp5_reference_%j.out
#SBATCH -e bp5_reference_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# BP5 V2 — SINGLE-SHOT reference run (dev queue, 2h, 8N x 400r).
# =======================================================================
# Produces the gold-standard 0 → 3 yr trajectory against which the
# restart pair (`bp5_restart_test_v2_dev_2hr.sbatch`'s Phase A 0 →
# 0.001 yr + Phase B 0.001 → 3 yr at the t=0.001 yr seam) is compared.
#
# Submit this sbatch in PARALLEL with the pair sbatch:
#   sbatch jobs/bp5/bp5_restart_test_v2_dev_2hr.sbatch
#   sbatch jobs/bp5/bp5_reference_test_v2_dev_2hr.sbatch
# Both fit independently in dev queue 2 h.  Once both finish, run
# the compare script on the login node:
#   scripts/bp5_restart_compare.sh \
#     bp5/results_restart_test_job<pair_jobid> \
#     bp5/results_reference_job<ref_jobid>
#
# Wall budget: BP5 adaptive TS dt is tiny in the 0 → ~0.001 yr
# initial transient (~12500 steps), then much larger in the post-
# transient interseismic regime.  Combined 0 → 3 yr crosses both
# regimes; ~16000-17000 PETSc TS steps total at ~140 steps/min on
# 8N x 400r → ~115 min wall.  Same parameters as the pair's Phase A
# and Phase B (except --tfinal 3 yr in one go and no --restart).
# Any parameter drift between this sbatch and the pair sbatch
# invalidates the comparison.

export LC_ALL=C
export LANG=C

set -euo pipefail

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

if [ -z "${SEAS_MFEM_ROOT:-}" ]; then
    _probe="${SLURM_SUBMIT_DIR:-$PWD}"
    while [ "${_probe}" != "/" ] && [ -n "${_probe}" ]; do
        if [ -f "${_probe}/build_frontera.sh" ]; then
            SEAS_MFEM_ROOT="${_probe}"
            break
        fi
        _probe="$(dirname "${_probe}")"
    done
    unset _probe
fi
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-/scratch2/10024/zhaochun/seas-project/seas-mfem}"
echo "  SEAS_MFEM_ROOT = ${SEAS_MFEM_ROOT}"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

HDF5_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/hdf5/install"
H5Z_ZFP_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install"

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${HDF5_BUILT_DIR}/lib:${LD_LIBRARY_PATH:-}"

if [ -d "${H5Z_ZFP_BUILT_DIR}/plugin" ]; then
   export HDF5_PLUGIN_PATH="${H5Z_ZFP_BUILT_DIR}/plugin"
fi

cd "${SEAS_MFEM_ROOT}"

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES'    config/config.mk || { echo "ERROR: no PETSc";   exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'    config/config.mk || { echo "ERROR: no MUMPS";   exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'     config/config.mk || { echo "ERROR: no HDF5";    exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES'  config/config.mk || { echo "ERROR: no H5Z-ZFP"; exit 1; }

cd miniapps/seas

make clean
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

# Per-job binary copy (concurrent-rebuild race protection — same
# Frontera jobid 7729945 mmap-hang root cause that hit TPV104).
DRIVER_BIN="./seas_bp5_full_${SLURM_JOB_ID:-$$}"
cp seas_bp5_full "${DRIVER_BIN}" \
   || { echo "Failed to copy binary to ${DRIVER_BIN}"; exit 1; }
chmod +x "${DRIVER_BIN}"
echo "Using per-job binary: ${DRIVER_BIN}"

# =======================================================================
# Reference single-shot output layout.  Distinct base dir from the
# pair sbatch's `bp5/results_restart_test_job<...>` so the compare
# script can locate both jobs' outputs by jobid.
# =======================================================================
REFERENCE_BASE="bp5/results_reference_job${SLURM_JOB_ID}"
RESULT_DIR_R="${REFERENCE_BASE}/single_shot"
OUTPUT_PREFIX_R="bp5_ref_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR_R}"

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Base dir:      ${REFERENCE_BASE}"
echo "Reference dir: ${RESULT_DIR_R}  (single_shot, 0 → 3 yr)"
echo "Footprint:     8N x 400r, dev queue, 2h total"
echo "Reference tfinal: 9.46728e7 s ≈ 3 yr (~115 min wall — crosses"
echo "                  the 0-to-0.001-yr small-dt transient AND the"
echo "                  post-transient interseismic regime)"
echo "PAIR COMPARISON: after this job AND the pair sbatch finish, run"
echo "  scripts/bp5_restart_compare.sh <pair_base> ${REFERENCE_BASE}"
echo "========================================"

# =======================================================================
# Reference run — single-shot 0 → 3 yr.  Same driver flags as the
# pair sbatch's Phase A (and Phase B except --restart and tfinal).
# =======================================================================
echo ""
echo "========================================"
echo "Reference: single-shot BP5 V2 run, 0 → 3 yr (~115 min wall)"
echo "========================================"
ibrun "${DRIVER_BIN}" \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --ref-dir bp5/benchmark_data \
      --output-dir "${RESULT_DIR_R}" \
      --output-prefix "${OUTPUT_PREFIX_R}" \
      --checkpoint-interval 5000 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --paraview-adaptive \
      --no-volume-pv \
      --paraview-fault-zfp-tol 1e-12 \
      --paraview-dt-co     1.0 \
      --paraview-dt-nu     43200 \
      --paraview-dt-inter-yr 1.0 \
      --tfinal 94672800 || REF_RC=$?
REF_RC=${REF_RC:-0}
echo "Reference exit code: ${REF_RC}"
# BP5 driver's post-simulation regression block returns 1 on missing
# per-jobid reference files (informational mode, expected).  PASS if
# the reference run wrote a fault.vtkhdf — that's all we need for the
# external compare script.
if [ -s "${RESULT_DIR_R}/fault.vtkhdf" ]; then
   if [ "${REF_RC}" -ne 0 ]; then
      echo "WARN: Reference exited ${REF_RC} but fault.vtkhdf exists; continuing."
      echo "  (informational regression check returned non-zero — expected"
      echo "   when per-jobid reference files are missing in --ref-dir)"
   fi
else
   echo "FAIL: Reference did not produce fault.vtkhdf."
   exit 1
fi

# =======================================================================
# Validation block — confirms the reference run actually completed.
# State-correctness comparison vs the pair sbatch happens in the
# external compare script.
#
# Acceptance checks (exit codes 2-4):
#   #1 Reference fault.vtkhdf landed (>1 KiB)
#   #2 Reference produced *_fltst_*.txt files
#   #3 Reference last station t ≈ tfinal = 9.46728e7 s (= 3 yr)
# =======================================================================
echo ""
echo "========================================"
echo "Validation summary (BP5 reference single-shot)"
echo "========================================"

# #1 — fault.vtkhdf
echo ""
echo "Validation #1: reference fault.vtkhdf landed"
echo "----------------------------------------"
FAULT_R="${RESULT_DIR_R}/fault.vtkhdf"
FAULT_R_SIZE=$(stat -c %s "${FAULT_R}" 2>/dev/null \
            || stat -f %z "${FAULT_R}" 2>/dev/null \
            || echo 0)
echo "  Reference fault.vtkhdf size: ${FAULT_R_SIZE} bytes"
if [ "${FAULT_R_SIZE}" -gt 1024 ]; then
   echo "  PASS — reference wrote fault.vtkhdf (> 1 KiB)"
else
   echo "  FAIL — reference fault.vtkhdf missing or empty"
   exit 2
fi

# #2 — *_fltst_*.txt station files
echo ""
echo "Validation #2: reference produced *_fltst_*.txt files"
echo "----------------------------------------"
R_STATIONS=$(ls "${RESULT_DIR_R}"/*_fltst_*.txt 2>/dev/null | wc -l | tr -d ' ')
echo "  Reference station files: ${R_STATIONS}"
if [ "${R_STATIONS}" -gt 0 ]; then
   echo "  PASS — reference produced ${R_STATIONS} station files"
else
   echo "  FAIL — reference produced no *_fltst_*.txt files"
   exit 3
fi

# #3 — last station t close to tfinal
echo ""
echo "Validation #3: reference last station t ≈ tfinal = 9.46728e7 s"
echo "----------------------------------------"
R_STATIONS_FIRST=$(ls "${RESULT_DIR_R}"/*_fltst_*.txt 2>/dev/null | head -1)
R_LAST_T=$(awk '/^[[:space:]]*[0-9-]/ { last=$1 } END { print last }' \
           "${R_STATIONS_FIRST}" 2>/dev/null || echo "")
echo "  Last t in ${R_STATIONS_FIRST}: ${R_LAST_T:-(unparseable)}"
# Tolerance: within 5% of tfinal (BP5 output cadence is adaptive,
# so last row may be a bit before tfinal).
R_LAST_T_OK=$(awk -v t="${R_LAST_T:-0}" \
              'BEGIN { print (t >= 9.0e7) ? "1" : "0" }')
if [ "${R_LAST_T_OK}" = "1" ]; then
   echo "  PASS — reference reached t=${R_LAST_T} s (>= 9.0e7, near tfinal)"
else
   echo "  FAIL — reference last t=${R_LAST_T} s, expected >= 9.0e7."
   echo "         The time loop terminated before reaching tfinal."
   exit 4
fi

# =======================================================================
# Summary + handoff to the compare script.
# =======================================================================
echo ""
echo "========================================"
echo "Reference run acceptance checks PASSED"
echo "========================================"

SUMMARY_FILE="${RESULT_DIR_R}/REFERENCE_SUMMARY.txt"
{
   echo "BP5 V2 reference single-shot — SUMMARY"
   echo "Git commit:  $(git rev-parse HEAD)"
   echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
   echo "Slurm job:   ${SLURM_JOB_ID}"
   echo "Footprint:   8N x 400r, dev queue, 2h reservation"
   echo ""
   echo "Output: ${RESULT_DIR_R}"
   echo "  fault.vtkhdf:  ${FAULT_R_SIZE} bytes"
   echo "  station files: ${R_STATIONS}"
   echo "  last t (s):    ${R_LAST_T}"
   echo ""
   echo "Acceptance checks:"
   echo "  #1 fault.vtkhdf landed:           PASS (${FAULT_R_SIZE} bytes)"
   echo "  #2 station files present:         PASS (${R_STATIONS})"
   echo "  #3 reached tfinal ≈ 3 yr:         PASS (t=${R_LAST_T} s)"
   echo ""
   echo "NEXT STEP: run the compare script after the pair sbatch finishes:"
   echo "  scripts/bp5_restart_compare.sh \\"
   echo "    <pair_base> ${REFERENCE_BASE}"
   echo "  where <pair_base> is the pair sbatch's"
   echo "  bp5/results_restart_test_job<pair_jobid> directory."
   echo ""
   echo "du -sh: $(du -sh ${RESULT_DIR_R} | cut -f1)"
} > "${SUMMARY_FILE}"
cat "${SUMMARY_FILE}"

# Clean up the per-job binary copy (build dir bloat protection).
rm -f "${DRIVER_BIN}" || true
