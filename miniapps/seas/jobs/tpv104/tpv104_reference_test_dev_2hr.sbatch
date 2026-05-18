#!/bin/bash
#SBATCH -J tpv104_ref_5s
#SBATCH -o tpv104_reference_%j.out
#SBATCH -e tpv104_reference_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 — SINGLE-SHOT reference run (dev queue, 2h, 8N x 400r).
# =======================================================================
# Produces the gold-standard 0 → 5.0 s trajectory against which the
# restart pair (`tpv104_restart_test_v1_dev_2hr.sbatch`'s Phase A +
# Phase B at the t=2.5 seam) is compared.
#
# Submit this sbatch in PARALLEL with the pair sbatch:
#   sbatch jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch
#   sbatch jobs/tpv104/tpv104_reference_test_dev_2hr.sbatch
# Both fit independently in dev queue 2 h; both can run concurrently
# on Frontera if the allocation has the cycles.  Once both finish,
# run the compare script on the login node:
#   scripts/tpv104_restart_compare.sh \
#     tpv104/results_pair_test_job<pair_jobid> \
#     tpv104/results_reference_job<ref_jobid>
# See debug_document/paraview_output_debug_document/
#   tpv104_restart_state_verification_2026-05-17.md
# for the verification design (Phase C reference + SEAM + REFERENCE
# tolerances + comparison column reference).
#
# Same driver flags / mesh / parameters as the pair's Phase A and
# Phase B (except --tfinal 5.0 in one go and NO --restart).  Any
# parameter drift between this sbatch and the pair sbatch invalidates
# the comparison.
# =======================================================================

export LC_ALL=C
export LANG=C

# Fail fast on module-load / build / pipe failures.
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
# NOTE: do NOT load phdf5 on Frontera.  build_frontera.sh builds HDF5
# 1.14.x from source into ${SEAS_MFEM_ROOT}/extern/hdf5/install.

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

# Mesh regen (TPV104 200m mesh is too large to commit; .geo is tracked).
MESH_SRC="tpv104/mesh/tpv104_200m.geo"
MESH_OUT="tpv104/mesh/tpv104_200m.msh"
GMSH_BIN="${SEAS_MFEM_ROOT}/extern/gmsh/bin/gmsh"
if [ ! -s "${MESH_OUT}" ] || [ "${FORCE_MESHGEN:-0}" = "1" ]; then
    if [ ! -x "${GMSH_BIN}" ]; then
        echo "ERROR: ${GMSH_BIN} not found.  Run build_frontera.sh first."
        exit 1
    fi
    if [ ! -s "${MESH_SRC}" ]; then
        echo "ERROR: ${MESH_SRC} missing from the checkout."
        exit 1
    fi
    echo "=== Regenerating ${MESH_OUT} from ${MESH_SRC} ==="
    "${GMSH_BIN}" -3 "${MESH_SRC}" -o "${MESH_OUT}" \
        || { echo "MESH GENERATION FAILED"; exit 1; }
fi

rm -f seas_tpv104_driver drivers/tpv104_driver.o
make seas_tpv104_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

# Per-job binary copy.  Two sbatches submitted concurrently would
# otherwise race on the single ./seas_tpv104_driver binary file:
# one sbatch's `rm` + `make` (truncate-and-write via gcc -o) clobbers
# the other sbatch's already-mapped binary mid-execution, producing
# the "Error mmapping the file: Bad file descriptor" + hard hang we
# saw at jobid 7729945.  Copying immediately after make to a
# per-jobid name decouples the run from any subsequent rebuilds.
DRIVER_BIN="./seas_tpv104_driver_${SLURM_JOB_ID:-$$}"
cp seas_tpv104_driver "${DRIVER_BIN}" \
   || { echo "Failed to copy binary to ${DRIVER_BIN}"; exit 1; }
chmod +x "${DRIVER_BIN}"
echo "Using per-job binary: ${DRIVER_BIN}"

# =======================================================================
# Reference single-shot output layout.  Distinct base dir from the
# pair sbatch's `results_pair_test_job<...>` so the compare script
# can locate both jobs' outputs by jobid.
# =======================================================================
REFERENCE_BASE="tpv104/results_reference_job${SLURM_JOB_ID}"
RESULT_DIR_R="${REFERENCE_BASE}/single_shot"
OUTPUT_PREFIX_R="tpv104_ref_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR_R}"

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Base dir:      ${REFERENCE_BASE}"
echo "Reference dir: ${RESULT_DIR_R}  (single_shot, 0 → 5.0 s)"
echo "Footprint:     8N x 400r, dev queue, 2h total"
echo "Reference tfinal: 5.0 s  (~2 hr wall, single-shot, no --restart)"
echo "PAIR COMPARISON: after this job AND the pair sbatch finish, run"
echo "  scripts/tpv104_restart_compare.sh <pair_base> ${REFERENCE_BASE}"
echo "========================================"

# =======================================================================
# Reference run — single-shot 0 → 5.0 s.  Same driver flags as the
# pair sbatch's Phase A (and Phase B except --restart and tfinal).
# Any drift in parameters between this sbatch and the pair sbatch
# invalidates the comparison.
# =======================================================================
echo ""
echo "========================================"
echo "Reference: single-shot TPV104 run, 0 → 5.0 s  (~2 hr wall)"
echo "========================================"
ibrun "${DRIVER_BIN}" \
      --mesh tpv104/mesh/tpv104_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 5.0 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law aging \
      --fault-iterator one-shot \
      --mixed-flux adjacent \
      --output-dir "${RESULT_DIR_R}" \
      --output-prefix "${OUTPUT_PREFIX_R}" \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.1 \
      --paraview-bulk-dt 0.1 \
      --paraview-volume-zfp-tol 1e-3 \
      --paraview-bulk-zfp-tol 1e-3 \
      --paraview-fault-zfp-tol 1e-12 \
      --paraview-max-snapshots 50 \
      --checkpoint-interval 100 \
      || REF_RC=$?
REF_RC=${REF_RC:-0}
echo "Reference exit code: ${REF_RC}"
if [ "${REF_RC}" -ne 0 ]; then
   echo "FAIL: Reference run did not exit cleanly."
   exit 1
fi

# =======================================================================
# Validation block — confirms the reference run actually completed.
# State-correctness comparison vs the pair sbatch happens in the
# external compare script.
#
# Acceptance checks (exit codes 2-4):
#   #1 Reference fault.vtkhdf landed (>1 KiB)
#   #2 Reference produced *_station_*.dat files
#   #3 Reference last station t ≈ tfinal = 5.0 s
# =======================================================================
echo ""
echo "========================================"
echo "Validation summary (reference single-shot)"
echo "========================================"

# Acceptance check #1 — reference fault.vtkhdf landed.
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

# Acceptance check #2 — reference produced station files.
echo ""
echo "Validation #2: reference produced *_station_*.dat files"
echo "----------------------------------------"
R_STATIONS=$(ls "${RESULT_DIR_R}"/*_station_*.dat 2>/dev/null | wc -l | tr -d ' ')
echo "  Reference station files: ${R_STATIONS}"
if [ "${R_STATIONS}" -gt 0 ]; then
   echo "  PASS — reference produced ${R_STATIONS} station files"
else
   echo "  FAIL — reference produced no *_station_*.dat files"
   exit 3
fi

# Acceptance check #3 — reference last station t ≈ tfinal.
echo ""
echo "Validation #3: reference last station t ≈ tfinal = 5.0 s"
echo "----------------------------------------"
R_STATIONS_FIRST=$(ls "${RESULT_DIR_R}"/*_station_*.dat 2>/dev/null | head -1)
R_LAST_T=$(awk '/^[[:space:]]*[0-9-]/ { last=$1 } END { print last }' \
           "${R_STATIONS_FIRST}" 2>/dev/null || echo "")
echo "  Last t in ${R_STATIONS_FIRST}: ${R_LAST_T:-(unparseable)}"
R_LAST_T_OK=$(awk -v t="${R_LAST_T:-0}" \
              'BEGIN { print (t >= 4.9) ? "1" : "0" }')
if [ "${R_LAST_T_OK}" = "1" ]; then
   echo "  PASS — reference reached t=${R_LAST_T} (>= 4.9, near tfinal=5.0)"
else
   echo "  FAIL — reference last t=${R_LAST_T}, expected >= 4.9."
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
   echo "TPV104 reference single-shot — SUMMARY"
   echo "Git commit:  $(git rev-parse HEAD)"
   echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
   echo "Slurm job:   ${SLURM_JOB_ID}"
   echo "Footprint:   8N x 400r, dev queue, 2h reservation"
   echo ""
   echo "Output: ${RESULT_DIR_R}"
   echo "  fault.vtkhdf: ${FAULT_R_SIZE} bytes"
   echo "  station files: ${R_STATIONS}"
   echo "  last t:       ${R_LAST_T}"
   echo ""
   echo "Acceptance checks:"
   echo "  #1 fault.vtkhdf landed:           PASS (${FAULT_R_SIZE} bytes)"
   echo "  #2 station files present:         PASS (${R_STATIONS})"
   echo "  #3 reached tfinal = 5.0 s:        PASS (t=${R_LAST_T})"
   echo ""
   echo "NEXT STEP: run the compare script after the pair sbatch finishes:"
   echo "  scripts/tpv104_restart_compare.sh \\"
   echo "    <pair_base> ${REFERENCE_BASE}"
   echo "  where <pair_base> is the pair sbatch's"
   echo "  tpv104/results_pair_test_job<jobid> directory."
   echo ""
   echo "du -sh: $(du -sh ${RESULT_DIR_R} | cut -f1)"
} > "${SUMMARY_FILE}"
cat "${SUMMARY_FILE}"

# Clean up the per-job binary copy (build dir bloat protection).
# Failure here is non-fatal — the SUMMARY has already been written.
rm -f "${DRIVER_BIN}" || true
