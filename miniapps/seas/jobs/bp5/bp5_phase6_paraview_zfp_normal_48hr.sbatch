#!/bin/bash
#SBATCH -J bp5_p6_zfp
#SBATCH -o bp5_p6_zfp_%j.out
#SBATCH -e bp5_p6_zfp_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# BP5 — Phase 6 ParaView compaction reference run (normal queue, 48h).
# =======================================================================
# Production-scale BP5 quasi-dynamic run aimed at REFERENCE_RUNS.md
# row #1 (Plan §Phase 7.5 acceptance criterion #1 — fault-only ZFP).
#
# Plan §7.5 #1 originally specified 800 ranks; this sbatch uses
# 8 nodes × 400 ranks (50 ranks/node on a 56-core CLX node) per the
# user's compute-footprint request.  REFERENCE_RUNS.md row #1 records
# the 400-rank run; the 800-rank variant is a future-calibration row.
#
# R-507 fix: `--no-volume-pv` is now passed (plan's "fault-only"
# criterion).  Volume displacement output suppressed; only fault
# surface (slip / slip-rate / traction / state-variable / normal-stress)
# is emitted, at ZFP tol 1e-12.
#
# Phase 6 ParaView flag pair in use:
#   --paraview-fault-zfp-tol 1e-12  (slip-rate floor preserved)
#   --no-volume-pv                  (volume displacement suppressed)
#
# BP5 has NO secondary collection so --paraview-bulk-* is not used.
#
# Snapshot cap (--paraview-max-snapshots 5000) bounds the writer count
# at the production default.  The Python schedule integrator
# (miniapps/seas/scripts/_io_size_schedule.py) estimates ~86k uncapped
# writes over 250 yr at n_events=1 (rounded from 250 yr / 240 yr
# recurrence per CLAUDE.md "What Constitutes a Regression"), dominated
# by the BP5 nucleation budget (1 event × NUCLEATION_DURATION_S["bp5"]
# = 86 400 s, at dt_nucleation = 1.0 s ⇒ 86 400 writes).  The 5000 cap
# therefore engages comfortably below the projected ceiling and bounds
# the interseismic write count via paraview_output.hpp:
# RecomputeIntervalForCap.  Cross-check:
#   python3 miniapps/seas/scripts/estimate_from_sbatch.py \
#       miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch
#
# R-701: --tfinal is the integer literal 7_889_400_000 (= 250 *
# BP5Params::seconds_per_year = 250 * 31_557_600), NOT a $(python3 ...)
# subshell.  Earlier revisions used a subshell that silently produced
# t_final = 0 if python3 was missing on the compute node (atof of the
# next flag), wasting the full 48 h reservation.
#
# Wall-time budget:
#   8 nodes x 400 ranks x 48h.  Production cycle target: 250 yr.

export LC_ALL=C
export LANG=C

# R-510: fail fast on module-load / build / pipe failures.
set -euo pipefail

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

# Project checkout root.  Auto-detect by walking up from $SLURM_SUBMIT_DIR
# (or $PWD if not running under Slurm) until we find build_frontera.sh.
# Override at submit time when running from a different clone:
#   SEAS_MFEM_ROOT=$(pwd) sbatch ...
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
# Last-resort fallback if the probe walked off the top of the filesystem.
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-/scratch2/10024/zhaochun/seas-project/seas-mfem}"
echo "  SEAS_MFEM_ROOT = ${SEAS_MFEM_ROOT}"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
# NOTE: do NOT load phdf5 on Frontera.  The phdf5/1.14.x module
# requires intel/23.1.0 + impi/21.9.0 (incompatible with the intel/19
# hypre/mumps/petsc stack we use), and older modules (1.10/1.12) lack
# symbols MFEM mesh/vtkhdf.cpp needs (H5S_BLOCK, the unsigned-long-long
# hsize_t typedef).  build_frontera.sh builds HDF5 1.14.x from source
# into ${SEAS_MFEM_ROOT}/extern/hdf5/install instead; we link/run
# against that.

HDF5_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/hdf5/install"
H5Z_ZFP_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install"

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${HDF5_BUILT_DIR}/lib:${LD_LIBRARY_PATH:-}"

if [ -d "${H5Z_ZFP_BUILT_DIR}/plugin" ]; then
   export HDF5_PLUGIN_PATH="${H5Z_ZFP_BUILT_DIR}/plugin"
else
   echo "WARNING: ${H5Z_ZFP_BUILT_DIR}/plugin not found — run build_frontera.sh first."
fi

cd "${SEAS_MFEM_ROOT}"

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES'    config/config.mk || { echo "ERROR: no PETSc";   exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'    config/config.mk || { echo "ERROR: no MUMPS";   exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'     config/config.mk || { echo "ERROR: no HDF5";    exit 1; }
# R-404: fail before ibrun if H5Z-ZFP is missing — the driver aborts at
# parse time on --paraview-fault-zfp-tol otherwise, wasting ~30 SU.
grep -q '^MFEM_USE_H5Z_ZFP *= YES'  config/config.mk || { echo "ERROR: no H5Z-ZFP (--paraview-fault-zfp-tol requires it)"; exit 1; }

cd miniapps/seas
make clean
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="bp5/results_phase6_zfp_normal_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_p6_zfp_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          bp5_tandem_exact.msh (1000 m, P=1 IP)"
echo "Volume PV:     SUPPRESSED (--no-volume-pv per Plan §7.5 #1)"
echo "Fault PV:      HDF5,           ZFP tol 1e-12 (slip-rate floor)"
echo "Snapshot cap:  --paraview-max-snapshots 5000"
echo "Queue:         normal (8N x 400r, 48h, tfinal 250 yr)"
echo "========================================"

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --ref-dir bp5/benchmark_data \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --checkpoint-interval 5000 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --paraview-adaptive \
      --no-volume-pv \
      --paraview-fault-zfp-tol 1e-12 \
      --paraview-max-snapshots 5000 \
      --tfinal 7889400000 \
      --verify

RESULT_FILE="${RESULT_DIR}/SIZE.txt"
{
   echo "BP5 Phase 6 ParaView ZFP — size summary"
   echo "Mesh:   bp5_tandem_exact.msh (1000 m)"
   echo "Output: ${RESULT_DIR}"
   echo
   echo "du -sb (total):"
   du -sb "${RESULT_DIR}"
   echo
   echo "du -sh per top-level entry:"
   du -sh "${RESULT_DIR}"/* 2>/dev/null | sort -h
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
