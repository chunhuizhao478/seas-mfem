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
# Phase 6 ParaView flag set in use:
#   --paraview-fault-zfp-tol 1e-12  (slip-rate floor preserved;
#                                    must stay <= 1e-12 because the
#                                    interseismic V_p floor is 1e-9 m/s
#                                    and ZFP absolute-error tol 1e-12
#                                    keeps ~3 decades of precision at
#                                    the floor; raising tol erases the
#                                    floor).
#   --no-volume-pv                  (volume displacement suppressed)
#   --paraview-adaptive             (V_max-driven schedule; see below)
#   --paraview-dt-co     1.0        (coseismic cadence,   1 s)
#   --paraview-dt-nu     43200      (nucleation cadence,  12 hr = 0.5 day)
#   --paraview-dt-inter-yr 1.0      (interseismic cadence, 1 yr)
#
# BP5 has NO secondary collection so --paraview-bulk-* is not used.
#
# --------------------------------------------------------------------
# Per-regime cadence and expected frame count
# --------------------------------------------------------------------
# Empirical phase durations (from results_phase7_production_job7650438
# over 1775.58 simulated years, ~7 events — see
# debug_document/paraview_output_debug_document/regime_budget_2026-05-16.md):
#
#     Phase        | duration / event | total / 1775 yr run
#     -------------+------------------+--------------------
#     interseismic | ~250 yr          | 1775 yr
#     nucleation   | ~19 days         | 134 days
#     coseismic    | ~7 min           | 51 min = 3,080 s
#
# Projected frame count at the cadences set above (1775 yr × 7 events):
#
#     Phase        | cadence  | frames / 1775 yr run
#     -------------+----------+---------------------
#     interseismic | 1 yr     | ~1,775
#     nucleation   | 12 hr    | ~268    (134 days / 12 hr)
#     coseismic    | 1 s      | ~3,080  (3,080 s / 1 s)
#     ----------------------------------+--------------
#                                       | ~5,123 total
#
# Same cadences scaled to the ~634 yr / ~2-3 event production target
# (--tfinal 2e10 s; steady-state recurrence ~240 yr):
#     interseismic ~634 + (2-3) × (nucleation ~38 + coseismic ~420)
#                       = ~1,550 to ~2,000 total frames
# At ~1.4 MB/frame (12 fields × ~6000 fault cells × ZFP 1e-12):
#     fault.vtkhdf ≈ 2-3 GB.
#
# These are the actual numbers of fault.vtkhdf timesteps the writer
# will produce.  No snapshot cap is applied.
#
# Why no --paraview-max-snapshots?
# An earlier revision tried `--paraview-max-snapshots 5000` as a
# defensive ceiling.  REVIEW.md round-2 (2026-05-16) showed the
# hard-cap implementation (R-002) silently sub-samples coseismic and
# nucleation to ~T/K = 18 days when K=5000, erasing all earthquake
# detail.  Rather than rebuild the cap to be smarter, this sbatch just
# trusts the cadence projection above (~1,500-2,000 frames for the
# ~634 yr / 2-3 event target, ~5,123 for an 1775-yr / 7-event extended
# run).  Both fit well within any reasonable Frontera scratch budget at
# 1.4 MB/frame fault.vtkhdf with ZFP 1e-12 — under 10 GB worst case.
#
# IMPORTANT CAVEAT.  --paraview-dt-co/-nu/-inter-yr only affect ParaView
# fault/volume snapshots.  The station probe files (*_fltst_*.txt the
# SCEC benchmark consumes) go through BP5BenchmarkOutput::OutputInterval
# which has its OWN hard-coded thresholds and ignores these CLI flags
# (miniapps/seas/io/bp5_benchmark_output.hpp, the
# `BP5BenchmarkOutput::OutputInterval` static method).  The bench_out
# probes therefore still write at:
#     coseismic    V > 1e-3 m/s         every 0.1 s    (SCEC spec)
#     nucleation   1e-6 < V <= 1e-3     every 0.1 s
#     interseismic V <= 1e-6 m/s        every 0.1 yr   (SCEC spec)
# regardless of the --paraview-dt-* flags.  Only the ParaView .vtkhdf
# size is reduced by those flags.
#
# Cross-check disk budget:
#   python3 miniapps/seas/scripts/estimate_from_sbatch.py \
#       miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch
#
# R-701: --tfinal is the integer literal 20_000_000_000 (= 2e10 s ≈
# 633.78 yr at BP5Params::seconds_per_year = 31_557_600 s/yr), NOT a
# $(python3 ...) subshell.  Earlier revisions used a subshell that
# silently produced t_final = 0 if python3 was missing on the compute
# node (atof of the next flag), wasting the full 48 h reservation.
#
# Why 633 yr instead of the original 250 yr?  BP5 quasi-dynamic
# events repeat at roughly the same period after the startup
# transient — running 250 yr captures ~1 event; running ~634 yr
# captures 2-3 events (at the SCEC-spec ~240 yr recurrence) and
# confirms the cycle is genuinely repetitive.  Beyond ~3 events the
# trajectory is just more of the same and produces little new
# diagnostic value.
#
# Wall-time budget:
#   8 nodes x 400 ranks x 48h.  Production cycle target: ~634 yr,
#   ~2-3 events (steady-state recurrence ~240 yr).

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
echo "Cadence:       coseismic=1 s, nucleation=43200 s (12 hr), interseismic=1 yr"
echo "Frames:        ~1,500-2,000 projected at ~634 yr / 2-3 events"
echo "               (steady-state recurrence ~240 yr; no snapshot cap)"
echo "Queue:         normal (8N x 400r, 48h, tfinal 2e10 s ≈ 634 yr)"
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
      --paraview-dt-co     1.0 \
      --paraview-dt-nu     43200 \
      --paraview-dt-inter-yr 1.0 \
      --tfinal 20000000000 \
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
