#!/bin/bash
#SBATCH -J bp5_suf500
#SBATCH -o bp5_system_update_500yr_%j.out
#SBATCH -e bp5_system_update_500yr_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 24:00:00
#SBATCH -A EAR20006

# BP5 system_update production run, 500 yr, slip_rate_strike-only VTU output.
#
# Differences from bp5_system_update_full_production.sbatch:
#   1. --tfinal 1.57788e10 s  (== 500 yr at 365.25 d/yr).
#   2. Wall-time 24 h (vs 48 h) — 500 yr is ~2 events + setup ≈ 500/1800 of
#      the full production, plus startup overhead.
#   3. Post-run VTU trimming: strip every fault-surface VTU down to the
#      slip_rate_strike field only.  Reduces per-file size ~12x but
#      preserves the inode count (one VTU per rank per cycle).
#   4. Optional per-cycle tar bundling (commented out by default;
#      uncomment if inode pressure returns).
#
# Before submitting, on the login node:
#   cd /scratch2/10024/zhaochun/seas-project/seas-mfem
#   git fetch origin && git checkout system_update && git pull
#   cd miniapps/seas && sbatch jobs/bp5/bp5_system_update_500yr_strike_only.sbatch

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${LD_LIBRARY_PATH}"

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
make clean
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="bp5/results_system_update_500yr_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_system_update_500yr_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "500-yr BP5 production: 400 ranks, 24hr, strike-slip-rate-only VTU"

# --tfinal is in SECONDS per bp5_verification_full.cpp:29,766.
# 500 yr at 365.25 days/yr = 500 * 31557600 s = 1.57788e10 s.
#
# ParaView cadence inherited from asymmetric production:
#   dt_co      = 5.0 s   -> ~20 frames per ~100-s rupture
#   dt_nu      = 100.0 s -> ~1 frame per nucleation phase
#   dt_inter   = 50 yr   -> 10 interseismic frames over 500 yr
# Expected frames: 10 (inter) + ~3 events * 20 (co) + ~3 (nu) ~= 73.
# Expected files:  73 * 401 (400 ranks + 1 PVTU) + 1 PVD ~= 29k.
#
# --paraview-fields emits ONLY the named CellData arrays in each VTU.
# Empty / omitted flag ⇒ all 12 standard fields (backward-compat).
# Here we keep slip_rate_strike only ⇒ ~12x per-file size reduction.
# Inode count unchanged by this flag — see the tar bundling block
# below if inode pressure is a concern.
ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --tfinal 1.57788e10 \
      --ref-dir bp5/benchmark_data/golden_parallel_400r_p1_1000m \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --checkpoint-interval 5000 \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --regression-tolerance 1e-2 \
      --verify \
      --paraview \
      --paraview-fault-only \
      --paraview-adaptive \
      --paraview-v-co        1e-3 \
      --paraview-v-nu        1e-6 \
      --paraview-hyst        10.0 \
      --paraview-dt-co       5.0 \
      --paraview-dt-nu       100.0 \
      --paraview-dt-inter-yr 50.0 \
      --paraview-fields      slip_rate_strike

# ---------------------------------------------------------------
# OPTIONAL: per-cycle tar bundle to reduce inode count ~400x.
# Uncomment if `lfs quota -u $USER $SCRATCH` shows inode pressure.
# PVD + PVTU files stay readable; only the per-rank VTUs are bundled.
# To use a bundled file later, extract the matching tar first:
#   tar -xf fault_surface_c<N>.tar
# ---------------------------------------------------------------
# FAULT_DIR="${RESULT_DIR}/ParaView/FaultSurface"
# if [ -d "${FAULT_DIR}" ]; then
#    echo "[POST] Bundling per-rank VTUs into per-cycle tar archives..."
#    cd "${FAULT_DIR}"
#    for pvtu in fault_surface_c*.pvtu; do
#       cyc="${pvtu#fault_surface_c}"; cyc="${cyc%.pvtu}"
#       tar -cf "fault_surface_c${cyc}.tar" fault_surface_r*_c${cyc}.vtu \
#          && rm fault_surface_r*_c${cyc}.vtu
#    done
#    cd -
#    echo "[POST] Bundling done"
# fi

echo "Done."
