#!/bin/bash
#SBATCH -J tpv104_p2_partfile
#SBATCH -o tpv104_p2_partfile_%j.out
#SBATCH -e tpv104_p2_partfile_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 P=2 production run — t = 2 s, dx = 200 m, ADER-DG p=2
# ============================================================
#
# Purpose
# -------
# Production verification of the precision-floor hypothesis for the
# TPV104 dip-direction pollution: at P=1 the per-side Pelties Riemann
# has loop gain > 1 in the V1 channel during rupture, amplifying 1-ULP
# FP-rounding seed to cm-scale slip_dip in ~1 s.  Higher polynomial
# order reduces the per-DOF FP-noise floor; local tests showed:
#   P=1 dx=500m: slip_dip = -1.008e-02 m (1 cm)   reference baseline
#   P=2 dx=500m: slip_dip = -5.80e-05 m (58 µm)   174× smaller
#   P=3 dx=500m: slip_dip ≈ 1e-4 m at t=1.81 s    (run stopped at 88%)
#
# This sbatch runs the same configuration as
# tpv104_symmirror_200m_partfile.sbatch but with --order 2 instead of
# --order 1.  All other settings preserved for apples-to-apples
# comparison.
#
# Why --order 2 (not 3)
# ---------------------
# P=2 doubles per-tet DOFs (4→10), CFL halves (factor 1/15 vs 1/9),
# and gives ~10× lower per-DOF noise floor.  Local result reduced
# slip_dip by 174× — well below mm-scale and likely below SCEC trace
# comparison resolution.  P=3 adds another 2× per-DOF reduction at
# ~5× compute cost; deferred to tpv104_symmirror_200m_p3_partfile.sbatch.
#
# Wall time / resources
# ---------------------
# P=2 increases per-step cost (~6× P=1) and step count (~1.7× P=1)
# vs the P=1 baseline (job 7677864 ran in ~2 h on 8N×400r).  Same
# 8N×400r / 2h development-queue allocation as baseline; if the run
# does not complete in 2h, the partial output up to wall-clock cut
# is still informative for the slip_dip-vs-time trajectory.  Resubmit
# with -p normal -t 12:00:00 if a complete tfinal=2.0s run is needed.
#
# Dispatch parity with P=1 baseline
# ---------------------------------
# Same: --ader-order 2, --bc-mode absorbing, --cfl 0.5, --tfinal 2.0,
# --paraview --paraview-dt 0.25 --paraview-bulk-dt 0.5, --verify-dispatch.
# Different: --order 1 → --order 2.
#
# Pre-requisites
# --------------
#   1. tpv104/mesh/tpv104_symmirror_200m.msh exists OR generated in-job.
#   2. tpv104_symmirror_200m.msh.np800.partition exists OR generated in-job
#      (note: partition is per-tet so partition file is order-independent;
#      the np400 file from baseline cannot be reused at np=800).
#   3. f_w = 0.2 (post-fix) in config/tpv104_params.hpp.

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

if ! grep -q '\-\-partition-file' drivers/tpv104_driver.cpp; then
   echo "ERROR: drivers/tpv104_driver.cpp does not carry --partition-file flag"
   exit 1
fi
if [ ! -f dynamic/fault_locality_partition.hpp ]; then
   echo "ERROR: dynamic/fault_locality_partition.hpp missing"
   exit 1
fi
echo "Source carries G1 partition support (commit $(git rev-parse --short HEAD))"

F_W_LINE=$(grep -E '^\s*static\s+constexpr\s+real_t\s+f_w' config/tpv104_params.hpp | head -1)
if ! echo "${F_W_LINE}" | grep -q '0\.2'; then
   echo "ERROR: f_w is NOT 0.2"; exit 1
fi
echo "f_w confirmed: ${F_W_LINE}"

MESH_PATH=tpv104/mesh/tpv104_symmirror_200m.msh
PART_PATH=tpv104/mesh/tpv104_symmirror_200m.msh.np${SLURM_NPROCS:-400}.partition

if [ ! -s "${MESH_PATH}" ]; then
   echo "Building symmetric 200 m mesh ..."
   cd tpv104/mesh
   python3 build_symmirror_mesh.py --dx 200 --lx 16000 --ly 4000 \
                                    --lz 16000 \
                                    --out tpv104_symmirror_200m.msh
   cd ../..
fi

if [ ! -s "${PART_PATH}" ]; then
   echo "Generating partition file for np=${SLURM_NPROCS:-400} ..."
   cd tpv104/mesh
   python3 build_symmirror_mesh.py --dx 200 --lx 16000 --ly 4000 \
                                    --lz 16000 \
                                    --out tpv104_symmirror_200m.msh \
                                    --emit-partition ${SLURM_NPROCS:-400}
   cd ../..
fi

if [ ! -s "${MESH_PATH}" ] || [ ! -s "${PART_PATH}" ]; then
   echo "ERROR: missing mesh or partition after generation attempt"
   ls -la "${MESH_PATH}" "${PART_PATH}" 2>&1
   exit 1
fi
echo "Mesh:      $(ls -lh ${MESH_PATH})"
echo "Partition: $(ls -lh ${PART_PATH})"

RESULT_DIR="tpv104/results_p2_partfile_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o dynamic/wave_operator.o
make seas_tpv104_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

unset SEAS_TPV104_FREEZE_SIGMA_N

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          tpv104_symmirror_200m.msh (768k tets, dx=200m)"
echo "Partition:     ${PART_PATH}"
echo "Polynomial:    --order 2 (DG p=2; ~6× cost vs P=1, 1.7× steps)"
echo "Target:        tfinal = 2.0 s at 8N×400r"
echo "Goal:          slip_dip ≪ 1 mm (precision-floor reduction)"
echo "Local ref:     P=2 dx=500m → slip_dip = 5.8e-5 m (vs 1e-2 P=1)"
echo "========================================"

ibrun ./seas_tpv104_driver \
      --mesh "${MESH_PATH}" \
      --partition-file "${PART_PATH}" \
      --mesh-scale 1 \
      --order 2 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 2.0 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law slip-srw \
      --fault-iterator substep \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv104_p2_partfile \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.25 \
      --paraview-bulk-dt 0.5 \
      --verify-dispatch

RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv104_p2_partfile_station_x2_0_x3_7.5.dat"

PEAK_V_HYPO="N/A"
PEAK_T_HYPO="N/A"
SIGMA_N_PERTURB_HYPO="N/A"
SLIP_DIP_END="N/A"
if [ -s "${HYPO_FILE}" ]; then
   read PEAK_V_HYPO PEAK_T_HYPO < <(awk '!/^#/ && NF>=3 {
      v=$3; if (v<0) v=-v;
      if (v>m) { m=v; tm=$1 }
   } END { printf "%g %g\n", m, tm }' "${HYPO_FILE}")
   SIGMA_N_PERTURB_HYPO=$(awk '!/^#/ && NF>=8 && $1+0==$1 {
      d=$8 - 120e6; if (d<0) d=-d;
      if (d>m) m=d
   } END { printf "%.4e", m+0 }' "${HYPO_FILE}")
   SLIP_DIP_END=$(awk '!/^#/ && NF>=5 && $1+0==$1 {a=$5} END{printf "%.4e", a+0}' "${HYPO_FILE}")
fi

declare -A SEISSOL_PEAKS=(
   ["x2_0_x3_7.5"]="7.39"
   ["x2_-12_x3_3"]="9.93"   ["x2_12_x3_3"]="9.52"
   ["x2_-12_x3_12"]="8.37"  ["x2_12_x3_12"]="8.28"
   ["x2_-9_x3_7.5"]="6.29"  ["x2_9_x3_7.5"]="6.33"
   ["x2_0_x3_3"]="5.14"     ["x2_0_x3_12"]="11.40"
)
OFF_HYPO_TABLE=""
for tag in x2_0_x3_7.5 x2_0_x3_3 x2_0_x3_12 \
           x2_-9_x3_7.5 x2_9_x3_7.5 \
           x2_-12_x3_3 x2_12_x3_3 \
           x2_-12_x3_12 x2_12_x3_12; do
   SF="${RESULT_DIR}/tpv104_p2_partfile_station_${tag}.dat"
   if [ -s "${SF}" ]; then
      MV=$(awk '!/^#/ && NF>=3 && $1+0==$1 {v=$3;if(v<0)v=-v;if(v>m)m=v} END{print m+0}' "${SF}")
      SV="${SEISSOL_PEAKS[$tag]}"
      RAT=$(awk "BEGIN{printf \"%.2f\", $MV/$SV}")
      SD=$(awk '!/^#/ && NF>=5 && $1+0==$1 {a=$5} END{printf "%.3e", a+0}' "${SF}")
      OFF_HYPO_TABLE+=$(printf "  %-18s  V=%6.3f (SS=%5.2f, ratio=%s)  slip_dip(end)=%s m\n" "$tag" "$MV" "$SV" "$RAT" "$SD")
      OFF_HYPO_TABLE+=$'\n'
   fi
done

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv104_p2_partfile_${SLURM_JOB_ID}.out"
DISPATCH_LINES=$(grep '\[dispatch\] rank=0' "${DISPATCH_LOG}" 2>/dev/null \
                 || echo "MISSING")

SIGMA_N_VERDICT="UNKNOWN"
if [ "${SIGMA_N_PERTURB_HYPO}" != "N/A" ]; then
   if awk "BEGIN { exit !(${SIGMA_N_PERTURB_HYPO} <= 1e3) }"; then
      SIGMA_N_VERDICT="PASS — σ_n flat at hypocenter"
   else
      SIGMA_N_VERDICT="FAIL — σ_n perturbation > 1 kPa"
   fi
fi

SLIP_DIP_VERDICT="UNKNOWN"
if [ "${SLIP_DIP_END}" != "N/A" ]; then
   ABS_SD=$(awk "BEGIN{a=${SLIP_DIP_END};if(a<0)a=-a;print a}")
   if awk "BEGIN { exit !(${ABS_SD} <= 1e-5) }"; then
      SLIP_DIP_VERDICT="PASS — slip_dip ≤ 10 µm (precision-floor fix verified)"
   elif awk "BEGIN { exit !(${ABS_SD} <= 1e-4) }"; then
      SLIP_DIP_VERDICT="PARTIAL — slip_dip ≤ 0.1 mm (>100× improvement)"
   elif awk "BEGIN { exit !(${ABS_SD} <= 1e-3) }"; then
      SLIP_DIP_VERDICT="WEAK — slip_dip ≤ 1 mm (10× improvement)"
   else
      SLIP_DIP_VERDICT="FAIL — slip_dip > 1 mm (precision-floor not enough)"
   fi
fi

{
   echo "TPV104 P=2 Cartesian-partition production run"
   echo "Mesh:                    tpv104_symmirror_200m (768k tets)"
   echo "Polynomial order:        --order 2 (DG p=2)"
   echo "Partition:               Cartesian X+Z (full Y per rank), np=${SLURM_NPROCS}"
   echo "tfinal:                  2.0 s"
   echo "f_w (compiled):          $(echo "${F_W_LINE}" | grep -oE '[01]\.[0-9]+' | head -1)"
   echo ""
   echo "peak |V_strike| at hypo  : ${PEAK_V_HYPO} m/s at t = ${PEAK_T_HYPO} s"
   echo "SeisSol reference (hypo) : 7.39 m/s at t = 1.07 s"
   echo ""
   echo "peak |σ_n − 120 MPa| at hypo : ${SIGMA_N_PERTURB_HYPO} Pa"
   echo "σ_n verdict                  : ${SIGMA_N_VERDICT}"
   echo ""
   echo "slip_dip at hypocenter (t=2s)  : ${SLIP_DIP_END} m"
   echo "  P=1 baseline (job 7677864) : ~3 mm (Cartesian partition, P=1)"
   echo "  P=1 local dx=500m          : -1.008e-02 m"
   echo "  P=2 local dx=500m          : -5.80e-05 m (174× reduction)"
   echo "slip_dip verdict             : ${SLIP_DIP_VERDICT}"
   echo ""
   echo "Per-station summary (V_strike + final slip_dip):"
   echo "${OFF_HYPO_TABLE}"
   echo "[dispatch] rank=0 lines  : ${DISPATCH_LINES}"
   echo ""
   if echo "${SIGMA_N_VERDICT}" | grep -q "^PASS" \
      && echo "${SLIP_DIP_VERDICT}" | grep -qE "^PASS|^PARTIAL" \
      && [ "${PEAK_V_HYPO}" != "N/A" ]; then
      echo "STATUS: PASS"
   else
      echo "STATUS: FAIL"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
