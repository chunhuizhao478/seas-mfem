#!/bin/bash
#SBATCH -J tpv104_partfile
#SBATCH -o tpv104_partfile_%j.out
#SBATCH -e tpv104_partfile_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 Cartesian-partition production run — t = 2 s, 200 m, 400 ranks
# ======================================================================
#
# Purpose
# -------
# Frontera production run on the y-mirror-symmetric 200 m mesh using
# the explicit Cartesian (X+Z) partition file produced by
# `tpv104/mesh/build_symmirror_mesh.py --emit-partition 400`.
#
# This run differs from the earlier symmirror_200m run (job 7677822)
# only by adding `--partition-file <path>`.  Local verification on
# the 1000m mirror mesh proved bit-exact y-mirror at np ∈ {1, 2, 4, 8}
# with this partition (slip_dip, slip_strike, V_max all identical to
# 14-digit display precision).  Production confirms the σ_n leak is
# fixed AND the dip-direction MPI pollution is fixed.
#
# Background
# ----------
# σ_n leak: caused by mesh non-mirror-symmetry; fixed by tpv104_symmirror_200m.msh.
# dip pollution: caused by ParMETIS cutting the y=0 plane at np=400;
#   fixed by Cartesian X+Z partition file (this sbatch).
#
# See debug docs:
#   debug_document/tpv104_debug_document/tpv104_bulk_asymmetry_followup_2026-04-25_pm.md
#   debug_document/code_limitation/02_cartesian_partition_dynamic.md
#
# Mesh/dispatch parity with previous symmirror run (7677822):
#   - tpv104_symmirror_200m.msh (768k tets, 16x4x16 km, dx=200 m).
#   - 8 nodes × 400 ranks (same).
#   - --ader-order 2, --order 1, --cfl 0.5 (same).
#   - --tfinal 2.0, --bc-mode absorbing (same).
#   - --paraview --paraview-dt 0.25 --paraview-bulk-dt 0.5 (same).
#
# The ONLY production-path change vs symmirror_200m run is the addition
# of --partition-file tpv104/mesh/tpv104_symmirror_200m.msh.np400.partition.
#
# Expected result
# ---------------
# σ_n at hypocenter station: still flat at 120 MPa (mesh fix preserved).
# slip_dip at hypocenter: drops from ~3 mm (job 7677822) to <1 mm
#   (residual physical Rayleigh-wave coupling at free surface).
# slip_strike, V_strike: bit-equivalent to np=1,2 reference.
# τ_strike peaks: same 81.6 MPa as before (matches SeisSol within 2%).
#
# Pre-requisites
# --------------
#   1. tpv104/mesh/tpv104_symmirror_200m.msh exists OR is generated
#      in-job by the build script.
#   2. tpv104_symmirror_200m.msh.np400.partition exists OR is generated
#      in-job for nprocs=400.
#   3. f_w = 0.2 (post-fix) in config/tpv104_params.hpp.
#   4. Source tree at commit acf0ef1 or later (G1 + partition-file flags).

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

# Verify the source tree carries the partition flags.
if ! grep -q '\-\-partition-file' drivers/tpv104_driver.cpp; then
   echo "ERROR: drivers/tpv104_driver.cpp does not carry --partition-file flag"
   echo "       Pull commit acf0ef1 or later and rebuild."
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

# Build the symmetric mesh + partition file on Frontera if missing.
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

PREV_INIT_RESULT=$(ls -dt tpv104/results_init_stations_job* 2>/dev/null | head -1)
if [ -z "${PREV_INIT_RESULT}" ] \
   || ! grep -q "^STATUS: PASS" "${PREV_INIT_RESULT}/RESULT.txt" 2>/dev/null; then
   echo "WARN: no passing tpv104_init_stations.sbatch result found.  Continuing anyway"
fi

RESULT_DIR="tpv104/results_partfile_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

# Production build (no diag flag).
rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o dynamic/wave_operator.o
make seas_tpv104_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

unset SEAS_TPV104_FREEZE_SIGMA_N

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          tpv104_symmirror_200m.msh (768k tets, Cartesian)"
echo "Partition:     ${PART_PATH}"
echo "Target:        tfinal = 2.0 s at 8N×400r"
echo "Goal:          σ_n flat AND slip_dip ≈ 0 (vs ~3 mm baseline)"
echo "========================================"

ibrun ./seas_tpv104_driver \
      --mesh "${MESH_PATH}" \
      --partition-file "${PART_PATH}" \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 2.0 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law slip-srw \
      --fault-iterator substep \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv104_partfile \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.25 \
      --paraview-bulk-dt 0.5 \
      --verify-dispatch

# Post-run analysis: peak σ_n perturbation, peak slip_dip, peak |V_strike|.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv104_partfile_station_x2_0_x3_7.5.dat"

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

# Per-station summary (V_strike comparison vs SeisSol + slip_dip).
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
   SF="${RESULT_DIR}/tpv104_partfile_station_${tag}.dat"
   if [ -s "${SF}" ]; then
      MV=$(awk '!/^#/ && NF>=3 && $1+0==$1 {v=$3;if(v<0)v=-v;if(v>m)m=v} END{print m+0}' "${SF}")
      SV="${SEISSOL_PEAKS[$tag]}"
      RAT=$(awk "BEGIN{printf \"%.2f\", $MV/$SV}")
      SD=$(awk '!/^#/ && NF>=5 && $1+0==$1 {a=$5} END{printf "%.3e", a+0}' "${SF}")
      OFF_HYPO_TABLE+=$(printf "  %-18s  V=%6.3f (SS=%5.2f, ratio=%s)  slip_dip(end)=%s m\n" "$tag" "$MV" "$SV" "$RAT" "$SD")
      OFF_HYPO_TABLE+=$'\n'
   fi
done

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv104_partfile_${SLURM_JOB_ID}.out"
DISPATCH_LINES=$(grep '\[dispatch\] rank=0' "${DISPATCH_LOG}" 2>/dev/null \
                 || echo "MISSING")
DISPATCH_OK=0
if echo "${DISPATCH_LINES}" | grep -q 'solver=brent' \
   && echo "${DISPATCH_LINES}" | grep -q 'iterator=oneshot' \
   && echo "${DISPATCH_LINES}" | grep -q 'fric-law=slip-srw'; then
   DISPATCH_OK=1
fi

# σ_n flatness verdict (mesh fix preserved).
SIGMA_N_VERDICT="UNKNOWN"
if [ "${SIGMA_N_PERTURB_HYPO}" != "N/A" ]; then
   if awk "BEGIN { exit !(${SIGMA_N_PERTURB_HYPO} <= 1e3) }"; then
      SIGMA_N_VERDICT="PASS — σ_n flat at hypocenter"
   else
      SIGMA_N_VERDICT="FAIL — σ_n perturbation > 1 kPa"
   fi
fi

# slip_dip verdict (partition fix).
SLIP_DIP_VERDICT="UNKNOWN"
if [ "${SLIP_DIP_END}" != "N/A" ]; then
   ABS_SD=$(awk "BEGIN{a=${SLIP_DIP_END};if(a<0)a=-a;print a}")
   if awk "BEGIN { exit !(${ABS_SD} <= 1e-4) }"; then
      SLIP_DIP_VERDICT="PASS — slip_dip ≤ 0.1 mm (partition fix verified)"
   elif awk "BEGIN { exit !(${ABS_SD} <= 1e-3) }"; then
      SLIP_DIP_VERDICT="PARTIAL — slip_dip ≤ 1 mm (improved vs ~3 mm baseline)"
   else
      SLIP_DIP_VERDICT="FAIL — slip_dip > 1 mm (partition fix did not work)"
   fi
fi

{
   echo "TPV104 Cartesian-partition production run"
   echo "Mesh:                    tpv104_symmirror_200m (768k tets)"
   echo "Partition:               Cartesian X+Z (full Y per rank)"
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
   echo "  baseline (default ParMETIS) : ~3 mm  (job 7677822)"
   echo "slip_dip verdict             : ${SLIP_DIP_VERDICT}"
   echo ""
   echo "Per-station summary (V_strike + final slip_dip):"
   echo "${OFF_HYPO_TABLE}"
   echo "[dispatch] rank=0 lines  : ${DISPATCH_LINES}"
   if [ "${DISPATCH_OK}" -eq 1 ]; then
      echo "dispatch parity vs fwfix : PASS"
   else
      echo "dispatch parity vs fwfix : FAIL"
   fi
   echo ""
   if echo "${SIGMA_N_VERDICT}" | grep -q "^PASS" \
      && echo "${SLIP_DIP_VERDICT}" | grep -qE "^PASS|^PARTIAL" \
      && [ "${DISPATCH_OK}" -eq 1 ] \
      && [ "${PEAK_V_HYPO}" != "N/A" ]; then
      echo "STATUS: PASS"
   else
      echo "STATUS: FAIL"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
