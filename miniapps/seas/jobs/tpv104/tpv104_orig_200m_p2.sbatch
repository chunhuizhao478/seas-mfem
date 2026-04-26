#!/bin/bash
#SBATCH -J tpv104_orig_p2
#SBATCH -o tpv104_orig_p2_%j.out
#SBATCH -e tpv104_orig_p2_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 P=2 production run on ORIGINAL Gmsh tet mesh
# ===================================================
#
# Purpose
# -------
# Test whether the precision-floor hypothesis (P=2 reduces slip_dip
# leak by ~10×) holds on the ORIGINAL non-symmetric Gmsh-built tet
# mesh `tpv104_200m.msh` (NOT the symmirror mesh).
#
# Why this run matters
# --------------------
# The symmirror mesh was built specifically to eliminate two coupled
# leaks observed earlier:
#   (1) σ_n perturbation across the fault (mesh-driven non-mirror at
#       sub-cell level)
#   (2) slip_dip exponential growth from FP-noise (loop-gain instability
#       in Pelties V1 channel; precision-floor reduction via higher P
#       reduces the seed)
#
# With the symmirror mesh, leak (1) is gone (mesh-side fix) and leak
# (2) is amplified by the loop gain regardless of mesh.  Local tests
# at P=2 dx=500m on symmirror reduced slip_dip by 174×.
#
# THIS run uses the original Gmsh tet mesh (asymmetric) at P=2.  Two
# possible outcomes:
#   A) slip_dip ≈ symmirror P=2 result (~µm scale) → leak (1) was the
#      mesh-symmetry issue ONLY; P=2 stabilizes (2) regardless of mesh.
#      Symmirror mesh is no longer needed for TPV104 at P≥2.
#   B) slip_dip is much larger than symmirror P=2 (mm-cm scale, sigma_n
#      perturbation > 0.5 MPa) → mesh asymmetry produces a separate
#      mesh-driven seed at the fault that the loop also amplifies, and
#      P=2's precision-floor reduction is insufficient.  Symmirror mesh
#      remains required.
#
# Either outcome is informative: result A would simplify the production
# path (drop the symmirror mesh requirement); result B confirms the
# mesh-symmetry constraint persists at higher polynomial order.
#
# Mesh
# ----
# tpv104/mesh/tpv104_200m.msh — structured Gmsh (built-in kernel),
# 200 m fault resolution, NOT mirror-symmetric.  ~215k tets vs 768k
# in the symmirror mesh.  No partition file (uses default ParMETIS).
#
# Wall time / resources
# ---------------------
# Original mesh has ~3.6× FEWER tets than symmirror (215k vs 768k),
# so per-step cost is correspondingly smaller.  P=2 on this mesh
# should fit comfortably in the 2h dev queue allocation.
#
# Dispatch parity with previous P=1 baselines
# -------------------------------------------
# Same: --ader-order 2, --bc-mode absorbing, --cfl 0.5, --tfinal 2.0.
# Different: --order 1 → --order 2.  Different: no --partition-file
# (the original mesh hasn't had a sidecar partition built; ParMETIS
# default is used).

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

F_W_LINE=$(grep -E '^\s*static\s+constexpr\s+real_t\s+f_w' config/tpv104_params.hpp | head -1)
if ! echo "${F_W_LINE}" | grep -q '0\.2'; then
   echo "ERROR: f_w is NOT 0.2"; exit 1
fi
echo "f_w confirmed: ${F_W_LINE}"

MESH_PATH=tpv104/mesh/tpv104_200m.msh

if [ ! -s "${MESH_PATH}" ]; then
   echo "Building original Gmsh tet mesh from .geo ..."
   cd tpv104/mesh
   gmsh -3 tpv104_200m.geo -o tpv104_200m.msh -format msh22
   cd ../..
fi

if [ ! -s "${MESH_PATH}" ]; then
   echo "ERROR: missing mesh after generation attempt"
   ls -la "${MESH_PATH}" 2>&1
   exit 1
fi
echo "Mesh: $(ls -lh ${MESH_PATH})"

RESULT_DIR="tpv104/results_orig_p2_job${SLURM_JOB_ID}"
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
echo "Mesh:          tpv104_200m.msh (original Gmsh, NOT mirror-symmetric)"
echo "Partition:     default ParMETIS (no sidecar partition file)"
echo "Polynomial:    --order 2 (DG p=2)"
echo "Target:        tfinal = 2.0 s at 8N×400r"
echo "Goal:          quantify slip_dip + sigma_n leak on non-symmetric mesh"
echo "Reference:     symmirror P=2 local dx=500m → slip_dip = 5.8e-5 m"
echo "========================================"

ibrun ./seas_tpv104_driver \
      --mesh "${MESH_PATH}" \
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
      --output-prefix tpv104_orig_p2 \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.25 \
      --paraview-bulk-dt 0.5 \
      --verify-dispatch

RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv104_orig_p2_station_x2_0_x3_7.5.dat"

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
   SF="${RESULT_DIR}/tpv104_orig_p2_station_${tag}.dat"
   if [ -s "${SF}" ]; then
      MV=$(awk '!/^#/ && NF>=3 && $1+0==$1 {v=$3;if(v<0)v=-v;if(v>m)m=v} END{print m+0}' "${SF}")
      SV="${SEISSOL_PEAKS[$tag]}"
      RAT=$(awk "BEGIN{printf \"%.2f\", $MV/$SV}")
      SD=$(awk '!/^#/ && NF>=5 && $1+0==$1 {a=$5} END{printf "%.3e", a+0}' "${SF}")
      OFF_HYPO_TABLE+=$(printf "  %-18s  V=%6.3f (SS=%5.2f, ratio=%s)  slip_dip(end)=%s m\n" "$tag" "$MV" "$SV" "$RAT" "$SD")
      OFF_HYPO_TABLE+=$'\n'
   fi
done

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv104_orig_p2_${SLURM_JOB_ID}.out"
DISPATCH_LINES=$(grep '\[dispatch\] rank=0' "${DISPATCH_LOG}" 2>/dev/null \
                 || echo "MISSING")

# Verdict bands tuned for the original-mesh expectation: pre-fix
# (P=1 on this mesh) showed σ_n perturbation ~0.5 MPa and slip_dip
# ~1 cm.  P=2 should cut both substantially if the mesh-symmetry
# was orthogonal to the loop-gain mechanism.
SIGMA_N_VERDICT="UNKNOWN"
if [ "${SIGMA_N_PERTURB_HYPO}" != "N/A" ]; then
   if awk "BEGIN { exit !(${SIGMA_N_PERTURB_HYPO} <= 1e3) }"; then
      SIGMA_N_VERDICT="PASS — σ_n flat at hypocenter (P=2 mesh-symmetry-free)"
   elif awk "BEGIN { exit !(${SIGMA_N_PERTURB_HYPO} <= 1e5) }"; then
      SIGMA_N_VERDICT="PARTIAL — σ_n perturbation ≤ 100 kPa"
   else
      SIGMA_N_VERDICT="FAIL — σ_n perturbation > 100 kPa (mesh symmetry still drives leak)"
   fi
fi

SLIP_DIP_VERDICT="UNKNOWN"
if [ "${SLIP_DIP_END}" != "N/A" ]; then
   ABS_SD=$(awk "BEGIN{a=${SLIP_DIP_END};if(a<0)a=-a;print a}")
   if awk "BEGIN { exit !(${ABS_SD} <= 1e-4) }"; then
      SLIP_DIP_VERDICT="PASS — slip_dip ≤ 0.1 mm (P=2 stabilizes regardless of mesh)"
   elif awk "BEGIN { exit !(${ABS_SD} <= 1e-3) }"; then
      SLIP_DIP_VERDICT="PARTIAL — slip_dip ≤ 1 mm"
   else
      SLIP_DIP_VERDICT="FAIL — slip_dip > 1 mm (mesh asymmetry persists at P=2)"
   fi
fi

{
   echo "TPV104 P=2 ORIGINAL-mesh production run"
   echo "Mesh:                    tpv104_200m (original Gmsh, NON-mirror)"
   echo "Polynomial order:        --order 2 (DG p=2)"
   echo "Partition:               default ParMETIS"
   echo "tfinal:                  2.0 s"
   echo "f_w (compiled):          $(echo "${F_W_LINE}" | grep -oE '[01]\.[0-9]+' | head -1)"
   echo ""
   echo "peak |V_strike| at hypo  : ${PEAK_V_HYPO} m/s at t = ${PEAK_T_HYPO} s"
   echo "SeisSol reference (hypo) : 7.39 m/s at t = 1.07 s"
   echo ""
   echo "peak |σ_n − 120 MPa| at hypo : ${SIGMA_N_PERTURB_HYPO} Pa"
   echo "  P=1 reference on this mesh : ~0.46 MPa (mesh-driven, pre-symmirror)"
   echo "σ_n verdict                  : ${SIGMA_N_VERDICT}"
   echo ""
   echo "slip_dip at hypocenter (t=2s)  : ${SLIP_DIP_END} m"
   echo "  P=1 reference on this mesh   : ~1 cm (loop gain × mesh seed)"
   echo "  P=2 symmirror local dx=500m  : -5.80e-05 m"
   echo "slip_dip verdict             : ${SLIP_DIP_VERDICT}"
   echo ""
   echo "Per-station summary (V_strike + final slip_dip):"
   echo "${OFF_HYPO_TABLE}"
   echo "[dispatch] rank=0 lines  : ${DISPATCH_LINES}"
   echo ""
   if echo "${SIGMA_N_VERDICT}" | grep -qE "^PASS|^PARTIAL" \
      && echo "${SLIP_DIP_VERDICT}" | grep -qE "^PASS|^PARTIAL" \
      && [ "${PEAK_V_HYPO}" != "N/A" ]; then
      echo "STATUS: PASS"
   else
      echo "STATUS: FAIL"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
