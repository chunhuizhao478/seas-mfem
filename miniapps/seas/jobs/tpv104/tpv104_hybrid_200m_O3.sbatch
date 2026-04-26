#!/bin/bash
#SBATCH -J tpv104_hybrid_O3
#SBATCH -o tpv104_hybrid_O3_%j.out
#SBATCH -e tpv104_hybrid_O3_%j.err
#SBATCH -p normal
#SBATCH -N 16
#SBATCH -n 800
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# TPV104 𝒪3 (P=2 + ADER-O3) production run on HYBRID symmetric-strip mesh
# =======================================================================
#
# Purpose
# -------
# First production run on the hybrid mesh that combines:
#   - a structured y-mirror prism strip near the fault
#     (5 layers of 200 m per side, ±1 km in y, Pelties Fig 4(a) Option A)
#   - graded unstructured tets in the outer 60 km half-space
#     (200 m at strip boundary -> 5000 m at far field, matching the
#      original tpv104_200m.geo Distance + MathEval + Threshold scheme).
#
# The hybrid mesh resolves three earlier concerns at once:
#   (1) y-mirror near fault (eliminates the partition-driven sigma_n
#       perturbation seed).
#   (2) larger outer domain than the symmirror mesh (lx=lz=60 km vs 16 km),
#       so absorbing-BC reflections do not contaminate the trace at
#       hypocenter station before t = 2 s.
#   (3) graded outer mesh (saves ~6× tets vs uniform 200 m) -> 3.4 M
#       tets total, vs an estimated ~20 M tets for a uniform 60 km dx=200
#       mesh.
#
# Why this run matters
# --------------------
# Combines the precision-floor reduction (P=2 + 𝒪3 ADER) with the only
# mesh that simultaneously satisfies near-fault y-mirror AND a large
# outer domain.  Earlier symmirror runs at lx=ly=4 km hit reflected
# S waves at t ~ 1.15 s; the hybrid puts reflections beyond t = 2 s.
#
# Mesh
# ----
# tpv104/mesh/tpv104_200m_hybrid.msh — built from
#   tpv104/mesh/tpv104_200m_hybrid.geo
# via Gmsh OCC (Extrude+Layers strip + BooleanFragments outer).  ~533 k
# nodes, ~3.24 M tets.  No sidecar partition file — uses default
# ParMETIS (acceptable: the structured strip preserves y-mirror at the
# element level regardless of MPI partition; partition-driven asymmetry
# would only enter if the partition cuts the strip face mesh, which
# ParMETIS avoids in practice for thin strips).
#
# Wall time / resources
# ---------------------
# Hybrid is 4.4× larger than symmirror (3.24 M vs 768 k tets).  Symmirror
# 𝒪3 at 8N×400r needed ~3 h estimated.  Hybrid at 16N×800r should fit
# in ~6-8 h; 48 h normal-queue allocation gives generous headroom.
#
# Pelties order
# -------------
# --order 2 (DG p=2) + --ader-order 3 = canonical Pelties 𝒪3.
# Local results (symmirror dx=500m): P=1 gave slip_dip = -1.0e-2 m;
# P=2 gave -5.8e-5 m (174× reduction).  Hybrid 𝒪3 expected to reproduce
# the symmirror P=2 reduction without the small-domain reflection
# contamination.

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

MESH_PATH=tpv104/mesh/tpv104_200m_hybrid.msh
GEO_PATH=tpv104/mesh/tpv104_200m_hybrid.geo

if [ ! -s "${MESH_PATH}" ]; then
   echo "Building hybrid mesh from .geo (Gmsh OCC, Extrude+Layers strip + outer) ..."
   if [ ! -s "${GEO_PATH}" ]; then
      echo "ERROR: ${GEO_PATH} missing — cannot build mesh"
      exit 1
   fi
   cd tpv104/mesh
   gmsh -3 tpv104_200m_hybrid.geo -o tpv104_200m_hybrid.msh -format msh22
   cd ../..
fi

if [ ! -s "${MESH_PATH}" ]; then
   echo "ERROR: missing mesh after generation attempt"
   ls -la "${MESH_PATH}" 2>&1
   exit 1
fi
echo "Mesh: $(ls -lh ${MESH_PATH})"

RESULT_DIR="tpv104/results_hybrid_O3_job${SLURM_JOB_ID}"
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
echo "Mesh:          tpv104_200m_hybrid.msh (533 k nodes, 3.24 M tets)"
echo "               y-mirror prism strip + graded unstructured outer"
echo "Partition:     default ParMETIS (no sidecar partition file)"
echo "Pelties order: 𝒪3 (--order 2 + --ader-order 3)"
echo "Iterator:      --fault-iterator substep (ADER sub-step time integration)"
echo "Target:        tfinal = 2.0 s at 16N×800r"
echo "Goal:          slip_dip ≪ 1 mm with no reflection contamination"
echo "Reference:     symmirror P=2 dx=500m (small domain) -> -5.80e-5 m"
echo "========================================"

ibrun ./seas_tpv104_driver \
      --mesh "${MESH_PATH}" \
      --mesh-scale 1 \
      --order 2 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 2.0 \
      --ader-order 3 \
      --friction-solver newton-stable \
      --fric-law slip-srw \
      --fault-iterator substep \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv104_hybrid_O3 \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.5 \
      --paraview-bulk-dt 0.5 \
      --verify-dispatch

RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv104_hybrid_O3_station_x2_0_x3_7.5.dat"

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
   SF="${RESULT_DIR}/tpv104_hybrid_O3_station_${tag}.dat"
   if [ -s "${SF}" ]; then
      MV=$(awk '!/^#/ && NF>=3 && $1+0==$1 {v=$3;if(v<0)v=-v;if(v>m)m=v} END{print m+0}' "${SF}")
      SV="${SEISSOL_PEAKS[$tag]}"
      RAT=$(awk "BEGIN{printf \"%.2f\", $MV/$SV}")
      SD=$(awk '!/^#/ && NF>=5 && $1+0==$1 {a=$5} END{printf "%.3e", a+0}' "${SF}")
      OFF_HYPO_TABLE+=$(printf "  %-18s  V=%6.3f (SS=%5.2f, ratio=%s)  slip_dip(end)=%s m\n" "$tag" "$MV" "$SV" "$RAT" "$SD")
      OFF_HYPO_TABLE+=$'\n'
   fi
done

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv104_hybrid_O3_${SLURM_JOB_ID}.out"
DISPATCH_LINES=$(grep '\[dispatch\] rank=0' "${DISPATCH_LOG}" 2>/dev/null \
                 || echo "MISSING")

SIGMA_N_VERDICT="UNKNOWN"
if [ "${SIGMA_N_PERTURB_HYPO}" != "N/A" ]; then
   if awk "BEGIN { exit !(${SIGMA_N_PERTURB_HYPO} <= 1e3) }"; then
      SIGMA_N_VERDICT="PASS — σ_n flat at hypocenter (< 1 kPa perturbation)"
   elif awk "BEGIN { exit !(${SIGMA_N_PERTURB_HYPO} <= 1e5) }"; then
      SIGMA_N_VERDICT="PARTIAL — σ_n perturbation ≤ 100 kPa"
   else
      SIGMA_N_VERDICT="FAIL — σ_n perturbation > 100 kPa"
   fi
fi

SLIP_DIP_VERDICT="UNKNOWN"
if [ "${SLIP_DIP_END}" != "N/A" ]; then
   ABS_SD=$(awk "BEGIN{a=${SLIP_DIP_END};if(a<0)a=-a;print a}")
   if awk "BEGIN { exit !(${ABS_SD} <= 1e-5) }"; then
      SLIP_DIP_VERDICT="PASS — slip_dip ≤ 10 µm (precision-floor + symmetric strip)"
   elif awk "BEGIN { exit !(${ABS_SD} <= 1e-4) }"; then
      SLIP_DIP_VERDICT="PARTIAL — slip_dip ≤ 0.1 mm (>100× improvement)"
   elif awk "BEGIN { exit !(${ABS_SD} <= 1e-3) }"; then
      SLIP_DIP_VERDICT="WEAK — slip_dip ≤ 1 mm (10× improvement)"
   else
      SLIP_DIP_VERDICT="FAIL — slip_dip > 1 mm"
   fi
fi

{
   echo "TPV104 𝒪3 (P=2 + ADER-O3) HYBRID-mesh production run"
   echo "Mesh:                    tpv104_200m_hybrid (3.24 M tets, structured strip + graded outer)"
   echo "Polynomial order:        --order 2 + --ader-order 3 (canonical Pelties 𝒪3)"
   echo "Iterator:                substep (ADER sub-step time integration)"
   echo "Partition:               default ParMETIS"
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
   echo "  symmirror P=2 dx=500m local  : -5.80e-05 m (174× reduction vs P=1)"
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
