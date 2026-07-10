#!/bin/bash
#SBATCH -J tpv205uv_p1
#SBATCH -o tpv205uv_p1_%j.out
#SBATCH -e tpv205uv_p1_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=12:00:00
#SBATCH --mem=240G
#SBATCH -A ddp408

# =========================================================================
# TPV205 unified-substep VERIFICATION + CALIPER SPEED, P=1 (ADER-O2)
# *** EXPANSE (SDSC) job — 1N x 128r / 12:00:00 ***
#
# Unify-plan Phase 6 steps 1+3 (PLAN_unify_interior_shared_fault_substep_
# 2026-07-09.md) on the PRODUCTION TPV205 mesh, branch
# system/substep-unified-tpv26-tpv27:
#
#   Phase A (VERIFY):  np sweep {TPV_NP_SWEEP} of the unified
#     `--fault-iterator substep` path; station traces must agree across np
#     to print precision (the shared-fault path is bit-neutral per the
#     local symmirror evidence; the station tie-break fix makes the
#     readout deterministic — np4_attractor_root_cause_2026-07-10.md).
#   Phase B (COMPARE): jobs/unify_phase6_expanse/compare_np_sweep.py
#     prints max |X_np - X_ref| per station into the log + a RESULT file.
#   Phase C (CALIPER): one profiled run (runtime-report + spot) at
#     {TPV_CALI_RANKS} ranks for the speed-vs-polynomial-order graph
#     (plot with jobs/unify_phase6_expanse/plot_caliper_speed.py after the
#     P1/P2/P3 jobs finish).  Skipped (with a warning) if the binary is
#     not Caliper-linked.
#
# BUILD (once, on Expanse):
#   USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh
#   cd miniapps/seas && make seas_tpv205_driver -j 16
# MESH (gitignored; generate with gmsh v2.2 ASCII or scp):
#   gmsh tpv205/mesh/tpv2053d_200m.geo -3 -format msh22 -o tpv205/mesh/tpv2053d_200m.msh
#
# Env knobs: TPV_NP_SWEEP (default "16 32 64 128"; use "1 2 4 8" only on a
# coarser TPV_MESH — np<=8 on the 200 m mesh will not fit the wall time),
# TPV_MESH, TPV_TFINAL_VERIFY (default 0.5 s), TPV_TFINAL_CALI (default
# 1.2 s), TPV_CALI_RANKS (default 128), TPV_CFL (default 0.5),
# TPV_SKIP_VERIFY / TPV_SKIP_CALI (=1 to skip a phase), TPV_CALI_CONFIG
# (full CALI_CONFIG override), TPV_OUTDIR.
# =========================================================================

export LC_ALL=C
export LANG=C

# --- Expanse modules (MUST match build_expanse.sh: cpu/0.17.3b gcc openmpi) ---
module purge
module load slurm 2>/dev/null || true
module load cpu/0.17.3b
module load gcc/10.2.0
for _mpi in openmpi/4.1.5 openmpi/4.1.3 openmpi/4.1.1 openmpi; do
    if module load "${_mpi}" 2>/dev/null && module is-loaded "${_mpi}" 2>/dev/null; then
        echo "OpenMPI module: ${_mpi}"; break
    fi
done

# Resolve seas-mfem root (walk UP to miniapps/seas/Makefile).
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR}}"
_walk_up_dir="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5; do
    if [[ -f "${_walk_up_dir}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"; break
    fi
    _walk_up_dir="${_walk_up_dir}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}."
    echo "  Set SEAS_MFEM_ROOT via sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi

# Runtime shared libs (libcaliper first; harmless when not Caliper-linked).
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"

cd "${SEAS_MFEM_ROOT}/miniapps/seas"
if [[ ! -x ./seas_tpv205_driver ]]; then
    echo "ERROR: ./seas_tpv205_driver not found.  Build it FIRST on Expanse:"
    echo "  USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh && cd miniapps/seas && make seas_tpv205_driver -j 16"
    exit 1
fi

MESH="${TPV_MESH:-tpv205/mesh/tpv2053d_200m.msh}"
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: mesh not found: $(pwd)/${MESH}  (gitignored)."
    echo "  Generate:  gmsh tpv205/mesh/tpv2053d_200m.geo -3 -format msh22 -o tpv205/mesh/tpv2053d_200m.msh"
    exit 1
fi

NP_SWEEP="${TPV_NP_SWEEP:-16 32 64 128}"
TF_VERIFY="${TPV_TFINAL_VERIFY:-0.5}"
TF_CALI="${TPV_TFINAL_CALI:-1.2}"
CALI_RANKS="${TPV_CALI_RANKS:-128}"
CFL="${TPV_CFL:-0.5}"
OUT_BASE="${TPV_OUTDIR:-tpv205/results/unify_phase6_tpv205_p1_${SLURM_JOB_ID}}"
mkdir -p "${OUT_BASE}"

echo "Git commit:  $(git -C "${SEAS_MFEM_ROOT}" rev-parse HEAD)"
echo "Branch:      $(git -C "${SEAS_MFEM_ROOT}" branch --show-current)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Problem:     TPV205  P=1  ADER-O2  fault-iterator=substep  LSW (newton-stable, mixed-flux none)"
echo "Mesh:        ${MESH}"
echo "Phase A np:  ${NP_SWEEP}  tfinal=${TF_VERIFY} s   Phase C: ${CALI_RANKS} ranks tfinal=${TF_CALI} s"
echo "Output:      ${OUT_BASE}"

RC_ALL=0

# ---------------- Phase A: np-sweep verification (substep path) -------------
if [[ "${TPV_SKIP_VERIFY:-0}" != "1" ]]; then
  for NP in ${NP_SWEEP}; do
    OUT="${OUT_BASE}/np${NP}"
    mkdir -p "${OUT}"
    echo; echo "=== Phase A: np=${NP} ==="
    srun --mpi=pmi2 -n "${NP}" stdbuf -oL -eL ./seas_tpv205_driver \
        --mesh "${MESH}" \
        --mesh-scale 1 \
        --order 1 \
        --ader-order 2 \
        --bc-mode absorbing \
        --cfl "${CFL}" \
        --tfinal "${TF_VERIFY}" \
        --friction-solver newton-stable \
        --fric-law lsw \
        --mixed-flux none \
        --fault-iterator substep \
        --output-dir "${OUT}" \
        --output-prefix tpv205_p1_np${NP} \
        --output-dt 0.05 \
        > "${OUT}/run.log" 2>&1
    rc=$?
    echo "np=${NP} exit=${rc}  $(grep -m1 'V_max' "${OUT}/run.log" | tail -1)"
    tail -2 "${OUT}/run.log"
    [[ ${rc} -ne 0 ]] && RC_ALL=1
  done

  # ---------------- Phase B: cross-np station comparison --------------------
  echo; echo "=== Phase B: cross-np station comparison ==="
  set -- ${NP_SWEEP}
  REF_NP=$1
  CMP_ARGS="${OUT_BASE}/np${REF_NP}"
  for NP in ${NP_SWEEP}; do
    [[ "${NP}" == "${REF_NP}" ]] && continue
    CMP_ARGS="${CMP_ARGS} ${OUT_BASE}/np${NP}"
  done
  python3 "${SEAS_MFEM_ROOT}/miniapps/seas/jobs/unify_phase6_expanse/compare_np_sweep.py" \
      ${CMP_ARGS} | tee "${OUT_BASE}/RESULT_np_sweep.txt"
  [[ ${PIPESTATUS[0]} -ne 0 ]] && RC_ALL=1
fi

# ---------------- Phase C: Caliper speed profile -----------------------------
if [[ "${TPV_SKIP_CALI:-0}" != "1" ]]; then
  echo; echo "=== Phase C: Caliper profile (P=1, ${CALI_RANKS} ranks) ==="
  if ! ldd ./seas_tpv205_driver 2>/dev/null | grep -qi 'libcaliper'; then
    echo "WARNING: ./seas_tpv205_driver is NOT Caliper-linked (MFEM_USE_CALIPER=NO) —"
    echo "  the MFEM_PERF_SCOPE regions are compiled out.  Rebuild with"
    echo "  USE_CALIPER=YES bash build_expanse.sh; SKIPPING Phase C."
  else
    OUT="${OUT_BASE}/cali"
    mkdir -p "${OUT}"
    LOGBASE="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/unify_phase6_expanse/tpv205_p1_cali_${SLURM_JOB_ID}"
    export CALI_CONFIG="${TPV_CALI_CONFIG:-runtime-report(calc.inclusive,region.count,profile.mpi,max_column_width=60,output=${LOGBASE}.cali-region-report.txt),spot(output=${LOGBASE}.cali)}"
    echo "CALI_CONFIG: ${CALI_CONFIG}"
    srun --mpi=pmi2 -n "${CALI_RANKS}" stdbuf -oL -eL ./seas_tpv205_driver \
        --mesh "${MESH}" \
        --mesh-scale 1 \
        --order 1 \
        --ader-order 2 \
        --bc-mode absorbing \
        --cfl "${CFL}" \
        --tfinal "${TF_CALI}" \
        --friction-solver newton-stable \
        --fric-law lsw \
        --mixed-flux none \
        --fault-iterator substep \
        --output-dir "${OUT}" \
        --output-prefix tpv205_p1_cali \
        --output-dt 0.05 \
        > "${OUT}/run.log" 2>&1
    rc=$?
    unset CALI_CONFIG
    echo "cali exit=${rc}  steps: $(grep -m1 'steps' "${OUT}/run.log" | tail -1)"
    [[ ${rc} -ne 0 ]] && RC_ALL=1
    if [[ -f "${LOGBASE}.cali-region-report.txt" ]]; then
      echo "--- Caliper region report (head) ---"
      head -30 "${LOGBASE}.cali-region-report.txt"
      echo "(full: ${LOGBASE}.cali-region-report.txt; Spot: ${LOGBASE}.cali)"
    fi
  fi
fi

echo; echo "DONE  RC_ALL=${RC_ALL}"
exit ${RC_ALL}
