#!/bin/bash
#SBATCH -J h_ucerf3d_cudss
#SBATCH -o h_ucerf3d_cudss_%j.out
#SBATCH -e h_ucerf3d_cudss_%j.err
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 08:00:00
#SBATCH -A EAR26024

# UCERF 3-D, p = 1: size and speed of the direct QD solve on the P0 meshes
# (ucerf3d/README.md), GPU (cuDSS) against CPU (MUMPS on the node's 144
# Grace cores), on one GB200 node.  A measurement, not a gate.
#
#   sbatch jobs/horizon/h_ucerf3d_cudss.sbatch   (from miniapps/seas of the
#                                                 Horizon 3-D clone)
#
# Uses tests/parallel/test_cudss_vs_mumps (GPU plan G3a) with --backends.  It
# assembles K with BP5's boundary set: fault attr 3, Dirichlet attr 5 (the
# UCERF side walls), natural attr 1 (top).  The UCERF bottom (2) and creep
# faces (6) get no special treatment (bottom free, creep faces ordinary
# interior DG faces), and the material is BP5's.  K has the size and
# sparsity of the UCERF QD matrix, which is what the factor size and the
# time per solve depend on; it is not the UCERF physics.  The meshes are in
# metres, so --mesh-scale 1.
#
#   1. target (1 km)  cuDSS analysis alone (--estimate-only), 1 and 4 GPUs
#   2. pilot  (2 km)  GPU: cuDSS alone (--backends cudss), 1 GPU
#   3. pilot          CPU: MUMPS alone (--backends mumps --device cpu),
#                     144 ranks
#   4. target         GPU: cuDSS alone, 4 GPUs (one rank per GPU)
#   5. target         CPU: MUMPS alone, 144 ranks
# Runs 2-5: analysis + factorization once, then 10 right-hand sides; each
# prints its setup time and the median time per solve, and checks its own
# relative residual (< 1e-10).  A run that fails (e.g. out of memory) is
# reported and the job goes on to the next one.
#
# Meshes: MESH_DIR (default $SCRATCH/ucerf3d-meshes) holds
# ucerf3d_pilot_2km.msh and ucerf3d_target_1km.msh, built from this commit
# with the README commands; the job checks their sha256 (README values
# unless PILOT_SHA256 / TARGET_SHA256 are set) and runs
# check_ucerf3d_mesh.py (M1-M5) on each, and stops on a mismatch.
#
# Binary: MFEM and its deps are the clone's last CUDA build
# (BUILD_COMMIT_CUDA, h_d2_build.sbatch); the job stops if anything outside
# miniapps/seas changed from that commit to HEAD, then builds
# seas_test_cudss_vs_mumps
# (build_horizon.sh STAGE=seas).
# Output: $SCRATCH/ucerf3d_cudss/job<id>/.

set -uo pipefail
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ROOT=$(cd "${SLURM_SUBMIT_DIR}/../.." && pwd)        # submitted from miniapps/seas

# ---- GitHub guard (CLAUDE.md) ----
cd "${ROOT}" || exit 2
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in ${ROOT}"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
[ "$(git rev-parse HEAD)" = "$(git rev-parse '@{u}' 2>/dev/null)" ] \
   || { echo "ERROR: HEAD is not at its upstream commit; pull first"; exit 2; }
SHA=$(git rev-parse HEAD)
echo "commit ${SHA}"

# shellcheck source=h_dyn_env.sh
source "${ROOT}/miniapps/seas/jobs/horizon/h_dyn_env.sh" || exit 2
OUT="${SCRATCH}/ucerf3d_cudss/job${SLURM_JOB_ID}"
mkdir -p "${OUT}"
T0=$(date +%s)
{ echo "job ${SLURM_JOB_ID} on $(hostname), $(date -u +%FT%TZ), commit ${SHA} (${ROOT})"
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader; } \
   | tee "${OUT}/RUN_INFO.txt"

# ---- Binary ----
BC="${ROOT}/miniapps/seas/BUILD_COMMIT_CUDA"
[ -f "${BC}" ] || { echo "ERROR: no ${BC}; run h_d2_build.sbatch first"; exit 2; }
BUILT=$(cat "${BC}")
# seas_test_cudss_vs_mumps is one object linked against libmfem only, and
# it is recompiled below, so miniapps/seas changes (headers included) are
# picked up; anything outside miniapps/seas would need a new libmfem.
OTHER=$(git diff --name-only "${BUILT}" HEAD -- . ':(exclude)miniapps/seas')
[ -z "${OTHER}" ] || { echo "ERROR: changes since the CUDA build ${BUILT} outside miniapps/seas:"; \
                       echo "${OTHER}"; echo "rebuild with h_d2_build.sbatch first"; exit 2; }
echo "MFEM: CUDA build at ${BUILT} (no change outside miniapps/seas since)" | tee -a "${OUT}/RUN_INFO.txt"
rm -f miniapps/seas/tests/parallel/test_cudss_vs_mumps.o miniapps/seas/seas_test_cudss_vs_mumps
if ! MODE=cuda GPU_ARCH="${GPU_ARCH}" STAGE="seas" WITH_HDF5=1 \
      DEPS_ROOT="${DEPS_ROOT}" HDF5_ROOT="${HDF5_ROOT}" \
      PROVENANCE_FILE="${OUT}/build_provenance_cuda.txt" \
      SEAS_TARGETS="seas_test_cudss_vs_mumps" JOBS=72 \
      bash build_horizon.sh > "${OUT}/build.log" 2>&1; then
   echo "BUILD: FAIL (${OUT}/build.log)"; tail -40 "${OUT}/build.log"; exit 1
fi
echo "BUILD: seas_test_cudss_vs_mumps at ${SHA} in $(( $(date +%s) - T0 )) s"

# ---- Meshes ----
cd "${ROOT}/miniapps/seas" || exit 2
MESH_DIR="${MESH_DIR:-${SCRATCH}/ucerf3d-meshes}"
PILOT="${MESH_DIR}/ucerf3d_pilot_2km.msh"
TARGET="${MESH_DIR}/ucerf3d_target_1km.msh"
PILOT_SHA256="${PILOT_SHA256:-91073aab7a80232472505a958082dac0f02ff7dbbf879ed90801faa538bd2ab1}"
TARGET_SHA256="${TARGET_SHA256:-8714bfb88899e01baf7ddd81b3687900a122abbe87abc85c9e7b63741be1bdf2}"
for pair in "${PILOT}:${PILOT_SHA256}" "${TARGET}:${TARGET_SHA256}"; do
   m=${pair%%:*}; want=${pair##*:}
   [ -f "${m}" ] || { echo "ERROR: no ${m}"; exit 2; }
   got=$(sha256sum "${m}" | cut -d' ' -f1)
   echo "${got}  ${m}" | tee -a "${OUT}/RUN_INFO.txt"
   [ "${got}" = "${want}" ] || { echo "ERROR: ${m} sha256 is not ${want}"; exit 2; }
   python3 ucerf3d/check_ucerf3d_mesh.py "${m}" > "${OUT}/$(basename "${m}" .msh).check.log" 2>&1 \
      || { echo "ERROR: check_ucerf3d_mesh.py failed on ${m}"; tail -20 "${OUT}/$(basename "${m}" .msh).check.log"; exit 2; }
   echo "MESH: $(basename "${m}") sha256 and M1-M5 OK"
done

# ---- Runs ----
module purge
# shellcheck disable=SC2086   # a list of modules (h_env.sh)
module load ${G0_MODULES} || { echo "ERROR: module load ${G0_MODULES} failed"; exit 2; }
module -t list 2>&1 | tr '\n' ' ' | sed 's/^/modules: /' | tee -a "${OUT}/RUN_INFO.txt"; echo
overall=0
run () {   # NAME NP ARGS...
   local name=$1 np=$2 t rc
   shift 2
   t=$(date +%s)
   echo "=== ${name} on ${np} GPU(s): start $(date -u +%T) ==="
   ibrun -n "${np}" -o 0 ./seas_test_cudss_vs_mumps --mesh-scale 1 --order 1 \
         "$@" > "${OUT}/${name}.log" 2>&1
   rc=$?
   grep -a -E "^G3a|QD-SOLVER|cuDSS analysis|^  rank|^  [0-9]|solve time|PASS|FAIL" "${OUT}/${name}.log"
   echo "${name}: exit ${rc} in $(( $(date +%s) - t )) s"
   [ "${rc}" = 0 ] || overall=1
}
GPU=(--device cuda --backends cudss)
CPU=(--device cpu --backends mumps)
run target_p1_r1_estimate   1   --mesh "${TARGET}" --device cuda --estimate-only
run target_p1_r4_estimate   4   --mesh "${TARGET}" --device cuda --estimate-only
run pilot_p1_gpu_r1         1   --mesh "${PILOT}"  "${GPU[@]}" --nrhs 10
run pilot_p1_cpu_r144       144 --mesh "${PILOT}"  "${CPU[@]}" --nrhs 10
run target_p1_gpu_r4        4   --mesh "${TARGET}" "${GPU[@]}" --nrhs 10
run target_p1_cpu_r144      144 --mesh "${TARGET}" "${CPU[@]}" --nrhs 10

echo "UCERF 3-D cuDSS p=1: $([ "${overall}" = 0 ] && echo "all runs exit 0" || echo "a run FAILED");" \
     "commit ${SHA}; logs in ${OUT}; total $(( $(date +%s) - T0 )) s"
exit "${overall}"
