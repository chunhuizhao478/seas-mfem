# g1_cuda_env.sh — environment for the GPU plan's G1 jobs on Vista's
# Grace-Hopper (gh) nodes, in the CUDA clone.  Sourced by jobs/vista/g1_*
# GPU jobs; not a job.
#
# Same as g0_env.sh (threads, v63 flags, meshes, provenance), with the
# CUDA module stack of build_horizon.sh MODE=cuda, and a guard that the
# checkout is clean and at its upstream commit (GitHub is the version
# control: clusters build and run only pushed commits).

G0_MODULES="${G0_MODULES:-gcc/13.2.0 ${CUDA_MODULE:-cuda/13.0} openmpi/5.0.11 nvpl/26.9 cmake/3.31.5}"
G0_DEPS_PROVENANCE="${G0_DEPS_PROVENANCE:-${WORK:-$HOME}/seas-gpu-deps/build_provenance_cuda.txt}"
# shellcheck source=g0_env.sh
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"

# G0's goldens live in the CPU checkout.
G1_GOLDEN_ROOT="${G1_GOLDEN_ROOT:-${WORK:-$HOME}/seas-mfem-gpu/miniapps/seas/bp5/goldens/g0_vista}"

g1_check_checkout() {
   if [ -n "$(git -C "${SEAS_MFEM_ROOT}" status --porcelain --untracked-files=no)" ]; then
      echo "ERROR: ${SEAS_MFEM_ROOT} has uncommitted changes; run only pushed commits." >&2
      git -C "${SEAS_MFEM_ROOT}" status --short --untracked-files=no >&2
      exit 1
   fi
   if [ "$(git -C "${SEAS_MFEM_ROOT}" rev-parse HEAD)" != \
        "$(git -C "${SEAS_MFEM_ROOT}" rev-parse '@{u}' 2>/dev/null)" ]; then
      echo "ERROR: HEAD is not at its upstream commit; pull first." >&2
      exit 1
   fi
   grep -q '^MFEM_USE_CUDA *= YES' "${SEAS_MFEM_ROOT}/config/config.mk" || {
      echo "ERROR: ${SEAS_MFEM_ROOT} is not a CUDA build (build_horizon.sh MODE=cuda)." >&2
      exit 1
   }
   echo "Checkout: ${SEAS_MFEM_ROOT} at $(git -C "${SEAS_MFEM_ROOT}" rev-parse HEAD)"
}
