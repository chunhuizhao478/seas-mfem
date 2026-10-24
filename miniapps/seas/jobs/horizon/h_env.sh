# h_env.sh — Horizon settings for the GPU plan jobs.  Sourced by
# jobs/horizon/h_*.sbatch, which then run the jobs/vista script of the same
# phase: those scripts read every setting below from the environment, so
# the Vista and Horizon runs share one script per phase.
#
# Horizon differs from Vista's gh nodes in what these settings cover:
#   - 4 GB200 GPUs per node (sm_100, about 186 GB HBM each), 144 Grace
#     cores; partitions debug and debug-shared only (Oct 2026);
#   - no gcc 13.2, cmake 3.31 or cuda 13.0 modules, so the build uses
#     gcc 15.3, cmake 4.4 and cuda 13.1 with the same pinned PETSc 3.24.6,
#     MUMPS, hypre and cuDSS 0.7.1.6 (built from source, build_horizon.sh).
#     Not cuda 13.2 or 13.3: their Thrust (CCCL 3.2 and later) no longer
#     brings thrust::pair, and in 13.3 thrust::tuple and make_tuple, in
#     through the headers hypre 3.0.0 includes, so hypre's CUDA build stops
#     at src/utilities/device_utils.c (job 65658, cuda/13.3);
#   - no $WORK for early users: the dependencies, the goldens and the
#     Green's-function cache live under $SCRATCH (purged when unused).
#
# The G0 goldens are Vista's (gcc 13.2), copied to G1_GOLDEN_ROOT; the
# first Horizon gates measure the compiler and machine change against them.
#
# Order of the chain (jobs/horizon/h_submit.sh submits it with afterok):
#   h_cuda_build -> h_gpu_checks -> h_green_cache -> h_gpu_gates -> h_full_cycle

: "${SCRATCH:?h_env.sh: SCRATCH is not set}"
if ! command -v ibrun >/dev/null 2>&1; then
   echo "ERROR: no ibrun on this machine; the jobs/vista scripts launch with it." >&2
   return 1
fi

export CUDA_MODULE="${CUDA_MODULE:-cuda/13.1}"
# Read by g0_env.sh (module load) and by build_horizon.sh (MODULES).
# openmpi/5.0.11 links libucc.so.1 but its module does not load ucc: with
# cuda/13.1 nothing else does, and ompi_info fails to start without it.
export G0_MODULES="${G0_MODULES:-gcc/15.3.0 ${CUDA_MODULE} ucc/1.9.0 openmpi/5.0.11 nvpl/26.9 cmake/4.4.0}"
export MODULES="${G0_MODULES}"
export GPU_ARCH=100
# MFEM and seas: device code optimized for compute_90 (as on Vista), then
# assembled to sm_100 machine code at build time.  nvcc 13.0 and 13.1
# optimizing for compute_100 take hours on some MFEM kernel files (job
# 67627 hit its 6 h limit in fem/integ/bilininteg_convection_ea.cpp,
# which takes under a minute for sm_90).  PETSc's hypre stays sm_100.
export MFEM_CUDA_ARCH="compute_90 -code=sm_100"
export DEPS_ROOT="${DEPS_ROOT:-${SCRATCH}/seas-gpu-deps}"
export G0_DEPS_PROVENANCE="${DEPS_ROOT}/build_provenance_cuda.txt"
export G1_GOLDEN_ROOT="${G1_GOLDEN_ROOT:-${SCRATCH}/seas-gpu-goldens/g0_vista}"
export GREEN_CACHE="${GREEN_CACHE:-${SCRATCH}/seas-qd-green/bp5_prod_p1.bin}"
echo "Horizon settings: ${G0_MODULES}; sm_${GPU_ARCH}; deps ${DEPS_ROOT};" \
     "goldens ${G1_GOLDEN_ROOT}; G cache ${GREEN_CACHE}"
