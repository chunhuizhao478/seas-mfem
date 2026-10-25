# shellcheck shell=bash
# h_dyn_env.sh — Horizon settings for the GPU dynamic-rupture plan's D2 jobs
# (jobs/horizon/h_d2_build.sbatch, h_d2_gate.sbatch, h_d2_submit.sh).
# Sourced, not a job.  It sources h_env.sh, the GPU plan's Horizon toolchain
# (the modules, sm_100 through MFEM_CUDA_ARCH, the CUDA PETSc and cuDSS of
# DEPS_ROOT, NSYS), so both GPU branches build with the same settings, then
# adds what the Nitsche driver's jobs need.  The jobs/vista D2 scripts read
# all of it from the environment.
#
# DYN_HOME is the directory that holds this clone, unless set.  In it:
#   seas-gpu-dyn-deps/   HDF5_ROOT: the parallel HDF5 of fault.vtkhdf, built
#                        by the first h_d2_build job from the tarball fetched
#                        on a login node (h_d2_submit.sh says how)
#   seas-gpu-dyn-refs/   D1_REF_ROOT, copied from Vista: the D1 gg checks at
#                        1, 2 and 4 ranks (d1_gg_check_r<NP>_<jobid>/, their
#                        RUN_INFO.txt and runs/<case>/ref/) and, in meshes/,
#                        D0a's two 1000 m meshes (MESH_SRC, at their paths
#                        under miniapps/seas; every job checks their sha256)
# DEPS_ROOT (h_env.sh: $SCRATCH/seas-gpu-deps) belongs to the GPU
# quasi-dynamic clone; these jobs only read it.

: "${SCRATCH:?h_dyn_env.sh: SCRATCH is not set}"
_h_dyn_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=h_env.sh
source "${_h_dyn_dir}/h_env.sh" || return 1
DYN_HOME="${DYN_HOME:-$(cd "${_h_dyn_dir}/../../../../.." && pwd)}"
unset _h_dyn_dir
export DYN_HOME
export HDF5_ROOT="${HDF5_ROOT:-${DYN_HOME}/seas-gpu-dyn-deps}"
export D1_REF_ROOT="${D1_REF_ROOT:-${DYN_HOME}/seas-gpu-dyn-refs}"
export MESH_SRC="${MESH_SRC:-${D1_REF_ROOT}/meshes}"
echo "Horizon D2 settings: HDF5 ${HDF5_ROOT}; D1 references ${D1_REF_ROOT}; meshes ${MESH_SRC}"
