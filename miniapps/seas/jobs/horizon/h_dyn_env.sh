# shellcheck shell=bash
# h_dyn_env.sh — Horizon settings for the GPU dynamic-rupture plan's jobs
# (jobs/horizon/h_d2_build.sbatch, h_d2_gate.sbatch, h_d4_200m.sbatch,
# h_d2_submit.sh).
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
#                        and, for h_d4_200m.sbatch, the 200 m meshes of the
#                        Vista CPU golds (D4_MESH_SHA256)
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

# The 200 m meshes of the Vista CPU golds (made on Vista by
# jobs/vista_nitsche/build_verify_vista.sbatch, gmsh 4.15.0; the sha256 in
# each gold's README)
D4_MESH_SHA256="30c63b923781177830330a36aebcc4159cfc8043374ea4dcd883a9b2eec13538  tpv205/mesh/tpv2053d_200m.msh
10cecfd73d0c7cb93cfa188aac08e4a99a6ec81d439a3363bd7175e3a8ae4aa0  tpv104/mesh/tpv104_200m.msh"
# d4_check_mesh DIR MESH: DIR/MESH (MESH as in D4_MESH_SHA256) has the
# gold's sha256; prints it.
d4_check_mesh () {
   local s
   s=$(cd "$1" 2>/dev/null && sha256sum "$2" 2>/dev/null)
   echo "${s:-no $1/$2}"
   grep -qxF "$s" <<< "${D4_MESH_SHA256}" || { echo "MESH: FAIL ($1/$2 is not the gold's mesh)"; return 1; }
   echo "MESH: the gold's ($1/$2)"
}
