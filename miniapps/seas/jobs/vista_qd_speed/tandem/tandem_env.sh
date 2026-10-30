# tandem_env.sh: common environment for the Tandem vs MFEM-SEAS BP5 QD speed
# test on Vista's gg nodes.  Sourced by jobs/vista_qd_speed/tandem/*.sbatch;
# not a job.
#
# Tandem is the upstream code as its authors ship it, at one recorded commit
# (TANDEM_COMMIT), in its own Vista folder (TANDEM_ROOT, Andy 2026-10-30
# "approve the Tandem clone").  It builds against the PETSc our CPU runs use
# (build_horizon.sh MODE=cpu: PETSc 3.24.6 with MUMPS, ParMETIS, nvpl BLAS),
# so both codes factor with the same MUMPS on the same BLAS.

source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"

TANDEM_ROOT="${TANDEM_ROOT:-${WORK}/seas-mfem-tandem}"
TANDEM_SRC="${TANDEM_ROOT}/tandem"
TANDEM_COMMIT="b75f66692d299673bf85632ed71e2a7da71ff2e0"   # upstream main, 2026-06-29
TANDEM_PETSC="${TANDEM_PETSC:-${WORK}/seas-gpu-deps/petsc-3.24.6-cpu}"
TANDEM_MESH_SHA256="f413b957cf1c6d02d2c84b65fc05d1f99232301ed8083cb7d8533c9fa442c3b4"

# numpy for yateto (Tandem's kernel generator): Vista's default python3 has
# none; the python3/3.11.8 module does (qd_green_compress.sbatch).
module load python3/3.11.8 eigen/3.4.0 2>/dev/null \
   || { echo "ERROR: module load python3/3.11.8 eigen/3.4.0 failed" >&2; exit 1; }
export PKG_CONFIG_PATH="${TANDEM_PETSC}/lib/pkgconfig${PKG_CONFIG_PATH:+:${PKG_CONFIG_PATH}}"
export CMAKE_PREFIX_PATH="${TANDEM_PETSC}${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}"

# tandem_guard: the Tandem checkout is clean and at TANDEM_COMMIT.
tandem_guard() {
   [ -d "${TANDEM_SRC}/.git" ] \
      || { echo "ERROR: no Tandem checkout at ${TANDEM_SRC}" >&2; exit 2; }
   local head
   head="$(git -C "${TANDEM_SRC}" rev-parse HEAD)"
   [ "${head}" = "${TANDEM_COMMIT}" ] \
      || { echo "ERROR: Tandem at ${head}, expected ${TANDEM_COMMIT}" >&2; exit 2; }
   [ -z "$(git -C "${TANDEM_SRC}" status --porcelain --untracked-files=no --ignore-submodules=none)" ] \
      || { echo "ERROR: uncommitted changes in ${TANDEM_SRC}" >&2;
           git -C "${TANDEM_SRC}" status --short >&2; exit 2; }
   echo "Tandem: ${TANDEM_SRC} at ${head}"
   git -C "${TANDEM_SRC}" submodule status
}

# seas_guard: the CLAUDE.md guard on our checkout (these scripts live there).
seas_guard() {
   cd "${SEAS_MFEM_ROOT}" || exit 2
   git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
   [ -z "$(git status --porcelain --untracked-files=no)" ] \
      || { echo "ERROR: uncommitted changes in ${SEAS_MFEM_ROOT}"; git status --short; exit 2; }
   [ -n "$(git branch -r --contains HEAD)" ] \
      || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
   echo "commit $(git rev-parse HEAD)"
   cd "${SEAS_DIR}" || exit 2
}

tandem_bin() { echo "${TANDEM_ROOT}/build-d3-p$1/app/tandem"; }
