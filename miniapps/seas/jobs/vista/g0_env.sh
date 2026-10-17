# g0_env.sh — common environment for the GPU plan's G0 jobs on Vista's
# Grace-Grace (gg) nodes.  Sourced by jobs/vista/g0_*.sbatch; not a job.
#
# Sets the module stack build_horizon.sh MODE=cpu used, pins every thread
# pool to one thread (one MPI rank per core; threaded BLAS under MPI
# oversubscribes cores and makes results run-dependent), locates the
# checkout, and defines the v63 production flags and a provenance writer.

export LC_ALL=C
export LANG=C

# Checkout: the directory holding build_horizon.sh, found from the submit
# directory upwards unless SEAS_MFEM_ROOT is set.
if [ -z "${SEAS_MFEM_ROOT:-}" ]; then
   _probe="${SLURM_SUBMIT_DIR:-$PWD}"
   while [ "${_probe}" != "/" ] && [ -n "${_probe}" ]; do
      if [ -f "${_probe}/build_horizon.sh" ]; then
         SEAS_MFEM_ROOT="${_probe}"
         break
      fi
      _probe="$(dirname "${_probe}")"
   done
   unset _probe
fi
if [ -z "${SEAS_MFEM_ROOT:-}" ]; then
   echo "ERROR: cannot find the checkout (no build_horizon.sh above ${SLURM_SUBMIT_DIR:-$PWD}); set SEAS_MFEM_ROOT." >&2
   exit 1
fi
export SEAS_MFEM_ROOT
SEAS_DIR="${SEAS_MFEM_ROOT}/miniapps/seas"
export SEAS_DIR

# Same modules as build_horizon.sh MODE=cpu.
G0_MODULES="${G0_MODULES:-gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5}"
# shellcheck disable=SC2086
module load ${G0_MODULES}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NVPL_BLAS_NUM_THREADS=1
export NVPL_LAPACK_NUM_THREADS=1

# The v63 production deck (jobs/bp5/bp5_v63_recompute_sliprate_48hr.sbatch)
# without its output and diagnostic flags.  Each job adds --mesh/--order.
G0_V63_FLAGS=(
   --mesh-scale 1000
   --dg-method IP
   --tandem-time-stepping
   --solver mumps
   --petsc-ts
   --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg
)

# Meshes (relative to miniapps/seas).
G0_MESH_PROD="bp5/mesh/reference/bp5_tandem_exact.msh"   # 62,874 tets
# 500 m / 250 m: the same Tandem bp5.geo at res_f 0.5 / 0.25 (fault y=0,
# tags 1/3/5, km).  bp5/mesh/bp5_{500m,250m}.msh use the old x=0 convention.
G0_MESH_500M="bp5/mesh/reference/bp5_tandem_500m.msh"    # 233,783 tets
G0_MESH_250M="bp5/mesh/reference/bp5_tandem_250m.msh"    # 935,485 tets

g0_check_binary() {
   if [ ! -x "${SEAS_DIR}/seas_bp5_full" ]; then
      echo "ERROR: ${SEAS_DIR}/seas_bp5_full missing; run jobs/vista/g0_build.sbatch first." >&2
      exit 1
   fi
   grep -q '^MFEM_USE_PETSC *= YES' "${SEAS_MFEM_ROOT}/config/config.mk" \
      || { echo "ERROR: MFEM was not built with PETSc." >&2; exit 1; }
   grep -q '^MFEM_USE_MUMPS *= YES' "${SEAS_MFEM_ROOT}/config/config.mk" \
      || { echo "ERROR: MFEM was not built with MUMPS." >&2; exit 1; }
}

# gold_readme_job DIR writes DIR/README.md: the commit, job, script, flags,
# meshes (sha256), binary, build and the commands that reproduce DIR.
# shellcheck source=gold_readme.sh
source "${SEAS_DIR}/jobs/vista/gold_readme.sh"

# g0_provenance DIR: what produced the files in DIR (provenance.txt and
# README.md).
g0_provenance() {
   local dir="$1"
   mkdir -p "${dir}"
   {
      echo "job: ${SLURM_JOB_ID:-none} ${SLURM_JOB_NAME:-} on $(hostname), $(date -u +%FT%TZ)"
      echo "nodes: ${SLURM_JOB_NUM_NODES:-?}, tasks: ${SLURM_NTASKS:-?}, partition: ${SLURM_JOB_PARTITION:-?}"
      echo "source: $(git -C "${SEAS_MFEM_ROOT}" rev-parse HEAD) ($(git -C "${SEAS_MFEM_ROOT}" rev-parse --abbrev-ref HEAD)), dirty files: $(git -C "${SEAS_MFEM_ROOT}" status --porcelain --untracked-files=no | wc -l)"
      echo "binary: $(ls -l --time-style=+%FT%T "${SEAS_DIR}/seas_bp5_full" | awk '{print $6, $5}') sha256 $(sha256sum "${SEAS_DIR}/seas_bp5_full" | cut -c1-16)"
      echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
      echo "threads: OMP=${OMP_NUM_THREADS} NVPL_BLAS=${NVPL_BLAS_NUM_THREADS}"
      if [ -f "${G0_DEPS_PROVENANCE:-${WORK:-$HOME}/seas-gpu-deps/build_provenance_cpu.txt}" ]; then
         echo "--- build provenance ---"
         cat "${G0_DEPS_PROVENANCE:-${WORK:-$HOME}/seas-gpu-deps/build_provenance_cpu.txt}"
      fi
   } > "${dir}/provenance.txt"
   gold_readme_job "${dir}"
}
