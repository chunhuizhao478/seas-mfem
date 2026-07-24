#!/bin/bash
#SBATCH --job-name=seas_driver_build
#SBATCH --account=lbl107
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48000M
#SBATCH -t 00:30:00
#SBATCH -o seas_driver_build.%j.out
#SBATCH -e seas_driver_build.%j.out
#SBATCH --export=ALL

# =============================================================================
# CORRECT Expanse build of ./miniapps/seas/seas_spatial_dyn_driver.
#
# Replaces jobs/safs/spatial_dyn_build_safs.sbatch, which is STALE: it loads
# non-existent modules (petsc/3.15, fftw3/3.3.8) and never puts a compiler on
# PATH, so `mpicxx: Command not found` (job 52415017 failed this way, and its
# CLEAN step had already removed the driver objects).  This job uses the SAME
# module + LD_LIBRARY_PATH block as the run sbatch (cpu/0.17.3b + gcc/10.2.0 +
# openmpi/4.1.x), so the compiler resolves and the binary links libcaliper from
# the already-built extern/caliper.  Recompiles the driver TUs against the
# existing MFEM (config.mk MFEM_USE_CALIPER=YES) and relinks — it does NOT
# rebuild MFEM or need internet (no gmsh fetch).
#
# Override: BUILD_JOBS (parallel make jobs; default = cpus-per-task).
# =============================================================================

set -o pipefail
export LC_ALL=C
export LANG=C

module purge
module load slurm 2>/dev/null || true
module load cpu/0.17.3b
module load gcc/10.2.0
_mpi_ok=0
for _mpi in openmpi/4.1.5 openmpi/4.1.3 openmpi/4.1.1 openmpi; do
    if module load "${_mpi}" 2>/dev/null && module is-loaded "${_mpi}" 2>/dev/null; then
        echo "OpenMPI module: ${_mpi}"; _mpi_ok=1; break
    fi
done
[[ "${_mpi_ok}" == "1" ]] || { echo "ERROR: no openmpi module loaded"; exit 1; }
command -v mpicxx >/dev/null 2>&1 || { echo "ERROR: mpicxx not on PATH after module load"; exit 1; }
echo "mpicxx:      $(command -v mpicxx)  ($(mpicxx --version 2>/dev/null | head -1))"

# --- resolve the seas-mfem root (walk UP to miniapps/seas/Makefile) ----------
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SUBMIT_DIR}}"
_walk="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5 6 7; do
    if [[ -f "${_walk}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk}" && pwd)"; break
    fi
    _walk="${_walk}/.."
done
[[ -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]] || { echo "ERROR: cannot locate seas-mfem root from ${SUBMIT_DIR}"; exit 1; }

export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"

cd "${SEAS_MFEM_ROOT}/miniapps/seas"
BUILD_JOBS="${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-16}}"
echo "Git commit:  $(git -C "${SEAS_MFEM_ROOT}" rev-parse HEAD 2>/dev/null)"
echo "Build:       make seas_spatial_dyn_driver -j${BUILD_JOBS}  (incremental against existing MFEM)"

make seas_spatial_dyn_driver -j"${BUILD_JOBS}" || { echo "BUILD FAILED"; exit 1; }

if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: make returned 0 but ./seas_spatial_dyn_driver is missing/not executable."
    exit 1
fi
# The run job REQUIRES Caliper; verify it linked, else A0 has no MPI_Waitall.
if ! ldd ./seas_spatial_dyn_driver 2>/dev/null | grep -qi 'libcaliper'; then
    echo "ERROR: built binary does NOT link libcaliper — A0 cannot measure MPI_Waitall."
    exit 2
fi
echo "BUILD OK: $(pwd)/seas_spatial_dyn_driver"
ldd ./seas_spatial_dyn_driver 2>/dev/null | grep -i libcaliper | sed 's/^/  /'
ls -la --time-style=long-iso ./seas_spatial_dyn_driver | awk '{print "  "$6,$7,$8}'
