#!/bin/bash
#SBATCH -J safs_partition_250m
#SBATCH -o safs_partition_250m_%j.out
#SBATCH -e safs_partition_250m_%j.err
#SBATCH -p flex
#SBATCH -N 40
#SBATCH -n 400
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# PRE-SPLIT the 250 m triqsubdiv mesh into per-rank ParPrint files (R-OOM).
# =========================================================================
# The production run job
#   spatial_dyn_ratestate_depthprofile_normalcap_pureupwind_triq_8N_400r_dev_2hr_safs.sbatch
# (8 nodes x 300-400 ranks) aborts with std::bad_alloc during mesh load: the
# default serial read builds the FULL 16 M-element serial mesh on EVERY rank
# before partitioning, so per-node startup memory ~ (ranks_per_node x elements)
# exhausts node RAM.
#
# This job runs ONCE to do that replicated build at the SAME task count the
# production run will use (-n 400) but spread over ENOUGH nodes that the
# one-time peak fits (here 40 nodes => 10 ranks/node; ~4x lower per-node memory
# than the failing 8N x ~38/node).  It writes one ParPrint file per rank +
# a manifest:
#   <prefix>.000000 ... <prefix>.000399
#   <prefix>.np   (== 400)
# Production then reads its own partition with --mesh-mode presplit on just
# 8 nodes (each rank loads ~1/400 of the mesh; no startup peak).
#
# IMPORTANT: the per-rank files are valid ONLY at the np they were written for.
# This job's -n MUST equal the production run job's -n.  If you change the
# production rank count, re-run this job at the new -n.
#
# If THIS job still OOMs at 10 ranks/node, move it to the normal/large queue
# with more nodes (fewer ranks/node) -- production node count is unaffected.
#
# Overrides:
#   SAFS_MESH        input serial .msh (default: the 250 m triqsubdiv mesh)
#   SAFS_MESH_PREFIX output prefix     (default: .../presplit/safs_250m_triqsubdiv_np<N>)

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

# Resolve seas-mfem root by walking UP until miniapps/seas/Makefile is found.
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR}}"
_walk_up_dir="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5; do
    if [[ -f "${_walk_up_dir}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"
        break
    fi
    _walk_up_dir="${_walk_up_dir}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR="
    echo "  ${SLURM_SUBMIT_DIR}.  Set SEAS_MFEM_ROOT explicitly via"
    echo "  sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi
cd "${SEAS_MFEM_ROOT}/miniapps/seas"

# This job does NOT build (avoid the relink race; mirror the run job).  The
# pre-split tool must be pre-built by spatial_dyn_build_safs.sbatch.
if [[ ! -x ./seas_partition_mesh ]]; then
    echo "ERROR: ./seas_partition_mesh not found/executable in $(pwd)."
    echo "       Build it FIRST:"
    echo "         sbatch miniapps/seas/jobs/safs/spatial_dyn_build_safs.sbatch"
    echo "       or on a login/idev node:  cd miniapps/seas && make seas_partition_mesh"
    exit 1
fi

SAFS_ROOT=safs/project_7.0_alternative

# Input serial mesh (default: the 250 m triqsubdiv mesh used by the run job).
TRIQ_MESH="${SAFS_ROOT}/experimental_mesh_refinement/safs_fault_box_nwcut_250m_lcfar3000_z0embed_triqsubdiv.msh"
MESH="${SAFS_MESH:-${TRIQ_MESH}}"
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: input mesh not found: $(pwd)/${MESH}"
    echo "       The .msh is gitignored; scp it to the path above first."
    exit 1
fi

# Output prefix keyed by task count (must match the production run's -n).  This
# is the SAME default the run job resolves for SAFS_MESH_MODE=presplit, so the
# two stay in lockstep.
PRESPLIT_PREFIX="${SAFS_MESH_PREFIX:-${SAFS_ROOT}/experimental_mesh_refinement/presplit/safs_250m_triqsubdiv_np${SLURM_NTASKS}}"
mkdir -p "$(dirname "${PRESPLIT_PREFIX}")"

echo "Git commit:  $(git rev-parse HEAD)"
echo "Queue/size:  development, ${SLURM_NNODES}N x ${SLURM_NTASKS}r, 30 min"
echo "Input mesh:  ${MESH}"
echo "Out prefix:  ${PRESPLIT_PREFIX}  (-> <prefix>.NNNNNN + <prefix>.np)"
echo "NOTE: run production at EXACTLY -n ${SLURM_NTASKS} with SAFS_MESH_MODE=presplit."

ibrun ./seas_partition_mesh --mesh "${MESH}" --out "${PRESPLIT_PREFIX}"
RC=$?
if [[ ${RC} -ne 0 ]]; then
    echo "partition_rc=${RC} — see the .err.  If std::bad_alloc, this job needs"
    echo "  MORE nodes (fewer ranks/node): resubmit on the normal/large queue"
    echo "  with a larger -N (same -n ${SLURM_NTASKS})."
    exit ${RC}
fi

# Post-check: manifest present, equals the task count, and the full file set.
MANIFEST="${PRESPLIT_PREFIX}.np"
if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: manifest ${MANIFEST} was not written."; exit 1
fi
NP_WRITTEN="$(cat "${MANIFEST}")"
N_FILES="$(ls "${PRESPLIT_PREFIX}".[0-9]* 2>/dev/null | grep -c '\.[0-9]\{6\}$')"
echo "manifest np=${NP_WRITTEN}; per-rank files found=${N_FILES} (expected ${SLURM_NTASKS})"
if [[ "${NP_WRITTEN}" != "${SLURM_NTASKS}" || "${N_FILES}" != "${SLURM_NTASKS}" ]]; then
    echo "ERROR: incomplete partition set (np=${NP_WRITTEN}, files=${N_FILES}, expected ${SLURM_NTASKS})."
    exit 1
fi
echo "OK: pre-split complete.  Submit production with:"
echo "  sbatch --export=ALL,SAFS_MESH_MODE=presplit \\"
echo "    miniapps/seas/jobs/safs/spatial_dyn_ratestate_depthprofile_normalcap_pureupwind_triq_8N_400r_dev_2hr_safs.sbatch"
