#!/bin/bash
#SBATCH -J v_ucerf3d_pcg
#SBATCH -o v_ucerf3d_pcg_%j.out
#SBATCH -e v_ucerf3d_pcg_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 06:00:00
#SBATCH -A EAR20006

# UCERF 3-D, p = 1: CPU speed of the iterative QD solve (Q1.1, PCG) on the
# P0 pilot (2 km) and target (1 km) meshes, on Vista gg nodes (2 x 72
# Neoverse V2 cores).  The CPU half of the CPU vs GPU measurement; the GPU
# half (cuDSS) is jobs/horizon/h_ucerf3d_cudss.sbatch.  A measurement, not
# a gate.  K is the SIPG stiffness with UCERF's attributes as
# tests/ucerf3d/test_qd_pcg.cpp sets them (walls 5 Dirichlet, top 1
# natural, bottom 2 and creep faces 6 unlisted, fault 3) and BP5's
# material: UCERF size and sparsity, not UCERF physics.
#
# One-time, on the login node (internet), from the clone root
# ($WORK/seas-mfem-fully-dynamic-fault-opening-3d, cloned with
# --recurse-submodules -b fully-dynamic-fault-opening-3d):
#   STAGE=fetch bash miniapps/seas/jobs/fd_cycle/build_vista.sh
# Then from miniapps/seas:
#   sbatch --export=ALL,PILOT_SHA256=<sha>,TARGET_SHA256=<sha> \
#          jobs/vista_ucerf3d/v_ucerf3d_pcg.sbatch
#   (PART=build|run|all, default all; -N 4 -n 576 with PART=run for the
#    4-node runs on the binary a build job made at this commit)
#
# Steps:
#   1. (build, all) build_vista.sh STAGE="deps mfem": PETSc 3.24.6 (MUMPS,
#      hypre, HDF5) and MFEM in <clone>-deps and in place, then
#      seas_test_qd_pcg (tests/ucerf3d/ucerf3d.mk); BUILD_COMMIT_VISTA
#   2. (run, all) MESH_DIR (default $SCRATCH/ucerf3d-meshes) holds
#      ucerf3d_pilot_2km.msh and ucerf3d_target_1km.msh; their sha256 must
#      be PILOT_SHA256 / TARGET_SHA256 (required: the meshes the Horizon
#      job 106454 used, built with gmsh 4.13.1 on macOS, which differ from
#      the README bytes by platform)
#   3. (run, all) seas_test_qd_pcg --reference none on each mesh, all
#      ranks of the job, both preconditioners (amg, aux_h1), NRHS (default
#      5) stage right-hand sides, to 1e-10 from x = 0: AMG setup, iterations
#      and seconds per solve, true residual <= 10 x tol checked.  A run
#      that fails (e.g. out of memory) is reported and the job goes on.
# Code in $WORK, results in $SCRATCH/ucerf3d_pcg/job<id>/.

export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
export MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
ROOT=$(cd "${SLURM_SUBMIT_DIR}/../.." && pwd)        # submitted from miniapps/seas
DEPS_ROOT=${DEPS_ROOT:-${ROOT}-deps}
PETSC_PREFIX=${DEPS_ROOT}/petsc-${PETSC_VERSION:-3.24.6}
PART=${PART:-all}
NRHS=${NRHS:-5}
case "${PART}" in all|build|run) ;; *) echo "ERROR: PART '${PART}' (all, build or run)"; exit 2 ;; esac

# ---- GitHub guard (CLAUDE.md) ----
cd "${ROOT}" || exit 2
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in ${ROOT}"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
SHA=$(git rev-parse HEAD)
echo "commit ${SHA}"

OUT="${SCRATCH:?}/ucerf3d_pcg/job${SLURM_JOB_ID}"
mkdir -p "${OUT}" || exit 2
T0=$(date +%s)
echo "job ${SLURM_JOB_ID} on $(hostname) (${SLURM_NNODES} node(s), ${SLURM_NTASKS} ranks), $(date -u +%FT%TZ), commit ${SHA}, part ${PART}" \
   | tee "${OUT}/RUN_INFO.txt"

SEAS="${ROOT}/miniapps/seas"
MK="make -f tests/ucerf3d/ucerf3d.mk SEAS_FS_LIB= HDF5_DIR=${PETSC_PREFIX}"

# ---- 1. Build ----
if [ "${PART}" != run ]; then
   # PETSc is skipped when installed; MFEM is rebuilt only when libmfem.a is
   # missing or MFEM's own sources changed since the build that made it.
   ST="deps mfem"
   if [ -f "${ROOT}/libmfem.a" ] && [ -f "${SEAS}/BUILD_COMMIT_VISTA" ] \
      && [ -z "$(git diff --name-only "$(cat "${SEAS}/BUILD_COMMIT_VISTA")" HEAD -- . ':(exclude)miniapps')" ]; then
      ST="deps"; echo "MFEM: libmfem.a from the build at $(cat "${SEAS}/BUILD_COMMIT_VISTA") (no MFEM change since)"
   fi
   JOBS=64 STAGE="${ST}" bash "${SEAS}/jobs/fd_cycle/build_vista.sh" > "${OUT}/build.log" 2>&1
   rc=$?; echo "build_vista.sh ${ST}: exit ${rc} after $(( $(date +%s) - T0 )) s"
   [ ${rc} -eq 0 ] || { tail -40 "${OUT}/build.log"; exit ${rc}; }
fi
for m in ${MODULES}; do module load "${m}" || { echo "ERROR: module load ${m}"; exit 2; }; done
cd "${SEAS}" || exit 2
if [ "${PART}" != run ]; then
   ${MK} -j 64 seas_test_qd_pcg >> "${OUT}/build.log" 2>&1 \
      || { echo "BUILD seas_test_qd_pcg: FAIL"; tail -40 "${OUT}/build.log"; exit 1; }
   echo "${SHA}" > BUILD_COMMIT_VISTA_UCERF3D
   echo "BUILD: seas_test_qd_pcg at ${SHA}"
   [ "${PART}" = build ] && exit 0
else
   BUILT=$(cat BUILD_COMMIT_VISTA_UCERF3D 2>/dev/null)
   CODE=$(git -C "${ROOT}" diff --name-only "${BUILT}" HEAD 2>/dev/null | grep -v -E '^miniapps/seas/jobs/')
   { [ -n "${BUILT}" ] && [ -z "${CODE}" ]; } \
      || { echo "ERROR: seas_test_qd_pcg not built at this code (built '${BUILT}'); PART=build first"; echo "${CODE}"; exit 2; }
   ${MK} -q seas_test_qd_pcg || { echo "ERROR: make would rebuild seas_test_qd_pcg"; exit 2; }
   echo "binary from the build job at ${BUILT}"
fi
echo "binary: $(sha256sum seas_test_qd_pcg | cut -c1-16)" | tee -a "${OUT}/RUN_INFO.txt"

# ---- 2. Meshes ----
MESH_DIR="${MESH_DIR:-${SCRATCH}/ucerf3d-meshes}"
PILOT="${MESH_DIR}/ucerf3d_pilot_2km.msh"
TARGET="${MESH_DIR}/ucerf3d_target_1km.msh"
for pair in "${PILOT}:${PILOT_SHA256:?set PILOT_SHA256}" "${TARGET}:${TARGET_SHA256:?set TARGET_SHA256}"; do
   m=${pair%%:*}; want=${pair##*:}
   [ -f "${m}" ] || { echo "ERROR: no ${m}"; exit 2; }
   got=$(sha256sum "${m}" | cut -d' ' -f1)
   echo "${got}  ${m}" | tee -a "${OUT}/RUN_INFO.txt"
   [ "${got}" = "${want}" ] || { echo "ERROR: ${m} sha256 is not ${want}"; exit 2; }
done

# ---- 3. Runs ----
overall=0
run () {   # NAME MESH
   local name=$1 mesh=$2 t rc
   t=$(date +%s)
   echo "=== ${name} on ${SLURM_NTASKS} ranks: start $(date -u +%T) ==="
   ibrun ./seas_test_qd_pcg --mesh "${mesh}" --order 1 --reference none \
         --nrhs "${NRHS}" --precs "amg aux_h1" > "${OUT}/${name}.log" 2>&1
   rc=$?
   grep -a -E "^PCG timing|QD-SOLVER|rhs |mean iterations|PASS|FAIL|CHECKS|abort|rror" "${OUT}/${name}.log" | head -40
   echo "${name}: exit ${rc} in $(( $(date +%s) - t )) s"
   [ "${rc}" = 0 ] || overall=1
}
run "pilot_p1_r${SLURM_NTASKS}"  "${PILOT}"
run "target_p1_r${SLURM_NTASKS}" "${TARGET}"
echo "UCERF 3-D PCG p=1: $([ "${overall}" = 0 ] && echo "all runs exit 0" || echo "a run FAILED");" \
     "commit ${SHA}; logs in ${OUT}; total $(( $(date +%s) - T0 )) s"
exit "${overall}"
