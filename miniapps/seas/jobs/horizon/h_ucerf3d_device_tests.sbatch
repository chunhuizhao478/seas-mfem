#!/bin/bash
#SBATCH -J h_ucerf3d_device_tests
#SBATCH -o h_ucerf3d_device_tests_%j.out
#SBATCH -e h_ucerf3d_device_tests_%j.err
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -A EAR26024

# UCERF 3-D plan rev5, P1.12, check V-G1: the device tests U-G1 to U-G5 of
# P1.7 to P1.11 on real GPUs (CUDA), against the host, on one GB200 node.
#
#   sbatch jobs/horizon/h_ucerf3d_device_tests.sbatch   (from miniapps/seas
#                                                        of the Horizon 3-D clone)
#
# Every test is a pair of runs of one binary: --device cpu writes the host
# result, --device cuda compares with it.  Two lengths, on 1, 2 and 4 ranks
# (all ranks share GPU 0: the test meshes are small):
#   one   one step.  The gate: CUDA within 1e-12 of the host, relative
#         (plan rev5: "CUDA against host: <= 1e-12 relative for one step").
#   full  the length the tests have on the debug device (U-G1, U-G2, U-G5
#         a-c: 300 steps; U-G3: 200 and 289; U-G4 and U-G5 W: 4,264 to 14,440).
#         Reported, with the limit FULL_TOL (default 1e-9; the round-off of
#         CUDA against the host grows with the steps: the D3 gates of the GPU
#         dynamic plan measured up to 2.6e-10 after hundreds of steps).
# A pair that fails is reported and the job goes on; the job's exit code is
# 0 only if every pair passed.
#
# Tests (tests/ucerf3d/ucerf3d.mk):
#   U-G1  seas_test_device_cycle_3d          the cycle's dynamic phase, T-plane
#   U-G2  seas_test_device_cycle_3d          + quench, relaxation r, floor (b, c)
#   U-G3  seas_test_contact_3d               contact: block, T-bend (a-d)
#   U-G4  seas_test_pml_lts_3d --only G4K,G4W  the PML band
#   U-G5  seas_test_pml_lts_3d --only G5W, seas_test_device_cycle_3d (a-c)
#                                            held walls, creep faces
#
# Binary: MFEM and its deps are the clone's last CUDA build
# (BUILD_COMMIT_CUDA, h_d2_build.sbatch); the job stops if anything outside
# miniapps/seas changed from that commit to HEAD, then builds the three test
# binaries from this commit with nvcc (the first CUDA build of the cycle's
# sources: a compile error stops the job with the log's tail).
# Output: $SCRATCH/ucerf3d_device_tests/job<id>/.

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
OUT="${SCRATCH}/ucerf3d_device_tests/job${SLURM_JOB_ID}"
mkdir -p "${OUT}/ref" || exit 2
T0=$(date +%s)
FULL_TOL="${FULL_TOL:-1e-9}"
RANKS="${RANKS:-1 2 4}"
{ echo "job ${SLURM_JOB_ID} on $(hostname), $(date -u +%FT%TZ), commit ${SHA} (${ROOT}); ranks ${RANKS}; full-length limit ${FULL_TOL}"
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader; } \
   | tee "${OUT}/RUN_INFO.txt"

# ---- Binaries ----
BC="${ROOT}/miniapps/seas/BUILD_COMMIT_CUDA"
[ -f "${BC}" ] || { echo "ERROR: no ${BC}; run h_d2_build.sbatch first"; exit 2; }
BUILT=$(cat "${BC}")
OTHER=$(git diff --name-only "${BUILT}" HEAD -- . ':(exclude)miniapps/seas')
[ -z "${OTHER}" ] || { echo "ERROR: changes since the CUDA build ${BUILT} outside miniapps/seas:"; \
                       echo "${OTHER}"; echo "rebuild with h_d2_build.sbatch first"; exit 2; }
grep -q '^MFEM_USE_CUDA *= YES' "${ROOT}/config/config.mk" \
   || { echo "ERROR: ${ROOT}/config/config.mk is not a CUDA build"; exit 2; }
echo "MFEM: CUDA build at ${BUILT} (no change outside miniapps/seas since)" | tee -a "${OUT}/RUN_INFO.txt"
module purge
# shellcheck disable=SC2086   # a list of modules (h_env.sh)
module load ${G0_MODULES} || { echo "ERROR: module load ${G0_MODULES} failed"; exit 2; }
module -t list 2>&1 | tr '\n' ' ' | sed 's/^/modules: /' | tee -a "${OUT}/RUN_INFO.txt"; echo
cd "${ROOT}/miniapps/seas" || exit 2
BINS="seas_test_device_cycle_3d seas_test_contact_3d seas_test_pml_lts_3d"
HDF5_PREFIX="${HDF5_PREFIX:-${HDF5_ROOT}/hdf5-1.14.6-cuda}"
# the makefiles have no complete header dependencies: fresh objects (as
# build_horizon.sh's seas stage)
find . -name '*.o' -not -path './extern/*' -delete
# shellcheck disable=SC2086
if ! make -f tests/ucerf3d/ucerf3d.mk -j 72 SEAS_FS_LIB= HDF5_DIR="${HDF5_PREFIX}" ${BINS} \
      > "${OUT}/build.log" 2>&1; then
   echo "BUILD: FAIL (${OUT}/build.log)"; grep -a -n -m 20 -E "error" "${OUT}/build.log"
   tail -40 "${OUT}/build.log"; exit 1
fi
for b in ${BINS}; do
   [ -x "${b}" ] || { echo "ERROR: ${b} missing after make"; exit 1; }
   echo "binary ${b}: $(sha256sum "${b}" | cut -c1-16)" | tee -a "${OUT}/RUN_INFO.txt"
done
echo "BUILD: ${BINS} at ${SHA} in $(( $(date +%s) - T0 )) s"

# ---- Runs ----
npass=0; nfail=0
SUMMARY="${OUT}/summary.txt"; : > "${SUMMARY}"
pair () {   # NAME NP: runs CMD + HOST, then CMD + DEV (arrays set by the caller)
   local name=$1 np=$2; shift 2
   local t rc1 rc2 log="${OUT}/${name}_np${np}"
   t=$(date +%s)
   ibrun -n "${np}" -o 0 "${CMD[@]}" "${HOST[@]}" > "${log}.host.log" 2>&1; rc1=$?
   ibrun -n "${np}" -o 0 "${CMD[@]}" "${DEV[@]}"  > "${log}.cuda.log" 2>&1; rc2=$?
   local res="PASS"
   if [ "${rc1}" != 0 ] || [ "${rc2}" != 0 ]; then res="FAIL"; nfail=$((nfail + 1)); else npass=$((npass + 1)); fi
   echo "${res} ${name} np ${np}: host exit ${rc1}, cuda exit ${rc2}, $(( $(date +%s) - t )) s" | tee -a "${SUMMARY}"
   grep -a -E "Device configuration|PASSED|FAILED|MFEM abort|ERROR" "${log}.cuda.log" | cut -c1-400 | tee -a "${SUMMARY}"
}
A=tests/ucerf3d/configs/T-bend_loading_A.toml
for np in ${RANKS}; do
   for len in one full; do
      R="${OUT}/ref/${len}_np${np}"
      if [ "${len}" = one ]; then TOL=1e-12; N=(--steps 1); GN=(--g4-steps 1)
      else TOL="${FULL_TOL}"; N=(); GN=(); fi
      # U-G1, U-G2, U-G5 a-c: the cycle's dynamic phase
      cyc () {   # NAME ARGS...
         local name=$1; shift
         CMD=(./seas_test_device_cycle_3d --name "${name}" --tol "${TOL}" "${N[@]}" "$@")
         HOST=(--device cpu --write "${R}.${name}.bin"); DEV=(--device cuda --compare "${R}.${name}.bin")
         pair "${name}_${len}" "${np}"
      }
      # one step: no event yet, so no peak slip rate to ask for
      MP=(); [ "${len}" = one ] && MP=(--min-peak -1)
      cyc U-G1 "${MP[@]}"
      # U-G2 as test-device-quench-3d: a (one step), b and c (300 steps)
      if [ "${len}" = one ]; then
         cyc U-G2a --quench-from 0 --relax --floor 24.5e6 "${MP[@]}"
      else
         cyc U-G2b --quench-from 100 --quench-to 200 --relax --floor 24.5e6
         cyc U-G2c --quench-from 0 --relax
      fi
      cyc U-G5a --deck "${A}" --t-s 1e6 --min-peak -1 $([ "${len}" = full ] && echo "--reenter 50")
      cyc U-G5b --deck "${A}" --t-s 1e6 --min-peak -1 --creep 0
      cyc U-G5c --absorbing off --t-s 1e6 "${MP[@]}"
      # U-G3: contact
      CMD=(./seas_test_contact_3d --tol "${TOL}" --tol-v "$([ "${len}" = one ] && echo 1e-12 || echo "${FULL_TOL}")" "${N[@]}")
      HOST=(--device cpu --write "${R}.g3"); DEV=(--device cuda --compare "${R}.g3")
      pair "U-G3_${len}" "${np}"
      # U-G4 K (no host run of its own: one process compares the two products)
      if [ "${len}" = one ]; then
         CMD=(./seas_test_pml_lts_3d --only G4K); HOST=(--device cpu); DEV=(--device cuda)
         pair "U-G4K" "${np}"
      fi
      # U-G4 W, U-G5 W: the band, held walls
      for g in G4W G5W; do
         CMD=(./seas_test_pml_lts_3d --only "${g}" --g4-tol "${TOL}" \
              --g4-tol-v "$([ "${len}" = one ] && echo 1e-12 || echo "${FULL_TOL}")" "${GN[@]}")
         HOST=(--device cpu --write "${R}.${g}"); DEV=(--device cuda --compare "${R}.${g}")
         pair "U-${g}_${len}" "${np}"
      done
   done
done

echo "UCERF 3-D V-G1 (device tests on CUDA): ${npass} pairs passed, ${nfail} failed;" \
     "commit ${SHA}; logs and summary.txt in ${OUT}; total $(( $(date +%s) - T0 )) s"
[ "${nfail}" = 0 ]
