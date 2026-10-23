#!/bin/bash
#SBATCH -J cpu_clone_build
#SBATCH -o cpu_clone_build_%j.out
#SBATCH -e cpu_clone_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 72
#SBATCH -t 01:30:00
#SBATCH -A EAR20006

# CPU build of MFEM and seas in a second clone, so that CPU gates can run
# while a long job still uses the first checkout ($WORK/seas-mfem-gpu, where
# the full-cycle job runs).  PETSc and its packages are G0's, in
# $WORK/seas-gpu-deps: the deps stage is not run, so the libraries a running
# job links are never rebuilt.  G0's goldens are not in git; they are copied
# into this clone from G0_GOLDEN_SRC if the clone has none.
#
#   git clone -b gpu-qd-bp5 https://github.com/chunhuizhao478/seas-mfem.git $WORK/seas-mfem-gpu-g2
#   cd $WORK/seas-mfem-gpu-g2/miniapps/seas
#   jid=$(sbatch --parsable jobs/vista/cpu_clone_build.sbatch)
#   sbatch --dependency=afterok:${jid} jobs/vista/g2_cpu_gate.sbatch
#
# GitHub is the version control: the clone must be clean and at the commit
# its upstream branch points to, or the job stops.

set -euo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"

cd "${SEAS_MFEM_ROOT}"
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
   echo "ERROR: checkout has uncommitted changes; build only pushed commits." >&2
   git status --short --untracked-files=no >&2
   exit 1
fi
if [ "$(git rev-parse HEAD)" != "$(git rev-parse '@{u}' 2>/dev/null)" ]; then
   echo "ERROR: HEAD $(git rev-parse --short HEAD) is not its upstream" \
        "($(git rev-parse --short '@{u}' 2>/dev/null || echo none)); pull first." >&2
   exit 1
fi
echo "Checkout: ${SEAS_MFEM_ROOT} at $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD))"

G0_GOLDEN_SRC="${G0_GOLDEN_SRC:-${WORK:-$HOME}/seas-mfem-gpu/miniapps/seas/bp5/goldens/g0_vista}"
GOLD="${SEAS_DIR}/bp5/goldens/g0_vista"
if [ "$(realpath -m "${GOLD}")" = "$(realpath -m "${G0_GOLDEN_SRC}")" ]; then
   echo "ERROR: ${SEAS_MFEM_ROOT} holds the goldens itself; this job is for a second clone." >&2
   exit 1
fi
if [ ! -d "${GOLD}" ]; then
   [ -d "${G0_GOLDEN_SRC}" ] || { echo "ERROR: no goldens in ${G0_GOLDEN_SRC}." >&2; exit 1; }
   mkdir -p "$(dirname "${GOLD}")"
   cp -a "${G0_GOLDEN_SRC}" "${GOLD}"
   echo "Copied G0 goldens from ${G0_GOLDEN_SRC}"
fi
for g in p1_r1 p1_r2 p1_r4 p1_r8 p2_r1; do
   if [ ! -f "${GOLD}/${g}/floor.json" ] || [ ! -d "${GOLD}/${g}/run_a" ]; then
      echo "ERROR: golden ${GOLD}/${g} is incomplete (floor.json, run_a)." >&2
      exit 1
   fi
done

MODE=cpu STAGE="mfem seas" JOBS=64 MODULES=none bash build_horizon.sh
g0_check_binary
echo "CPU clone build done: ${SEAS_DIR}/seas_bp5_full"
