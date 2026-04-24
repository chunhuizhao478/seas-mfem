#!/bin/bash
#SBATCH -J tpv104_mesh
#SBATCH -o tpv104_mesh_%j.out
#SBATCH -e tpv104_mesh_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# TPV104 Phase-2 Step-10 .msh file generation on Frontera.
#
# Purpose: run `gmsh -3` on the three .geo templates to produce
#   tpv104_1000m.msh  (coarse debug — T_TPV104_MESH_2 wave-speed check)
#   tpv104_500m.msh   (intermediate Frontera runs)
#   tpv104_200m.msh   (production — canonical Phase-3 comparison mesh)
#
# Pre-requisite for tpv104_init_stations.sbatch and
# tpv104_t3_mesh_coupled.sbatch below.
#
# Requires Frontera's `gmsh` module — verify availability with
# `module avail gmsh` before submitting.  If gmsh is not available as
# a module, the script falls back to the conda `pythonenv` path (per
# CLAUDE.md: gmsh under pythonenv).

export LC_ALL=C
export LANG=C

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas/tpv104/mesh

echo "========================================"
echo "Building TPV104 meshes (Gmsh → .msh)"
echo "========================================"
echo "Git commit: $(git -C /scratch2/10024/zhaochun/seas-project/seas-mfem rev-parse HEAD)"
echo "PWD: $(pwd)"
echo "========================================"

# Try Frontera module first, then conda fallback.
GMSH_BIN=""
if module load gmsh 2>/dev/null; then
   GMSH_BIN=$(which gmsh)
elif command -v gmsh >/dev/null 2>&1; then
   GMSH_BIN=$(which gmsh)
else
   source "${HOME}/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || true
   source "${HOME}/miniforge/etc/profile.d/conda.sh"  2>/dev/null || true
   conda activate pythonenv 2>/dev/null || true
   if command -v gmsh >/dev/null 2>&1; then
      GMSH_BIN=$(which gmsh)
   fi
fi

if [ -z "${GMSH_BIN}" ]; then
   echo "ERROR: gmsh binary not found (tried module gmsh, PATH, conda pythonenv)."
   echo "       Install Gmsh or `module load gmsh` to proceed."
   exit 1
fi
echo "Using gmsh: ${GMSH_BIN}"
${GMSH_BIN} --version 2>&1 | head -1
echo "========================================"

for res in 1000m 500m 200m; do
   GEO="tpv104_${res}.geo"
   MSH="tpv104_${res}.msh"
   if [ ! -f "${GEO}" ]; then
      echo "ERROR: ${GEO} not found"
      exit 1
   fi
   echo "[${res}] ${GMSH_BIN} -3 ${GEO} -o ${MSH}"
   ${GMSH_BIN} -3 "${GEO}" -o "${MSH}" 2>&1 | tail -5
   if [ ! -s "${MSH}" ]; then
      echo "FAIL: ${MSH} was not produced"
      exit 1
   fi
   ls -lh "${MSH}"
done

echo "========================================"
echo "All three .msh files built successfully."
ls -lh tpv104_*.msh
echo "========================================"
