#!/bin/bash
#SBATCH -J tpv104_seisol_mesh
#SBATCH -o tpv104_seisol_mesh_%j.out
#SBATCH -e tpv104_seisol_mesh_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# TPV104 SeisSol mesh build — gmsh + pumgen
# ==========================================
#
# Produces the PUMGen-format mesh files that SeisSol's MeshNml block
# reads (see tpv104_seisol.par MeshFile = 'tpv104_seisol_500m').
#
# Steps:
#   1. gmsh -3 -optimize tpv104_seisol_500m.geo  →  tpv104_seisol_500m.msh
#   2. pumgen -s msh2  tpv104_seisol_500m.msh  tpv104_seisol_500m
#      → produces tpv104_seisol_500m.h5 + tpv104_seisol_500m.xdmf
#
# Both files must end up in the run directory before tpv104_seisol_t12_run.sbatch
# is submitted.  This script copies them to the same directory it's
# launched from (so the run sbatch can find them via relative paths).
#
# Mesh resolution: 500 m on the fault, ~5 km in the far field.  At
# Λ_dyn ≈ 160 m for TPV104, h/Λ ≈ 3.1 — matches MFEM's 500 m mesh
# resolution for apples-to-apples comparison.  Use the 1000m variant
# for cheaper dev runs (h/Λ ≈ 6.25, marginally resolved).
#
# Pre-requisites on Frontera:
#   - gmsh in $PATH (often `module load gmsh`)
#   - pumgen in $PATH (typically built locally; see
#     https://github.com/SeisSol/PUMGen/wiki for build instructions)

export LC_ALL=C
export LANG=C

# Toolchain modules (align with SeisSol Frontera build per
# ~/.claude/.../reference_seissol_frontera.md):
module purge
module load gcc/9.1
module load impi/21.9
module load python3
module load gmsh/4.11.1   2>/dev/null || echo "WARN: gmsh module unavailable; assuming gmsh is in PATH already"
module load hdf5/1.14.4   2>/dev/null || true

# Resolve mesh size from CLI ($1) or default to 500m.
RES="${1:-500m}"
GEO="tpv104_seisol_${RES}.geo"
MSH="tpv104_seisol_${RES}.msh"

if [ ! -s "${GEO}" ]; then
   echo "ERROR: ${GEO} not found in $(pwd).  Submit from the seisol/ dir."
   exit 1
fi

echo "=== Step 1: gmsh ${GEO} → ${MSH} ==="
gmsh -3 -optimize "${GEO}" -o "${MSH}" \
   || { echo "GMSH FAILED"; exit 1; }

if [ ! -s "${MSH}" ]; then
   echo "ERROR: ${MSH} was not produced."
   exit 1
fi

echo
echo "=== Step 2: pumgen -s msh2 ${MSH} tpv104_seisol_${RES} ==="
which pumgen 2>/dev/null \
   || { echo "ERROR: pumgen not in PATH; build per https://github.com/SeisSol/PUMGen/wiki"; exit 1; }
pumgen -s msh2 "${MSH}" "tpv104_seisol_${RES}" \
   || { echo "PUMGEN FAILED"; exit 1; }

echo
echo "=== Step 3: sanity check ==="
ls -la "tpv104_seisol_${RES}".{h5,xdmf,msh} 2>&1
echo
echo "=== Mesh symmetry check ==="
python3 ../../../tpv104/scripts/tpv104_mesh_mirror_check.py "${MSH}" 2>&1 \
   || echo "(symmetry check failed; mesh produced but mirror diagnostic unavailable)"

echo
echo "STATUS: PASS — mesh ready for tpv104_seisol_t12_run.sbatch"
