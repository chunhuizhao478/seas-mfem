#!/usr/bin/env bash
# =============================================================================
# Build the OFFICIAL SeisSol PUMGen on Linux and convert the TPV104 gmsh mesh
# (tpv104.msh) -> SeisSol PUML (tpv104.puml.h5).  Run this ON A LINUX MACHINE.
#
# Why Linux: PUMGen mixes uint64_t and std::size_t for connectivity, which only
# compiles where they are the same type (Linux LP64).  It does NOT build on
# macOS arm64.  No source changes are needed on Linux.
#
# NO RETAG NEEDED.  tpv104.msh already carries SeisSol's boundary-condition
# physical tags directly:  Physical Surface 1 = free surface,
# 3 = dynamic rupture, 5 = absorbing  (Physical Volume 1 = bulk).  This matches
# the working QuakeWorx tpv13 mesh exactly (verified: tpv13_training.puml.h5 and
# safs_mesh.puml.h5 both decode to per-face BC codes {0,1,3,5}).  PUMGen maps
# these tags straight through, so do NOT remap them to 101/103/105.
#
# Usage:
#   bash build_puml_mesh.sh [/abs/path/to/tpv104.msh]
#   (defaults to ./tpv104.msh next to this script)
# Produces:  tpv104.puml.h5  (+ tpv104.xdmf for visualization), next to the .msh.
# =============================================================================
set -o pipefail   # NB: not `set -u` — `conda activate` trips on unbound vars

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MSH="${1:-$SCRIPT_DIR/tpv104.msh}"
[[ -f "$MSH" ]] || { echo "ERROR: mesh not found: $MSH"; exit 1; }
WORK="$(cd "$(dirname "$MSH")" && pwd)"
NPROC="$(nproc 2>/dev/null || echo 4)"

# ---- 1. Dependencies: HDF5 (parallel) + MPI + OpenMP + CMake + C++ -----------
# Option A (recommended, reproducible) — conda-forge:
if command -v conda >/dev/null 2>&1; then
  if ! conda env list | grep -q '^pumgen-build '; then
    conda create -y -n pumgen-build -c conda-forge \
      cxx-compiler c-compiler cmake make openmpi "hdf5=*=mpi_openmpi*"
  fi
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate pumgen-build
  PREFIX="$CONDA_PREFIX"
  export CC=mpicc CXX=mpicxx
else
  # Option B — system packages (Debian/Ubuntu):
  #   sudo apt-get install -y cmake g++ libopenmpi-dev libhdf5-openmpi-dev
  echo "conda not found — assuming system HDF5(parallel)+MPI+cmake are installed"
  PREFIX="/usr"
  export CC=mpicc CXX=mpicxx
fi

# ---- 2. Build PUMGen (official commands; PUMGen wiki) -------------------------
cd "$WORK"
if [[ ! -d PUMGen ]]; then
  git clone --depth=1 --recursive https://github.com/SeisSol/PUMGen.git PUMGen
fi
mkdir -p PUMGen/build && cd PUMGen/build
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_PREFIX_PATH="$PREFIX" \
         -DCMAKE_INSTALL_PREFIX="$PREFIX"
cmake --build . -j "$NPROC"
PUMGEN="$(pwd)/pumgen"
[[ -x "$PUMGEN" ]] || { echo "ERROR: pumgen did not build"; exit 1; }

# ---- 3. Convert: gmsh msh2 -> PUML .puml.h5 ----------------------------------
cd "$WORK"
# -s msh2  : input is a gmsh v2.2 ASCII .msh
# default output name = <input without ext>.puml.h5  => tpv104.puml.h5
mpirun -np "$NPROC" "$PUMGEN" -s msh2 "$MSH"
echo
echo "DONE. Produced: ${MSH%.msh}.puml.h5  (+ .xdmf)"
echo "It must sit next to parameters.par (MeshFile='tpv104.puml.h5')."
