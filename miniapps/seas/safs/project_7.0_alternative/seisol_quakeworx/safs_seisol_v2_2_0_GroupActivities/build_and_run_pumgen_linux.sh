#!/usr/bin/env bash
# =============================================================================
# Build the OFFICIAL SeisSol PUMGen on Linux and convert the (already retagged)
# SAFS gmsh mesh -> SeisSol PUML (.puml.h5).  Run this ON THE LINUX MACHINE.
#
# Why Linux: PUMGen mixes uint64_t and std::size_t for connectivity, which only
# compiles where they are the same type (Linux LP64).  It does NOT build on
# macOS arm64.  No source changes are needed on Linux.
#
# Prerequisite file (copy from your Mac):
#   safs_mesh.msh   <- already retagged to SeisSol's gmsh convention
#                      (101=free surface, 103=dynamic rupture, 105=absorbing;
#                       see retag_msh_to_seissol.py / docs/gmsh.rst). DO NOT
#                       re-run pumgen on the original MFEM-tagged .msh.
#
# Usage:
#   bash build_and_run_pumgen_linux.sh /abs/path/to/safs_mesh.msh
# =============================================================================
set -euo pipefail

MSH="${1:?usage: build_and_run_pumgen_linux.sh /abs/path/to/safs_mesh.msh}"
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
  # Then HDF5/MPI are found automatically; PREFIX can stay empty.
  echo "conda not found — assuming system HDF5(parallel)+MPI+cmake are installed"
  echo "  (Debian/Ubuntu: sudo apt-get install cmake g++ libopenmpi-dev libhdf5-openmpi-dev)"
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
# default output name = <input without ext>.puml.h5  => safs_mesh.puml.h5
mpirun -np "$NPROC" "$PUMGEN" -s msh2 "$MSH"
echo
echo "DONE. Produced: ${MSH%.msh}.puml.h5  (+ .xdmf)"
echo "Copy safs_mesh.puml.h5 next to parameters.par (MeshFile='safs_mesh.puml.h5')."
