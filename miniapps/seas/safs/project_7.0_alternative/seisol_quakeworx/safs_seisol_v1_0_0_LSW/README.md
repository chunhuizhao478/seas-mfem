# Build `safs_mesh.puml.h5` on a Linux machine — step by step

SeisSol reads the PUML format (`.puml.h5`), not gmsh `.msh`. The official converter
is **pumgen**, which builds cleanly on Linux (it does not compile on macOS due to a
`uint64_t`/`size_t` type clash in its own source). The mesh is already **retagged** to
SeisSol's gmsh boundary convention (`safs_mesh.msh`), so all you do on Linux is build
pumgen and run the conversion.

> Files you need on Linux: **`safs_mesh.msh`** (the retagged mesh — NOT the original
> MFEM-tagged `.msh`) and **`build_and_run_pumgen_linux.sh`**.

---

## Step 0 — copy two files from the Mac → Linux

On the Mac (replace `USER@LINUX` and the target dir):

```bash
cd /Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seissol
scp safs_mesh.msh build_and_run_pumgen_linux.sh USER@LINUX:~/safs/
```

## Step 1 — build pumgen + convert (one command)

On the Linux machine:

```bash
cd ~/safs
bash build_and_run_pumgen_linux.sh ~/safs/safs_mesh.msh
```

This installs dependencies, clones + builds the official pumgen, and runs the
conversion. Result: **`~/safs/safs_mesh.puml.h5`** (plus `safs_mesh.xdmf`).

### Manual alternative (if you prefer explicit steps)

```bash
# (a) dependencies — pick ONE:
#   Debian/Ubuntu:
sudo apt-get install -y cmake g++ git libopenmpi-dev libhdf5-openmpi-dev
#   or conda:
#   conda create -y -n pumgen-build -c conda-forge cxx-compiler cmake make openmpi "hdf5=*=mpi_openmpi*"
#   conda activate pumgen-build

# (b) build official pumgen
git clone --depth=1 --recursive https://github.com/SeisSol/PUMGen.git
cd PUMGen && mkdir build && cd build
CC=mpicc CXX=mpicxx cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j

# (c) convert (run from the directory containing safs_mesh.msh)
cd ~/safs
mpirun -np 4 ~/safs/PUMGen/build/pumgen -s msh2 safs_mesh.msh
# -> produces safs_mesh.puml.h5
```

For GMSH conversion, HDF5 is the only required dependency — **no** SCOREC/PUMI or
Simmetrix (per the PUMGen wiki). This mesh (~1 M tets) is well within the HDF5-only path.

## Step 2 — verify the result (recommended)

```bash
python3 - <<'EOF'
import h5py, numpy as np
f = h5py.File('safs_mesh.puml.h5','r'); b = f['boundary'][:].astype(np.int64)
faces = {}
for k in range(4):
    for c,n in zip(*np.unique((b>>(8*k))&0xFF, return_counts=True)):
        faces[int(c)] = faces.get(int(c),0) + int(n)
print('elements', f['connect'].shape[0], 'nodes', f['geometry'].shape[0])
print('boundary codes', faces)
EOF
```

**PASS** = boundary codes are only `{0, 1, 3, 5}` with
`3`≈60,658 (fault / dynamic rupture), `1`≈34,096 (free surface), `5`≈37,398 (absorbing).

## Step 3 — place it next to the parameter file

```bash
scp ~/safs/safs_mesh.puml.h5 USER@MAC:/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seissol/
```

`parameters.par` already references `MeshFile='safs_mesh.puml.h5'`, so the case is
ready to run once the SeisSol solver is built (with your chosen `-DORDER=`).

---

## Background (how the mesh got here)

- **Retag (already done on the Mac):** the MFEM tags (101 fault / 102 top / 103 bottom /
  104 sides) clash with SeisSol's gmsh convention (`docs/gmsh.rst`: 101 = free surface,
  103 = dynamic rupture, 105 = absorbing; pumgen subtracts 100 → BC 1/3/5). `retag_msh_to_seissol.py`
  remapped them: fault 101→103, top 102→101, bottom 103→105, sides 104→105 (rock volume stays 1),
  producing `safs_mesh.msh`. Geometry/connectivity unchanged.
- **Why Linux:** pumgen mixes `uint64_t` and `std::size_t` for connectivity, which only
  compiles where they are the same type (Linux LP64), not macOS arm64.

## Related files
- `safs_mesh.msh` — retagged gmsh mesh (pumgen input)
- `build_and_run_pumgen_linux.sh` — official pumgen build + convert (run on Linux)
- `retag_msh_to_seissol.py` — the retagger (tag remap only; reproducible)
- `MESH_BUILD_README.md` — longer notes; `README_port.md` — the SeisSol run guide
