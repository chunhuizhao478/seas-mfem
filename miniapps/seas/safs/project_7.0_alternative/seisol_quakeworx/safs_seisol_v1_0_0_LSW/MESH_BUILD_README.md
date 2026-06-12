# Building `safs_mesh.puml.h5` (SeisSol PUML) — Linux runbook

SeisSol reads PUML (`.puml.h5`), not gmsh `.msh`. The official converter is
**pumgen**. pumgen does not compile on macOS (a `uint64_t`/`size_t` clash in its
own source); on Linux LP64 it builds unmodified. So: prep the mesh on either
machine, build+run pumgen on **Linux**.

## Step 1 — retag the mesh (DONE, on the Mac)

The MFEM mesh tags (101 fault / 102 top / 103 bottom / 104 sides) do **not**
match SeisSol's gmsh convention (`docs/gmsh.rst`: 101=free surface, 103=dynamic
rupture, 105=absorbing). `retag_msh_to_seissol.py` remapped them:

| MFEM | → SeisSol tag | pumgen BC | meaning | triangles |
|---|---|---|---|---|
| 101 fault | 103 | 3 | dynamic rupture | 60,658 |
| 102 top | 101 | 1 | free surface | 34,096 |
| 103 bottom + 104 sides | 105 | 5 | absorbing | 37,398 |
| 1 rock (vol) | 1 | — | material group | — |

Output: **`safs_mesh.msh`** (geometry/connectivity byte-preserved; only physical
tags changed). Use THIS file with pumgen, never the original MFEM-tagged `.msh`.

To reproduce: `python retag_msh_to_seissol.py <orig>.msh safs_mesh.msh`.

## Step 2 — copy to Linux & build/run pumgen

Copy `safs_mesh.msh` and `build_and_run_pumgen_linux.sh` to the Linux machine, then:

```bash
bash build_and_run_pumgen_linux.sh /abs/path/to/safs_mesh.msh
```

The script (official PUMGen wiki procedure):
1. gets deps — HDF5(parallel)+MPI+OpenMP+cmake+C++ (conda-forge, or `apt-get install
   cmake g++ libopenmpi-dev libhdf5-openmpi-dev`);
2. `git clone --depth=1 --recursive https://github.com/SeisSol/PUMGen.git`; cmake + make;
3. `mpirun -np N pumgen -s msh2 safs_mesh.msh` → **`safs_mesh.puml.h5`** (+ `.xdmf`).

For GMSH conversion, HDF5 alone is required — **no** SCOREC/PUMI or Simmetrix
(per the PUMGen wiki). This mesh (~1M tets) is well within the HDF5-only path.

## Step 3 — verify the PUML (on Linux, or copy back to the Mac)

```python
import h5py, numpy as np
f = h5py.File('safs_mesh.puml.h5','r')
b = f['boundary'][:].astype(np.int64)
faces = {}
for k in range(4):
    for c,n in zip(*np.unique((b>>(8*k))&0xFF, return_counts=True)):
        faces[int(c)] = faces.get(int(c),0)+int(n)
print('elements', f['connect'].shape[0], 'nodes', f['geometry'].shape[0])
print('boundary codes:', faces)   # expect only {0,1,3,5}
# sanity: code 3 (DR) face count should ~match the 60,658 retagged fault triangles;
#         code 1 ~34,096 (free surface); code 5 ~37,398 (absorbing).
```
PASS = boundary codes are a subset of {0 interior, 1 free surface, 3 dynamic
rupture, 5 absorbing}, and the per-code face counts match the Step-1 triangle
counts. Then place `safs_mesh.puml.h5` next to `parameters.par`
(`MeshFile='safs_mesh.puml.h5'`).

## Files
- `safs_mesh.msh` — retagged gmsh mesh (input to pumgen)
- `retag_msh_to_seissol.py` — the retagger (tag remap only; not a converter)
- `build_and_run_pumgen_linux.sh` — official pumgen build + convert on Linux
