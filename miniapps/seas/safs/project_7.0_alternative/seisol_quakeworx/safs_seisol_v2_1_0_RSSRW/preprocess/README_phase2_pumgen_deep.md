# Phase 2 runbook — convert `safs_mesh_deep.msh` → `safs_mesh_deep.puml.h5`

Plan: `PLAN_safv4_deep_mesh_port_2026-06-15.md` Phase 2. SeisSol reads PUML, not
gmsh; `pumgen` does the conversion and **only builds/runs on Linux** (uint64_t/
size_t clash on macOS arm64). Run on Linux/Expanse/Frontera, then copy the PUML
back next to `parameters.par`.

The retagged mesh is already built and verified locally (Phase 1):
`safs_mesh_deep.msh` — single DR group (faults 101/102/103 → 103), free surface
201 → 101, absorbing 202/203 → 105.

> Layout note: all Phase-2 tooling lives in `safs_seisol_v2_1_0_RSSRW/preprocess/`.
> The produced `safs_mesh_deep.puml.h5` is a final-submission file and must be copied
> back to the PARENT folder `safs_seisol_v2_1_0_RSSRW/` (next to `parameters.par`),
> NOT into `preprocess/`.

## Files to copy TO the Linux machine

From this folder (`safs_seisol_v2_1_0_RSSRW/preprocess/`):
- `safs_mesh_deep.msh`            (88 MB, the retagged input)
- `build_and_run_pumgen_linux.sh` (builds official PUMGen + runs the conversion)
- `verify_puml_deep.py`           (boundary-code check; optional on Linux)

```bash
# example (edit USER@LINUX and the target dir):
scp safs_mesh_deep.msh build_and_run_pumgen_linux.sh verify_puml_deep.py \
    USER@LINUX:~/safs_deep/
```

## Step 1 — build pumgen + convert (one command, on Linux)

```bash
cd ~/safs_deep
bash build_and_run_pumgen_linux.sh ~/safs_deep/safs_mesh_deep.msh
# -> produces safs_mesh_deep.puml.h5  (+ safs_mesh_deep.xdmf)
```

The script creates/uses a `pumgen-build` conda env (cxx/c compiler, cmake, openmpi,
parallel hdf5), clones+builds `https://github.com/SeisSol/PUMGen.git`, then runs
`mpirun -np <nproc> pumgen -s msh2 safs_mesh_deep.msh`. GMSH conversion needs HDF5
only (no SCOREC/PUMI/Simmetrix). Output name = input stem + `.puml.h5`.

If pumgen is already built on the cluster, you can skip the build and run directly:
```bash
mpirun -np <N> /path/to/pumgen -s msh2 safs_mesh_deep.msh
```

Watch the log for negative-Jacobian / flipped-element aborts. The tetgen mesh
loaded cleanly in MFEM (preferred-project RUNLOG), so this is expected clean —
if pumgen reports a non-negligible flipped count, STOP and report it.

## Step 2 — verify the PUML (Linux or back on the Mac)

```bash
python3 verify_puml_deep.py safs_mesh_deep.puml.h5
```

**PASS requires** (these are exact, derived from the verified `.msh` triangle counts):

| quantity | expected |
|---|---|
| `connect` (tets) | 1,180,159 |
| `geometry` (nodes) | 279,448 |
| BC 1 free surface (former 201) | **175,033** (exact) |
| BC 5 absorbing (former 202+203) | **122,593** (exact) |
| BC 3 dynamic rupture (former 101+102+103) | **132,333** (1×) or **264,666** (2×) |
| boundary codes present | ⊆ {0, 1, 3, 5} |

The DR count may be 1× or 2× the 132,333 fault triangles depending on the pumgen
version (DR is an internal surface; some versions tag both adjacent tet faces).
Either is accepted; the verifier prints which. Free-surface and absorbing are true
boundary faces (one tet each) and must match exactly.

Manual h5py one-liner equivalent (if you prefer not to copy the script):
```python
import h5py, numpy as np
f = h5py.File('safs_mesh_deep.puml.h5','r'); b = f['boundary'][:].astype(np.int64)
faces = {}
for k in range(4):
    for c,n in zip(*np.unique((b>>(8*k))&0xFF, return_counts=True)):
        faces[int(c)] = faces.get(int(c),0)+int(n)
print('elements', f['connect'].shape[0], 'nodes', f['geometry'].shape[0])
print('boundary codes', faces)   # expect {0:.., 1:175033, 3:132333|264666, 5:122593}
```

## Step 3 — copy the PUML back

Copy `safs_mesh_deep.puml.h5` back into the PARENT submission folder
`safs_seisol_v2_1_0_RSSRW/` (next to `parameters.par` — NOT into `preprocess/`):

```bash
scp USER@LINUX:~/safs_deep/safs_mesh_deep.puml.h5 \
    /Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative/seisol_quakeworx/safs_seisol_v2_1_0_RSSRW/
```

(The `.xdmf` is optional — only needed to view the mesh in ParaView.)

Then verify from `preprocess/` (its default now points at the parent folder):
```bash
cd preprocess && python3 verify_puml_deep.py    # reads ../safs_mesh_deep.puml.h5
```

## After this

Phase 3 (case assembly) consumes `safs_mesh_deep.puml.h5`: point
`parameters.par` at `MeshFile = 'safs_mesh_deep.puml.h5'`. Note (Phase 0 finding):
nucleation seeds nearest the Garnet Hill strand, and the faults dip ~64° — the
per-strand strike/dip (`XRef`) check in Phase 5 is the place to validate traction
signs across all three strands.
