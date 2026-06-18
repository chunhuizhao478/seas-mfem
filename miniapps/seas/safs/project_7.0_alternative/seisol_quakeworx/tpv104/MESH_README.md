# TPV104 mesh — provenance and how to build `tpv104.puml.h5`

## What's here

| File             | Tracked? | Role |
|------------------|----------|------|
| `tpv104.msh`     | no (`*.msh` is git-ignored) | gmsh v2.2 ASCII mesh, **input to pumgen** |
| `tpv104.puml.h5` | yes      | SeisSol PUML mesh, **what `parameters.par` reads** (`MeshFile = 'tpv104.puml.h5'`) |

`tpv104.puml.h5` is the deliverable mesh (exactly as `tpv13_training.puml.h5` is for
tpv13).  It must be produced by **pumgen on a Linux machine** — see below.

## Provenance of `tpv104.msh`

`tpv104.msh` is a byte-for-byte copy of the MFEM TPV104 200 m mesh:

    miniapps/seas/tpv104/mesh/tpv104_200m.msh

It is a **full** (not half/mirrored) tetrahedral mesh of the SCEC TPV104 domain:

- Domain: x, y in [-60, 60] km, z in [-60, 0] km (free surface at z = 0).
- Fault: vertical strike-slip plane at **y = 0**, sliding region
  x in [-18, 18] km, z in [-18, 0] km (the central 30 x 15 km VW core plus the
  3 km SCEC transition margin on each edge / bottom).
- Fault-local resolution ~200 m (SCEC-adequate: process zone Lambda_dyn ~160 m
  is resolved), coarsening to ~10 km in the far field.
- Coordinate frame matches the SeisSol TPV104 example: x = along-strike,
  y = fault-normal, z = up.

### Boundary-condition tags — NO retag needed

`tpv104.msh` already uses SeisSol's BC physical tags **directly**:

    Physical Surface 1 = free surface
    Physical Surface 3 = dynamic rupture (the y = 0 fault)
    Physical Surface 5 = absorbing (4 outer side walls + bottom)
    Physical Volume  1 = bulk

This is the **same convention the working QuakeWorx tpv13 mesh uses**.  Verified by
decoding the per-face BC codes stored in two known-good PUML meshes — both contain
exactly `{0, 1, 3, 5}`:

    tpv13_training.puml.h5  -> {0,1,3,5}
    safs_mesh.puml.h5       -> {0,1,3,5}

pumgen maps tags 1/3/5 straight through to those BC codes.  **Do NOT** remap them to
101/103/105 — that scheme is only for meshes authored with the +100 gmsh convention
(e.g. the SeisSol Examples `tpv104_half.geo`); feeding 101/103/105 through a
direct-mapping pumgen would yield unrecognized BCs.

## Build `tpv104.puml.h5` (Linux only)

pumgen does **not** compile on macOS arm64, so the conversion must run on Linux
(a login/compute node on the cluster behind QuakeWorx, or any Linux box):

    bash build_puml_mesh.sh            # uses ./tpv104.msh, writes ./tpv104.puml.h5

`build_puml_mesh.sh` builds the official SeisSol PUMGen (conda-forge toolchain) and
runs `pumgen -s msh2 tpv104.msh`.  The result, `tpv104.puml.h5`, must sit next to
`parameters.par`.

If you already have pumgen on PATH, the single relevant command is simply:

    pumgen -s msh2 tpv104.msh          # -> tpv104.puml.h5 (+ tpv104.xdmf)
