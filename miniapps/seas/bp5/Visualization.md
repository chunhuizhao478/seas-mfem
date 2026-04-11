# BP5 Mesh Visualization

## Generating VTK from Gmsh meshes

```bash
gmsh bp5_1000m.msh -0 -o bp5_1000m.vtk
```

The `-0` flag tells Gmsh to read and convert without re-meshing. Physical group tags are preserved as `CellEntityIds` in the VTK file.

## Physical group tags

| Tag | Group | Description |
|-----|-------|-------------|
| 1 | xm | x = -Lx (far-field, natural BC) |
| 2 | xp | x = +Lx (far-field, natural BC) |
| 3 | yp | y = +Ly (Dirichlet: +Vp/2) |
| 4 | ym | y = -Ly (Dirichlet: -Vp/2) |
| 5 | ztop | z = 0 (free surface, natural BC) |
| 6 | zbot | z = Lz (deep boundary, natural BC) |
| 10 | domain | Volume elements (tets) |
| 100 | fault | Fault plane (internal interface at x = 0) |

## Viewing in ParaView

The VTK file contains both 3D tetrahedra (tag 10) and 2D boundary triangles (tags 1–6, 100) in a single file. The 3D tets obscure the 2D faces by default.

### Threshold filter (recommended)

1. Open `bp5_1000m.vtk`
2. **Filters > Threshold**
3. Set Scalars to `CellEntityIds`
4. Set range to isolate groups:
   - `1` to `6` — all boundary faces
   - `100` to `100` — fault surface only
   - `1` to `100`, then uncheck "All Scalars" or exclude 10 — boundaries + fault without volume

### Extract cells by type

1. **Filters > Extract Cells By Type**
2. Uncheck tetrahedra, keep only triangles
3. Color by `CellEntityIds` — shows all boundary and fault faces

### Opacity trick

1. Set the full mesh opacity to 0.1
2. Apply a Threshold for tag 100 on top — fault shows through the transparent volume

## How the simulation code uses these tags

- **Dirichlet BCs**: `elasticity_operator.hpp` checks boundary attributes 3 and 4 for plate loading
- **Fault detection**: geometry-based (centroid at x ≈ 0 within fault bounds), not tag-based
- **Natural BCs** (tags 1, 2, 5, 6): automatic in DG (zero traction, no special handling)
- **Volume attribute 10**: not used by the code (assembles over all elements regardless)
