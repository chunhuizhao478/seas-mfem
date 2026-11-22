# cmame_revision 2D mesh

- `2d_mesh.geo`, `2d_mesh.msh`: byte copies of farms@cdms 46ad2eb
  `pulsepower/cmame_revision/2d_mesh/` (msh 4.1, the MOOSE input's mesh).
- `2d_mesh_v22.msh`: `gmsh 2d_mesh.msh -0 -format msh22 -o 2d_mesh_v22.msh`
  (gmsh 4.15.2). MFEM's reader rejects the 4.1 file ("vertices indices are
  not unique"); same 14418 nodes and 28376 triangles.
