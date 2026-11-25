# spe_hftc meshes

`make_spe_section_2d.py` is a byte copy of farms@cdms
`pulsepower/spe_hftc/mesh/make_spe_section_2d.py`.

| file | command | triangles | nodes |
| --- | --- | --- | --- |
| `spe_section_h0p3.msh` | `python3 make_spe_section_2d.py --h 0.3e-3 --out spe_section_h0p3 --msh-version 2.2 --threads 1` (gmsh 4.15.2) | 66 945 | 33 861 |

The 2D section mesh is committed (3.7 MB) so that pff and the MOOSE
reference runs (M-C1, M-C2a) use the same file. The farms runs used their own
copy made with the same command (msh 4.1, 8 threads), which can differ in
node placement; compare against a MOOSE run made on this file.
Physical names: curves `borehole`, `xmin`, `xmax`, `ymin`, `ymax`; surface
`acrylic`. Borehole radius 4.15 mm (8.3 mm, as in code_2d).

## 3D cube (cube_d7p94_mid)

`make_spe_cube_mesh.py` is a byte copy of farms@cdms 46ad2eb
`pulsepower/spe_hftc/mesh/make_spe_cube_mesh.py` (git blob 9be31f13).
The cube meshes are not committed (850 MB); regenerate them:

| file | command | tets | nodes |
| --- | --- | --- | --- |
| `spe_cube_d7p94_h0p3.msh` | `python3 make_spe_cube_mesh.py --rb 3.97e-3 --h-fine 0.3e-3 --h-coarse 1.0e-3 --r-fine 28e-3 --out spe_cube_d7p94_h0p3 --msh-version 2.2 --threads 4` (gmsh 4.15.2, md5 0773e5e6) | 14 618 527 | 2 386 745 |
| `spe_cube_d7p94_smoke.msh` | `... --rb 3.97e-3 --h-fine 1.5e-3 --h-coarse 1.5e-3 --out spe_cube_d7p94_smoke --msh-version 2.2 --threads 4` | 151 695 | 28 760 |

The Vista MOOSE run cube_d7p94_mid_eta0p45_absorb used the same command in
msh 4.1 (md5 9109e05c, 2.39 M nodes, 14.6 M tets); node placement can differ.
Large 3D meshes: split once with `pff_partition -m <msh> -np <ranks> -o <dir>`
and run `pff_driver -pm <dir>` (every rank otherwise loads the whole mesh,
> 13 GB per rank for the production cube).
