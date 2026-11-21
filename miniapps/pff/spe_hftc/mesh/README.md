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
