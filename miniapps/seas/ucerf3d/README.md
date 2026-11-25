# UCERF 3-D fault-opening model: geometry and meshes (plan rev5, phase P0)

Builder, checks and small test meshes for the 3-D San Andreas / San Jacinto
model of branch `fully-dynamic-fault-opening-3d`.  The plan is
`/mnt/project-files/ucerf-3d/plan-rev5/ucerf3d-plan-rev5.pdf`, sections 4, 9.4
and 10.2.

## Files

| file | what |
|---|---|
| `build_ucerf3d_mesh.py` | geometry and gmsh mesh: `ucerf` (the model) and `test` (T-plane, T-dip, T-bend, T-pml) |
| `check_ucerf3d_mesh.py` | checks M1-M5 of plan 10.2 on a built mesh (numpy only) |
| `plot_ucerf3d_mesh.py` | figures: fault and creep faces, the San Gorgonio Pass cut-back, the gaps against depth |
| `mesh/T-*.msh`, `.json` | the four test meshes and the builder's record of each |
| `mesh/T-*.check.json` | the checker's report of each test mesh |

## Geometry

- Source: the five R 70 km chains of the 2.5D mesh
  `ucerf2d/mesh/ucerf2d_saf_sjf_500m_tilted_D415km_R70.msh` (sha256
  `d3df54b3...`, checked by the builder).  They are the 10 km depth slice of
  UCERF3 FM3.1.
- saf1, saf3, sjf1, sjf2 are vertical: the chain extruded straight down.
- saf2 (San Gorgonio Pass - Garnet Hill) dips 58 deg toward azimuth 19.8 deg:
  at depth d it is the chain moved along the dip direction by
  (d - 10 km) / tan 58 deg, so it reaches the surface 6.25 km SSW of the chain.
- Each surface is a gmsh spline through the top slice, extruded by
  translation: friction surface 0 to 25 km, option A creep surface 25 to
  50 km (E4).
- D6 (three separate surfaces): the ends of saf2 are cut back until saf2 is at
  least 2 km from saf1 and saf3 at every depth from 0 to 50 km, sampled every
  0.5 km.  Result: 1.85 km at each end (55.9 km -> 52.2 km).  Before the cut
  the gaps were 0.18 km (saf3, at the surface) and 0.21 km (saf1, at 5 km).
- Where two surfaces are closer than 2 h_f, the mesh size there is capped at
  0.4 times their local distance, so every gap holds at least 2 elements
  (only in the 2 km pilot: saf1-saf2, saf3-saf2, saf1-sjf1 2.35 km, sjf1-sjf2
  3.49 km).

## Box, sizes and attributes

- Box: the 2.5D box (tilted to azimuth 123 deg, walls 415 km from every fault
  point, 1394 x 894 km), 415 km deep.  The PML band is inside it, as in 2.5D:
  the cycle code builds the band from the wall faces and `[cycle.pml]
  thickness_km`.
- Size: h = h_f within 5 km of a friction surface, then h_f + 0.2 (d - 5 km),
  capped at 20 km.
- gmsh built-in kernel (Vista's gmsh has no OpenCASCADE), 2-D Frontal-Delaunay,
  3-D HXT, one thread, optimise threshold 0.45, gmsh 2.2 ASCII.  One thread
  makes a rebuild byte-identical (checked: T-dip and the pilot rebuilt to the
  same sha256 with gmsh 4.13.1).

| attribute | faces |
|---|---|
| 3 | friction fault faces (interior) |
| 6 | option A creep faces (interior) |
| 1 | top, z = 0 (free surface) |
| 2 | bottom |
| 5 | side walls |
| volume 1 | everything |

## Meshes

| mesh | command (from `miniapps/seas/ucerf3d`) | tets | sha256 |
|---|---|---|---|
| pilot, 2 km, faults to 20 km, creep to 45 km | `python3 build_ucerf3d_mesh.py ucerf --hf 2000 --fault_depth 20e3 --creep_depth 45e3 --sampling 500 --out ucerf3d_pilot_2km.msh` | 504,801 | `91073aab7a80232472505a958082dac0f02ff7dbbf879ed90801faa538bd2ab1` |
| target, 1 km, faults to 25 km, creep to 50 km | `python3 build_ucerf3d_mesh.py ucerf --hf 1000 --fault_depth 25e3 --creep_depth 50e3 --sampling 1000 --out ucerf3d_target_1km.msh` | 1,767,660 | `8714bfb88899e01baf7ddd81b3687900a122abbe87abc85c9e7b63741be1bdf2` |
| T-plane, vertical 40 x 20 km fault, 120 km box | `python3 build_ucerf3d_mesh.py test --name T-plane --hf 1500 --hmax 15000 --out mesh/T-plane.msh` | 35,753 | `ae3e1d6bd186952cbbfd5dffe77b10f289991306df74c0a2bd9fbecd5f25b7fa` |
| T-dip, the same fault dipping 58 deg toward +y | `... --name T-dip ...` | 39,642 | `7e15ff44b25ac6669e829770c7143a491703ddc2861e71e3e08ff6d55571e422` |
| T-bend, two 20 km segments with a 20 deg bend | `... --name T-bend ...` | 35,814 | `7cd5cc4b65668102267f67511c992110ee276d6b4716e60293aa9fa92624d377` |
| T-pml, T-plane with a 30 km band outside the 120 km box | `... --name T-pml ...` | 39,072 | `a11274318ac32b10a27d0bc68ab46cb8d71aa190b04730491d6190daf9f28a60` |

The test meshes have a friction fault from 0 to 20 km and a creep surface from
20 to 30 km.  The pilot (27 MB) and target (98 MB) files are not committed:
build them from the commit with the command above and compare the sha256.

Check a mesh: `python3 check_ucerf3d_mesh.py MESH.msh` (exit 0 when M1-M5 pass).
