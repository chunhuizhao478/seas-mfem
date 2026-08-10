# ALT small mesh, footprint extended to the ShakeOut 2008 box

`results/safalt_small_shakeoutbox.puml.h5` — 1,770,259 tets / 347,757 verts.

Parent: `~/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_small_.../mesh_alt.puml.h5`
(1,319,294 tets), byte-identical to `meshing/results/msh/safalt_3fault_daylight_2M.puml.h5`.

## What changed

| | parent | new |
|---|---|---|
| footprint | 253.2 x 420.4 km rotated 30 deg, **106,476 km2** | axis-aligned **706.619 x 464.225 km = 328,030 km2** |
| depth | -22,100.9 m | unchanged |
| tets | 1,319,294 | 1,770,259 (+450,965, +34.2 %) |
| fault | 320,560 facets | **identical, area delta exactly 0.000e+00 m2** |

The new box strictly contains the ShakeOut v1 bbox
(E 74,850.144–781,468.938, N 3,542,641.466–3,993,303.385 m, UTM 11N),
so the PGV comparison is no longer windowed to 28 % of the product's footprint.

## Why the box is 464 km tall and not 451

ShakeOut's bbox is 706.6 x 450.7 km, but the ALT footprint is a rectangle
**rotated 30 deg to the San Andreas strike**, and its north corner sits at
N 3,996,866.9 — already 3.56 km *beyond* ShakeOut's north edge. Snapping the
new box to ShakeOut's north edge would mean cutting into verified mesh.
Keeping it exactly at the parent's own extreme is no good either: the collar
then pinches to a single point at that corner, the annulus stops being an
annulus, and gmsh reports `2 intersections in the 1D mesh` and emits **no
elements at all**. `--margin 10000` gives the collar a floor width there.

## How it is built (2 steps, ~10 s total)

The parent is FROZEN. Its vertical absorbing side wall becomes the inner
boundary of a new collar, reproduced vertex-for-vertex by tetgen `-Y`, so the
collar welds on conformally and nothing inside the old footprint moves. This is
the lateral analogue of `meshing_deep40km/code/{build_deep_slab,merge_slab}.py`.

```bash
P=~/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_small_plast_phi30_40_gradedfw_k1p40_attenuation/mesh_alt.puml.h5
CVM=~/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc

python code/build_collar.py --parent $P --out build_tmp/collar_small.npz \
    --h-mode const --h-const 4800 --h-max 6000

python code/merge_collar.py --parent $P --collar build_tmp/collar_small.npz \
    --out results/safalt_small_shakeoutbox.puml.h5

python code/verify_extension.py --parent $P \
    --new results/safalt_small_shakeoutbox.puml.h5 --cvm $CVM --gate 0.6667
```

### Why `--h-mode const 4800` and not the 0.5 Hz gate

This tier has never carried a volume frequency gate. Measured on the parent
itself: worst Vs/dx = **0.0782**, with 166,142 of its 1,319,294 cells below
0.6667, and a far-field element size of dx median 4,748 m (0–20 km inside the
edge) / 4,627 m (20–60 km). The collar is sized to that same rule, and comes
out **better** than the mesh it extends: collar worst Vs/dx **0.1314**.

Sizing the collar to 3,000 m instead was tried first and rejected — it produced
3,734,144 tets, 4.5x the parent's own far-field cell density per unit volume,
for resolution the parent does not have anywhere.

## Acceptance (12/12 PASS, `code/verify_extension.py`)

| check | result |
|---|---|
| P1 parent block bit-identical | geometry + connect identical; only change is 5,460 tets whose absorbing side face became interior |
| P2 fault triangle multiset | 320,560 -> 320,560, identical |
| P3 fault area delta | **exactly 0.000e+00 m2** (13,073.313420 km2) |
| P4 free surface flat | z spread 7.4e-14 m at z = 0 |
| P4 free surface area | 328,030.4 km2 == box area |
| P5 BC round-trip | tagged 497,666 == hull 177,106 + fault x2 320,560 |
| P6 inverted tets | 0 |
| P7 covers ShakeOut v1 bbox | yes |
| P8 collar vs parent gate | collar worst 0.1314 >= parent worst 0.0782 |
| P9 collar eta<0.05 | **0** (collar eta_min 0.0790 > parent 0.0367) |
| P9 collar min edge | 1,092.9 m vs parent 9.34 m — no new dt floor |

Quality decomposes exactly as parent + collar for every count, which is the
positive evidence that the parent's own cells were untouched.

## Notes for the deck

* the lid is still **exactly flat at z = 0**, so the free-surface definition is
  unchanged and receivers at z = -1 m remain below their local top by
  construction;
* the new eastern edge reaches E 781.5 km while the CVM hull stops at
  E 705 km, so ~76 km of the eastern collar has no CVM column. ASAGI
  edge-clamps and `safs_material_cvm.yaml`'s `!ConstantMap` (Vs 3462 m/s,
  Qs 173.1) is the fallback. Ground motion there is on a constant medium and
  should not be read as a basin prediction. Run `code/stage_f_check.py` for the
  exact per-grid coverage numbers.
