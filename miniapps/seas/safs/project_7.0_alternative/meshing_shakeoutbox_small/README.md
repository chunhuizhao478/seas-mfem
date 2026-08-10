# ALT small mesh, footprint extended to the ShakeOut 2008 box

`results/safalt_small_shakeoutbox.puml.h5` — 1,804,173 tets / 354,233 verts,
80,672,608 B, md5 `f25ffec82ca6343532cb8edaec05b86d`.

Parent: `~/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_small_.../mesh_alt.puml.h5`
(1,319,294 tets), byte-identical to `meshing/results/msh/safalt_3fault_daylight_2M.puml.h5`.

Domain box and its derivation: **`../DOMAIN_SHAKEOUTBOX.md`** (shared by all
three tiers, and by construction ALT n PREFERRED = the PREFERRED box).

| | parent | new |
|---|---|---|
| footprint | 253.2 x 420.4 km rotated 30 deg, **106,476 km2** | axis-aligned **714.000 x 483.000 km = 344,862 km2** |
| depth | -22,100.9 m | unchanged |
| tets | 1,319,294 | 1,804,173 (+484,879, +36.8 %) |
| fault | 320,560 facets (160,280 tri) | **identical, area delta exactly 0.000e+00 m2** |
| free surface | 62,387 facets | 87,320, still exactly flat at z = 0 |

## Reproduce

```bash
P=~/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_small_plast_phi30_40_gradedfw_k1p40_attenuation/mesh_alt.puml.h5
CVM=~/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc

python code/build_collar.py --parent $P --out build_tmp/collar_small_v2.npz \
    --h-mode const --h-const 4800 --h-max 6000

python code/merge_collar.py --parent $P --collar build_tmp/collar_small_v2.npz \
    --out results/safalt_small_shakeoutbox.puml.h5

python code/verify_extension.py --parent $P \
    --new results/safalt_small_shakeoutbox.puml.h5 --cvm $CVM --gate 0.6667
```

### Why `--h-mode const 4800` and not the 0.5 Hz gate

**This tier has never carried a volume frequency gate**, so the collar is sized
to the parent's own far-field rule instead. Measured on the parent itself:
worst Vs/dx **0.0782**, with 166,142 of its 1,319,294 cells below 0.6667, and a
far-field element size of dx median 4,748 m (0–20 km inside the edge) /
4,627 m (20–60 km).

The collar comes out **better** than the mesh it extends: collar worst Vs/dx
**0.1314** against the parent's 0.0782. Whole-mesh gate failures are 457,885
(25.4 %) — that number is inherited from the parent, not introduced here.

Sizing the collar to 3,000 m was tried first and rejected: it produced
3,734,144 tets, 4.5x the parent's own far-field cell density per unit volume,
for resolution the parent does not have anywhere.

## Acceptance — 12/12 PASS (`code/verify_extension.py`)

| check | result |
|---|---|
| P1 parent block bit-identical | geometry + connect identical; only change is 5,460 tets whose absorbing side face became interior |
| P2 fault triangle multiset | 320,560 -> 320,560 |
| P3 fault area delta | **exactly 0.000e+00 m2** (13,073.313420 km2) |
| P4 free surface flat | z spread 7.4e-14 m at z = 0 |
| P4 free surface area | 344,862.0 km2 == box area |
| P5 BC round-trip | tagged 501,196 == hull 180,636 + fault x2 320,560 |
| P6 inverted tets | 0 |
| P7 covers ShakeOut v1 bbox | yes, with 2.6–16.4 km margin |
| P8 collar vs parent gate | collar worst 0.1314 >= parent worst 0.0782 |
| P9 collar eta<0.05 | **0** (collar eta_min 0.0953 > parent 0.0367) |
| P9 collar min edge | 1,098.1 m vs parent 9.34 m |

Every quality COUNT decomposes exactly as parent + collar, which is the
positive evidence that the parent's own cells were untouched.

## Mesh facts

```
tets       1,804,173        verts   354,233
x          72,000.0 .. 786,000.0        (714.000 km)
y          3,524,000.0 .. 4,007,000.0   (483.000 km)
z          -22,100.9 .. 0.0             (22.101 km)
BC faces   fault 320,560 (= 2 x 160,280 tri)  free 87,320  absorbing 93,316
fault area 13,073.313420 km2      lid area 344,862.0 km2
edges      min 9.34   med 1,835   max 8,097 m
eta        min 0.0367  med 0.799   <0.05 15   <0.1 298
gate       worst 0.0782 (0.0587 Hz at p3)   below 0.6667: 457,885 (25.379 %)
```

`eta_min` 0.0367, `eta<0.05` 15 and min edge 9.34 m are all the PARENT's own
values, carried over unchanged; the collar's are 0.0953 / 0 / 1,098 m.

## Notes for the deck

* the lid is still **exactly flat at z = 0**, so the free-surface definition is
  unchanged and receivers at z = -1 m remain below their local top by
  construction;
* the domain now reaches E 786 km while the CVM hull stops at E 705 km, so part
  of the eastern collar has no CVM column and falls to ASAGI's edge-clamp and
  `safs_material_cvm.yaml`'s `!ConstantMap` (Vs 3462 m/s, Qs 173.1). Ground
  motion there is on a constant medium and should not be read as a basin
  prediction. Run `code/stage_f_check.py` for exact per-grid coverage.
