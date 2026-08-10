# ALT heavy mesh, footprint extended to the ShakeOut 2008 box

`results/safalt_fb200_deep40km_refine2_shakeoutbox.puml.h5`
— 144,428,706 tets / 25,459,993 verts, 6,388,192,168 B,
md5 `d8f3176d5c0a3e3908eb4070fc4d9954`.

Parent: `meshing_deep40km_refine2/results/safalt_fb200_deep40km_refine2.puml.h5`
(122,162,105 tets), the mesh behind SeisSol run `7885821`.

Domain box and its derivation: **`../DOMAIN_SHAKEOUTBOX.md`**.

| | parent | new |
|---|---|---|
| footprint | 253.2 x 420.4 km rotated 30 deg, **106,476 km2** | axis-aligned **714.000 x 483.000 km = 344,862 km2** |
| depth | -40,000 m | unchanged |
| tets | 122,162,105 | 144,428,706 (+22,266,601, **+18.2 %**) |
| fault | 2,564,480 triangles | **identical, area delta exactly 0.000e+00 m2** |
| free surface | flat at z = 0 | 5,424,501 facets, still exactly flat |

## The collar is SHARED with the intermediate mesh

`../meshing_shakeoutbox_intermediate/build_tmp/collar_deep_v8.npz` is welded
onto BOTH deep parents. That is sound because they carry a **bit-identical side
wall**: 26,704 vertical absorbing triangles / 15,205 vertices, verified as an
exact coordinate multiset. `refine2` only ever touched the fault band, so the
far field it inherited from `safalt_0d5Hz_p3_deep40km` is untouched — measured
independently as 10,690,093 vs 10,695,022 cells in the >40 km-from-fault
columns, with an identical per-Vs-bin table.

`merge_collar.py` relocates the wall on whichever parent it is handed, by exact
coordinate key, never by the vertex or tet indices stored in the collar npz —
the two parents number their vertices completely differently (20.9M vs 3.2M).

```bash
BASE=$PWD/../..
python code/merge_collar.py \
  --parent $BASE/meshing_deep40km_refine2/results/safalt_fb200_deep40km_refine2.puml.h5 \
  --collar $BASE/meshing_shakeoutbox_intermediate/build_tmp/collar_deep_v8.npz \
  --out results/safalt_fb200_deep40km_refine2_shakeoutbox.puml.h5
```

## Weld and BC census (from the merge)

```
WELD OK: 26,704 interface faces identical on both sides
hull identity: tagged 10,630,428 == hull 5,501,468 + fault 5,128,960 -> OK
collar signed volume: 0 negative
merged: 144,428,706 tets, 25,459,993 verts
```

## Acceptance — 9/10 PASS (`code/verify_extension_lean.py`)

| check | result |
|---|---|
| P1 parent block bit-identical | geometry + connect streamed, identical |
| P2 fault triangle multiset | 5,128,960 -> 5,128,960 |
| P3 fault area delta | **exactly 0.000e+00 m2** |
| P4 free surface flat | z spread 1.7e-10 m at z = 0 |
| P4 free surface area | 344,862.0 km2 == box area |
| P6 inverted tets | 0 |
| P7 covers ShakeOut v1 bbox | yes, 2.6–16.4 km margin |
| P9 collar adds no short edge | collar min edge **61.31 m** vs parent 2.34 m |
| P9 collar has no eta < 0.05 sliver | **0** (collar eta_min 0.0663, parent 0.0763) |
| **P8 gate** | **FAIL — 1,226,028 of the collar's 22,266,601 cells (5.506 %)** |

`verify_extension_lean.py` streams every pass. The non-streaming verifier holds
both meshes' `connect` as int64 at once and peaked at **15.9 GB with 75 MB of
free RAM** on this mesh before being killed. P5's independent hull census is
absent for the same reason (it sorts 4x nt faces, ~12 GB here); merge_collar
checks the equivalent identity in O(1), cross-validated against the full census
on the small and intermediate meshes.

## The gate residual

**5.506 % of the collar**, worst 0.2615 (0.196 Hz at p3). Over the whole mesh
that is 0.849 %. **98.3 % of it is in the top 500 m.** The parent's own
122,162,105 cells are bit-identical and keep their compliance: measured
**0 failures, worst exactly 0.6667**.

Full sizing sweep, and the measured reason a lower residual was NOT taken (it
cost a 135x collapse of the timestep floor), is in
`../meshing_shakeoutbox_intermediate/README.md`.

## Mesh facts

```
tets       144,428,706      verts   25,459,993
x          72,000.0 .. 786,000.0        (714.000 km)
y          3,524,000.0 .. 4,007,000.0   (483.000 km)
z          -40,000.0 .. 0.0             (40.000 km)
BC faces   fault 5,128,960 (= 2 x 2,564,480 tri)  free 5,424,501  absorbing 76,967
fault area 13,073.313420 km2      lid area 344,862.0 km2
edges      min 2.34   med 116   max 7,110 m      (min edge is the PARENT's own)
eta        min 0.0663  med 0.860   <0.05 0   <0.1 366
gate       worst 0.2615 (0.1961 Hz at p3)   below 0.6667: 1,226,028 (0.849 %)
```

## Cost — quote LTS, not tets

On the intermediate tier the same collar measured **LTS cost 2.144x for 2.393x
the cells, with dt_min UNCHANGED** and only one extra cluster level. The heavy
parent is 5.5x larger and its dt_min is set by the fault band, which is
untouched, and the collar's min r_insphere (1.330 m) is well above the heavy
parent's own floor (2.34 m min edge / its own insphere), so the collar adds
**no cluster level below the parent's floor** here either.

+18.2 % tets should therefore cost well under +18.2 % run time. **Measure the
LTS cost on this mesh before converting to node-hours.**
