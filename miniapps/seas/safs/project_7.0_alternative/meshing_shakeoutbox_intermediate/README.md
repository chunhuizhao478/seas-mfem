# ALT intermediate mesh, footprint extended to the ShakeOut 2008 box

`results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5` — 38,246,504 tets / 7,749,471 verts.

Parent: `meshing_deep40km/results/safalt_0d5Hz_p3_deep40km.puml.h5` (15,979,903 tets).
Domain box and its derivation: `../DOMAIN_SHAKEOUTBOX.md`.

| | parent | new |
|---|---|---|
| footprint | 253.2 x 420.4 km rotated 30 deg, **106,476 km2** | axis-aligned **714.000 x 483.000 km = 344,862 km2** |
| depth | -40,000 m | unchanged |
| tets | 15,979,903 | 38,246,504 (+22,266,601, +139.3 %) |
| fault | 320,560 facets | **identical, area delta exactly 0.000e+00 m2** |

The same collar (`build_tmp/collar_deep_v8.npz`) is welded onto the HEAVY mesh
too — the two parents share a bit-identical 26,704-triangle side wall.

## Reproduce

```bash
BASE=$PWD/../..
CVM=~/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
P=$BASE/meshing_deep40km/results/safalt_0d5Hz_p3_deep40km.puml.h5

python code/build_collar.py --parent $P --out build_tmp/collar_deep_v8.npz \
    --h-mode gate --cvm $CVM --gate 0.6667 --h-min 180 --h-max 5000 \
    --h-safety 0.75 --h-safety-lid 0.68 --seed-top 4000 \
    --seed-clear 0.55 --mindihedral 20

python code/merge_collar.py --parent $P --collar build_tmp/collar_deep_v8.npz \
    --out results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5

python code/verify_extension_lean.py --parent $P \
    --new results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5 --cvm $CVM --gate 0.6667
```

## Acceptance: 9 of 10 PASS

| check | result |
|---|---|
| P1 parent block bit-identical | geometry + connect streamed, identical |
| P2 fault triangle multiset | 320,560 -> 320,560 |
| P3 fault area delta | **exactly 0.000e+00 m2** |
| P4 free surface flat | z spread 1.5e-10 m at z = 0 |
| P4 free surface area | 344,862.0 km2 == box area |
| P6 inverted tets | 0 |
| P7 covers ShakeOut v1 bbox | yes, with 2.6-16.4 km margin |
| P9 collar adds no short edge | collar min edge 61.31 m vs parent 9.34 m |
| P9 collar has no eta < 0.05 sliver | **0** (collar eta_min 0.0663 vs parent 0.0636) |
| **P8 gate** | **FAIL — see below** |

## Cost: quote the LTS number, not the tet count

| | cells | dt_min | clusters | **LTS cost** |
|---|---|---|---|---|
| parent alone | 15,979,903 | 2.895e-04 | 9 | 1.000 |
| **this mesh** | 38,246,504 | **2.895e-04 (UNCHANGED)** | 10 | **2.144** |

cost = sum_e 1/dt_cluster(e) under clustered rate-2 LTS, dt_e ~ r_insphere/Vp.
The collar adds 2.393x the cells but only 2.144x the cost, because the new
cells are shallow far-field and join COARSE clusters — and, critically, it does
**not move dt_min**, so no cluster level is added below the parent's floor.

## P8 — the collar does not fully hold the 0.5 Hz gate

**1,226,028 of the collar's 22,266,601 cells (5.506 %)** are below
Vs/dx = 0.6667; worst 0.2615 (0.196 Hz at p3). Over the whole mesh that is
3.206 %. **98.3 % of the residual is in the top 500 m** (1,205,181 cells,
median barycentre depth −90 m); below 500 m the collar is essentially clean
(70 cells in −2,000..−500 m, 9 in −4,000..−2,000 m).

Interpret PGV in the new area accordingly: resolved to 0.5 Hz away from the
immediate free surface, not at it.

### Why this residual and not a smaller one — the trade that was measured

An earlier collar reached **1.357 %** on the smaller box, but only by seeding
so aggressively that tetgen produced flat slivers:

| collar | box | tets | gate fail | eta_min | r_insphere min | dt_min | clusters | LTS cost |
|---|---|---|---|---|---|---|---|---|
| v6 | old 706.6x464.2 | 26.80M | **1.357 %** | 0.0198 | 0.103 m | 9.80e-05 (**3x worse**) | 11 | 2.766 |
| v7 | new 714x483 | 28.38M | — | 0.0007 | **0.014 m** | 2.15e-06 (**135x worse**) | 17 | 2.935 |
| **v8 SHIPPED** | new 714x483 | 22.27M | 5.506 % | **0.0663** | **1.330 m** | **2.90e-04 (unchanged)** | 10 | **2.144** |

**eta and min-edge both hid this.** v7's min EDGE was a perfectly ordinary
14 m, but a flat tet keeps an ordinary edge while its INSPHERE collapses — and
dt is proportional to insphere, not to edge length. Only measuring
`min r_insphere / Vp` exposed a 135x cut in the timestep floor caused by 18
cells, which cost +37 % LTS for a worse mesh.

v8 buys back the timestep and the quality by keeping every seed at least
0.55 h from any PLC point (`--seed-clear`), at the price of 4x more
near-surface gate cells. That is the right way round: the gate residual is a
bounded, documented, near-surface accuracy limit, while a collapsed dt is an
unbounded cost on every run.

**Not attempted: mmg.** Closing the residual properly wants a true 3-D
metric-driven remesh (`meshing_0d5Hz/code/mmg_refine_sizemap.py`, absolute
metric, wall frozen via RequiredTriangles).

## Notes for the deck

* lid still **exactly flat at z = 0**; receivers at z = -1 m are below their
  local top by construction;
* the domain now reaches E 786 km while the CVM hull stops at E 705 km, so part
  of the eastern collar falls to ASAGI's edge-clamp and
  `safs_material_cvm.yaml`'s `!ConstantMap` (Vs 3462 m/s, Qs 173.1). Run
  `code/stage_f_check.py` for exact per-grid coverage.
