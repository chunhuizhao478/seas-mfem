# ALT intermediate mesh, footprint extended to the ShakeOut 2008 box

`results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5` — 42,783,975 tets / 8,407,615 verts.

Parent: `meshing_deep40km/results/safalt_0d5Hz_p3_deep40km.puml.h5` (15,979,903 tets).

| | parent | new |
|---|---|---|
| footprint | 253.2 x 420.4 km rotated 30 deg, **106,476 km2** | axis-aligned **706.619 x 464.225 km = 328,030 km2** |
| depth | -40,000 m | unchanged |
| tets | 15,979,903 | 42,783,975 (+26,804,072, **+167.7 %**) |
| fault | 320,560 facets | **identical, area delta exactly 0.000e+00 m2** |
| free surface | 62,387 facets | 5,117,184 facets, still exactly flat at z = 0 |

The same collar (`build_tmp/collar_deep_v6.npz`) is welded onto the HEAVY mesh
too — the two parents share a bit-identical 26,704-triangle side wall.

## Reproduce

```bash
BASE=$PWD/../..
CVM=~/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
P=$BASE/meshing_deep40km/results/safalt_0d5Hz_p3_deep40km.puml.h5

python code/build_collar.py --parent $P --out build_tmp/collar_deep_v6.npz \
    --h-mode gate --cvm $CVM --gate 0.6667 \
    --h-min 180 --h-max 5000 --h-safety 0.75 --h-safety-lid 0.68 --seed-top 4000

python code/merge_collar.py --parent $P --collar build_tmp/collar_deep_v6.npz \
    --out results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5

python code/verify_extension.py --parent $P \
    --new results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5 --cvm $CVM --gate 0.6667
```

## Acceptance: 9 of 12 PASS

**Everything structural passes.** The parent is carried over verbatim and the
domain is correct:

| check | result |
|---|---|
| P1 parent block bit-identical | geometry + connect identical; only change is 26,704 tets whose absorbing side face became interior |
| P2 fault triangle multiset | 320,560 -> 320,560, identical |
| P3 fault area delta | **exactly 0.000e+00 m2** |
| P4 free surface flat | z spread 1.5e-10 m at z = 0 |
| P4 free surface area | 328,030.4 km2 == box area |
| P5 BC round-trip | tagged 5,512,556 == hull 5,191,996 + fault x2 320,560 |
| P6 inverted tets | 0 |
| P7 covers ShakeOut v1 bbox | yes |
| P9 collar adds no short edge | collar min edge **16.08 m** vs parent 9.34 m — **no new dt floor** |

**Three fail, and they are real, not bookkeeping.**

### P8 — the collar does not fully hold the 0.5 Hz gate

| | worst Vs/dx | cells below 0.6667 |
|---|---|---|
| parent block | 0.5772 | 2 of 15,979,903 (pre-existing) |
| collar block | **0.3163** | **363,608 of 26,804,072 (1.357 %)** |

**89.7 % of the residual is in the top 500 m** (326,169 cells; median barycentre
depth −102 m), i.e. the layer whose Vs comes from the CVM's slowest bin
(186–560 m/s). Worst resolves 0.237 Hz at p3 instead of 0.5 Hz. Below 500 m the
collar is essentially clean: 215 failures in −2,000..−500 m, 39 in
−4,000..−2,000 m, and 37,185 deep cells (0.14 %).

**Interpret PGV in the new area accordingly**: away from the immediate free
surface it is resolved to 0.5 Hz; the very-near-surface layer over part of the
collar is not.

### P9 — the collar contains slivers the parent does not

`eta_min` 0.0198 (parent 0.0636), 12 cells below 0.05, 1,042 below 0.1, out of
26.8M. Introduced by the interior seed cloud. The dt-relevant metric is
unaffected (min edge 16.08 m > the parent's 9.34 m).

Every quality COUNT decomposes exactly as parent + collar, which is the
positive evidence that the parent's own cells were untouched.

## How the sizing was arrived at (measured, not guessed)

tetgen has no 3-D size field — only boundary point density and one global
volume cap — so the collar's interior sizing had to be built, not repaired.

| variant | lid safety | h_min | collar tets | gate fail | eta_min |
|---|---|---|---|---|---|
| v1 no seed | — | — | 8,795,083 | 13.27 % | 0.0674 |
| v1 + LEB bisection | — | — | 28,929,609 | **14.9 % (WORSE)** | 0.0469 |
| v2 | 1.00 | 250 | 14,763,187 | 14.08 % | 0.0211 |
| v3 | 0.75 | 200 | 25,975,600 | 2.34 % | 0.0053 |
| v4 | 0.60 | 200 | — | tetgen OOM | — |
| v5 | 0.68 | 250 | 27,488,268 | 3.52 % | 0.0127 |
| **v6 SHIPPED** | **0.68** | **180** | **26,804,072** | **1.357 %** | **0.0198** |

Three findings worth keeping:

1. **Longest-edge bisection cannot fix this.** It fixes element SIZE, not
   GRADING. Applied to v1 it tripled the collar to 28.9M tets, left MORE gate
   failures than it started with (1,166,859 -> 1,309,413), and drove the
   minimum edge from 65 m to 8.0 m by hammering one small fine region near the
   wall rim. The marked count peaked at round 5 and then fell 0.15 %/round
   while the mesh grew 0.4-1.2M cells per round.
2. **The size field needs a safety factor, and it is not slack.** It sets
   TRIANGLE edge lengths, but the gate is judged on each TET's MAXIMUM edge,
   which is always longer; gmsh also returns edges ~1.0-1.3x its field. At
   safety 1.0, 97.5 % of failures sat in the top 500 m — exactly the layer
   sized at the limit.
3. **h_min, not the lid factor, binds on the worst columns.** v5 (lid 0.68,
   floor 250 m) came out WORSE than v3 (lid 0.75, floor 200 m) despite more
   tets, because the slowest columns need ~190 m and were clipped to 250.

**Not attempted: mmg.** Closing the last 1.36 % properly wants a true 3-D
metric-driven remesh (`meshing_0d5Hz/code/mmg_refine_sizemap.py`, absolute
metric, wall frozen via RequiredTriangles). That is the designed tool for this
and was not run here.

## Notes for the deck

* the lid is still **exactly flat at z = 0**; receivers at z = -1 m are below
  their local top by construction;
* the new domain reaches E 781.5 km while the CVM hull stops at E 705 km, so
  part of the eastern collar has no CVM column and falls to ASAGI's edge-clamp
  and `safs_material_cvm.yaml`'s `!ConstantMap` (Vs 3462 m/s, Qs 173.1).
  Run `code/stage_f_check.py` for the exact per-grid coverage.
