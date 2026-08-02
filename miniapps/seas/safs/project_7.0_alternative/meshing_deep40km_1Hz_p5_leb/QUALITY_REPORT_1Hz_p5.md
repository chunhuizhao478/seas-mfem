# 1 Hz @ p5 mesh — `safalt_fb200_deep40km_1Hz_p5.puml.h5`

Built 2026-08-01 from `safalt_fb200_deep40km_refine2.puml.h5` by conforming 3-D
longest-edge (Rivara) bisection, driven by a z-pooled Vs target. **Gate PASS: 0 of
133,699,789 cells below the target; the mesh resolves 1.0000 Hz at p5 everywhere.**

## 1. Result

| | parent (refine2) | **this mesh** |
|---|---:|---:|
| tets | 122,162,105 | **133,699,789** (+11,537,684, **+9.44%**) |
| vertices | 20,920,528 | 23,123,324 |
| cells below the gate (Vs/dx < 0.8) | **274,299** | **0** |
| worst Vs/dx | 0.6667 | **0.8000** |
| **resolved everywhere @ p5** | 0.8334 Hz | **1.0000 Hz** |
| **resolved everywhere @ p3** | 0.5000 Hz | **0.6000 Hz** |

Convention (locked, unchanged): `f = Vs/dx`, `dx` = element MAX edge, `Vs =
sqrt(mu/rho)` NEAREST-GRID at the element BARYCENTER from the deck's
`safs_material_cvm.nc`; resolved `f` at order p is `(p/4)·Vs/dx`, so 1 Hz at p5
⇔ `Vs/dx ≥ 0.8`. Verified by `census_1Hz_p5.py` run independently on the
finished file, not by the builder's own check.

## 2. Quality — no regression that matters

| metric | parent | new | note |
|---|---:|---:|---|
| eta_min | 0.0762612 | 0.0762612 | **identical** |
| eta < 0.05 | 0 | 0 | **identical — no new slivers** |
| eta < 0.1 | 127 | 127 | **identical** |
| eta < 0.3 | 23,549 | 23,612 | +63 (+0.27%) |
| eta median | 0.858903 | 0.851297 | **−0.9% — the only real regression** |
| min edge | 2.33589 m | 2.33589 m | **identical** |
| inverted tets | 0 | 0 | identical |
| edges < 50 m | 7,885,386 | 7,885,515 | +129 |
| edges < 100 m | 79,824,898 | 79,969,415 | +144,517 |
| edges < 150 m | 96,111,757 | 97,022,596 | +910,839 |

The extra short edges are the intended product. The worst-case shape metrics
(`eta_min`, `eta<0.05`, `eta<0.1`, `min edge`) are **bit-identical to the parent** —
bisection never produced anything worse than what the parent already contained.

## 3. Fault identity and PUML contract — all PASS

`check_fault_identity.py` against the parent:

```
[PASS] F1 fault triangle multiset identical   5,128,960 vs 5,128,960
[PASS] F2 fault area identical  13073.313420396 km^2  (delta 0.000e+00 m^2)
[PASS] F3 fault faces interior (x2)           2,564,480 unique fault faces
[PASS] F3 boundary faces appear once
[PASS] F4 free surface exactly flat  z in [-0.000000000, 0.000000000], 1,114,474 verts
[PASS] F5 no inverted tets
```

Fault edges were **frozen** throughout (never bisected), so every fault-referenced
deck input — stress nc, friction nc, the `[rs_muw]` LuaMap, nucleation, pickpoints —
transfers **unchanged**. The fault area delta is exactly zero.

## 4. Deck compatibility

| check | result |
|---|---|
| DR facets | 5,128,960 — **identical** to the parent |
| PGV receivers inside the footprint | **97,711 / 97,711**, 0 outside |
| free-surface facets | 1,985,976 → 2,225,273 |
| absorbing facets | 57,790 → 77,127 |
| domain extent | unchanged → `OutputRegionBounds` transfers |

## 5. Run cost — essentially free

| | parent | new | ratio |
|---|---:|---:|---:|
| dt_min | 1.26476e−05 s | 1.26476e−05 s | **1.000** |
| LTS clusters | 11 | 11 | 1.000 |
| LTS updates / simulated s | 3.86559e+11 | 3.93242e+11 | **1.017** |

**+9.4% tets costs only +1.7% run time.** The refinement added shallow far-field
cells that join COARSE LTS clusters, and never touched the dt-limiting cells — so
`dt_min` is unchanged to the last digit. Versus the intermediate deck mesh (job
7885060) this mesh is 57.3× rather than 56.4×, so the 512-node / 48 h sizing of
`..._fbrefine2_...` transfers essentially unchanged.

## 6. How it was built, and three things that did not work

Two passes of `leb_refine_1Hz.py` (pass 1 hops=3, 60 rounds, 75 min, +11,165,152
tets; pass 2 hops=8, 46 rounds, 12 min, +372,532 tets). Peak RSS 20.0 GB.

**Why bisection and not mmg.** Two prior mmg size-map regrades of this mesh
*regressed* the gate (worst 0.6667 → 0.4104 and → 0.4137) because mmg's ±41%
acceptance band licenses COARSENING: compliant cells drift over tolerance and the
pass mints failures as fast as it fixes them. Bisection only splits, so an unsplit
cell's barycenter never moves and cannot newly fail.

**Failure 1 — nearest-grid target (treadmill).** Targeting the raw barycenter Vs
stalled at ~65,000 failures, falling 1.8%/round while adding 157,000 tets/round.
Cause: the CVM is nearest-grid on a 250 m vertical lattice and Vs steps ~4× at the
z = −125 m bin edge, so bisecting a straddling cell throws one child into the SLOW
bin where it needs a 4× smaller dx than its parent did.

**Failure 2 — full 3-D min-pooling (over-refinement).** Pooling Vs over the whole
cell converged but cost **+55.8%** tets on the small ALT mesh vs +0.04% for the raw
rule — it takes the slowest Vs anywhere inside a multi-km cell. **Fix: pool in z
only.** The CVM is 1500 m laterally but 250 m vertically while these cells are
≤ ~500 m across, so bin changes are essentially always vertical. The shipped target
is `min Vs over the cell's own vertical extent at the barycenter's (x,y)` — a lower
bound on any child's measured Vs, hence refinement-stable, while a cell wholly
inside one bin is not penalised at all.

**Failure 3 — LEPP cycling (fixed, but was not the bottleneck).** Ties in
`E.argmax(1)` broke by local slot. This mesh is built by exact red refinement, whose
children are *similar* to their parent, so exactly-equal edge lengths are everywhere
and adjacent tets could name each other's edge, closing the LEPP chain into a cycle.
Ordering ties by the globally unique edge key makes "longest edge" a strict total
order. This cut rim-freezing from 16,536 to 566 terminal edges — but the marked-count
trajectory was unchanged (880,330 vs 880,342 at round 4), so **cycling was not the
cause of the plateau**; the z-pooled target was. Both fixes are kept.

## 7. Residual / caveats

* **None on the gate** — 0 cells fail, verified independently.
* The builder's internal *pooled* target still showed ~17,700 cells at the end of
  pass 1; that set is 5.17× larger than the gate set by construction and is a
  deliberately conservative proxy, not a defect. Pass 2 cleared it to the point where
  the gate is exactly satisfied.
* Pass 1 terminated on rim-frozen terminal edges (patch hops=3), which is why pass 2
  exists. Pass 2 ran with `froze rim 0` throughout.
* Median eta is 0.9% lower than the parent's. Every worst-case shape metric is
  unchanged.

## 8. Reproduce

```bash
CVM=<deck>/safs_material_cvm.nc
python code/census_1Hz_p5.py  --mesh <parent>.puml.h5 --cvm $CVM --gate 0.8 --out results/census_round0
python code/leb_refine_1Hz.py --mesh <parent>.puml.h5 --cvm $CVM --out pass1.puml.h5 \
       --gate 0.8 --hops 3 --max-rounds 80 --stats results/leb_stats.json
python code/leb_refine_1Hz.py --mesh pass1.puml.h5    --cvm $CVM --out final.puml.h5 \
       --gate 0.8 --hops 8 --max-rounds 120 --stats results/leb_stats_p2.json
python code/census_1Hz_p5.py       --mesh final.puml.h5 --cvm $CVM --gate 0.8 --out results/census_final
python code/check_fault_identity.py --parent <parent>.puml.h5 --new final.puml.h5
```

Knobs: `--gate` (0.8 = 1 Hz @ p5; 0.6667 = 0.5 Hz @ p3), `--hops` (patch radius —
raise if the log reports terminal edges frozen at the rim), `--max-rounds`.
