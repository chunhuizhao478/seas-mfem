# PREFERRED 1 Hz @ p5 mesh — `safpref_fb200_refine2_gatefix_1Hz_p5.puml.h5`

Built 2026-08-01 from `safpref_fb200_refine2_gatefix.puml.h5` by conforming 3-D
longest-edge (Rivara) bisection with a z-pooled Vs target — skill architecture 4,
the same recipe as the ALT `meshing_deep40km_1Hz_p5_leb` build.
**Gate PASS: 0 of 124,232,100 cells below the target; resolves 1.0000 Hz at p5
everywhere.**

## 1. Result

| | parent (fb200 ref2 gatefix) | **this mesh** |
|---|---:|---:|
| tets | 110,736,175 | **124,232,100** (+13,495,925, **+12.19 %**) |
| vertices | 18,924,562 | 21,504,494 |
| cells below the gate (Vs/dx < 0.8) | **428,060** | **0** |
| worst Vs/dx | 0.6667 | **0.8000** |
| **resolved everywhere @ p5** | 0.8334 Hz | **1.0000 Hz** |
| **resolved everywhere @ p3** | 0.5000 Hz | **0.6000 Hz** |

    5,485,395,952 bytes   md5 963e0a2e2669543f5b94b46154383a10

Convention (locked): `f = Vs/dx`, `dx` = element MAX edge, `Vs` NEAREST-GRID at the
element BARYCENTER; resolved `f` at order p is `(p/4)·Vs/dx`, so 1 Hz @ p5 ⇔
`Vs/dx ≥ 0.8`. Verified by `census_1Hz_p5.py` run independently on the finished
file, not by the builder's own check.

## 2. Quality — no regression that matters

| metric | parent | new | note |
|---|---:|---:|---|
| eta_min | 0.0306865 | 0.0306865 | **identical** |
| eta < 0.05 | 8 | 8 | **identical — no new slivers** |
| eta < 0.1 | 41 | 41 | **identical** |
| eta < 0.3 | 1,988 | 2,083 | +95 (+4.78 %) |
| eta median | 0.853546 | 0.843698 | **−1.15 % — the only real regression** |
| min edge | 10.5805 m | 10.5805 m | **identical** |
| inverted tets | 0 | 0 | identical |
| edges < 50 m | 2,014,425 | 2,014,735 | +310 |
| edges < 100 m | 63,847,874 | 64,093,465 | +245,591 |
| edges < 150 m | 84,701,328 | 86,205,155 | +1,503,827 |

Every worst-case shape metric is **bit-identical to the parent**. The extra short
edges are the intended product.

## 3. Fault identity, PUML contract and the TOPOGRAPHIC lid — all PASS

```
[PASS] F1 fault triangle multiset identical   5,522,976 vs 5,522,976
[PASS] F2 fault area identical  15925.102707625 km^2  (delta 0.000e+00 m^2)
[PASS] F3 fault faces interior (x2)           2,761,488 unique fault faces
[PASS] F3 boundary faces appear once
[PASS] F4a free-surface area unchanged  101193.090767938 km^2 (rel 1.508e-16)
[PASS] F4b free-surface z-range unchanged  z in [-49.923296, 25.818675]
       (TOPOGRAPHIC lid), 975,156 verts
[PASS] F5 no inverted tets
```

**The PREFERRED lid is topographic** (z up to +25.8 m), unlike the flat ALT lid, so
"free surface exactly flat" is the wrong invariant here and `check_fault_identity.py`
was extended. LEB inserts edge MIDPOINTS, and the midpoint of an edge of a planar
triangle lies ON that triangle — so refining the lid must leave its total AREA and
z-range bit-unchanged. Area is the sharp test (any motion off the surface changes
it) and it agrees to **1.5e-16 relative**, i.e. machine epsilon.

Fault edges were frozen throughout, so every fault-referenced deck input transfers
unchanged; the fault area delta is exactly zero.

## 4. Deck compatibility (Stage F)

Against `..._PREFERRED_THERMAL_CASE2_fb200ref2gatefix_plast_phi30_40_fw15pregate_nw08_k1p75_mw783_attenuation`:

| check | result |
|---|---|
| DR facets | 5,522,976 — **identical** to the parent |
| PGV receivers located | **93,700 / 93,700** |
| receivers below their LOCAL lid | **93,700 / 93,700**, min clearance **0.290 m** |
| pickpoints on the fault | 9 / 9, max plane distance 6.5e-05 m |
| domain extent | unchanged → `OutputRegionBounds` transfers |

Receiver containment used barycentric containment under the **local** top triangle
(`check_receivers_local_top.py`), not a 2-D footprint test — mandatory on a
topographic lid, because SeisSol v1.1.3 silently DROPS receivers above their local
free surface. The min clearance of **0.290 m reproduces the parent's documented
value exactly**, which is independent confirmation that the lid did not move.

## 5. Run cost — essentially free

| | parent | new | ratio |
|---|---:|---:|---:|
| dt_min | 1.61293e−05 s | 1.61293e−05 s | **1.000** |
| LTS clusters | 11 | 11 | 1.000 |
| LTS updates / simulated s | 3.04104e+11 | 3.12404e+11 | **1.0273** |

**+12.2 % tets costs only +2.7 % run time** — the added cells are shallow far-field
and join COARSE LTS clusters, and `dt_min` is unchanged to the last digit. Judge this
campaign on the LTS delta, not the tet delta.

## 6. ⚠ NODE SIZING AT p5 — AN ARITHMETIC INCONSISTENCY TO RESOLVE BEFORE RUNNING

A sizing note in circulation reasons that *"ORDER 5 has 35 basis functions/cell vs
ORDER 4's 20, so memory per cell is 1.75×"* and calls 35 the "p5-equivalent",
concluding 1024 nodes on `-p large`. **35 basis functions is polynomial degree 4, not
degree 5.** In 3-D the count is `(p+1)(p+2)(p+3)/6`:

| SeisSol ORDER | degree p | basis fns | memory vs ORDER 4 |
|---|---|---:|---:|
| 4 | p3 | 20 | 1.00× |
| 5 | p4 | 35 | 1.75× |
| **6** | **p5** | **56** | **2.80×** |

These decks state "ORDER 4 = p3", so a genuine **p5** run is **ORDER 6** and costs
**2.80×** per cell, not 1.75×. If the 1024-node figure was derived from the 1.75×
factor it is **~1.6× too optimistic** for a true p5 run. This mesh exists to be run at
p5, so resolve which order is actually intended before committing SU — the two
blockers already noted (no ORDER-5/6 ParMETIS binary, and unmeasured cost) apply
either way. This is a flag, not a measurement: confirm with a short probe run.

## 7. How it was built

Two passes of `leb_refine_1Hz.py`:

| pass | hops | rounds | wall | tets added | gate after |
|---|---|---|---|---|---|
| 1 | 3 | 73 (stopped: all terminal edges rim-frozen) | 96 min | +12,869,257 | 9,548 |
| 2 | 8 | 50 | 14 min | +626,668 | **0** |

Peak RSS 16.6 GB (pass 1), 12.1 GB (pass 2).

The census classified the 428,060 failures as `fault-edge exempt 0 |
fault-vertex pinned 9 | free 428,051`. The 9 fault-vertex-pinned cells are **not** a
structural floor under this method: the refiner freezes fault EDGES, not fault
VERTICES, and bisecting a non-fault edge that merely touches a fault vertex does not
move the fault surface — all 9 were refined and the gate closed to zero.

Everything else about the method — why bisection rather than mmg, why the refinement
target must not be the gate, why pooling is z-only, and the strict edge-key LEPP
tie-break — is documented in the ALT build's report
(`project_7.0_alternative/meshing_deep40km_1Hz_p5_leb/QUALITY_REPORT_1Hz_p5.md`) and
in skill architecture 4. Nothing new failed on this lineage.

## 8. Reproduce

```bash
CVM=<deck>/safs_material_cvm.nc
P=<...>/meshing_faultband_refine2/results/safpref_fb200_refine2_gatefix.puml.h5

python code/census_1Hz_p5.py  --mesh $P --cvm $CVM --gate 0.8 --out results/census_round0
python code/leb_refine_1Hz.py --mesh $P --cvm $CVM --out results/pass1.puml.h5 \
       --gate 0.8 --hops 3 --max-rounds 80  --stats results/leb_stats_p1.json
python code/leb_refine_1Hz.py --mesh results/pass1.puml.h5 --cvm $CVM \
       --out results/safpref_fb200_refine2_gatefix_1Hz_p5.puml.h5 \
       --gate 0.8 --hops 8 --max-rounds 140 --stats results/leb_stats_p2.json
python code/census_1Hz_p5.py           --mesh results/...1Hz_p5.puml.h5 --cvm $CVM --gate 0.8 --out results/census_final
python code/check_fault_identity.py    --parent $P --new results/...1Hz_p5.puml.h5
python code/check_receivers_local_top.py --mesh results/...1Hz_p5.puml.h5 \
       --receivers <deck>/safs_pgv_receivers_preferred_100k.dat
```
